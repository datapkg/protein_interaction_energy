import datetime
import logging
import numpy as np
import pandas as pd
import sqlalchemy as sa
from IPython.display import display, HTML
import kmtools as ascommon


logger = logging.getLogger(__name__)

# #################################################################################################
# Blast
def annotate_blast_results(result_df, domain_start, domain_sequence_length):
    """Add additional information to blast results.

    .. note::

        This code is used for parsing Profs libraries, and may have to be adapted
        for for libraries following a different format.
    """
    result_df = result_df.copy()

    if result_df.empty:
        return result_df

    if result_df['subject_id'].values[0].count('|') == 2:
        HEADER_TYPE = 1
    elif result_df['subject_id'].values[0].count('|') == 5:
        HEADER_TYPE = 2
    else:
        raise Exception

    # Parse blast results
    result_df['pdb_id'] = result_df['subject_id'].apply(lambda x: x.split('|')[0].split('_')[0])
    result_df['pdb_chain'] = result_df['subject_id'].apply(lambda x: x.split('|')[0].split('_')[1])
    result_df['pdb_pdbfam_name'] = result_df['subject_id'].apply(lambda x: x.split('|')[1])
    result_df['pdb_pdbfam_idx'] = result_df['subject_id'].apply(lambda x: int(x.split('|')[2]))
    if HEADER_TYPE == 2:
        result_df['pdb_pfam_clan'] = result_df['subject_id'].apply(lambda x: x.split('|')[3])
        result_df['pdb_domain_def'] = result_df['subject_id'].apply(lambda x: x.split('|')[4])
        result_df['pdb_cath_id'] = result_df['subject_id'].apply(lambda x: x.split('|')[5])

    # Score sequence identity
    alpha = 0.95
    result_df['alignment_identity'] = result_df['pc_identity']
    result_df['alignment_coverage'] = (
        (result_df['q_end'] - result_df['q_start'] + 1) / float(domain_sequence_length) * 100.0
    )
    result_df['alignment_score'] = (
        alpha * (result_df['alignment_identity'].values / 100) *
        (result_df['alignment_coverage'].values / 100) +
        (1 - alpha) * (result_df['alignment_coverage'].values / 100)
    )

    # Domain definitions
    result_df['domain_start_new'] = domain_start + result_df['q_start'].astype(int) - 1
    result_df['domain_end_new'] = domain_start + result_df['q_end'].astype(int) - 1
    result_df['domain_def_new'] = (
        result_df['domain_start_new'].astype(str) + ':' +
        result_df['domain_end_new'].astype(str)
    )
    result_df['t_date_modified'] = datetime.datetime.now()
    return result_df


# #################################################################################################
# Taipale
def load_taipale_mmc3(table_name):
    """Load study data from excel files."""
    mmc3_df = pd.read_excel('../downloads/taipale/mmc3.xlsx', table_name)
    print(table_name, ':', mmc3_df.shape)

    # Get `refseq_base_id` for entries with mutations
    (mmc3_df['refseq_base_id'],
     mmc3_df['refseq_mutation'],
     mmc3_df['refseq_mutation_pos']) = (
        zip(*mmc3_df['Mutation_RefSeq_AA'].apply(ascommon.sequence_tools.format_hgvs_mutation))
    )

    return mmc3_df


def get_conversion_tables(mmc3_table1_df):
    """.

    .. note::
        Obsolete. No longer used!

    Returns
    -------
    gene_symbol_to_refseq : DataFrame
        `Symbol` (gene name) column mapped to RefSeq ID (using `Mutation_RefSeq_AA` column).
    gene_id_to_uniprot : DataFrame
        `Entrez_Gene_ID` mapped to UniProt ACC using the UniProt KB identifier table.
    """
    """Create a table for mapping `Symbol` -> `refseq_base_id` for entries without mutations."""
    gene_symbol_to_refseq = (
        mmc3_table1_df[['Symbol', 'refseq_base_id']]
        .rename(columns={'refseq_base_id': 'refseq_base'})
        .dropna()
        .drop_duplicates()
    )
    gene_symbol_to_refseq.head()

    # Create a table for mapping `Interactor_Gene_ID` -> `uniprot_id`
    # for entries without mutations (interaction partners)
    engine = sa.create_engine('mysql://strokach:@192.168.6.19:3306/uniprot_kb')
    sql_query_template = """
select identifier_id, uniprot_id
from uniprot_kb.uniprot_identifier
where identifier_id in ('{}')
and identifier_type = 'GeneID';
"""
    sql_query = sql_query_template.format("', '".join(
        set(str(x) for x in mmc3_table1_df['Entrez_Gene_ID'].values) |
        set(str(x) for x in mmc3_table1_df['Interactor_Gene_ID'].values)
    ))
    gene_id_to_uniprot = pd.read_sql_query(sql_query, engine)
    gene_id_to_uniprot['identifier_id'] = gene_id_to_uniprot['identifier_id'].astype(int)

    return gene_symbol_to_refseq, gene_id_to_uniprot


def convert_interacting_partner(mmc3_df, gene_symbol_to_refseq, gene_id_to_uniprot):
    """Convert refseq id to uniprot id."""
    # Get `refseq_base_id_1` and `refseq_base_id_2` for all entries
    mmc3_df = (
        mmc3_df
        .merge(
            gene_symbol_to_refseq.rename(
                columns={'refseq_base': 'refseq_base_id_1'}),
            on=['Symbol'],
            how='left')
        # .merge(
        #     gene_id_to_uniprot.rename(
        #         columns={'identifier_id': 'Entrez_Gene_ID', 'uniprot_id': 'uniprot_id_1'}),
        #     on=['Entrez_Gene_ID'],
        #     how='left')
        .merge(
            gene_id_to_uniprot.rename(
                columns={'identifier_id': 'Interactor_Gene_ID', 'uniprot_id': 'uniprot_id_2'}),
            on=['Interactor_Gene_ID'],
            how='left')
    )
    _print_stats('After merging:', mmc3_df)

    # Remove `refseq_base_id_1` which does not match the mutation, when a mutation is present
    mmc3_df = mmc3_df.drop(
        mmc3_df[
            mmc3_df[['refseq_base_id', 'refseq_base_id_1']]
            .apply(lambda x: pd.notnull(x[0]) and x[0] != x[1], axis=1)
        ].index
    )

    # Validate the results
    if mmc3_df['refseq_base_id_1'].isnull().any():
        missing_refseq_base_id = (
            (mmc3_df['refseq_base_id_1'].isnull()) &
            (~mmc3_df['Category'].isin(['PRS', 'RRS']))
        )
        error_message = """\
We have {} rows from the Positive Reference Set (PRS) and Random Reference Set (RRS) \
without a `refseq_base_id_1` entry.\
""".format(missing_refseq_base_id.sum())
        assert not missing_refseq_base_id.any(), error_message
    assert (
        mmc3_df[['refseq_base_id', 'refseq_base_id_1']]
        .apply(lambda x: pd.isnull(x[0]) or x[0] == x[1], axis=1)
        .all()
    )

    _print_stats('After processing:', mmc3_df)
    display(mmc3_df.head(4))
    return mmc3_df


# === Helper Functions ===
def _print_stats(timepoint, mmc3_df):
    """Helper function printing statistics about the data being analyzed."""
    unique_cols = ['Symbol', 'Mutation_RefSeq_AA', 'Interactor_Gene_ID']
    display(HTML("<h4>" + timepoint))
    ascommon.df_tools.print2(
        "All rows:",
        mmc3_df.shape)
    ascommon.df_tools.print2(
        "Unique interfaces / mutations:",
        mmc3_df[unique_cols].dropna(how='any').drop_duplicates().shape, '<--')
    ascommon.df_tools.print2(
        "Unique interfaces / mutations mapped to uniprot:",
        mmc3_df[
            mmc3_df[['Entrez_Gene_ID', 'Interactor_Gene_ID']]
            .notnull().all(axis=1)
        ][unique_cols].dropna(how='any').drop_duplicates().shape)
    ascommon.df_tools.print2(
        "Unique interfaces / mutations mapped to uniprot and refseq:",
        mmc3_df[
            mmc3_df[['refseq_base_id_1', 'Entrez_Gene_ID', 'Interactor_Gene_ID']]
            .notnull().all(axis=1)
        ][unique_cols].dropna(how='any').drop_duplicates().shape)


# === SIFTS ===

def get_pfam_info(sifts_df, pdb_chain, pdb_mutation):
    pdb_aa = pdb_mutation[0]
    pdb_resnum = pdb_mutation[1:-1]
    sifts_df_select = (
        sifts_df[
            (sifts_df['pdb_chain'] == pdb_chain) &
            (sifts_df['pdb_aa'] == pdb_aa) &
            (sifts_df['resnum'] == pdb_resnum)
        ]
    )
    try:
        uniprot_id = sifts_df_select.iloc[0]['uniprot_id']
    except (KeyError, IndexError):
        logger.info("No 'uniprot_id' data for mutation '{}'!".format(pdb_mutation))
        uniprot_id = np.nan
    try:
        pfam_id = sifts_df_select.iloc[0]['pfam_id']
    except (KeyError, IndexError):
        logger.info("No 'pfam_id' data for mutation '{}'!".format(pdb_mutation))
        pfam_id = np.nan
    return uniprot_id, pfam_id
