# '{{source_table}}' should contain 'uniprot_id' and 'uniprot_id_2' columns

# Populate tables
insert into elaspic_training_interface.uniprot_domain
select * from elaspic.uniprot_domain
where uniprot_id in (select uniprot_id from {{source_table}})
or uniprot_id in (select uniprot_id_2 from {{source_table}});


insert into elaspic_training_interface.uniprot_domain_template
select * from elaspic.uniprot_domain_template
where uniprot_domain_id in (
    select uniprot_domain_id from elaspic_training_interface.uniprot_domain);


insert into elaspic_training_interface.uniprot_domain_model
select * from elaspic.uniprot_domain_model
where uniprot_domain_id in (
    select uniprot_domain_id from elaspic_training_interface.uniprot_domain);


insert into elaspic_training_interface.uniprot_domain_pair (uniprot_domain_id_1, uniprot_domain_id_2, uniprot_id_1, uniprot_id_2)
select distinct
ud1.uniprot_domain_id uniprot_domain_id_1,
ud2.uniprot_domain_id uniprot_domain_id_2,
ud1.uniprot_id uniprot_id_1,
ud2.uniprot_id uniprot_id_2
from elaspic_training_interface.uniprot_domain ud1
join elaspic_training_interface.uniprot_domain ud2
join {{source_table}} i ON (
    i.uniprot_id = ud1.uniprot_id and i.uniprot_id_2 = ud2.uniprot_id)
where ud1.uniprot_domain_id <= ud2.uniprot_domain_id
union
select distinct
ud2.uniprot_domain_id uniprot_domain_id_1,
ud1.uniprot_domain_id uniprot_domain_id_2,
ud2.uniprot_id uniprot_id_1,
ud1.uniprot_id uniprot_id_2
from elaspic_training_interface.uniprot_domain ud1
join elaspic_training_interface.uniprot_domain ud2
join {{source_table}} i ON (
    i.uniprot_id = ud1.uniprot_id and i.uniprot_id_2 = ud2.uniprot_id)
where ud1.uniprot_domain_id > ud2.uniprot_domain_id
;
