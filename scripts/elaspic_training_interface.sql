drop database if exists elaspic_training_interface;
create database elaspic_training_interface;
use elaspic_training_interface;

# Create views
create view elaspic_training_interface.domain as select * from elaspic.domain;
create view elaspic_training_interface.domain_contact as select * from elaspic.domain_contact;
create view elaspic_training_interface.provean as select * from elaspic.provean;


# Create tables
create table elaspic_training_interface.uniprot_domain like
    elaspic.uniprot_domain;
create table elaspic_training_interface.uniprot_domain_template like
    elaspic.uniprot_domain_template;
create table elaspic_training_interface.uniprot_domain_model like
    elaspic.uniprot_domain_model;
create table elaspic_training_interface.uniprot_domain_mutation like
    elaspic.uniprot_domain_mutation;

create table elaspic_training_interface.uniprot_domain_pair like
    elaspic.uniprot_domain_pair;
create table elaspic_training_interface.uniprot_domain_pair_template like
    elaspic.uniprot_domain_pair_template;
create table elaspic_training_interface.uniprot_domain_pair_model like
    elaspic.uniprot_domain_pair_model;
create table elaspic_training_interface.uniprot_domain_pair_mutation like
    elaspic.uniprot_domain_pair_mutation;
