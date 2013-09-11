MuSiC2
===========
Mutational Significance in Cancer (Cancer Mutation Analysis) version 2.

Usage
-----

    Program:     music2 - Mutational Significance in Cancer (Cancer Mutation Analysis) version 2.
    Version:     V0.1
    Author:      Beifang Niu && Cyriac Kandoth

    Usage:  music2 <command> [options]

Key commands:

    bmr                    ...  Calculate gene coverages and background mutation rates.
    clinical-correlation        Correlate phenotypic traits against mutated genes, or       
                                 against individual variants.
    cosmic                      Match a list of variants to those in COSMIC, and highlight  
                                 druggable targets.
    cosmic-omim                 Compare the amino acid changes of supplied mutations to
                                 COSMIC and OMIM databases.
    create-visualizations       no description!!!: define doc in the class definition for
                                 Genome::Model::Tools::Music::CreateVisualizations
    data                   ...  Parse, standardize, and index various third-party datasets  
                                 used by the tools in MuSiC.
    mutation-relation           Identify relationships of mutation concurrency or mutual    
                                 exclusivity in genes across cases.
    path-scan                   Find signifcantly mutated pathways in a cohort given a list 
                                 of somatic mutations.
    pfam                        Add Pfam annotation to a MAF file.
    play                        Run the full suite of MuSiC tools sequentially.
    plot                   ...  Generate relevant plots and visualizations for MuSiC2.
    proximity                   Perform a proximity analysis on a list of mutations.
    proximity-window            Perform a sliding window proximity analysis on a list of mutations.
    smg                         Identify significantly mutated genes.
    survival                    Create survival plots and P-values for clinical and mutational phenotypes.      

SUPPORT

For user support please mail bniu@genome.wustl.edu


Install
-------

        Just run: cpanm MuSiC2-0.1.tar.gz

xxx

