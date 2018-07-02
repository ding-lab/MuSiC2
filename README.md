MuSiC2
===========
Mutational Significance in Cancer (Cancer Mutation Analysis) version 2.

Usage
-----

    Program:     music2 - Mutational Significance in Cancer (Cancer Mutation Analysis) version 2.
    Version:     V0.2
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


Install (Ubuntu 14.04.01)
-------

Prerequisites:

Install samtools ( Download the samtools-0.1.19 from SOURCEFORGE (http://sourceforge.net/projects/samtools/files/samtools/0.1.19) )

        tar jxf samtools-0.1.19.tar.bz2
        cd samtools-0.1.19
        make
        export SAMDIR=$PWD
        sudo mv samtools /usr/local/bin/
      
Install joinx 

        git clone --recursive https://github.com/genome/joinx.git
        sudo apt-get install build-essential cmake libbz2-dev libgtest-dev
        cd joinx
        mkdir build
        cd build
        cmake .. -DCMAKE_BUILD_TYPE=release
        make deps
        make
        sudo make install

Install calcRoiCovg 

        sudo apt-get install git libbam-dev zlib1g-dev
        git clone https://github.com/Beifang/calcRoiCovg.git
        cd calc-roi-covg
        make
        sudo mv calcRoiCovg /usr/local/bin/

Install bedtools 

        curl http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz > BEDTools.tar.gz
        tar -zxvf BEDTools.tar.gz
        cd BEDTools-v2.17.0
        make
        sudo mv ./bin /usr/local/bin/

In order to install MuSiC2 package, we need CPANM program
(cpanm - get, unpack build and install modules from CPANM)

        sudo apt-get install cpanminus

Intall Perl5 local lib

        cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

Intall Test::Most module
        
        wget http://search.cpan.org/CPAN/authors/id/O/OV/OVID/Test-Most-0.34.tar.gz
        cpanm Test-Most-0.34.tar.gz

Intall Statistics::Descriptive module
        cpanm Statistics::Descriptive
        cpanm Statistics::Distributions

Install MuSiC2 package
        
        git clone https://github.com/ding-lab/MuSiC2
        cd MuSiC2
        cpanm MuSiC2-#.#.tar.gz

Notes: Python is needed to be installed if you run music2 dendrix & dendrix-permutation 

example
-------



SUPPORT
-------

If you have any questions, please contact one or more of the following folks:

Beifang Niu <bniu@sccas.cn>
Cyriac Kandoth <ckandoth@gmail.com>
Li Ding <lding@wustl.edu>

