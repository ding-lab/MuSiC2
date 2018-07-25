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


Install (Ubuntu & CentOS)
-------
Note: We provided binaries for joinx, samtools, calcRoiCovg and bedtools in /bin dir, and which were compiled on CentOS, and tested on CentOS/Ubuntu.

Prerequisites for Ubuntu:

        sudo apt-get install build-essential \
        git \
        cmake \
        curl \
        cpanminus
        libbz2-dev \
        libgtest-dev \
        libbam-dev \
        zlib1g-dev 

Prerequisites for CentOS:

        sudo yum install yum-utils
        sudo yum install curl
        sudo yum install git
        sudo yum install cmake
        sudo yum groupinstall "Development Tools"
        sudo yum update -y nss curl libcurl
        sudo yum install perl-devel
        sudo yum install perl-CPAN
        sudo yum install bzip2-libs
        sudo yum install zlib-devel
        sudo curl -L http://cpanmin.us | perl - --sudo App::cpanminus


Change C++11 compiler for CentOS (required for joinx installation)

   Reference 
> https://www.softwarecollections.org/en/scls/rhscl/devtoolset-3/ 

    1. Install a package with repository for your system:
    On CentOS, install package centos-release-scl available in CentOS repository:
        $ sudo yum install centos-release-scl
    On RHEL, enable RHSCL repository for you system:
        $ sudo yum-config-manager --enable rhel-server-rhscl-7-rpms
    2. Install the collection:
        $ sudo yum install devtoolset-3
    3. Start using software collections:
        $ scl enable devtoolset-3 bash
    Set env variables --optional
        CC=gcc CXX=g++ 

Install samtools ( Download the samtools-0.1.19 from SOURCEFORGE (http://sourceforge.net/projects/samtools/files/samtools/0.1.19) )

        tar jxf samtools-0.1.19.tar.bz2
        cd samtools-0.1.19
        make
        export SAMTOOLS_DIR=$PWD
        sudo mv samtools /usr/local/bin/

Install calcRoiCovg 

        git clone https://github.com/Beifang/calcRoiCovg.git
        cd calc-roi-covg
        make
        sudo mv calcRoiCovg /usr/local/bin/

Install bedtools 

        wget https://github.com/arq5x/bedtools2/archive/v2.27.1.tar.gz
        tar -zxvf v2.27.1.tar.gz
        cd bedtools2-2.27.1/
        make
        sudo mv ./bin /usr/local/bin/

Install joinx 

        git clone --recursive https://github.com/genome/joinx.git
        cd joinx
        mkdir build
        cd build
        cmake ..
        make deps
        make
        sudo make install

Fix joinx bugs

        StreamLineSource.cpp
        bool StreamLineSource::getline(std::string& line) {
            std::getline(_in, line);
            return true;
        }

Intall Perl modules

        sudo cpanm Test::Most 
        sudo cpanm Statistics::Descriptive
        sudo cpanm Statistics::Distributions
        sudo cpanm Bit::Vector

Install MuSiC2 package
        
        git clone https://github.com/ding-lab/MuSiC2
        cd MuSiC2
        sudo cpanm MuSiC2-#.#.tar.gz

Notes: Python is needed to be installed if you run music2 dendrix & dendrix-permutation 


example
-------



SUPPORT
-------

If you have any questions, please contact one or more of the following folks:

Beifang Niu <bniu@sccas.cn>
Cyriac Kandoth <ckandoth@gmail.com>
Li Ding <lding@wustl.edu>

