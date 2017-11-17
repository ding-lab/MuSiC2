package TGI::MuSiC2::DendrixPermutation;
##

use warnings;
use strict;

use IO::File;
use IPC::Open3;

# Using permutationTestDendrix.py from Dendrix v 0.3

sub new {
    my $class = shift;
    my $this = {};

    $this->{_MUTATIONS_FILE} = undef;
    $this->{_K} = 2;
    $this->{_MIN_FREQ_GENE} = 1;
    $this->{_NUMBER_ITERATIONS} = 1000000;
    $this->{_ANALYZED_GENES_FILE} = undef;
    $this->{_NUMBER_PERMUTATIONS} = 100;
    $this->{_VALUE_TESTED} = 48;
    $this->{_RANK} = 1;

    bless $this, $class;
    $this->process();

    return $this;
}


sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'mutations-file=s'           => \$this->{_MUTATIONS_FILE},
        'set-size=i'                 => \$this->{_K},
        'minimum-freq=i'             => \$this->{_MIN_FREQ_GENE},
        'number-interations=i'       => \$this->{_NUMBER_ITERATIONS},
        'analyzed-genes-file=s'      => \$this->{_ANALYZED_GENES_FILE},
        'number-permutationss=i'     => \$this->{_NUMBER_PERMUTATIONS},
        'value-tested=i'             => \$this->{_VALUE_TESTED},
        'rank=i'                     => \$this->{_RANK},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    #### processing ####
    # Check on all the input data before starting work
    print STDERR "Mutation matrix file not found or is empty: $this->{_MUTATIONS_FILE}\n" unless( -s $this->{_MUTATIONS_FILE} );
    print STDERR "Analyzed gene list file not found or is empty: $this->{_ANALYZED_GENES_FILE}\n" unless( -s $this->{_ANALYZED_GENES_FILE} );
    # Collect args
    my $args = $this->{_MUTATIONS_FILE} . ' ' . $this->{_K} . ' ' . $this->{_MIN_FREQ_GENE};
    $args = $args . $this->{_NUMBER_ITERATIONS} . ' ' . $this->{_ANALYZED_GENES_FILE} . ' ';
    $args = $args . $this->{_NUMBER_PERMUTATIONS} . ' ' . $this->{_VALUE_TESTED} . ' ' . $this->{_RANK};
    my ( $input, $output, $err );
    use Symbol 'gensym'; 
    $err = gensym;
    my $pid = open3( $input, $output, $err, 'permutationTestDendrix $args' ); 
    # puts output in the specified text files
    # python PermutationTestDendrix.py example/mutation_matrix 3 1 1000000 example/analyzed_genes 100 48 1
    waitpid( $pid, 0 );
    my $child_exit_status = $? >> 8;
    if ( $output ) { print "Output: \n"; while (<$output>) { print; }; print "\n"; };
    if ( $err ) { print "Errors: \n"; while (<$err>) { print; }; print "\n"; };
   
    return 1;

}


## usage
sub help_text {
    my $this = shift;
        return <<HELP

USAGE
 music2 dendrix-permutation --mutations-file=? --set-size=? minimum-freq=? 
        --number-interations=? --analyzed-genes-file=? 
        --number-permutations=? --value-tested=? --rank=? 

SYNOPSIS
 music2 dendrix-permutation \
      --mutations-file example/dendrix/mutation_matrix \
      --set-size 2 \
      --minimum-freq 1 \
      --number-interations 1000000 \
      --analyzed-genes-file example/dendrix/analyzed_genes \
      --number-permutations 100 \
      --value-tested 48
      --rank 1 \

(A "mutations-file" can be generated using the raw maf file.)

REQUIRED INPUTS
  mutations-file
    Input file with mutation matrix (see example/dendrix/mutation_matrix format)
  set-size
    Size of the sets to be sampled
  minimum-freq
    minimum frequency of mutation for a gene/group to be considered in the analysis
  number-interations
    number of iterations of the MCMC
  analyzed-genes-file
    file with list of analyzed genes, one per line
  number-permutations
    number of times the permuted datasets to consider in the permutation test
  value-tested
    value of the weight for which the permutation test will be run
  rank
    rank of the weight tested (value_tested) in the dendrix output: sets_weightOrder_experiment.txt file

DESCRIPTION
    This script runs the permutation test for Dendrix output. The input for this script
    is from Dendrix output file (see EXAMPLE) .

ARGUMENTS

    --mutations-file

      Input file with mutation matrix (see example/dendrix/mutation_matrix format)

    --set-size

      Size of the sets to be sampled 

    --minimum-freq

      Minimum frequency of mutation for a gene/group to be considered in the analysis

    --number-interations

      Number of iterations of the MCMC
      
    --analyzed-genes-file

      File with list of analyzed genes, one per line

    --number-permutations
      
      Number of times the permuted datasets to consider in the permutation test

    --value-tested
      
      Value of the weight for which the permutation test will be run

    --rank

      Rank of the weight tested (value_tested) in the dendrix output: sets_weightOrder_experiment.txt file

OUTPUT

    p_value_dendrix.txt: contains the p-value of the permutation test.

EXAMPLE

    music2 dendrix --mutations-file example/dendrix/mutation_matrix --set-size 3 --minimum-freq 1 \
      --number-interations 1000000 --analyzed-genes-file example/dendrix/analyzed_genes \
      --number-experiments 1 --step-length 1000

    Runs the MCMC for 1000000 iterations, sampling sets of size 3 every 1000
    iterations. Produces two files  (since 1 experiment is run):

    sets_frequencyOrder_experiment0.txt 
    and
    sets_weightOrder_experiment0.txt. 
    
    The first lines of set_frequencyOrder_experiment0.txt look like the following:

    13	FGFR2	PTEN	RB1	48
    10	FGFR2	MTAP	PTEN	47
    10	FGFR2	PTCH1	PTEN	47
    10	MTAP	PTEN	TRIM2	47
    10	MTAP	PDGFRA	PTEN	47

    In each line: the first number is the number of times the set has been
    sampled; the following K tokens are the genes in the set; the last number is
    the weight of the set.

    The first lines of set_wrightOrder_experiment0.txt file look like the following:

    48	FGFR2	PTEN	RB1	13
    47	FGFR2	MTAP	PTEN	10
    47	FGFR2	PTCH1	PTEN	10
    47	MTAP	PTEN	TRIM2	10
    47	MTAP	PDGFRA	PTEN	10

    The format is as the set_frequencyOrder_experiment0.txt file, but with the
    first and last token switched.

    If you want to compute the p-value for the first set having weight 48, you can run:
      
    music2 dendrix-permutation --mutations-file example/dendrix/mutation_matrix --set-size 2 --minimum-freq 1 \
      --number-interations 1000000 --analyzed-genes-file example/dendrix/analyzed_genes \
      --number-permutations 100 --value-tested 48 --rank 1

    If instead you want to compute the p-value for the second set having weight 47, you can run:
    
    music2 dendrix-permutation --mutations-file example/dendrix/mutation_matrix --set-size 2 --minimum-freq 1 \
      --number-interations 1000000 --analyzed-genes-file example/dendrix/analyzed_genes \
      --number-permutations 100 --value-tested 47 --rank 2

    The (empirical) p-value will be found in the file p_value_dendrix.txt.

NOTES:

    - what we call "genes" can be any group of mutations
    - the analyzed_genes_file can be used to remove genes from the analysis without
      the need of modifying the mutation matrix.
    - example/ contains the small example described above, and two mutations tables assembled from real data:
      
      a) GBM_mutationsAndCN: assembled from somatic mutations (single
      nucleotide mutations and small indels) and focal copy number aberrations
      described in: The Cancer Genome Atlas Network. (2008) Comprehensive genomic
      characterization defines human glioblastoma genes and core pathways.
      Nature 455, 1061-1068. (The 7 hypermutated samples have been removed from
      the mutation table.)
      
      b) Lung_mutations: assembled from somatic mutations (single nucleotide
      mutations and small indels) described in: L. Ding et al. (2008) Somatic
      mutations affect key pathways in lung adenocarcinoma. Nature 455, 1069-1075.


REFERENCES:

    If you use Dendrix in your research, please cite:

    F. Vandin, E. Upfal, and B.J. Raphael. (2012) De novo Discovery of Mutated
    Driver Pathways in Cancer. Genome Research, 22, 2:375-85.

WEBSITE:

    http://cs.brown.edu/people/braphael/software.html

HELP

}

1;

