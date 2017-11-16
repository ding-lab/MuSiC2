package TGI::MuSiC2::Dendrix;
##

use warnings;
use strict;

use IO::File;
use IPC::Open3;

# Using Dendrix.py from Dendrix v 0.3

sub new {
    my $class = shift;
    my $this = {};

    $this->{_MUTATIONS_FILE} = undef;
    $this->{_K} = 2;
    $this->{_MIN_FREQ_GENE} = 1;
    $this->{_NUMBER_ITERATIONS} = 1000000;
    $this->{_ANALYZED_GENES_FILE} = undef;
    $this->{_NUM_EXPER} = 1;
    $this->{_STEP_LENGTH} = 1000;

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
        'number-experiments=i'       => \$this->{_NUM_EXPER},
        'step-length=i'              => \$this->{_STEP_LENGTH},
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
    $args = $args . $this->{_NUM_EXPER} . ' ' . $this->{_STEP_LENGTH};
    my ( $input, $output, $err );
    use Symbol 'gensym'; 
    $err = gensym;
    my $pid = open3( $input, $output, $err, 'dendrix $args' ); 
    # puts output in the specified text files
    #mutations_file K minFreqGene number_iterations analyzed_genes_file num_exper step_length' )
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
 music2 dendrix --mutations-file=? --set-size=? minimum-freq=? 
        --number-interations=? --analyzed-genes-file=? 
        --number-experiments=? --step-length=? 

SYNOPSIS
 music2 dendrix \
      --mutations-file example/dendrix/mutation_matrix \
      --set-size 2 \
      --minimum-freq 1 \
      --number-interations 1000000 \
      --analyzed-genes-file example/dendrix/analyzed_genes \
      --number-experiments 1 \
      --step-length 1000 \

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
  number-experiments
    number of times the experiment is going to be run. For each experiment, a 
    random solution is used as initial state for the MCMC
  step-length
    number of iterations of the MCMC between two samples

DESCRIPTION
    This script runs Dendrix. Dendrix (De _no_vo _Dri_ver E_x_clusivity) is an algorithm for 
    discovery of mutated driver pathways in cancer using only mutation data. It finds sets of 
    genes, domains, or nucleotides whose mutations exhibit both high coverage and high 
    exclusivity in the analyzed samples.

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

    --number-experiments

      Number of times the experiment is going to be run. For each experiment, a 
      random solution is used as initial state for the MCMC

    --step-length
    
      Number of iterations of the MCMC between two samples

OUTPUT

     sets_frequencyOrder_experiment#j.txt: for each experiment, the 1000 sets sampled with highest frequency by 
     the MCMC. One set per line, and the corresponding weight is reported. #j denote the number 
     of the experiment (starting from 0).
     sets_weightOrder_experiment#j.txt: for each experiment, the 1000 sets of highest weight sampled by the MCMC. 
     One set per line, and the corresponding weight is reported. #j denote the number
     of the experiment (starting from 0).

EXAMPLE

    music2 dendrix --mutations-file example/dendrix/mutation_matrix --set-size 2 --minimum-freq 1 \
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
 
HELP

}

1;

