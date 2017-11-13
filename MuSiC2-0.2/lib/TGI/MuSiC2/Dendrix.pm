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
    $this->{_MIN_FREQ_GENE} = undef;
    $this->{_NUMBER_ITERATIONS} = undef;
    $this->{_ANALYZED_GENES_FILE} = undef;
    $this->{_NUM_EXPER} = undef;
    $this->{_STEP_LENGTH} = undef;

    bless $this, $class;
    $this->process();

    return $this;
}

sub process {
    my $this = shift;
    # do regular Dendrix
    my $args; 
    for my $arg (@ARGV) { $args .= $arg . ' '; }
    my ( $input, $output, $err );
    use Symbol 'gensym'; 
    $err = gensym;
    my $pid = open3( $input, $output, $err, 'Dendrix.py $args' ); # puts output in the specified text files
    #mutations_file K minFreqGene number_iterations analyzed_genes_file num_exper step_length' )
    waitpid( $pid, 0 );
    my $child_exit_status = $? >> 8;
    if ( $output ) { print 'Output: ' . $output . "\n"; }
    if ( $err ) { print 'Errors: ' . $err; }
}


sub help_text {
    my $this = shift;
    return <<HELP
Sub-commands for music2 dendrix:
p-test			Runs a permutation test for Dendrix and calculates the corresponding p-value.

HELP
}

#$pid2 = open3($input2, $output2, $err2, 'python permutationTestDendrix.py mutations_file K minFreqGene number_iterations analyzed_genes_file num_permutations value_tested rank' )



1;

