package TGI::MuSiC2::CalcWindowMaf;
##
# given gene based MAF file, calculate window based MAF file
##
use strict;
use warnings;
#
#
use IO::File;
use Getopt::Long;
use File::Temp qw/ tempfile /;

sub new {
    my $class = shift;
    my $this = {};

    $this->{_MAF_FILE} = undef;
    $this->{_OUTPUT_MAF_FILE}  = "window_based_maf";

    $this->{_WINDOW_SIZE} = 1000000;

    bless $this, $class;
    $this->process();

    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (

        'maf-file=s'               => \$this->{_MAF_FILE},
        'output-maf-file=s'        => \$this->{_OUTPUT_MAF_FILE},
        'window-size=i'            => \$this->{_WINDOW_SIZE},

        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    #### processing ####
    #
    # Check on all the input data
    print STDERR "MAF file not found or is empty: $this->{_MAF_FILE}\n" unless( -s $this->{_MAF_FILE} );
    return undef unless( -s $this->{_MAF_FILE} );
    unless ( $this->{_WINDOW_SIZE} =~ /\d+/ ) {  
        print STDERR "Window size format is not valid !\n"; 
        return undef; 
    };
    #
    #
    my $outputfh = IO::File->new( $this->{_OUTPUT_MAF_FILE}, ">" ) or die "Couldn't open $this->{_OUTPUT_MAF_FILE}. $!";
    my $maffh = IO::File->new( $this->{_MAF_FILE} ) or die "Couldn't open $this->{_MAF_FILE}. $!";
    my @outputs;
    # perl -an -F'\t' -e '$F[0]="chr$F[4]:".(int(($F[5]-1)/1000000)*1000000+1) unless($F[0]=~m/Hugo_Symbol/);
    # print join("\t",@F)' merged.cds.filtered.maf > merged_modded_gene_names_based_1M_regions.maf
    #
    while ( my $line = $maffh->getline ) {
        my @t = split( /\t/, $line );
        $t[0] = "CHR$t[4]:".(int(($t[5]-1)/$this->{_WINDOW_SIZE})*$this->{_WINDOW_SIZE}+1) unless($t[0]=~m/Hugo_Symbol/);
        push( @outputs, join("\t", @t) ); 
    }
    $maffh->close();
    $outputfh->print( @outputs );
    $outputfh->close();
    #
    return 1;
}

## usage
sub help_text {
    my $this = shift;
        return <<HELP
USAGE
 music2 bmr calc-window-maf --maf-file=? --output-maf-file=? [--window-size=?]

SYNOPSIS
 music2 bmr calc-window-maf \
    --maf-file input_dir/input.maf \
    --output-maf-file output_dir/window_based.maf \

 music2 bmr calc-window-maf \
    --maf-file input_dir/input.maf \
    --output-maf-file output_dir/window_based.maf \
    --window-size 1000000

REQUIRED INPUTS
  maf-file
    List of mutations using TCGA MAF specification v2.3 
  output-maf-file
    Output window based MAF file  

OPTIONAL INPUTS
  window-size
    The window size for local gene/regional specific BMR calculation 
    Default value '1000000' if not specified

REQUIRED OUTPUTS
  maf-file

DESCRIPTION
    Given a mutation list (MAF), this script calculates the window based MAF file based on given
    window size.


ARGUMENTS

    --window-size

      The window size for local gene/regional specific BMR calculation. The typical size is 1M and
      user can choose appropriate window size based on input gene/region ROIs.

HELP

}

1;

