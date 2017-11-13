package TGI::MuSiC2::CalcWindowRoi;
##
# given gene based ROIs, calculate window based ROI 
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

    $this->{_ROI_FILE} = undef;
    $this->{_OUTPUT_ROI_FILE}  = "window_based_roi";

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

        'roi-file=s'               => \$this->{_ROI_FILE},
        'output-roi-file=s'        => \$this->{_OUTPUT_ROI_FILE},
        'window-size=i'            => \$this->{_WINDOW_SIZE},

        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    #### processing ####
    #
    # Check on all the input data
    print STDERR "ROI file not found or is empty: $this->{_ROI_FILE}\n" unless( -s $this->{_ROI_FILE} );
    return undef unless( -s $this->{_ROI_FILE} );
    unless ( $this->{_WINDOW_SIZE} =~ /\d+/ ) {  
        print STDERR "Window size format is not valid !\n"; 
        return undef; 
    };
    #
    #
    my ( undef, $temp_win_bed ) = tempfile();
    my $outputfh = IO::File->new( $temp_win_bed, ">" ) or die "Temporary file could not be created. $!";
    my $roifh = IO::File->new( $this->{_ROI_FILE} ) or die "Couldn't open $this->{_ROI_FILE}. $!";
    while ( my $line = $roifh->getline ) {
        chomp( $line );
        my ( $chr, $start, $stop ) = split( /\t/, $line );
        ( $stop >= $start ) or die "Stop locus is less than start in:\n$line\n\n";
        my ( $start_window, $stop_window ) = (( int(($start-1)/$this->{_WINDOW_SIZE}) * $this->{_WINDOW_SIZE} + 1 ), ( int(($stop-1)/$this->{_WINDOW_SIZE}) * $this->{_WINDOW_SIZE} + 1 ));
        if ( $start_window == $stop_window ) {
            $outputfh->print( "$chr\t$start\t$stop\tCHR$chr:$start_window\n" );
        } elsif ( $start_window < $stop_window ) {
            $outputfh->print( "$chr\t$start\t",  $start_window + $this->{_WINDOW_SIZE} - 1, "\tCHR$chr:$start_window\n" );
            $start_window += $this->{_WINDOW_SIZE};
            while ( $start_window != $stop_window ) {
                $outputfh->print( "$chr\t$start_window\t", $start_window + $this->{_WINDOW_SIZE} - 1, "\tCHR$chr:$start_window\n" );
                $start_window += $this->{_WINDOW_SIZE};
            }
            $outputfh->print( "$chr\t$stop_window\t$stop\tCHR$chr:$stop_window\n" );
        } else { die "Unhandled exception! \n"; }
    }
    $roifh->close();
    $outputfh->close();
    #
    my ( undef, $temp_win_sorted_bed ) = tempfile();
    unless ( -e $temp_win_sorted_bed ) { die "Temporary file could not be created. $!" };
    system( "joinx sort -i $temp_win_bed -o $temp_win_sorted_bed" );
    #
    # write results
    #
    system( "joinx bed-merge -n -u -i $temp_win_sorted_bed -o $this->{_OUTPUT_ROI_FILE}" );
    #
    #
    return 1;

}

## usage
sub help_text {
    my $this = shift;
        return <<HELP
USAGE
 music2 bmr calc-window-roi --roi-file=? --output-roi-file=? [--window-size=?]

SYNOPSIS
 music2 bmr calc-window-roi \
    --roi-file input_dir/all_coding_exons.tsv \
    --output-roi-file output_dir/window_roi.tsv \

 music2 bmr calc-window-roi \
    --roi-file input_dir/all_coding_exons.tsv \
    --output-roi-file output_dir/window_roi.tsv \
    --window-size 1000000

REQUIRED INPUTS
  roi-file
    Tab delimited list of ROIs [chr start stop gene_name] (See DESCRIPTION) 
  output-roi-file
    Output window based ROIs file  

OPTIONAL INPUTS
  window-size
    The window size for local gene/regional specific BMR calculation 
    Default value '1000000' if not specified

REQUIRED OUTPUTS
  output-roi-file 

DESCRIPTION
    Given a tab delimited list of ROIs [chr start stop gene_name] and window size. 
    The script generates a file with window based ROIs.


ARGUMENTS

    --roi-file

      The regions of interest (ROIs) of each gene are typically regions targeted for sequencing or
      are merged exon loci (from multiple transcripts) of genes with 2-bp flanks (splice
      junctions). ROIs from the same chromosome must be listed adjacent to each other in this file.
      This allows the underlying C-based code to run much more efficiently and avoid re-counting
      bases seen in overlapping ROIs (for overall covered base counts). For per-gene base counts,
      an overlapping base will be counted each time it appears in an ROI of the same gene. To avoid
      this, be sure to merge together overlapping ROIs of the same gene. BEDtools' mergeBed can
      help if used per gene.

    --window-size

      The window size for local gene/regional specific BMR calculation. The typical size is 1M and
      user can choose appropriate window size based on input gene/region ROIs.

HELP

}

1;

