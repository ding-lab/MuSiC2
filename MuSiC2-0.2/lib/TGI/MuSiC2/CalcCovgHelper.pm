package TGI::MuSiC2::CalcCovgHelper;
##
# calculate coverage one pair
##
#
use strict;
use warnings;

use IO::File;
use Getopt::Long;

sub new {
    my $class = shift;
    my $this = {};

    $this->{_ROI_FILE} = undef;
    $this->{_REF_SEQ}  = undef;
    $this->{_NOR_TUM_PAIR} = undef;
    $this->{_OUTPUT_FILE} = undef;
    $this->{_BP_CLASS_TYPES} = 'AT,CpG,CG';
    $this->{_NOR_MIN_DEPTH} = 6;
    $this->{_TUM_MIN_DEPTH} = 8;
    $this->{_MIN_MAPQ} = 20;

    bless $this, $class;
    $this->process();

    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'roi-file=s'              => \$this->{_ROI_FILE},
        'reference-sequence=s'    => \$this->{_REF_SEQ},
        'normal-tumor-bam-pair=s' => \$this->{_NOR_TUM_PAIR},
        'output-file=s'           => \$this->{_OUTPUT_FILE},
        'bp-class-types=s'        => \$this->{_BP_CLASS_TYPES},
        'normal-min-depth=i'      => \$this->{_NOR_MIN_DEPTH},
        'tumor-min-depth=i'       => \$this->{_TUM_MIN_DEPTH},
        'min-mapq=i'              => \$this->{_MIN_MAPQ},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    
    my ( $sample_name, $normal_bam, $tumor_bam, ) = split /\t/, $this->{_NOR_TUM_PAIR};

    # Check on all the required input data 
    print STDERR "ROI file not found or is empty: $this->{_ROI_FILE}\n"    unless( -s $this->{_ROI_FILE} );
    print STDERR "Reference sequence file not found: $this->{_REF_SEQ}\n"  unless( -e $this->{_REF_SEQ} );
    print STDERR "Normal BAM file not found or is empty: $normal_bam\n"    unless( -s $normal_bam );
    print STDERR "Tumor BAM file not found or is empty: $tumor_bam\n"      unless( -s $tumor_bam );
    return undef unless( -s $this->{_ROI_FILE} && -e $this->{_REF_SEQ} && -s $normal_bam && -s $tumor_bam );

    #### processing ####
    #
    # Check whether the annotated regions of interest 
    # are clumped together by chromosome
    my $roiFh = IO::File->new( $this->{_ROI_FILE} ) or die "ROI file could not be opened. $!\n";
    my @chroms = ( "" );
    # Emulate Unix's uniq command on the chromosome column
    while ( my $line = $roiFh->getline ) {
        my ( $chrom ) = ( $line =~ m/^(\S+)/ );
        push( @chroms, $chrom ) if ( $chrom ne $chroms[-1] );
    }
    $roiFh->close;
    # Get the actual number of unique chromosomes
    my %chroms = map { $_ => 1 } @chroms; 
    if ( scalar( @chroms ) != scalar( keys %chroms ) ) {
        print STDERR "ROIs from the same chromosome must be listed adjacent to each other in file. ";
        print STDERR "If in UNIX, try:\nsort -k 1,1 $this->{_ROI_FILE}\n";
        return undef;
    }
    # If the reference sequence FASTA file hasn't been indexed, do it
    my $ref_seq_idx = "$this->{_REF_SEQ}.fai";
    system( "samtools faidx $this->{_REF_SEQ}" ) unless( -e $ref_seq_idx );
    $normal_bam = '' unless( defined $normal_bam );
    $tumor_bam  = '' unless( defined $tumor_bam );
    print STDERR "Normal BAM not found: \"$normal_bam\"\n" unless( -e $normal_bam );
    print STDERR "Tumor BAM not found: \"$tumor_bam\"\n" unless( -e $tumor_bam );
    next unless( -e $normal_bam && -e $tumor_bam );
    # Construct the command that calculates coverage per ROI
    #my $calcRoiCovg_cmd = "calcRoiCovg $normal_bam $tumor_bam $roi_file $ref_seq $output_file $normal_min_depth $tumor_min_depth $min_mapq";
    my $cprogram_paras = "calcRoiCovg -q $this->{_MIN_MAPQ} -n $this->{_NOR_MIN_DEPTH} -t $this->{_TUM_MIN_DEPTH} -c $this->{_BP_CLASS_TYPES} ";
    my $required_paras = "$normal_bam $tumor_bam $this->{_ROI_FILE} $this->{_REF_SEQ} $this->{_OUTPUT_FILE}";
    my $calcRoiCovg_cmd = $cprogram_paras . $required_paras;
    # If the calcRoiCovg output was already 
    # generated, then don't rerun it
    if ( -s $this->{_OUTPUT_FILE} ) {
        print "Output file $this->{_OUTPUT_FILE} found. Skipping re-calculation.\n";
    }
    # Run the calcRoiCovg command on this tumor-normal pair. This could take a while
    elsif( system( "$calcRoiCovg_cmd" ) != 0 ) {
        print STDERR "Failed to execute: $calcRoiCovg_cmd\n";
        return;
    } else {
        print "$this->{_OUTPUT_FILE} generated and stored.\n";
        return 1;
    }

}

## usage
sub help_text {
    my $this = shift;
        return <<HELP

USAGE
 music2 bmr calc-covg-helper --roi-file=? --reference-sequence=?
    --normal-tumor-bam-pair=? --output-file=? [--normal-min-depth=?]
    [--tumor-min-depth=?] [--min-mapq=?] [--bp-class-types=?]

SYNOPSIS
General usage:

 music2 bmr calc-covg-helper \
    --normal-tumor-bam-pair "sample-name path/to/normal_bam path/to/tumor_bam" \
    --reference-sequence input_dir/all_sequences.fa \
    --output-file output_file \
    --roi-file input_dir/all_coding_exons.tsv

REQUIRED INPUTS
  roi-file
    Tab delimited list of ROIs [chr start stop gene_name] (See Description) 
  reference-sequence
    Path to reference sequence in FASTA format 
  normal-tumor-bam-pair
    Tab delimited line with sample name, path to normal bam file, and path to tumor bam file (See
    Description) 
  output-file
    Output file path.  Specify either output-file or output-directory. 

OPTIONAL INPUTS
  normal-min-depth
    The minimum read depth to consider a Normal BAM base as covered 
    Default value '6' if not specified
  tumor-min-depth
    The minimum read depth to consider a Tumor BAM base as covered 
    Default value '8' if not specified
  min-mapq
    The minimum mapping quality of reads to consider towards read depth counts 
    Default value '20' if not specified
  bp-class-types
    Bp class types for coverage calculating,delimited by comma
    Default value 'AT,CpG,CG' if not specified

DESCRIPTION
    This script counts bases with sufficient coverage in the ROIs of each gene in the given pair of
    tumor-normal BAM files and categorizes them into - AT, CG (non-CpG), and CpG counts. It also
    adds up these base-counts across all ROIs of each gene in the sample, but covered bases that
    lie within overlapping ROIs are not counted more than once towards these total counts.


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

    --reference-sequence

      The reference sequence in FASTA format. If a reference sequence index is not found next to
      this file (a .fai file), it will be created.

    --normal-tumor-bam-pair

      "sample-name path/to/normal_bam path/to/tumor_bam"

    --output-file

      Specify an output file where the per-ROI covered base counts will be written

HELP

}

1;

