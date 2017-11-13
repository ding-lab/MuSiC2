package TGI::MuSiC2::CalcWigCovg;
##
# calculate coverage given multiple 
# wig list file
#
use strict;
use warnings;

use IO::File;
use Getopt::Long;
use File::Temp qw/ tempfile /;

sub new {
    my $class = shift;
    my $this = {};

    $this->{_ROI_FILE} = undef;
    $this->{_REF_SEQ}  = undef;
    $this->{_WIG_LIST} = undef;
    $this->{_OUTPUT_DIR} = undef;
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
        'wig-list=s'              => \$this->{_WIG_LIST},
        'output-dir=s'            => \$this->{_OUTPUT_DIR},
        'bp-class-types=s'        => \$this->{_BP_CLASS_TYPES},
        'normal-min-depth=i'      => \$this->{_NOR_MIN_DEPTH},
        'tumor-min-depth=i'       => \$this->{_TUM_MIN_DEPTH},
        'min-mapq=i'              => \$this->{_MIN_MAPQ},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }

    #### processing ####
    #
    # Check on all the input data
    print STDERR "ROI file not found or is empty: $this->{_ROI_FILE}\n"     unless( -s $this->{_ROI_FILE} );
    print STDERR "Reference sequence file not found: $this->{_REF_SEQ}\n"   unless( -e $this->{_REF_SEQ}  );
    print STDERR "List of WIGs not found or is empty: $this->{_WIG_LIST}\n" unless( -s $this->{_WIG_LIST} );
    print STDERR "Output directory not found: $this->{_OUTPUT_DIR}\n"       unless( -e $this->{_OUTPUT_DIR} );
    return undef unless( -s $this->{_ROI_FILE} && -e $this->{_REF_SEQ} && -s $this->{_WIG_LIST} && -e $this->{_OUTPUT_DIR} );

    # Remove trailing forward slashes if any
    $this->{_OUTPUT_DIR} =~ s/(\/)+$//; 
    # Stores output from calcRoiCovg per sample
    my $roi_covg_dir = "$this->{_OUTPUT_DIR}/roi_covgs";
    # Stores per-gene coverages per sample
    my $gene_covg_dir = "$this->{_OUTPUT_DIR}/gene_covgs"; 
    # Stores total coverages per sample
    my $tot_covg_file = "$this->{_OUTPUT_DIR}/total_covgs";

    # If the reference sequence FASTA file hasn't been indexed, do it
    my $ref_seq_idx = "$this->{_REF_SEQ}.fai";
    system( "samtools faidx $this->{_REF_SEQ}" ) unless( -e $ref_seq_idx );
    #
    # Create a temporary 0-based ROI BED-file that we can use with 
    # joinx, and also measure gene lengths
    #
    my %geneLen = ();
    #
    ## using file::temp module for temp file

    # didn't touch Cyriac's code in this part
    #
    my ( undef, $roi_bed ) = tempfile();
    my $roiBedFh = IO::File->new( $roi_bed, ">" ) or die "Temporary ROI BED file could not be created. $!\n";
    my $roiFh = IO::File->new( $this->{_ROI_FILE} ) or die "ROI file could not be opened. $!\n";
    while ( my $line = $roiFh->getline ) {
        chomp( $line );
        my ( $chr, $start, $stop, $gene ) = split( /\t/, $line );
        --$start;
        unless( $start >= 0 && $start < $stop ) {
            print STDERR "Invalid ROI: $line\nPlease use 1-based loci and ensure that start <= stop\n";
            return undef;
        }
        $geneLen{$gene} += ( $stop - $start );
        $roiBedFh->print( "$chr\t$start\t$stop\t$gene\n" );
    }
    $roiFh->close;
    $roiBedFh->close;
    #
    #
    # Also create a merged BED file where overlapping ROIs are joined 
    # together into contiguous regions
    # ::TODO:: Use joinx instead of mergeBed, because 
    # we'd rather add an in-house dependency
    #
    #my $merged_roi_bed = "$this->{_OUTPUT_DIR}/.merged_bed_file";
    my ( undef, $merged_roi_bed ) = tempfile();
    #
    system( "mergeBed -i $roi_bed | joinx1.7 sort -s - -o $merged_roi_bed" );
    # or die "Failed to run mergeBed or joinx!\n$roi_bed\n$merged_roi_bed\n $!\n";
    #
    # Create the output directories unless they already exist
    mkdir $roi_covg_dir  unless( -e $roi_covg_dir );
    mkdir $gene_covg_dir unless( -e $gene_covg_dir );
    # bp class types 
    my @bp_types_array = split /,/, $this->{_BP_CLASS_TYPES};
    #
    # This is a file that will report the overall 
    # non-overlapping coverages per WIG
    #
    my $totCovgFh = IO::File->new( $tot_covg_file, ">" );
    $totCovgFh->print( "#Sample\tCovered_Bases\tAT_Bases_Covered\tCG_Bases_Covered\tCpG_Bases_Covered\n" );
    # Parse through each pair of WIG files provided and run calcRoiCovg as necessary
    my $wigFh = IO::File->new( $this->{_WIG_LIST} );
    while( my $line = $wigFh->getline ) {
        next if ( $line =~ m/^#/ );
        chomp( $line );
        my ( $sample, $wig_file ) = split( /\t/, $line );
        $wig_file = '' unless( defined $wig_file );
        print STDERR "Wiggle track format file for $sample not found: \"$wig_file\"\n" unless( -e $wig_file );
        next unless( -e $wig_file );
        #
        # Use joinx to parse the WIG file and return per-ROI 
        # coverages of AT, CG (non-CpG), and CpG
        #
        system( "joinx1.7 wig2bed -Zc $wig_file | joinx1.7 sort -s | joinx1.7 intersect -F \"I A3\" $roi_bed - | joinx1.7 ref-stats - $this->{_REF_SEQ} | cut -f 1-7 > $roi_covg_dir/$sample.covg" );
        # or die "Failed to run joinx to calculate per-gene coverages in $sample! $!\n";
        #
        # Read the joinx formatted coverage file and count covered bases per gene
        #
        my %geneCovg = ();
        my $roiCovgFh = IO::File->new( "$roi_covg_dir/$sample.covg" );
        #
        # TODO:: need to develop MuSiCmate to do any XpX coverage information 
        # didn't change "at,cg,cpg" coverage here
        #
        while ( my $line = $roiCovgFh->getline ) {
            chomp( $line );
            if ( $line !~ m/^#/ ) {
                my ( undef, undef, undef, $gene, $at_covd, $cg_covd, $cpg_covd ) = split( /\t/, $line );
                $geneCovg{$gene}{covd} += ( $at_covd + $cg_covd + $cpg_covd );
                $geneCovg{$gene}{at} += $at_covd;
                $geneCovg{$gene}{cg} += $cg_covd;
                $geneCovg{$gene}{cpg} += $cpg_covd;
            }

        }
        $roiCovgFh->close;
        #
        # Write the per-gene coverages to a file named after this sample_name
        my $geneCovgFh = IO::File->new( "$gene_covg_dir/$sample.covg", ">" );
        $geneCovgFh->print( "#Gene\tLength\tCovered\tAT_covd\tCG_covd\tCpG_covd\n" );
        foreach my $gene ( sort keys %geneLen ) {
            if ( defined $geneCovg{$gene} ) {
                $geneCovgFh->print( join( "\t", $gene, $geneLen{$gene}, $geneCovg{$gene}{covd}, $geneCovg{$gene}{at}, $geneCovg{$gene}{cg}, $geneCovg{$gene}{cpg} ), "\n" );
            } else { $geneCovgFh->print( "$gene\t" . $geneLen{$gene} . "\t0\t0\t0\t0\n" ); }
        }
        $geneCovgFh->close;
        # Measure coverage stats on the merged ROI file, so that 
        # bps across the genome are not counted twice
        my ( undef, $merged_roi_bed_covg ) = tempfile();
        system( "joinx1.7 wig2bed -Zc $wig_file | joinx1.7 sort -s | joinx1.7 intersect $merged_roi_bed - | joinx1.7 ref-stats - $this->{_REF_SEQ} | cut -f 1-6 > $merged_roi_bed_covg" );
        # or die "Failed to run joinx to calculate overall coverages in $sample! $!\n";
        #
        # Read the joinx formatted coverage file and sum up the coverage stats per region
        my ( $tot_covd, $tot_at_covd, $tot_cg_covg, $tot_cpg_covd );
        my $totRoiCovgFh = IO::File->new( $merged_roi_bed_covg );
        while( my $line = $totRoiCovgFh->getline ) {
            chomp( $line );
            if ( $line !~ m/^#/ ) {
                my ( $chr, $start, $stop, $at_covd, $cg_covd, $cpg_covd ) = split( /\t/, $line );
                $tot_covd += ( $at_covd + $cg_covd + $cpg_covd );
                $tot_at_covd += $at_covd;
                $tot_cg_covg += $cg_covd;
                $tot_cpg_covd += $cpg_covd;
            }
        }
        $totRoiCovgFh->close;
        $totCovgFh->print( "$sample\t$tot_covd\t$tot_at_covd\t$tot_cg_covg\t$tot_cpg_covd\n" );
    }
    $wigFh->close;
    $totCovgFh->close;

    return 1;
}

## usage
sub help_text {
    my $this = shift;
        return <<HELP
USAGE
 music2 bmr calc-wig-covg --roi-file=? --reference-sequence=? --wig-list=?
    --output-dir=?

SYNOPSIS
General usage:

 music2 bmr calc-wig-covg \
    --wig-list input_dir/wig_list \
    --output-dir output_dir/ \
    --reference-sequence input_dir/all_sequences.fa \
    --roi-file input_dir/all_coding_exons.tsv

REQUIRED INPUTS
  roi-file
    Tab-delimited list of ROIs [chr start stop gene_name] (See Description) 
  reference-sequence
    Path to reference sequence in FASTA format 
  wig-list
    Tab-delimited list of WIG files [sample_name wig_file] (See Description) 
  output-dir
    Directory where output files and subdirectories will be written 

REQUIRED OUTPUTS
  output-dir
    Directory where output files and subdirectories will be written 

DESCRIPTION
    This script counts bases with sufficient coverage in the ROIs of each gene from given wiggle
    track format files, and categorizes them into - AT, CG (non-CpG), and CpG counts. It also adds
    up these base-counts across all ROIs of each gene for each sample, but covered bases that lie
    within overlapping ROIs are not counted more than once towards these total counts.


ARGUMENTS

    --roi-file

      The regions of interest (ROIs) of each gene are typically regions targeted for sequencing or
      are merged exon loci (from multiple transcripts) of genes with 2-bp flanks (splice
      junctions). For per-gene base counts, an overlapping base will be counted each time it
      appears in an ROI of the same gene. To avoid this, be sure to merge together overlapping ROIs
      of the same gene. BEDtools' mergeBed can help if used per gene.

    --reference-sequence

      The reference sequence in FASTA format. If a reference sequence index is not found next to
      this file (a .fai file), it will be created.

    --wig-list

      Provide a file containing sample names and the wiggle track format file locations for each.
      Use the tab-delimited format [sample_name wig_file] per line. Additional columns like
      clinical data are allowed, but ignored. The sample_name must be the same as the tumor sample
      names used in the MAF file (16th column, with the header Tumor_Sample_Barcode).

    --output-dir

      Specify an output directory where the following will be created/written: roi_covgs:
      Subdirectory containing per-ROI covered base counts for each sample. gene_covgs: Subdirectory
      containing per-gene covered base counts for each sample. total_covgs: File containing the
      overall non-overlapping coverages per sample.


HELP

}

1;

