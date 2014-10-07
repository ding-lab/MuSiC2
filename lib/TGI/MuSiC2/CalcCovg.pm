package TGI::MuSiC2::CalcCovg;
##
# calculate coverage given multiple 
# normal tumor bam pairs 
##
use strict;
use warnings;

use IO::File;
use Getopt::Long;

sub new {
    my $class = shift;
    my $this = {};
    $this->{_ROI_FILE} = undef;
    $this->{_REF_SEQ}  = undef;
    $this->{_BAM_LIST} = undef;
    $this->{_OUTPUT_DIR} = undef;
    $this->{_CMD_LIST_FILE} = undef;
    $this->{_CMD_PREFIX} = undef;
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
    my ($help, $options);
    unless(@ARGV) { die $this->help_text(); }
    $options = GetOptions (
        'roi-file=s'              => \$this->{_ROI_FILE},
        'reference-sequence=s'    => \$this->{_REF_SEQ},
        'bam-list=s'              => \$this->{_BAM_LIST},
        'output-dir=s'            => \$this->{_OUTPUT_DIR},
        'cmd-list-file=s'         => \$this->{_CMD_LIST_FILE},
        'cmd-prefix=s'            => \$this->{_CMD_PREFIX},
        'bp-class-types=s'        => \$this->{_BP_CLASS_TYPES},
        'normal-min-depth=i'      => \$this->{_NOR_MIN_DEPTH},
        'tumor-min-depth=i'       => \$this->{_TUM_MIN_DEPTH},
        'min-mapq=i'              => \$this->{_MIN_MAPQ},
        'help' => \$help,
    );
    if ($help) { print STDERR help_text(); exit 0; }
    unless($options) { die $this->help_text(); }
    #### processing ####
    #
    # Check on all the input data
    print STDERR "ROI file not found or is empty: $this->{_ROI_FILE}\n"     unless( -s $this->{_ROI_FILE} );
    print STDERR "Reference sequence file not found: $this->{_REF_SEQ}\n"   unless( -e $this->{_REF_SEQ} );
    print STDERR "List of BAMs not found or is empty: $this->{_BAM_LIST}\n" unless( -s $this->{_BAM_LIST} );
    print STDERR "Output directory not found: $this->{_OUTPUT_DIR}\n"       unless( -e $this->{_OUTPUT_DIR} );
    return undef unless( -s $this->{_ROI_FILE} && -e $this->{_REF_SEQ} && -s $this->{_BAM_LIST} && -e $this->{_OUTPUT_DIR} );
    # Outputs of this script will be written to these 
    # locations in the output directory
    #
    # Remove trailing forward slashes if any
    $this->{_OUTPUT_DIR} =~ s/(\/)+$//; 
    # Stores output from calcRoiCovg per sample
    my $roi_covg_dir = "$this->{_OUTPUT_DIR}/roi_covgs";
    # Stores per-gene coverages per sample
    my $gene_covg_dir = "$this->{_OUTPUT_DIR}/gene_covgs"; 
    # Stores total coverages per sample
    my $tot_covg_file = "$this->{_OUTPUT_DIR}/total_covgs";
    ## optional paras 
    my $optional_params = "--normal-min-depth=$this->{_NOR_MIN_DEPTH} --tumor-min-depth=$this->{_TUM_MIN_DEPTH} --min-mapq=$this->{_MIN_MAPQ} --bp-class-types=$this->{_BP_CLASS_TYPES}";
    # Check whether the annotated regions of interest are clumped together by chromosome
    my $roiFh = IO::File->new( $this->{_ROI_FILE} ) or die "ROI file could not be opened. $!\n";
    my @chroms = ( "" );
    # Emulate Unix's uniq command on the chromosome column
    while ( my $line = $roiFh->getline ) {
        my ( $chrom ) = ( $line =~ m/^(\S+)/ );
        push( @chroms, $chrom ) if( $chrom ne $chroms[-1] );
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

    # Create the output directories unless they already exist
    mkdir $roi_covg_dir unless( -e $roi_covg_dir );
    mkdir $gene_covg_dir unless( -e $gene_covg_dir );

    my ( $cmdFh, $totCovgFh );
    if ( defined $this->{_CMD_LIST_FILE} ) {
        $cmdFh = IO::File->new( $this->{_CMD_LIST_FILE}, ">" );
        print "Creating a list of parallelizable jobs at $this->{_CMD_LIST_FILE}.\n";
        print "After successfully running all the jobs in $this->{_CMD_LIST_FILE},\n",
        "be sure to run this script a second time (without defining the cmd-list-file argument) to merge results in roi_covgs.\n";
    } else {
        $totCovgFh = IO::File->new( $tot_covg_file, ">" );
        ## process bp class types here
        # instead of using fixed AT, CpG, CG types
        #
        #
        $totCovgFh->print( "#Sample\tCovered_Bases\t" );
        $totCovgFh->print( join("_Bases_Covered\t", split /,/, $this->{_BP_CLASS_TYPES}) ); 
        $totCovgFh->print( "_Bases_Covered\n" );
    }
    # Parse through each pair of BAM files provided 
    # and run calcRoiCovg as necessary
    my $bamFh = IO::File->new( $this->{_BAM_LIST} );
    while ( my $line = $bamFh->getline ) {
        next if( $line =~ m/^#/ );
        chomp( $line );
        my ( $sample, $normal_bam, $tumor_bam ) = split( /\t/, $line );
        $normal_bam = '' unless( defined $normal_bam );
        $tumor_bam = '' unless( defined $tumor_bam );
        print STDERR "Normal BAM for $sample not found: \"$normal_bam\"\n" unless( -e $normal_bam );
        print STDERR "Tumor BAM for $sample not found: \"$tumor_bam\"\n"   unless( -e $tumor_bam );
        next unless( -e $normal_bam && -e $tumor_bam );
        # Construct the command that calculates coverage per ROI
        #
        my $calcRoiCovg_cmd = "\'music2 bmr calc-covg-helper --normal-tumor-bam-pair=\"$line\" --roi-file=$this->{_ROI_FILE} --reference-sequence=$this->{_REF_SEQ} --output-file=$roi_covg_dir\/$sample.covg $optional_params\'";

        # If user only wants the calcRoiCovg commands, write 
        # them to file and skip running calcRoiCovg
        if ( defined $this->{_CMD_LIST_FILE} ) {
            $calcRoiCovg_cmd = $this->{_CMD_PREFIX} . " $calcRoiCovg_cmd" if ( defined $this->{_CMD_PREFIX} );
            $cmdFh->print( "$calcRoiCovg_cmd\n" );
            next;
        }

        # If the calcRoiCovg output was already 
        # generated, then don't rerun it
        if ( -s "$roi_covg_dir/$sample.covg" ) {
            print "$sample.covg found in $roi_covg_dir. Skipping re-calculation.\n";
        } else {
            print STDERR "$sample.covg not found in $roi_covg_dir. please make a command list file to run calcRoiCovg !\n";
            return undef; 
        }
        # Read the calcRoiCovg output and count 
        # covered bases per gene
        my %geneCovg = ();
        my @bp_types_array = split /,/, $this->{_BP_CLASS_TYPES};

        #my ( $tot_covd, $tot_at_covd, $tot_cg_covg, $tot_cpg_covd );
        my ( $tot_covd, );
        my %tot_covd_bp_types  = map{ ($_, 0) } @bp_types_array;
        my %gene_covd_bp_types = map{ ($_, 0) } @bp_types_array;

        my $roiCovgFh = IO::File->new( "$roi_covg_dir/$sample.covg" );
        while ( my $line = $roiCovgFh->getline ) {
            chomp( $line );
            if ( $line =~ m/^#NonOverlappingTotals/ ) {
                my @ta = split( /\t/, $line );
                $tot_covd = $ta[3]; my $i = 1;
                map { $tot_covd_bp_types{$_} = $ta[3+$i++]; } @bp_types_array;

                #( undef, undef, undef, $tot_covd, $tot_at_covd, $tot_cg_covg, $tot_cpg_covd ) = split( /\t/, $line );
                #
            } elsif ( $line !~ m/^#/ ) {
                # my ( $gene, undef, $length, $covd, $at_covd, $cg_covd, $cpg_covd ) = split( /\t/, $line );
                my @tag = split( /\t/, $line );
                 my ( $gene, undef, $length, $covd, ) = @tag;
                $geneCovg{$gene}{len} += $length;
                $geneCovg{$gene}{covd_len} += $covd;
                my $i = 1;
                map { $geneCovg{$gene}{$_} += $tag[3+$i++]; } @bp_types_array;

                #$geneCovg{$gene}{at} += $at_covd;
                #$geneCovg{$gene}{cg} += $cg_covd;
                #$geneCovg{$gene}{cpg} += $cpg_covd;
            }
        }

        $roiCovgFh->close;
        # Write the per-gene coverages to a file 
        # named after this sample_name
        #
        my $geneCovgFh = IO::File->new( "$gene_covg_dir/$sample.covg", ">" );
        #$geneCovgFh->print( "#Gene\tLength\tCovered\tAT_covd\tCG_covd\tCpG_covd\n" );
        $geneCovgFh->print( "#Gene\tLength\tCovered\t", join("_covd\t", @bp_types_array ), "_covd\n" ); 
        foreach my $gene ( sort keys %geneCovg ) {
            #$geneCovgFh->print( join( "\t", $gene, $geneCovg{$gene}{len}, $geneCovg{$gene}{covd_len},$geneCovg{$gene}{at}, $geneCovg{$gene}{cg}, $geneCovg{$gene}{cpg} ), "\n" );
            $geneCovgFh->print( join( "\t", $gene, $geneCovg{$gene}{len}, $geneCovg{$gene}{covd_len}, map { $geneCovg{$gene}{$_} } @bp_types_array ), "\n" );
        }
        $geneCovgFh->close;
        # Write total coverages for this sample to a file
        #$totCovgFh->print( "$sample\t$tot_covd\t$tot_at_covd\t$tot_cg_covg\t$tot_cpg_covd\n" );

        $totCovgFh->print( join( "\t", $sample, $tot_covd, map{ $tot_covd_bp_types{$_} } @bp_types_array ), "\n" );

    }
    $bamFh->close;
    $cmdFh->close if ( defined $this->{_CMD_LIST_FILE} );
    $totCovgFh->close unless ( defined $this->{_CMD_LIST_FILE} );

    return 1;
}

## usage
sub help_text {
    my $this = shift;
        return <<HELP
USAGE
 music2 bmr calc-covg --gene-covg-dir=? --roi-file=? --reference-sequence=? --bam-list=?
    --output-dir=? [--cmd-list-file=?] [--cmd-prefix=?] [--normal-min-depth=?]
    [--tumor-min-depth=?] [--min-mapq=?]

SYNOPSIS
General usage:

 music2 bmr calc-covg \
    --bam-list input_dir/bam_list \
    --output-dir output_dir/ \
    --reference-sequence input_dir/all_sequences.fa \
    --roi-file input_dir/all_coding_exons.tsv

To create a list of commands that will allow the processing of each tumor-normal pair in parallel
with an LSF job scheduler:

 music2 bmr calc-covg \
    --bam-list input_dir/bam_list \
    --output-dir output_dir/ \
    --reference-sequence input_dir/all_sequences.fa \
    --roi-file input_dir/all_coding_exons.tsv \
    --cmd_list_file parallelizable_commands \
    --cmd_prefix bsub

In the above case, the commands printed into the output file "parallelizable_commands" can be run
in parallel. After they complete, rerun this script as printed directly below (--cmd_list_file
and --cmd_prefix have been removed) to merge the parallelized calculations:

 music2 bmr calc-covg \
    --bam-list input_dir/bam_list \
    --output-dir output_dir/ \
    --reference-sequence input_dir/all_sequences.fa \
    --roi-file input_dir/all_coding_exons.tsv

REQUIRED INPUTS
  roi-file
    Tab delimited list of ROIs [chr start stop gene_name] (See Description) 
  reference-sequence
    Path to reference sequence in FASTA format 
  bam-list
    Tab delimited list of BAM files [sample_name normal_bam tumor_bam] (See Description) 
  output-dir
    Directory where output files and subdirectories will be written 

OPTIONAL INPUTS
  cmd-list-file
    A file to write calcRoiCovg commands to (See Description) 
  cmd-prefix
    A command that submits a job to your cluster (See Description) 
  normal-min-depth
    The minimum read depth to consider a Normal BAM base as covered 
  tumor-min-depth
    The minimum read depth to consider a Tumor BAM base as covered 
  min-mapq
    The minimum mapping quality of reads to consider towards read depth counts 
  bp-class-types
    Bp class types for coverage calculating,delimited by comma
    Default value 'AT,CpG,CG' if not specified

DESCRIPTION
    This script counts bases with sufficient coverage in the ROIs of each gene in the given pairs
    of tumor-normal BAM files and categorizes them into - AT, CG (non-CpG), and CpG counts. It also
    adds up these base-counts across all ROIs of each gene for each sample, but covered bases that
    lie within overlapping ROIs are not counted more than once towards these total counts.

    By default, this script runs a C-based tool named calcRoiCovg for each sample one after
    another, taking ~30 mins per sample to generate per-ROI covered base counts. If the results of
    calcRoiCovg for a sample already exists in the output subdirectory roi_covgs, re-calculation is
    skipped. This allows you to run your own calcRoiCovg jobs in parallel or on multiple machines
    (Keep reading).

    Speed things up by running calcRoiCovg jobs in parallel: If a compute cluster or multiple
    machines are available, run this script twice as follows:

      * Define cmd-list-file and cmd-prefix to generate a file with commands that can be submitted
      to a cluster or run manually. These jobs will write per-ROI base counts in a subdirectory
      roi_covgs.

      * After all the parallelized calcRoiCovg jobs are completed, run this script again to add
      them up and generate the final per-gene base counts in a subdirectory gene_covgs. Remember to
      remove the cmd-list-file and cmd-prefix arguments or you will just be re-creating a list of
      commands.


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

    --bam-list

      Provide a file containing sample names and normal/tumor BAM locations for each. Use the tab-
      delimited format [sample_name normal_bam tumor_bam] per line. Additional columns like
      clinical data are allowed, but ignored. The sample_name must be the same as the tumor sample
      names used in the MAF file (16th column, with the header Tumor_Sample_Barcode).

    --output-dir

      Specify an output directory where the following will be created/written: roi_covgs:
      Subdirectory containing per-ROI covered base counts for each sample. gene_covgs: Subdirectory
      containing per-gene covered base counts for each sample. total_covgs: File containing the
      overall non-overlapping coverages per sample.

    --cmd-list-file

      Specify a file into which a list of calcRoiCovg jobs will be written to. These can be
      scheduled in parallel, and will write per-ROI covered base-counts into the output
      subdirectory roi_covgs. If cmd-list-file is left unspecified, this script runs calcRoiCovg
      per sample one after another, taking ~30 mins per sample, but it skips samples whose output
      is already in roi_covgs.

    --cmd-prefix

      Specify a job submission command that will be prefixed to each command in cmd-list-file. This
      makes batch submission easier. Just run the cmd-list-file file as a shell script to submit
      jobs. cmd-prefix is "bsub" if your cluster uses the LSF job scheduler, or "qsub" in Torque.
      Add arguments as necessary. For example, "bsub -M 4GB" sets a soft memory limit of 4GB.



HELP

}

1;

