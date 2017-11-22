package TGI::MuSiC2::CalcBmr;
##
# calculate background mutation rate (BMR) 
##
use strict;
use warnings;
#
use Bit::Vector;
use List::Util qw( min sum );
use File::Temp qw/ tempfile /;
#
# These constants let us use space-efficient arrays instead 
# of hashes, while keeping the code fairly readable
#
use constant { AT_Transitions    => 0, 
               AT_Transversions  => 1, 
               CG_Transitions    => 2, 
               CG_Transversions  => 3,
               CpG_Transitions   => 4, 
               CpG_Transversions => 5, 
               Indels            => 6, 
               Truncations       => 7, 
               Overall           => 8,
               covd_bases        => 0, 
               mutations         => 1, 
               bmr => 2 };
#
#
use IO::File;
use Getopt::Long;

sub new {
    my $class = shift;
    my $this = {};

    $this->{_ROI_FILE} = undef;
    $this->{_REF_SEQ}  = undef;
    $this->{_BAM_LIST} = undef;
    $this->{_OUTPUT_DIR} = undef;
    $this->{_MAF_FILE} = undef;
    $this->{_BMR_GROUPS} = 1;
    $this->{_SHOW_SKIPPED} = 0;
    $this->{_NOSHOW_SKIPPED} = 1;
    $this->{_SEPERATE_TRUNCATIONS} = 0;
    $this->{_NOSEPERATE_TRUNCATIONS} = 1;
    $this->{_MERGE_CONCURRENT_MUTS} = 0;
    $this->{_NOMERGE_CONCURRENT_MUTS} = 1;
    $this->{_GENES_TO_IGNORE} = undef;
    $this->{_SKIP_NON_CODING} = 1;
    $this->{_NOSKIP_NON_CODING} = 0;
    $this->{_SKIP_SILENT} = 1;
    $this->{_NOSKIP_SILENT} = 0;
    $this->{_BMR_OUTPUT} = 0;
    $this->{_BP_CLASS_TYPES} = 'AT,CpG,GC';

    bless $this, $class;
    $this->process();

    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'roi-file=s'                => \$this->{_ROI_FILE},
        'reference-sequence=s'      => \$this->{_REF_SEQ},
        'bam-list=s'                => \$this->{_BAM_LIST},
        'output-dir=s'              => \$this->{_OUTPUT_DIR},
        'maf-file=s'                => \$this->{_MAF_FILE},
        'bmr-groups=i'              => \$this->{_BMR_GROUPS},
        'show-skipped'              => \$this->{_SHOW_SKIPPED},
        'noshow-skipped'            => \$this->{_NOSHOW_SKIPPED},
        'seperate-truncations'      => \$this->{_SEPERATE_TRUNCATIONS},
        'noseperate-truncations'    => \$this->{_NOSEPERATE_TRUNCATIONS},
        'merge-concurrent-muts'     => \$this->{_MERGE_CONCURRENT_MUTS},
        'nomerge-concurrent-muts'   => \$this->{_NOMERGE_CONCURRENT_MUTS},
        'genes-to-ignore=s'         => \$this->{_GENES_TO_IGNORE},
        'skip-non-coding'           => \$this->{_SKIP_NON_CODING},
        'noskip-non-coding'         => \$this->{_NOSKIP_NON_CODING},
        'skip-silent'               => \$this->{_SKIP_SILENT},
        'noskip-silent'             => \$this->{_NOSKIP_SILENT},
        'bp-class-types=s'          => \$this->{_BP_CLASS_TYPES},

        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    #### processing ####
    #
    # Check on all the input data
    print STDERR "ROI file not found or is empty: $this->{_ROI_FILE}\n"     unless( -s $this->{_ROI_FILE} );
    print STDERR "Reference sequence file not found: $this->{_REF_SEQ}\n"   unless( -e $this->{_REF_SEQ} );
    print STDERR "List of BAMs not found or is empty: $this->{_BAM_LIST}\n" unless( -s $this->{_BAM_LIST} );
    print STDERR "Output directory not found: $this->{_OUTPUT_DIR}\n"       unless( -e $this->{_OUTPUT_DIR} );
    print STDERR "MAF file not found or is empty: $this->{_MAF_FILE}\n"     unless( -s $this->{_MAF_FILE} );
    return undef unless( -s $this->{_ROI_FILE} && -e $this->{_REF_SEQ} && -s $this->{_BAM_LIST} && -e $this->{_OUTPUT_DIR} && -s $this->{_MAF_FILE} );
    #
    #
    ## Mutation clusters potentially from the same event
    #  should be treated together
    #
    my $ref_rep_muts = $this->mut_clustering( $this->{_MAF_FILE} );
    my ( undef, $temp_clustering_maf_file ) = tempfile();
    my $temp_clustering_maf_fileFh = IO::File->new( $temp_clustering_maf_file, ">" ) or die "Temporary file could not be created. $!";
    foreach my $rep_mut ( @$ref_rep_muts ) { $temp_clustering_maf_fileFh->print( $rep_mut ) };
    $temp_clustering_maf_fileFh->close;
    $this->{_MAF_FILE} = $temp_clustering_maf_file;
    #
    #
    # process boolean paras
    #
    $this->{_SKIP_NON_CODING} = 0 if ( $this->{_NOSKIP_NON_CODING} );
    $this->{_SKIP_SILENT} = 0 if ( $this->{_NOSKIP_SILENT} );
    #
    # Outputs of this script will be written to these 
    # locations in the output directory
    #
    # Remove trailing forward slashes if any
    $this->{_OUTPUT_DIR} =~ s/(\/)+$//; 
    # Stores per-gene coverages per sample
    my $gene_covg_dir = "$this->{_OUTPUT_DIR}/gene_covgs"; 
    # Stores total coverages per sample
    my $total_covgs_file = "$this->{_OUTPUT_DIR}/total_covgs";
    print STDERR "Directory with per-gene coverages not found: $gene_covg_dir\n"   unless( -e $gene_covg_dir );
    print STDERR "Total coverages file not found or is empty: $total_covgs_file\n" unless( -s $total_covgs_file );
    return undef unless( -e $gene_covg_dir && -s $total_covgs_file );
    #
    # Outputs of this script will be written to these 
    # locations in the output directory
    my $overall_bmr_file = "$this->{_OUTPUT_DIR}/overall_bmrs";
    my $gene_mr_file = "$this->{_OUTPUT_DIR}/gene_mrs";
    #
    # Build a hash to quickly lookup the genes to be 
    # ignored for overall BMRs
    #
    my %ignored_genes = ();
    if ( defined $this->{_GENES_TO_IGNORE} ) {
        %ignored_genes = map { $_ => 1 } split( /,/, $this->{_GENES_TO_IGNORE} );
    }
    #
    # Parse out the names of the samples which should match 
    # the names of the coverage files needed
    #
    my ( @all_sample_names, %sample_idx );
    my $idx = 0;
    my $sampleFh = IO::File->new( $this->{_BAM_LIST} ) or die "Couldn't open $this->{_BAM_LIST}. $!";
    while( my $line = $sampleFh->getline ) {
        next if ( $line =~ m/^#/ );
        chomp( $line );
        my ( $sample ) = split( /\t/, $line );
        push( @all_sample_names, $sample );
        $sample_idx{$sample} = $idx++;
    }
    $sampleFh->close;
    #
    # If the reference sequence FASTA file hasn't 
    # been indexed, do it
    #
    my $ref_seq_idx = "$this->{_REF_SEQ}.fai";
    system( "samtools faidx $this->{_REF_SEQ}" ) unless( -e $ref_seq_idx );
    #
    # Parse gene names and ROIs. Mutations outside 
    # these ROIs will be skipped
    #
    my ( @all_gene_names, %gene_idx );
    $idx = 0;
    my $roi_bitmask = $this->create_empty_genome_bitmask( $ref_seq_idx );
    #
    my $roiFh = IO::File->new( $this->{_ROI_FILE} ) or die "Couldn't open $this->{_ROI_FILE}. $!";
    while ( my $line = $roiFh->getline ) {
        next if ( $line =~ m/^#/ );
        chomp $line;
        my ( $chr, $start, $stop, $gene ) = split( /\t/, $line );
        if ( !$roi_bitmask->{$chr} or $start > $roi_bitmask->{$chr}->Size ) {
            print STDERR "Skipping invalid ROI bitmask $chr:$start-$stop\n";
            next;
        }
        $roi_bitmask->{$chr}->Interval_Fill( $start, $stop );
        unless( defined $gene_idx{$gene} ) {
            push( @all_gene_names, $gene );
            $gene_idx{$gene} = $idx++;
        }
    }
    $roiFh->close;
    #
    # These are the various categories that each 
    # mutation will be classified into
    #
    my @mut_classes = ( AT_Transitions, AT_Transversions, CG_Transitions, CG_Transversions, CpG_Transitions, CpG_Transversions, Indels );
    push( @mut_classes, Truncations ) if ( $this->{_SEPERATE_TRUNCATIONS} );
    #
    # Save the actual class names for reporting purposes, because the 
    # elements above are really just numerical constants
    #
    my @mut_class_names = qw( AT_Transitions AT_Transversions CG_Transitions CG_Transversions CpG_Transitions CpG_Transversions Indels );
    push( @mut_class_names, 'Truncations' ) if ( $this->{_SEPERATE_TRUNCATIONS} );
    #
    # Stores per sample covg 
    # and mutation information
    #
    my @sample_mr;
    foreach my $sample ( @all_sample_names ) {
        $sample_mr[$sample_idx{$sample}][$_][mutations] = 0 foreach( @mut_classes );
        $sample_mr[$sample_idx{$sample}][$_][covd_bases] = 0 foreach( @mut_classes );
    }
    #
    # Load the covered base-counts per sample from the 
    # output of "music bmr calc-covg"
    #
    print STDERR "Loading per-sample coverages stored in $total_covgs_file\n";
    my $sample_cnt_in_total_covgs_file = 0;
    my $totCovgFh = IO::File->new( $total_covgs_file ) or die "Couldn't open $total_covgs_file. $!";
    while( my $line = $totCovgFh->getline ) {
        next unless( $line =~ m/^\S+\t\d+\t\d+\t\d+\t\d+$/ and $line !~ m/^#/ );
        chomp( $line );
        ++$sample_cnt_in_total_covgs_file;
        my ( $sample, $covd_bases, $covd_at_bases, $covd_cg_bases, $covd_cpg_bases ) = split( /\t/, $line );
        $sample_mr[$sample_idx{$sample}][AT_Transitions][covd_bases] = $covd_at_bases;
        $sample_mr[$sample_idx{$sample}][AT_Transversions][covd_bases] = $covd_at_bases;
        $sample_mr[$sample_idx{$sample}][CG_Transitions][covd_bases] = $covd_cg_bases;
        $sample_mr[$sample_idx{$sample}][CG_Transversions][covd_bases] = $covd_cg_bases;
        $sample_mr[$sample_idx{$sample}][CpG_Transitions][covd_bases] = $covd_cpg_bases;
        $sample_mr[$sample_idx{$sample}][CpG_Transversions][covd_bases] = $covd_cpg_bases;
        $sample_mr[$sample_idx{$sample}][Indels][covd_bases] = $covd_bases;
        $sample_mr[$sample_idx{$sample}][Truncations][covd_bases] = $covd_bases if ( $this->{_SEPERATE_TRUNCATIONS} );
    }
    $totCovgFh->close;
    #
    #
    unless( $sample_cnt_in_total_covgs_file == scalar( @all_sample_names )) {
        print STDERR "Mismatching number of samples in $total_covgs_file and $this->{_BAM_LIST}\n";
        return undef;
    }
    #
    # Stores per gene covg and mutation information
    #
    my %gene_mr;
#    foreach my $gene ( @all_gene_names ) {
#        foreach my $sample ( @all_sample_names ) {
#            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$_][mutations] = 0 foreach( @mut_classes );
#            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$_][covd_bases] = 0 foreach( @mut_classes );
#        }
#    }
    #
    # Sum up the per-gene covered base-counts across samples from 
    # the output of "music bmr calc-covg"
    #
    print STDERR "Loading per-gene coverage files stored under $gene_covg_dir/\n";
    foreach my $sample ( @all_sample_names ) {
        my $sample_covg_file = "$gene_covg_dir/$sample.covg";
        my $sampleCovgFh = IO::File->new( $sample_covg_file ) or die "Couldn't open $sample_covg_file. $!";
        while( my $line = $sampleCovgFh->getline ) {
            next unless( $line =~ m/^\S+\t\d+\t\d+\t\d+\t\d+\t\d+$/ and $line !~ m/^#/ );
            chomp( $line );
            my ( $gene, undef, $covd_bases, $covd_at_bases, $covd_cg_bases, $covd_cpg_bases ) = split( /\t/, $line );
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[AT_Transitions][covd_bases] += $covd_at_bases;
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[AT_Transversions][covd_bases] += $covd_at_bases;
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[CG_Transitions][covd_bases] += $covd_cg_bases;
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[CG_Transversions][covd_bases] += $covd_cg_bases;
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[CpG_Transitions][covd_bases] += $covd_cpg_bases;
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[CpG_Transversions][covd_bases] += $covd_cpg_bases;
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[Indels][covd_bases] += $covd_bases;
            $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[Truncations][covd_bases] += $covd_bases if ( $this->{_SEPERATE_TRUNCATIONS} );
        }
        $sampleCovgFh->close;
    }
    #
    # Run "joinx ref-stats" to classify SNVs as being at 
    # AT, CG, or CpG sites in the reference
    #
    print STDERR "Running 'joinx ref-stats' to read reference FASTA and identify SNVs at AT, CG, CpG sites\n";
    #
    #my $maf_bed = "$this->{_OUTPUT_DIR}/.temp_maf_bed_file";
    my ( undef, $maf_bed ) = tempfile();
    my $mafBedFh = IO::File->new( $maf_bed, ">" ) or die "Temporary file could not be created. $!";
    my $mafFh = IO::File->new( $this->{_MAF_FILE} ) or die "Couldn't open $this->{_MAF_FILE}. $!";
    while( my $line = $mafFh->getline ) {
        next if ( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $chr, $start, $stop, $mutation_type, $ref, $var1, $var2 ) = @cols[4..6,9..12];
        if ( $mutation_type =~ m/^(SNP|DNP|ONP|TNP)$/ ) {
            $mafBedFh->print( "$chr\t" . ( $start - 2 ) . "\t" . ( $start + 1 ) . "\n" );
        }
    }
    $mafFh->close;
    $mafBedFh->close;
    #
    #my $refstats_file = "$this->{_OUTPUT_DIR}/.temp_refstats_file";
    #
    my ( undef, $refstats_file )  = tempfile();
    `joinx ref-stats --ref-bases --bed $maf_bed --fasta $this->{_REF_SEQ} --output $refstats_file`;
    #
    #
    # Parse through the ref-stats output and load it 
    # into hashes for quick lookup later
    #
    my ( %ref_base, %cpg_site );
    my $refStatsFh = IO::File->new( $refstats_file ) or die "Couldn't open $refstats_file. $!";
    while( my $line = $refStatsFh->getline ) {
        next if ( $line =~ m/^#/ );
        chomp $line;
        my ( $chr, undef, $pos, undef, undef, undef, $ref ) = split( /\t/, $line );
        my $locus = "$chr\t" . ( $pos - 1 );
        $ref_base{$locus} = substr( $ref, 1, 1 );
        $cpg_site{$locus} = 1 if ( $ref =~ m/CG/i );
    }
    $refStatsFh->close;
    #
    # Create a hash to help classify SNVs
    #
    my %classify;
    $classify{$_} = AT_Transitions foreach( qw( AG TC ));
    $classify{$_} = AT_Transversions foreach( qw( AC AT TA TG ));
    $classify{$_} = CG_Transitions foreach( qw( CT GA ));
    $classify{$_} = CG_Transversions foreach( qw( CA CG GC GT ));
    #
    # Parse through the MAF file and categorize 
    # each somatic mutation
    #
    print STDERR "Parsing MAF file to classify mutations\n";
    my %skip_cnts;
    $mafFh = IO::File->new( $this->{_MAF_FILE} ) or die "Couldn't open $this->{_MAF_FILE}. $!";
    while ( my $line = $mafFh->getline ) {
        next if ( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $chr, $start, $stop, $mutation_class, $mutation_type, $ref, $var1, $var2, $sample ) = @cols[0,4..6,8..12,15];
        # Skip mutations in samples that are not in the provided bam list
        unless( defined $sample_idx{$sample} ) {
            $skip_cnts{"belong to unrecognized samples"}++;
            print STDERR "Skipping unrecognized sample ($sample not in BAM list): $gene, $chr:$start-$stop\n" if ( $this->{_SHOW_SKIPPED} );
            next;
        }
        # If the mutation classification is odd, quit with error
        if ( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            print STDERR "Unrecognized Variant_Classification \"$mutation_class\" in MAF file: $gene, $chr:$start-$stop\n";
            print STDERR "Please use TCGA MAF Specification v2.3.\n";
            return undef;
        }
        #
        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if (( $this->{_SKIP_NON_CODING} && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
          ( $this->{_SKIP_SILENT} && $mutation_class =~ m/^Silent$/ )) {
            $skip_cnts{"are classified as $mutation_class"}++;
            print STDERR "Skipping $mutation_class mutation: $gene, $chr:$start-$stop\n" if ( $this->{_SHOW_SKIPPED} );
            next;
        }
        #
        # If the mutation type is odd, quit with error
        if ( $mutation_type !~ m/^(SNP|DNP|TNP|ONP|INS|DEL|Consolidated)$/ ) {
            print STDERR "Unrecognized Variant_Type \"$mutation_type\" in MAF file: $gene, $chr:$start-$stop\n";
            print STDERR "Please use TCGA MAF Specification v2.3.\n";
            return undef;
        }
        #
        # Skip mutations that were consolidated into others (E.g. SNP consolidated into a TNP)
        if ( $mutation_type =~ m/^Consolidated$/ ) {
            $skip_cnts{"are consolidated into another"}++;
            print STDERR "Skipping consolidated mutation: $gene, $chr:$start-$stop\n" if ( $this->{_SHOW_SKIPPED} );
            next;
        }
        #
        # Skip mutations that fall completely outside any of the provided regions of interest
        if ( $this->count_bits( $roi_bitmask->{$chr}, $start, $stop ) == 0 ) {
            $skip_cnts{"are outside any ROIs"}++;
            print STDERR "Skipping mutation that falls outside ROIs: $gene, $chr:$start-$stop\n" if ( $this->{_SHOW_SKIPPED} );
            next;
        }
        #
        # Skip mutations whose gene names don't match any of those in the ROI list
        unless( defined $gene_idx{$gene} ) {
            $skip_cnts{"have unrecognized gene names"}++;
            print STDERR "Skipping unrecognized gene name (not in ROI file): $gene, $chr:$start-$stop\n" if ( $this->{_SHOW_SKIPPED} );
            next;
        }
        #
        my $class = '';
        # Check if the mutation is the truncating type, if the user wanted a separate category of those
        if ( $this->{_SEPERATE_TRUNCATIONS} && $mutation_class =~ m/^(Nonsense_Mutation|Splice_Site|Frame_Shift_Del|Frame_Shift_Ins)/ ) {
            $class = Truncations;
        }
        # Classify the mutation as AT/CG/CpG Transition, AT/CG/CpG Transversion
        elsif ( $mutation_type =~ m/^(SNP|DNP|ONP|TNP)$/ ) {
            # ::TBD:: For DNPs and TNPs, we use only the first base for mutation classification
            $ref = substr( $ref, 0, 1 );
            $var1 = substr( $var1, 0, 1 );
            $var2 = substr( $var2, 0, 1 );
            #
            # If the alleles are anything but A, C, G, or T then quit with error
            if ( $ref !~ m/[ACGT]/ || $var1 !~ m/[ACGT]/ || $var2 !~ m/[ACGT]/ ) {
                print STDERR "Unrecognized allele in column Reference_Allele, Tumor_Seq_Allele1, or Tumor_Seq_Allele2: $gene, $chr:$start-$stop\n";
                print STDERR "Please use TCGA MAF Specification v2.3.\n";
                return undef;
            }
            #
            # Use the classify hash to find whether this SNV is an AT/CG Transition/Transversion
            #
            $class = $classify{ "$ref$var1" } if( defined $classify{ "$ref$var1" } );
            $class = $classify{ "$ref$var2" } if( defined $classify{ "$ref$var2" } );
            #
            # Check if the ref base in the MAF matched what we fetched from the ref-seq
            #
            my $locus = "$chr\t$start";
            if ( defined $ref_base{$locus} && $ref_base{$locus} ne $ref ) {
                print STDERR "Reference allele $ref for $gene variant at $chr:$start-$stop is " . $ref_base{$locus} . " in the FASTA. Using it anyway.\n";
            }
            #
            # Check if a C or G reference allele belongs to a CpG pair in the refseq
            #
            if (( $ref eq 'C' || $ref eq 'G' ) && defined $cpg_site{$locus} ) {
                $class = (( $class == CG_Transitions ) ? CpG_Transitions : CpG_Transversions );
            }
        }
        # Classify it as an indel (excludes splice-site and frame-shift if user wanted truncations separately)
        #
        elsif ( $mutation_type =~ m/^(INS|DEL)$/ ) {
            $class = Indels;
        }
        # The user's gene exclusion list affects only 
        # the overall BMR calculations
        $sample_mr[$sample_idx{$sample}][$class][mutations]++ unless( defined $ignored_genes{$gene} );
        $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$class][mutations]++;
    }
    #
    $mafFh->close;
    #
    # Display statistics related to parsing the MAF
    #
    print STDERR "Finished Parsing the MAF file to classify mutations\n";
    foreach my $skip_type ( sort {$skip_cnts{$b} <=> $skip_cnts{$a}} keys %skip_cnts ) {
        print STDERR "Skipped " . $skip_cnts{$skip_type} . " mutation(s) that $skip_type\n" if( defined $skip_cnts{$skip_type} );
    }
    #
    # If the user wants, merge together concurrent mutations of 
    # a gene in the same sample
    if ( $this->{_MERGE_CONCURRENT_MUTS} ) {
        foreach my $sample ( @all_sample_names ) {
            foreach my $gene ( @all_gene_names ) {
                next unless( defined $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}} );
                my $num_muts = 0;
                $num_muts += $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$_][mutations] foreach( @mut_classes );
                if ( $num_muts > 1 ) {
                    foreach my $class ( @mut_classes ) {
                        my $muts_in_class = $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$class][mutations]; # Num of muts of gene in this class
                        $sample_mr[$sample_idx{$sample}][$class][mutations] -= $muts_in_class; # Take it out of the sample total
                        $muts_in_class /= $num_muts; # Turn it into a fraction of the total number of muts in this gene
                        $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$class][mutations] = $muts_in_class; # Use the fraction as the num muts of gene in this class
                        $sample_mr[$sample_idx{$sample}][$class][mutations] += $muts_in_class; # Add the same fraction to the sample total
                    }
                }
            }
        }
    }
    #
    # Calculate per-sample BMRs, and also subtract out covered bases in genes the user wants ignored
    #
    foreach my $sample ( @all_sample_names ) {
        my $tot_muts = 0;
        foreach my $class ( @mut_classes ) {
            # Subtract the covered bases in this class that belong to the genes to be ignored
            # ::TBD:: Some of these bases may also belong to another gene (on the other strand maybe?), and those should not be subtracted
            foreach my $ignored_gene ( keys %ignored_genes ) {
                $sample_mr[$sample_idx{$sample}][$class][covd_bases] -= $gene_mr{$sample_idx{$sample}}{$gene_idx{$ignored_gene}}[$class][covd_bases] if ( defined $gene_mr{$sample_idx{$sample}}{$gene_idx{$ignored_gene}} );
            }
            $tot_muts += $sample_mr[$sample_idx{$sample}][$class][mutations];
        }
        $sample_mr[$sample_idx{$sample}][Overall][bmr] = $tot_muts / $sample_mr[$sample_idx{$sample}][Indels][covd_bases];
    }
    #
    # Cluster samples into bmr-groups using k-means clustering
    #
    my @sample_bmrs = map { $sample_mr[$sample_idx{$_}][Overall][bmr] } @all_sample_names;
    my @bmr_clusters = k_means( $this->{_BMR_GROUPS}, \@sample_bmrs );
    #
    # Calculate overall BMRs for each cluster of samples, 
    # and print them to file
    #
    my %cluster_bmr; # Stores per cluster categorized BMR
    my $totBmrFh = IO::File->new( $overall_bmr_file, ">" ) or die "Couldn't open $overall_bmr_file. $!";
    $totBmrFh->print( "#User-specified genes skipped in these calculations: $this->{_GENES_TO_IGNORE}\n" ) if ( defined $this->{_GENES_TO_IGNORE} );
    my ( $covered_bases_sum, $mutations_sum ) = ( 0, 0 );
    for ( my $i = 0; $i < scalar( @bmr_clusters ); ++$i ) {
        my @samples_in_cluster = map { $all_sample_names[$_] } @{$bmr_clusters[$i]};
        unless( $this->{_BMR_GROUPS} == 1 ) {
            $totBmrFh->print( "#BMR sub-group ", $i + 1, " (", scalar( @{$bmr_clusters[$i]} ), " samples)\n" );
            $totBmrFh->print( "#Samples: ", join( ",", @samples_in_cluster ), "\n" );
        }
        #$totBmrFh->print( "#Mutation_Class\tCovered_Bases\tMutations\tOverall_BMR\n" );
        $totBmrFh->print( "#Mutation_Class\tCovered_Bases\tMutations\tOverall_BMR\tCovered_Bases_per\n" );
        my ( $tot_covd_bases, $tot_muts ) = ( 0, 0 );
        foreach my $class ( @mut_classes ) {
            my ( $covd_bases, $mutations ) = ( 0, 0 );
            foreach my $sample ( @samples_in_cluster ) {
                $covd_bases += $sample_mr[$sample_idx{$sample}][$class][covd_bases];
                $mutations += $sample_mr[$sample_idx{$sample}][$class][mutations];
            }
            $tot_covd_bases = $covd_bases if ( $class == Indels ); # Save this to calculate overall BMR below
            # Calculate overall BMR for this mutation class and print it to file
            $cluster_bmr{$i}[$class][bmr] = ( $covd_bases == 0 ? 0 : ( $mutations / $covd_bases ));

            # added for qunyuan's requirement 
            #
            #
            $cluster_bmr{$i}[$class]["cov"] = $covd_bases;
            #$totBmrFh->print( join( "\t", $mut_class_names[$class], $covd_bases, $mutations, $cluster_bmr{$i}[$class][bmr] ), "\n" );
            $totBmrFh->print( join( "\t", $mut_class_names[$class], $covd_bases, $mutations, $cluster_bmr{$i}[$class][bmr], $cluster_bmr{$i}[$class]["cov"] ), "\n" );

            $tot_muts += $mutations;
        }
        #$totBmrFh->print( join( "\t", "Overall_BMR", $tot_covd_bases, $tot_muts, $tot_muts / $tot_covd_bases ), "\n\n" );
        $totBmrFh->print( join( "\t", "Overall_BMR", $tot_covd_bases, $tot_muts, $tot_muts / $tot_covd_bases, $tot_covd_bases ), "\n\n" );
        $covered_bases_sum += $tot_covd_bases;
        $mutations_sum += $tot_muts;
    }
    $totBmrFh->close;
    #
    #$this->{_BMR_OUTPUT}( $mutations_sum / $covered_bases_sum );
    $this->{_BMR_OUTPUT} = $mutations_sum / $covered_bases_sum ;
    #
    # Print out a file containing per-gene mutation counts and 
    # covered bases for use by "music smg"
    #
    my $geneBmrFh = IO::File->new( $gene_mr_file, ">" ) or die "Couldn't open $gene_mr_file. $!";
    $geneBmrFh->print( "#Gene\tMutation_Class\tCovered_Bases\tMutations\tBMR\tCovered_Bases_per\n" );
    foreach my $gene ( sort @all_gene_names ) {
        my ( $tot_covd_bases, $tot_muts ) = ( 0, 0 );
        for ( my $i = 0; $i < scalar( @bmr_clusters ); ++$i ) {
            my @samples_in_cluster = map { $all_sample_names[$_] } @{$bmr_clusters[$i]};
            foreach my $class ( @mut_classes ) {
                my ( $covd_bases, $mutations ) = ( 0, 0 );
                foreach my $sample( @samples_in_cluster ) {
                    if ( defined $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}} ) {
                        $covd_bases += $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$class][covd_bases];
                        $mutations += $gene_mr{$sample_idx{$sample}}{$gene_idx{$gene}}[$class][mutations];
                    }
                }
                my $rename_class = $mut_class_names[$class];
                $rename_class = ( $rename_class . "_SubGroup" . ( $i + 1 )) if (  $this->{_BMR_GROUPS} > 1 );
                $geneBmrFh->print( join( "\t", $gene, $rename_class, $covd_bases, $mutations, $cluster_bmr{$i}[$class][bmr], $cluster_bmr{$i}[$class]["cov"] ), "\n" );

                $tot_muts += $mutations;
                $tot_covd_bases += $covd_bases if( $class == Indels );
            }
        }
        #$geneBmrFh->print( join( "\t", $gene, "Overall", $tot_covd_bases, $tot_muts, $this->{_BMR_OUTPUT} ), "\n" );
        $geneBmrFh->print( join( "\t", $gene, "Overall", $tot_covd_bases, $tot_muts, $this->{_BMR_OUTPUT}, $covered_bases_sum ), "\n" );
    }
    $geneBmrFh->close;

    return 1;
}
#
# Creates an empty whole genome bitmask based on 
# the given reference sequence index
#
sub create_empty_genome_bitmask {
    my ( $this, $ref_seq_idx_file ) = @_;
    my %genome;
    my $refFh = IO::File->new( $ref_seq_idx_file ) or die "Couldn't open $ref_seq_idx_file. $!";
    while( my $line = $refFh->getline ) {
        my ( $chr, $length ) = split( /\t/, $line );
        # Adding a base for 1-based coordinates
        $genome{$chr} = Bit::Vector->new( $length + 1 ); 
    }
    $refFh->close;
    return \%genome;
}
#
# Given a list of numerical values, returns k 
# clusters based on k-means clustering
sub k_means {
    my ( $k, $list_ref ) = @_;
    my @vals = @{$list_ref};
    my $num_vals = scalar( @vals );

    # Start with the first k values as the centroids
    my @centroids = @vals[0..($k-1)];
    my @prev_centroids = map { 0 } @centroids;
    my @groups = ();

    my $diff_means = 1; # Arbitrary non-zero value
    # Repeat until the difference between these centroids and the previous ones, converges to zero
    while( $diff_means > 0 ) {
        @groups = ();
        # Group values into clusters based on closest centroid
        for( my $i = 0; $i < $num_vals; ++$i ) {
            my @distances = map { abs( $vals[$i] - $_ ) } @centroids;
            my $closest = min( @distances );
            for( my $j = 0; $j < $k; ++$j ) {
                if ( $distances[$j] == $closest ) { push( @{$groups[$j]}, $i ); last; }
            }
        }

        # Calculate means to be the new centroids, and the sum of differences
        $diff_means = 0;
        for( my $i = 0; $i < $k; ++$i ) {
            $centroids[$i] = sum( map {$vals[$_]} @{$groups[$i]} );
            $centroids[$i] /= scalar( @{$groups[$i]} );
            $diff_means += abs( $centroids[$i] - $prev_centroids[$i] );
        }

        # Save the current centroids for comparisons with those in the next iteration
        @prev_centroids = @centroids;
    }
    return @groups;
}

# Counts the number of bits that are set in 
# the given region of a Bit:Vector
sub count_bits {
    my ( $this, $vector, $start, $stop ) = @_;
    my $count = 0;
    for my $pos ( $start..$stop ) {
        ++$count if ( $vector->bit_test( $pos ));
    }
    return $count;
}

## mutation clustering | provide one new MAF file
#
# Mutation clusters potentially from the same event, should be treated together 
#  
# Step 1: Check for samples with > 2 mutations per gene.
# Step 2: Only retain two mutations per sample per gene in this order:

#    1.  Nonsense 
#    2.  Frameshift 
#    3.  Splice-Site (intronic positions 1-2bp from the exon)
#    4.  In_Frame
#    5.  Missense
#    6.  No_Stop,Nonstop,Readthrough

#Step 3. Write a file with the discarded mutations for users to check
#
#
sub mut_clustering {
    my ( $this, $p_maf_file, ) = @_;
    my ( %cluster_hash, @collecting, );
    #
    # keep maximal 2 mutations 
    # for one class
    #
    my $max_number_of_muts_keep = 2;
    my @total_class_types = qw( Nonsense Frameshift Splicesite Inframe Missense Nonstop Other );
    my $p_mafFh = IO::File->new( $p_maf_file ) or die "Couldn't open $p_maf_file. $!";
    while ( my $p_line = $p_mafFh->getline ) {
        if ( $p_line =~ m/^(#|Hugo_Symbol)/ ) { push( @collecting, $p_line ); next; };
        my @p_cols = split( /\t/, $p_line );
        my ( $gene, $v_class, $tumor_barcode ) = ( split /\t/, $p_line )[0,8,15];
        ## convert class type to upcase
        my $v_class_format = {
            'Nonsense'    => 'Nonsense',
            'Frame_Shift' => 'Frameshift',
            'Splice_Site' => 'Splicesite',
            'In_Frame'    => 'Inframe',
            'Missense'    => 'Missense', 
            'No_stop'     => 'Nonstop',
            'Nonstop'     => 'Nonstop',
            'Nonstop_Mutation'     => 'Nonstop',
            'Readthrough' => 'Nonstop',
        }->{$v_class} || "Other";
        $cluster_hash{$tumor_barcode}{$gene}{$v_class_format}{$p_line} = 1;
        $cluster_hash{$tumor_barcode}{$gene}{Number}++;
    }
    ## scanning 
    foreach my $tumor ( keys %cluster_hash ) {
        foreach my $gene ( keys %{$cluster_hash{$tumor}} ) {
            next unless ( defined $cluster_hash{$tumor}{$gene}{Number} );
            if ( ($cluster_hash{$tumor}{$gene}{Number} <= $max_number_of_muts_keep) && ($cluster_hash{$tumor}{$gene}{Number} >= 1) ) {
                foreach my $class ( @total_class_types ) {
                    foreach my $mut ( keys %{$cluster_hash{$tumor}{$gene}{$class}} ) { push( @collecting, $mut ); }
                }
            } else {
                ## select muts based on 
                ###  class priority 
                my $temp_mut_number = 0;
                foreach my $class ( @total_class_types ) {
                    foreach my $mut ( keys %{$cluster_hash{$tumor}{$gene}{$class}} ) {
                        last if ( $temp_mut_number >= $max_number_of_muts_keep );
                        push( @collecting, $mut );
                        $temp_mut_number++;
                    }
                    last if ( $temp_mut_number >= $max_number_of_muts_keep );
                }
            }
        }
    }
    ## return collected muts
    ####
    return \@collecting;
}

## usage
sub help_text {
    my $this = shift;
        return <<HELP
USAGE
 music2 bmr calc-bmr --bmr-output=? --roi-file=? --gene-mr-file=? --reference-sequence=?
    --bam-list=? --output-dir=? --maf-file=? [--skip-non-coding] [--skip-silent] [--bmr-groups=?]
    [--show-skipped] [--separate-truncations] [--merge-concurrent-muts] [--genes-to-ignore=?]

SYNOPSIS
 music2 bmr calc-bmr \
    --bam-list input_dir/bam_list \
    --maf-file input_dir/myMAF.tsv \
    --output-dir output_dir/ \
    --reference-sequence input_dir/all_sequences.fa \
    --roi-file input_dir/all_coding_exons.tsv

 music2 bmr calc-bmr \
    --bam-list input_dir/bam_list \
    --maf-file input_dir/myMAF.tsv \
    --output-dir output_dir/ \
    --reference-sequence input_dir/all_sequences.fa \
    --roi-file input_dir/all_coding_exons.tsv \
    --genes-to-ignore GENE1,GENE2

REQUIRED INPUTS
  roi-file
    Tab delimited list of ROIs [chr start stop gene_name] (See DESCRIPTION) 
  reference-sequence
    Path to reference sequence in FASTA format 
  bam-list
    Tab delimited list of BAM files [sample_name normal_bam tumor_bam] (See DESCRIPTION) 
  output-dir
    Directory where output files will be written (Use the same one used with calc-covg) 
  maf-file
    List of mutations using TCGA MAF specification v2.3 

OPTIONAL INPUTS
  skip-non-coding
    Skip non-coding mutations from the provided MAF file 
    Default value 'true' if not specified
  noskip-non-coding
    Make skip-non-coding 'false' 
  skip-silent
    Skip silent mutations from the provided MAF file 
    Default value 'true' if not specified
  noskip-silent
    Make skip-silent 'false' 
  bmr-groups
    Number of clusters of samples with comparable BMRs (See DESCRIPTION) 
    Default value '1' if not specified
  show-skipped
    Report each skipped mutation, not just how many 
    Default value 'false' (--noshow-skipped) if not specified
  noshow-skipped
    Make show-skipped 'false' 
  separate-truncations
    Group truncational mutations as a separate category 
    Default value 'false' (--noseparate-truncations) if not specified
  noseparate-truncations
    Make separate-truncations 'false' 
  merge-concurrent-muts
    Multiple mutations of a gene in the same sample are treated as 1 
    Default value 'false' (--nomerge-concurrent-muts) if not specified
  nomerge-concurrent-muts
    Make merge-concurrent-muts 'false' 
  genes-to-ignore
    Comma-delimited list of genes to ignore for background mutation rates 
  bp-class-types
    Comma-delimited list of bp types.
    Default value 'AT,CpG,CG'

REQUIRED OUTPUTS
  bmr-output
    TODO 
  gene-mr-file
    TODO 

DESCRIPTION
    Given a mutation list (MAF), and per-gene coverage data calculated using "music bmr
    calc-covg"), this script calculates overall Background Mutation Rate (BMR) and BMRs in the
    categories of AT/CG/CpG Transitions, AT/CG/CpG Transversions, and Indels. An optional category
    for truncational mutations can also be specified. The script generates a file with per-gene
    mutation rates that can be used with the tool that tests for significantly mutated genes (music
    smg).


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

    --bmr-groups

      Ideally, we want to test the mutation rate (MR) of a gene in a sample against the background
      mutation rate (BMR) across that sample. But if the BMRs of some samples are comparable, we
      can instead test the MR of a gene across a group of samples with comparable BMR, against the
      overall BMR of that group. This argument specifies how many such groups you want to cluster
      samples into. By default, it is assumed that all samples have comparable BMRs (bmr-groups =
      1).

    --output-dir

      This should be the same output directory used when running "music bmr calc-covg". The
      following outputs of this script will also be created/written: overall_bmrs: File containing
      categorized overall background mutation rates. gene_mrs: File containing categorized per-gene
      mutation rates.

    --genes-to-ignore

      A comma-delimited list of genes to ignore for overall BMR calculations. List genes that are
      known factors in this disease and whose mutations should not be classified as background.

HELP

}

1;

