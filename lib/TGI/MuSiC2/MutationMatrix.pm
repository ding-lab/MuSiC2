package TGI::MuSiC2::MutationMatrix;

use warnings;
use strict;
use Carp;
use IO::File;
use POSIX qw( WIFEXITED );
use Getopt::Long;
use File::Temp qw/ tempfile /;

# Extracting mutation matrix creation functionality from ClinicalCorrelation.pm


# Testing of Mutation Matrix creation: /gscuser/mwyczalk/projects/Virus/Virus_2013.9a/analysis/Music/Cordat/analysis/test
# implementation /gscuser/mwyczalk/projects/Virus/Virus_2013.9a/analysis/Music/Cordat/analysis/1_
# Note that this is part of regular processing, UnifiedVirus2/A_MAF_data/3_make_mutation_matrix.sh
#       - just use that for testing

sub new {
    my $class = shift;
    my $this = {};

# this updated
    $this->{BAM_LIST} = undef;
    $this->{MAF_FILE} = undef;
    $this->{OUTPUT_VARIANT_FILE} = undef;
    $this->{GENETIC_DATA_TYPE} = "gene";
    $this->{TRAIT_DATA_FILE} = undef;
    $this->{SKIP_NON_CODING} = 1;  
    $this->{SKIP_SILENT} = 1;

    bless $this, $class;
    $this->process();

    return $this;
}

# Standardizing nomenclature: 
#  * Variant (aka X, Mutation Matrix, Clinical Correlation Matrix)
#  * Trait (aka Y, Clinical Data)


sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }

    $options = GetOptions (  # http://search.cpan.org/~chips/perl5.004_05/lib/Getopt/Long.pm
        'bam-list=s'                                => \$this->{BAM_LIST},
        'maf-file=s'                                => \$this->{MAF_FILE},
        'output-variant-data-file=s'        => \$this->{OUTPUT_VARIANT_FILE},
        'genetic-data-type=s'                       => \$this->{GENETIC_DATA_TYPE},
        'trait-data-file=s'                  => \$this->{TRAIT_DATA_FILE},
        'skip-non-coding'                         => \$this->{SKIP_NON_CODING},
        'skip-silent'                             => \$this->{SKIP_SILENT},
        'help' => \$help,
    );

    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->usage_text(); }




### Implementation notes ###
# conversion GMS -> MuSiC2
# * convert $self to $this
# * error_message is {warn(""); return}
# * debug_message is {warn(""); next}
# * WIFEXITED is a wrapper around system calls (unchanged from GMS)
# * Errors (e.g. file open) trigger a die()

    # parse input arguments
    my $bam_list = $this->{BAM_LIST};
    my $genetic_data_type = $this->{GENETIC_DATA_TYPE};

    # check genetic data type
    unless( $genetic_data_type =~ /^gene|variant$/i ) {
        warn("Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter.");
        return;
    }

    # declarations
    my @all_sample_names; # names of all the samples, no matter if it's mutated or not

    # parse out the sample names from the bam-list which should match the names in the MAF file
    # bam_list content saved in @all_sample_names.  
    if ($bam_list) { 
        my $sampleFh = IO::File->new( $bam_list ) or die "Couldn't open $bam_list. $!\n";
        while( my $line = $sampleFh->getline ) {
            next if ( $line =~ m/^#/ );
            chomp( $line );
            my ( $sample ) = split( /\t/, $line );
            push( @all_sample_names, $sample );
        }
        $sampleFh->close;
    }

    # note that only one clinical data file supported


    #read through clinical data file to see which samples are represented and create input matrix for R

    my %samples;
    my $matrix_file;
    my $samples = \%samples;
    my $clin_fh = new IO::File $this->{TRAIT_DATA_FILE},"r";
    unless( $clin_fh ) {
        die "failed to open $this->{TRAIT_DATA_FILE} for reading: $!";
    }
    my $header = $clin_fh->getline;
    while( my $line = $clin_fh->getline ) {
        chomp $line;
        my ( $sample ) = split( /\t/, $line );
        $samples{$sample}++;
    }
    
    my $output_matrix = $this->{OUTPUT_VARIANT_FILE};  # this must be defined
    my $maf_file = $this->{MAF_FILE};
    my $skip_non_coding = $this->{SKIP_NON_CODING};
    my $skip_silent = $this->{SKIP_SILENT};

    if( $genetic_data_type =~ /^gene$/i ) {
        $matrix_file = $this->create_sample_gene_matrix_gene( $samples, @all_sample_names );
    }
    elsif( $genetic_data_type =~ /^variant$/i ) {
        $matrix_file = $this->create_sample_gene_matrix_variant( $samples, @all_sample_names );
    }
    else {
        warn( "Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter." );
        return;
    }

    return( 1 );
}

sub create_sample_gene_matrix_gene {

    my ( $this, $samples, @all_sample_names, $maf_file, $output_matrix, $skip_non_coding, $skip_silent ) = @_;

    #create a hash of mutations from the MAF file
    my ( %mutations, %all_genes, @all_genes );

    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $maf_file ) or die "Couldn't open MAF file!\n";
    while( my $line = $maf_fh->getline ){
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $mutation_class, $sample ) = @cols[0,8,15];

        #check that the mutation class is valid
        if( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            warn( "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\nPlease use TCGA MAF v2.3.\n" );
            return;
        }

        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if(( $skip_non_coding && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
           ( $skip_silent && $mutation_class =~ m/^Silent$/ )) {
#            print "Skipping $mutation_class mutation in gene $gene.\n";   # annoying - creates millions of useless lines
            next;
        }

        $all_genes{$gene}++;
        $mutations{$sample}{$gene}++;
    }
    $maf_fh->close;

    #sort @all_genes for consistency
    @all_genes = sort keys %all_genes;

    #write the input matrix for R code to a file
    my $matrix_fh = new IO::File $output_matrix,"w";
    unless ($matrix_fh) {
        die "Failed to create matrix file $output_matrix!: $!";
    }
    #print input matrix file header
    my $header = join("\t","Sample",@all_genes);
    $matrix_fh->print("$header\n");

    #print mutation relation input matrix
    for my $sample (sort @all_sample_names) {
        $matrix_fh->print($sample);
        for my $gene (@all_genes) {
            if (exists $mutations{$sample}{$gene}) {
                $matrix_fh->print("\t1");
            }
            else {
                $matrix_fh->print("\t0");
            }
        }
        $matrix_fh->print("\n");
    }

} # create_sample_gene_matrix_gene 

sub create_sample_gene_matrix_variant {

    my ( $this, $samples, @all_sample_names, $maf_file, $output_matrix, $skip_non_coding, $skip_silent ) = @_;

    #create hash of mutations from the MAF file
    my ( %variants_hash, %all_variants );

    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $maf_file ) or die "Couldn't open MAF file!\n";
    while( my $line = $maf_fh->getline ) {
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $chr, $start, $stop, $mutation_class, $mutation_type, $ref, $var1, $var2, $sample ) = @cols[0,4..6,8..12,15];

        #check that the mutation class is valid
        if( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            warn( "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\nPlease use TCGA MAF v2.3.\n" );
            return;
        }

        unless( exists $samples->{$sample} ) {
            warn "Sample Name: $sample from MAF file does not exist in Clinical Data File";
            next;
        }

        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if(( $skip_non_coding && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
           ( $skip_silent && $mutation_class =~ m/^Silent$/ )) {
            print "Skipping $mutation_class mutation in gene $gene.\n";
            next;
        }

        my $var;
        my $variant_name;
        if( $ref eq $var1 ) {
            $var = $var2;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        }
        elsif( $ref eq $var2 ) {
            $var = $var1;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        }
        elsif( $ref ne $var1 && $ref ne $var2 ) {
            $var = $var1;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
            $var = $var2;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        }
    }
    $maf_fh->close;

    #sort variants for consistency
    my @variant_names = sort keys %all_variants;

    #write the input matrix for R code to a file
    my $matrix_fh;
    $matrix_fh = new IO::File $output_matrix,"w";
    unless( $matrix_fh ) {
        die "Failed to create matrix file $output_matrix!: $!";
    }

    #print input matrix file header
    my $header = join( "\t", "Sample", @variant_names );
    $matrix_fh->print("$header\n");

    #print mutation relation input matrix
    for my $sample ( sort @all_sample_names ) {
        $matrix_fh->print( $sample );
        for my $variant ( @variant_names ) {
            if( exists $variants_hash{$sample}{$variant} ) {
                $matrix_fh->print("\t$variants_hash{$sample}{$variant}");
            }
            else {
                $matrix_fh->print("\t0");
            }
        }
        $matrix_fh->print("\n");
    }

    return $output_matrix;
}   # create_sample_gene_matrix_variant 

# MuSiC2 help
# taken from Smg.pm

sub usage_text {
    my $this = shift;
        return <<HELP

USAGE (OLD)
 music2 clinical-correlation --gene-mr-file=? [--max-fdr=?] [--skip-low-mr-genes]
    [--downsample-large-genes] [--bmr-modifier-file=?] [--processors=?]

SYNOPSIS 
 music2 clinical-correlation \\  # GLM test
        --maf-file /path/myMAF.tsv \\
        --bam-list /path/myBamList.tsv \\
        --trait-data-file /path/trait_data.tsv \\

SEE ALSO

 music2 clinical-correlation --help for details.

HELP
}

sub help_text {
    my $this = shift;
    my $usage = usage_text();

    return $usage . <<HELP

REQUIRED INPUTS 

OPTIONAL INPUTS 
    bam-list 
        Tab delimited list of BAM files [sample-name, normal-bam, tumor-bam]. (See DESCRIPTION)
    maf-file 
        List of mutations using TCGA MAF specification v2.3
    output-variant-data-file 
        Specify a file to store the sample-vs-gene matrix created during calculations
    genetic-data-type 
        Correlate clinical trait data to "gene" or "variant" level data.  Default "gene"
    trait-data-file 
        Clinical traits, mutational profiles, other traits.  May be numerical or categorical (See DESCRIPTION)
    skip-non-coding 
        Skip non-coding mutations from the provided MAF file. 
        Default value 'true' if not specified, add prefix "no" to negate
    skip-silent 
        Skip silent mutations from the provided MAF file.  
        Default value 'true' if not specified, add prefix "no" to negate

DESCRIPTION
    P-values indicate lower randomness, or likely true correlations.

ARGUMENTS

    --bam-list

    Provide a file containing sample names and normal/tumor BAM locations for each.
    Use the tab- delimited format [sample_name normal_bam tumor_bam] per line. This
    tool only needs sample_name, so all other columns can be skipped. The
    sample_name must be the same as the tumor sample names used in the MAF file
    (16th column, with the header Tumor_Sample_Barcode).

HELP

}

1;

