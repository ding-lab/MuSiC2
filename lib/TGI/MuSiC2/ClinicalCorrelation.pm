package TGI::MuSiC2::ClinicalCorrelation;

use warnings;
use strict;
use Carp;
use IO::File;
use POSIX qw( WIFEXITED );
use Getopt::Long;
use File::Temp qw/ tempfile /;

# TODO: 
#
# work through logic of defaults for use-maf-in-glm, skip-non-coding, skip-silent
# since =! does not seem to work...  make docs consistent.

# Removing non-glm functionality

# For GLM, workflow should be 1) create mutation matrix file from MAF and 2) perform correlation.
# These steps should be separate from the user's point of view; keeping them together leads to unnecessary
# parameters and complicated logic.
# Testing of Mutation Matrix creation: /gscuser/mwyczalk/projects/Virus/Virus_2013.9a/analysis/Music/Cordat/analysis/test
# implementation /gscuser/mwyczalk/projects/Virus/Virus_2013.9a/analysis/Music/Cordat/analysis/1_
# Note that this is part of regular processing, UnifiedVirus2/A_MAF_data/3_make_mutation_matrix.sh
#       - just use that for testing

# break up the logic into more manageable and testable pieces

# consider in clinical correlation only those participants which are found in trait and variant data
# see /Users/mwyczalk/Data/MuSiC2/MuSiC.ClinicalCorrelation.TCGA_Virus/src/VariantTraitMaker.R

# treat B and Q together. Make sure they have the same output format

# There seem to be random numerical errors during non-convergence.  Catch these cases and don't print them out
# Types of errors we observe which should be handled cleanly:
# * Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
# * Warning: glm.fit: algorithm did not converge
# * Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) :
#     contrasts can be applied only to factors with 2 or more levels
#       Error caught, continuing.  yi = clinical_M.is.M0  xi = BK.int  covi = Disease

# Use Smg.pm as model for MuSiC2 patterns

# working through mutation matrix creation
# Written to OUTPUT_VARIANT_FILE



sub new {
    my $class = shift;
    my $this = {};

# this updated
    $this->{BAM_LIST} = undef;
    $this->{MAF_FILE} = undef;
    $this->{OUTPUT_FILE} = undef;
    $this->{OUTPUT_VARIANT_FILE} = undef;
    $this->{INPUT_VARIANT_FILE} = undef;
    $this->{GENETIC_DATA_TYPE} = "gene";
    $this->{MODEL_FILE} = undef;
    $this->{TRAIT_DATA_FILE} = undef;
    $this->{USE_MAF_IN_GLM} = 0;
    $this->{SKIP_NON_CODING} = 1;  
    $this->{SKIP_SILENT} = 1;
    $this->{SKIP_CORRELATION} = 0;

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
        'output-file=s'                             => \$this->{OUTPUT_FILE},
        'output-variant-data-file=s'        => \$this->{OUTPUT_VARIANT_FILE},
        'variant-data-file=s'  => \$this->{INPUT_VARIANT_FILE},
        'genetic-data-type=s'                       => \$this->{GENETIC_DATA_TYPE},
        'model-file=s'                          => \$this->{MODEL_FILE},
        'trait-data-file=s'                  => \$this->{TRAIT_DATA_FILE},
        'use-maf-in-glm'                          => \$this->{USE_MAF_IN_GLM},
        'skip-non-coding'                         => \$this->{SKIP_NON_CODING},
        'skip-silent'                             => \$this->{SKIP_SILENT},
        #'use-maf-in-glm=!'                          => \$this->{USE_MAF_IN_GLM},
        #'skip-non-coding=!'                         => \$this->{SKIP_NON_CODING},
        #'skip-silent=!'                             => \$this->{SKIP_SILENT},
        'skip-correlation'                          => \$this->{SKIP_CORRELATION},
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
    my $output_file = $this->{OUTPUT_FILE};
    my $genetic_data_type = $this->{GENETIC_DATA_TYPE};
    my $glm_model = $this->{MODEL_FILE};

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

    my $full_output_filename;

    $full_output_filename = $output_file;  # User specifies filename exactly.

    #read through clinical data file to see which samples are represented and create input matrix for R
    # This section is not necessary if this->skip_correlation, and makes creating a mutation_matrix file annoying.

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
    
    # create correlation matrix unless it's glm analysis without using a maf file
    # contents of all_sample_names (and bam_list) ignored if GLM and not using MAF.  This is why 
    # bam_list argument should be optional. -- really?  
    # Empirially, content of bam_list changes nothing, but @all_sample_names is propagaged
    # Also, there are times when I want to create a matrix file and save it without doing processing.
    # This is currently not supported.

    # recall: unless(x) == if (!x)
    unless( !$this->{USE_MAF_IN_GLM} || $this->{INPUT_VARIANT_FILE} ) {

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
    }

    # TODO want to have the option to exit loop here, after mutation matrix has been constructed and written to disk.

    next if ( $this->{SKIP_CORRELATION} );  # we want this to work in future
    #warn("Skipping correlation, hard coded.");  # uncomment to stop after mutation matrix created
    #next;                                       # uncomment to stop after mutation matrix created


    if( $this->{INPUT_VARIANT_FILE} ) {
        $matrix_file = $this->{INPUT_VARIANT_FILE};
    }
    unless( defined $matrix_file ) { $matrix_file = "'*'"; }   
    ### What does matrix_file = '*' mean?
    # from R code, with x.file = matrix_file: 
    #   x.file="*" if x.file already merged into y.file

    #set up R command
    my $R_cmd = "R --slave --args < " . __FILE__ . ".R glm ";

    $R_cmd .= "$glm_model $this->{TRAIT_DATA_FILE} $matrix_file $full_output_filename";

    #run R command
    print "R_cmd:\n$R_cmd\n";
    WIFEXITED( system $R_cmd ) or croak "Couldn't run: $R_cmd ($?)";

    return( 1 );
}

sub create_sample_gene_matrix_gene {
# Writes to OUTPUT_VARIANT_FILE, or a temp file if not defined.
# Returns filename of file written to

    my ( $this, $samples, @all_sample_names ) = @_;
    my $output_matrix = $this->{OUTPUT_VARIANT_FILE};  # test this to make sure it works.  Did it ever?

    #create a hash of mutations from the MAF file
    my ( %mutations, %all_genes, @all_genes );

    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $this->{MAF_FILE} ) or die "Couldn't open MAF file!\n";
    while( my $line = $maf_fh->getline ) {
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $mutation_class, $sample ) = @cols[0,8,15];

        #check that the mutation class is valid
        if( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            warn( "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\nPlease use TCGA MAF v2.3.\n" );
            return;
        }

        #check if sample exists in clinical data
        # this is not helpful when creating a mutation matrix from MAF file.
        # presumably, user has curated MAF file appropriately.  This leads to problems when hoping to use one 
        # MAF file to create multiple clinical datasets. 
#        unless( defined $samples->{$sample} ) {
#            warn "Sample Name: $sample from MAF file does not exist in Clinical Data File";
#            next;
#        }

        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if(( $this->{SKIP_NON_CODING} && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
           ( $this->{SKIP_SILENT} && $mutation_class =~ m/^Silent$/ )) {
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
    my $matrix_fh;
    unless (defined $output_matrix) {
        $output_matrix = Genome::Sys->create_temp_file_path();  ### this is a problem, not in Genome anymore
    }
    $matrix_fh = new IO::File $output_matrix,"w";
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

    return $output_matrix;
} # create_sample_gene_matrix_gene 

sub create_sample_gene_matrix_variant {

    my ( $this, $samples, @all_sample_names ) = @_;
    my $output_matrix = $this->{OUTPUT_VARIANT_FILE};

    #create hash of mutations from the MAF file
    my ( %variants_hash, %all_variants );

    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $this->{MAF_FILE} ) or die "Couldn't open MAF file!\n";
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
        if(( $this->{SKIP_NON_CODING} && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
           ( $this->{SKIP_SILENT} && $mutation_class =~ m/^Silent$/ )) {
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
    unless( defined $output_matrix ) {
        $output_matrix = Genome::Sys->create_temp_file_path();
    }
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
 music2 clinical-correlation --gene-mr-file=? --output-file=? [--max-fdr=?] [--skip-low-mr-genes]
    [--downsample-large-genes] [--bmr-modifier-file=?] [--processors=?]

SYNOPSIS 
 music2 clinical-correlation \\  # GLM test
        --maf-file /path/myMAF.tsv \\
        --bam-list /path/myBamList.tsv \\
        --output-file /path/output_file \\
        --model-file /path/model.tsv \\
        --trait-data-file /path/trait_data.tsv \\
        --use-maf-in-glm

SEE ALSO

 music2 clinical-correlation --help for details.

HELP
}

sub help_text {
    my $this = shift;
    my $usage = usage_text();

    return $usage . <<HELP

REQUIRED INPUTS 
  output-file 
    Results of clinical-correlation tool. 

OPTIONAL INPUTS 
    bam-list 
        Tab delimited list of BAM files [sample-name, normal-bam, tumor-bam], optional if mutation matrix is supplied. (See Description)
    maf-file 
        List of mutations using TCGA MAF specification v2.3
    output-variant-data-file 
        Specify a file to store the sample-vs-gene matrix created during calculations
    variant-data-file 
        Instead of creating this from the MAF, input the sample-vs-gene matrix for calculations
        *** expand on this.  Typical flow is this is NOT a maf ***
    genetic-data-type 
        Correlate clinical trait data to "gene" or "variant" level data.  Default "gene"
    model-file 
        File outlining the type of model, response variable, covariants, etc. for the GLM analysis. (See DESCRIPTION)
    trait-data-file 
        Clinical traits, mutational profiles, other traits.  May be numerical or categorical (See DESCRIPTION)

# these are specific to creating variant file from MAF
    use-maf-in-glm 
        Create a variant matrix from the MAF file as variant input to GLM analysis.  
        Default value 'false' if not specified, add prefix "no" to negate
    skip-non-coding 
        Skip non-coding mutations from the provided MAF file. 
        Default value 'true' if not specified, add prefix "no" to negate
    skip-silent 
        Skip silent mutations from the provided MAF file.  
        Default value 'true' if not specified, add prefix "no" to negate

    skip_correlation 
        Do not perform correlation analysis, but exit after writing clinical_correlation_matrix_file.  
        Default value 'false' if not specified, add prefix "no" to negate

DESCRIPTION
    This command relates clinical traits and mutational data using a generalized linear 
    model (GLM) analysis.  Traditionally mutations come from a MAF and particular phenotypic 
    traits recorded in clinical data files. However, we support correlations between 
    any type of variant and trait, including numerical and categorical data for both.

    Internally, the input data is fed into an R script which calculates a
    P-value representing the probability that the correlation seen between the
    mutations in each gene (or variant) and each phenotype trait are random. Lower
    P-values indicate lower randomness, or likely true correlations.

    VARIANT -- TODO, update
    The analysis accepts a mixed numeric and categoric variant data file,
    defined by --variant-data-file, with these conventions: TODO confirm
        * Headers are required, columns are tab-separated, leading # indicates comment
        * Each file must include at least 1 sample_id column and 1 attribute
          column, with the format being [sample_id  trait_1  trait_2 ... ]
        * The sample ID must match the sample ID listed in the MAF under
          "Tumor_Sample_Barcode" for relating the mutations of this sample.
    variant data can be created from a MAF with the command, 
        music2 variant-from-maf ...

    TRAIT
    The GLM analysis accepts a mixed numeric and categoric trait data file,
    defined by --trait-data-file, with these conventions:
        * Headers are required, columns are tab-separated, leading # indicates comment
        * Each file must include at least 1 sample_id column and 1 attribute
          column, with the format being [sample_id  trait_1  trait_2 ... ]
        * The sample ID must match the sample ID listed in the MAF under
          "Tumor_Sample_Barcode" for relating the mutations of this sample.
    The header for each trait will appear in the output file to
    denote relationships with the variant data.

    MODEL
    GLM also requires the user to input a --model-file with the following conventions:
        * Headers are required, columns are tab-separated, leading # indicates comment
        * The columns are:
          [ analysis_type    trait_name    variant/gene_name   covariates  memo ]
        * The 'analysis_type' column must contain either "Q", indicating a
          quantative trait, or "B", indicating a binary trait will be examined.
        * The 'trait_name' is the name of a clinical data trait
          defined by a header in the --trait-data-file.

        *** simplify this ***  Why can this be either variant or trait name?
        * The 'variant/gene_name' can either be the name of one or more columns
          from the --trait-data-file, or the name of one or more mutated gene
          names from the MAF, separated by "|". If this column is left blank, or instead
          contains "NA", then each column from either the variant mutation matrix
          (--use-maf-in-glm) or alternatively the --trait-data-file is used
          consecutively as the variant column in independent analyses.

        * 'covariates' are the names of one or more columns from the
          --trait-data-file, separated by "+".
        * 'memo' is any note deemed useful to the user. It will be printed in the
          output data file for reference.

    GLM analysis may be performed using solely the data input into
    --trait-data-file, as described above, or alternatively, mutational data
    from the MAF may be included as variants in the GLM analysis, as also described
    above. Use the --use-maf-in-glm flag to include the mutation matrix derived
    from the maf as variant data.

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

