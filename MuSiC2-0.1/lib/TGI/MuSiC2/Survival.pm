package TGI::Music2::Survival;

use warnings;
use strict;

use Carp;
use IO::File;
use POSIX qw( WIFEXITED );
use File::Temp qw/ tempfile /;

our $VERSION = 'v0.1';

sub new {
    my $class = shift;
    my $this = {};

    $this->{'bam_list'} = undef;
    $this->{'output_dir'} = undef;
    $this->{'maf_file'} = undef;
    $this->{'numeric_clinical_data_file'} = undef;
    $this->{'categorical_clinical_data_file'} = undef;
    $this->{'glm_clinical_data_file'} = undef;
    $this->{'genetic_data_type'} = 'gene';
    $this->{'phenotypes_to_include'} = undef;
    $this->{'legend_placement'} = 'bottomleft';
    $this->{'skip_non_coding'} = 1;
    $this->{'noskip_non_coding'} = 0;
    $this->{'skip_silent'} = 1;
    $this->{'noskip_silent'} = 0;

    bless $this, $class;
    $this->process();

    return $this;
}

sub process {

    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'bam-list=s'                     => \$this->{'bam_list'},
        'output-dir=s'                   => \$this->{'output_dir'},
        'maf-file=s'                     => \$this->{'maf_file'},
        'numeric-clinical-data-file'     => \$this->{'numeric_clinical_data_file'},
        'categorical-clinical-data-file' => \$this->{'categorical_clinical_data_file'},
        'glm-clinical-data-file'         => \$this->{'glm_clinical_data_file'},
        'genetic-data-type'              => \$this->{'genetic_data_type'},
        'phenotypes-to-include'          => \$this->{'phenotypes_to_include'},
        'legend-placement'               => \$this->{'legend_placement'},
        'skip-non-coding'                => \$this->{'skip_non_coding'},
        'noskip-non-coding'              => \$this->{'noskip_non_coding'},
        'skip-silent'                    => \$this->{'skip_silent'},
        'noskip-silent'                  => \$this->{'noskip_silent'},

        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    #### processing ####
    # handle phenotype inclusions
    my ( @phenotypes_to_include, @clinical_phenotypes_to_include, @mutated_genes_to_include, );
    if ($this->{'phenotypes_to_include'}) { @phenotypes_to_include = split /,/, $this->{'phenotypes_to_include'} }
    # check genetic data type
    unless ($this->{'genetic_data_type'} =~ /^gene|variant$/i) {
        warn ("Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter.");
        return;
    }
    # load clinical data and analysis types
    my %clinical_data;
    if ($this->{'numeric_clinical_data_file'}} { $clinical_data{'numeric'} = $this->{'numeric_clinical_data_file'} }
    if ($this->{'categorical_clinical_data_file'}) { $clinical_data{'categ'} = $self->{'categorical_clinical_data_file'} }
    if ($this->{'glm_clinical_data_file'}) { $clinical_data{'glm'} = $this->{'glm_clinical_data_file'} }
    # create array of all sample names possibly included from clinical data and MAF
    my $sampleFh = IO::File->new( $this->{'bam_list'} ) or die "Couldn't open $this->{'bam_list'}. $!\n";
    my @all_sample_names = map{ unless(/^#/){ chomp; (split /\t/) } } $sampleFh->getlines;
    $sampleFh->close;
    # loop through clinical data files and assemble survival data hash (vital_status and days_to_last_followup required);
    my (%survival_data, $vital_status_flag, $days_to_last_follow_flag, );
    $vital_status_flag = $days_to_last_follow_flag = 0;
    map {
        my $clin_fh = new IO::File $clinical_data{$_}, "r";
        unless ($clin_fh) { warn( "Failed to open $clinical_data{$_} for reading: $!" ); return; }
        #initiate variables to hold column info
        my ( %phenotypes_to_print, $vital_status_col, $days_to_last_follow_col, $i, );
        $vital_status_col = $days_to_last_follow_col = $i = 0;
        #parse header and record column locations for needed data
        my @t_array = split /\t/, $clin_fh->getline;
        shift( @t_array );
        map {
            $i++;
            if (/vital_status|vitalstatus/i) { $vital_status_col = $i; $vital_status_flag++; }
            if (/days_to_last_(follow_up|followup)|daystolastfollowup/i) { $days_to_last_follow_col = $i; $days_to_last_follow_flag++; }
            if (scalar grep { /^$_$/i } @phenotypes_to_include ) { $phenotypes_to_print{$_} = $i; }
        } @t_array;
        while ( my $line = $clin_fh->getline ) {
            chomp $line;
            my @fields = split /\t/,$line;
            my $sample = $fields[0];
            unless (scalar grep { m/^$sample$/ } @all_sample_names) { warn( "Skipping sample $sample. (Sample is not in --bam-list)." ); next; }
            if ( $vital_status_col ) {
                my $vital_status;
                if ($fields[$vital_status_col] =~ /^(0|living)$/i) { $vital_status = 0; }
                elsif ($fields[$vital_status_col] =~ /^(1|deceased)$/i) { $vital_status = 1; }
                else { $vital_status = "NA"; }
                $survival_data{$sample}{'vital_status'} = $vital_status;
            }
            if ($days_to_last_follow_col) { $survival_data{$sample}{'days'} = $fields[$days_to_last_follow_col]; }
            for my $pheno (keys %phenotypes_to_print) { $survival_data{$sample}{$pheno} = $fields[$phenotypes_to_print{$pheno}]; }
        }
        $clin_fh->close;
        # record phenotypes included from clinical data
        push @clinical_phenotypes_to_include, keys %phenotypes_to_print;
    } keys %clinical_data;
    # check for necessary header fields
    unless ($vital_status_flag) { warn( 'Clinical data does not seem to contain a column labeled "vital_status".' ); return; }
    unless ($days_to_last_follow_flag) { warn( 'Clnical data does not seem to contain a column labeled "days_to_last_followup".' ); return; }
    # create temporary files for R command
    #
    my ( undef, $survival_data_file ) = tempfile();
    my ( undef, $mutation_matrix ) = tempfile();
    my $surv_fh = new IO::File $survival_data_file, "w" or die "Couldn't open survival data filehandle.";
    print $surv_fh join( "\t","Sample","Days_To_Last_Followup","Vital_Status" );
    if (@clinical_phenotypes_to_include) { print $surv_fh "\t" . join("\t", @clinical_phenotypes_to_include); }
    print $surv_fh "\n";
    map {
        my $sample = $_;
        unless (exists $survival_data{$sample}{'days'}) { $survival_data{$sample}{'days'} = "NA"; }
        unless (exists $survival_data{$sample}{'vital_status'}) { $survival_data{$sample}{'vital_status'} = "NA"; }
        print $surv_fh join("\t",$sample,$survival_data{$sample}{'days'},$survival_data{$sample}{'vital_status'});
        map {
            my $pheno = $_;
            unless (exists $survival_data{$sample}{$pheno}) { $survival_data{$sample}{$pheno} = "NA"; }
            print $surv_fh "\t" . $survival_data{$sample}{$pheno};
        } @clinical_phenotypes_to_include;
        print $surv_fh "\n";
    } keys %survival_data;
    $surv_fh->close;

    # find if any of the "phenotypes_to_include" are genes, and if so, limit the MAF mutation matrix to those genes
    my %clinical_pheno_to_include;
    @clinical_pheno_to_include{ @clinical_phenotypes_to_include } = ();
    map { push @mutated_genes_to_include, $_ unless exists $clinical_pheno_to_include{$_}; } @phenotypes_to_include;
    my $mutated_genes_to_include = \@mutated_genes_to_include;
    # create mutation matrix file
    if ( $genetic_data_type =~ /^gene$/i ) {
        $this->create_sample_gene_matrix_gene( $mutation_matrix, $mutated_genes_to_include, @all_sample_names );
    } elsif ( $genetic_data_type =~ /^variant$/i ) {
        $this->create_sample_gene_matrix_variant( $mutation_matrix, $mutated_genes_to_include, @all_sample_names );
    } else {
        warn( "Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter." );
        return;
    }
    # check and prepare output directory
    my $output_dir = $this->{'output_dir'} . "/";
    unless (-e $output_dir) {
        warn( "Creating output directory: $output_dir..." );
        unless(mkdir $output_dir) { warn( "Failed to create output directory: $!" ); return; }
    }
    # set up R command
    my $R_cmd = "R --slave --args < " . __FILE__ . ".R " . join( " ", $survival_data_file, $mutation_matrix, $legend_placement, $output_dir );
    print "R_cmd:\n$R_cmd\n";
    #run R command
    WIFEXITED( system $R_cmd ) or croak "Couldn't run: $R_cmd ($?)";

    return(1);
}


sub create_sample_gene_matrix_gene {
    my ( $this, $mutation_matrix, $mutated_genes_to_include, @all_sample_names ) = @_;
    #create hash of mutations from the MAF file
    my ( %mutations, %all_genes );
    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $this->{'maf_file'} ) or die "Couldn't open MAF file!\n";
    while ( my $line = $maf_fh->getline ) {
        next if ( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line; my @cols = split( /\t/, $line );
        my ( $gene, $mutation_class, $sample ) = @cols[0,8,15];
        #check that the mutation class is valid
        if ($mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            warn( "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\nPlease use TCGA MAF v2.3.\n" );
            return;
        }
        #check to see if this gene is on the list (if there is a list at all)
        if ( defined($mutated_genes_to_include) and @{$mutated_genes_to_include} ) { next unless( scalar grep { m/^$gene$/ } @{$mutated_genes_to_include} ) }
        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if ( ($this->{'skip_non_coding'} && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/) || ($this->{'skip_silent'} && $mutation_class =~ m/^Silent$/) ) {
            print "Skipping $mutation_class mutation in gene $gene.\n";
            next;
        }
        $all_genes{$gene}++;
        $mutations{$sample}{$gene}++;
    }
    $maf_fh->close;
    #sort @all_genes for consistency in header and loops
    my @all_genes = sort keys %all_genes;
    #write the input matrix for R code to a temp file
    my $matrix_fh = new IO::File $mutation_matrix,"w" or die "Failed to create matrix file $mutation_matrix!: $!";
    #print input matrix file header
    my $header = join( "\t", "Sample", @all_genes );
    $matrix_fh->print( "$header\n" );
    #print mutation relation input matrix
    map { 
        my $sample = $_; 
        $matrix_fh->print( $sample );
        map { if ( exists $mutations{$sample}{$gene} ) { $matrix_fh->print( "\tMutated" ); } else { $matrix_fh->print( "\tWildtype" ); } } @all_genes;
        $matrix_fh->print( "\n" );
    } sort @all_sample_names;

}

sub create_sample_gene_matrix_variant {

    my ( $this, $mutation_matrix, $mutated_genes_to_include, @all_sample_names ) = @_;
    #create hash of mutations from the MAF file
    my ( %variants_hash, %all_variants );
    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $this->{'maf_file'} ) or die "Couldn't open MAF file!\n";
    while( my $line = $maf_fh->getline ) {
        next if ( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $chr, $start, $stop, $mutation_class, $mutation_type, $ref, $var1, $var2, $sample ) = @cols[0,4..6,8..12,15];
        #check that the mutation class is valid
        if ( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            warn( "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\nPlease use TCGA MAF v2.3.\n" );
            return;
        }
        #check to see if this gene is on the list (if there is a list at all)
        if ( defined @{$mutated_genes_to_include} ) { next unless (scalar grep { m/^$gene$/ } @{$mutated_genes_to_include}); }
        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if (( $self->skip_non_coding && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
           ( $self->skip_silent && $mutation_class =~ m/^Silent$/ )) {
            print "Skipping $mutation_class mutation in gene $gene.\n";
            next;
        }
        my ( $var, $variant_name, );
        if ( $ref eq $var1 ) {
            $var = $var2;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        } elsif ( $ref eq $var2 ) {
            $var = $var1;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        } elsif ( $ref ne $var1 && $ref ne $var2 ) {
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
    my $matrix_fh = new IO::File $mutation_matrix,"w" or die "Failed to create matrix file $mutation_matrix!: $!";
    #print input matrix file header
    my $header = join( "\t", "Sample", @variant_names );
    $matrix_fh->print("$header\n");
    #print mutation relation input matrix
    for my $sample (sort @all_sample_names) {
        $matrix_fh->print( $sample );
        for my $variant (@variant_names) {
            if (exists $variants_hash{$sample}{$variant}) {
                $matrix_fh->print("\t$variants_hash{$sample}{$variant}");
            } else {
                $matrix_fh->print("\t0");
            }
        }
        $matrix_fh->print("\n");
    }

}

1;
