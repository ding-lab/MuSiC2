package TGI::MuSiC2::Smg;
##
# Significant Mutated Gene test (SMG)
#  
#
use strict;
use warnings;
#
use Statistics::Descriptive;
use Carp;
use POSIX qw( WIFEXITED );
use IO::File;
use File::Temp qw/ tempfile /;
use Getopt::Long;

sub new {
    my $class = shift;
    my $this = {};

    $this->{_GENE_MR_FILE} = undef;
    $this->{_OUTPUT_FILE}  = undef;
    $this->{_MAX_FDR} = 0.2;
    $this->{_SKIP_LOW_MR_GENES} = 1;
    $this->{NO_SKIP_LOW_MR_GENES} = 0;
    $this->{_DOWNSAMPLE_LARGE_GENES} = 0;
    $this->{NO_DOWNSAMPLE_LARGE_GENES} = 1;
    $this->{_SKIP_NON_EXPRESSED_GENES} = 1;
    $this->{NO_SKIP_NON_EXPRESSED_GENES} = 0;
    $this->{_SKIP_PSEUDOGENES} = 1;
    $this->{NO_SKIP_PSEUDOGENES} = 0;

    $this->{_BMR_MODIFIER_FILE} = undef;
    $this->{_BLACK_GENE_LIST_FILE} = undef;
    $this->{_PROCESSORS} = 1;
    $this->{_QQ_PLOT_FILE} = "smg_test_qq_plot.pdf";

    bless $this, $class;
    $this->process();

    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'gene-mr-file=s'             => \$this->{_GENE_MR_FILE},
        'output-file=s'              => \$this->{_OUTPUT_FILE},
        'max-fdr=f'                  => \$this->{_MAX_FDR},
        'skip-low-mr-genes'          => \$this->{_SKIP_LOW_MR_GENES},
        'noskip-low-mr-genes'        => \$this->{NO_SKIP_LOW_MR_GENES},
        'skip-non-expressed-genes'   => \$this->{_SKIP_NON_EXPRESSED_GENES},
        'noskip-non-expressed-genes' => \$this->{NO_SKIP_NON_EXPRESSED_GENES},
        'skip-pseudogenes'           => \$this->{_SKIP_},
        'noskip-pseudogenes'         => \$this->{NO_SKIP_NON_EXPRESSED_GENES},

        'downsample-large-genes'     => \$this->{_DOWNSAMPLE_LARGE_GENES},
        'nodownsample-large-genes'   => \$this->{NO_DOWNSAMPLE_LARGE_GENES},
        'bmr-modifier-file=s'        => \$this->{_BMR_MODIFIER_FILE},
        'black-gene-list-file=s'     => \$this->{_BLACK_GENE_LIST_FILE},
        'processors=i'               => \$this->{_PROCESSORS},
        'qq-plot-file=s'             => \$this->{_QQ_PLOT_FILE},

        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    #### processing ####
    # Check on all the input data
    my $output_file_detailed = $this->{_OUTPUT_FILE} . "_detailed";
    # Check on all the input data before starting work
    print STDERR "Gene mutation rate file not found or is empty: $this->{_GENE_MR_FILE}\n" unless( -s $this->{_GENE_MR_FILE} );
    return undef unless( -s $this->{_GENE_MR_FILE} && ( !defined $this->{_BMR_MODIFIER_FILE} || -s $this->{_BMR_MODIFIER_FILE} ));
    # Check boolean paras
    if ( $this->{NO_SKIP_LOW_MR_GENES} ) { $this->{_SKIP_LOW_MR_GENES} = 0; };
    if ( $this->{NO_SKIP_NON_EXPRESSED_GENES} ) { $this->{_SKIP_NON_EXPRESSED_GENES} = 0; }
    # Collect per-gene mutation rates for 
    # reporting in results later
    my ( %gene_muts, %gene_bps, %mut_categ_hash );
    my $inMrFh = IO::File->new( $this->{_GENE_MR_FILE} ) or die "Couldn't open $this->{_GENE_MR_FILE}. $!\n";
    while ( my $line = $inMrFh->getline ) {
        next if ( $line =~ m/^#/ );
        my ( $gene, $type, $covd_bps, $mut_cnt, undef ) = split( /\t/, $line );
        #
        # Warn user about cases where there could be fewer 
        # covered bps than mutations detected
        #
        if ( $type ne "Overall" and $covd_bps < $mut_cnt ) { warn "#More $type seen in $gene than there are bps with sufficient coverage!\n"; }
        SWITCH:{
            $type =~ m/^Overall/ && do { 
                $gene_muts{$gene}{Overall} = $mut_cnt; 
                $gene_bps{$gene} = $covd_bps;       
                last SWITCH; };
            $type =~ m/^Truncations/ && do { 
                $gene_muts{$gene}{Truncations} += $mut_cnt;
                $mut_categ_hash{Truncations} = 1;
                last SWITCH; };
            $type =~ m/^Indels/ && do {
                $gene_muts{$gene}{Indels} += $mut_cnt;
                $mut_categ_hash{Indels} = 1;
                last SWITCH; }; 
            $type =~ m/^(AT_|CG_|CpG_)(Transitions|Transversions)/ && do {
                $gene_muts{$gene}{SNVs} += $mut_cnt;
                $mut_categ_hash{SNVs} = 1;
                last SWITCH; };
            # default 
            die "Unrecognized mutation category in gene-mr-file. $!\n";
        }
    }
    $inMrFh->close;
    my @mut_categs = sort keys %mut_categ_hash;
    #
    # If requested, downsample large genes 
    # and create another gene_mr_file
    #
    if ( $this->{_DOWNSAMPLE_LARGE_GENES} ) {
        #
        # Collect covd_bps across genes per mutation category
        #
        my ( %gene_categ_bps );
        my $inMrFh = IO::File->new( $this->{_GENE_MR_FILE} ) or die "Couldn't open $this->{_GENE_MR_FILE}. $!\n";
        while ( my $line = $inMrFh->getline ) {
            next if ( $line =~ m/^#/ );
            my ( $gene, $type, $covd_bps, undef, undef ) = split( /\t/, $line );
            push( @{$gene_categ_bps{$type}}, $covd_bps );
        }
        $inMrFh->close;
        #
        # For each mutation category, find the Q1 and Q3 quartiles, and the median
        #
        my ( %q1, %q2, %q3 );
        foreach my $type ( keys %gene_categ_bps ) {
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data( @{$gene_categ_bps{$type}} );
            $q1{$type} = $stat->quantile( 1 );
            $q2{$type} = $stat->quantile( 2 );
            $q3{$type} = $stat->quantile( 3 );
        }
        #
        # For each mutation category, calculate the outlier cutoff
        my %cutoff;
        print "#Selecting large genes to downsample based on these cutoffs...\n";
        print "#Mutation_Type\tMedian\tIQR\tCutoff\n";
        foreach my $type ( sort keys %gene_categ_bps ) {
            unless( $type eq "Overall" ) {
                my $IQR = $q3{$type} - $q1{$type}; # Interquartile range
                $cutoff{$type} = sprintf( "%.0f", $q3{$type} + ( $IQR * 4.5 ));
                print "$type\t" . $q2{$type} . "\t$IQR\t" . $cutoff{$type} . "\n";
            }
        }
        # Create a new gene_mr_file with large genes downsampled 
        # in #covd_bps and #muts
        print "\n#Changed lines in $this->{_GENE_MR_FILE} as follows...\n";
        my $new_gene_mr_file = $this->{_GENE_MR_FILE} . "_downsampled_large_genes";
        my $outMrFh = IO::File->new( $new_gene_mr_file, ">" ) or die "Couldn't create $new_gene_mr_file. $!\n";
        $inMrFh = IO::File->new( $this->{_GENE_MR_FILE} ) or die "Couldn't open $this->{_GENE_MR_FILE}. $!\n";
        while ( my $line = $inMrFh->getline ) {
            if ( $line =~ m/^#/ ) {
                $outMrFh->print( $line );
                next;
            }
            my ( $gene, $type, $covd_bps, $mut_cnt, $bmr ) = split( /\t/, $line );
            #
            # If the #covd bps in this line is an outlier for its category, then downsample it
            # The category "Overall" is only for reporting purposes, so we can leave it unchanged
            #
            if ( $type ne "Overall" and $covd_bps > $cutoff{$type} ) {
                # If this line has zero mutations, then simply set $covd_bps to the cutoff
                if ( $mut_cnt == 0 ) {
                    $covd_bps = $cutoff{$type};
                }
                # Otherwise, iteratively reduce $mut_cnt, no lower than 1 mutation, bringing
                # $covd_bps as close to the cutoff as possible, but keeping the MR unchanged
                #
                else {
                    for ( my $i = $mut_cnt - 1; $i >= 0; --$i ) {
                        if ( ($covd_bps * $i / $mut_cnt) <= $cutoff{$type} || $i == 0 ) {
                            ++$i;
                            $covd_bps = sprintf( "%.0f", $covd_bps * $i / $mut_cnt );
                            $mut_cnt = $i;
                            last;
                        }
                    }
                }
                # Print out the changed #bps and #muts, for users' benefit
                if ( $line ne "$gene\t$type\t$covd_bps\t$mut_cnt\t$bmr" ) {
                    print "< $line";
                    print "> $gene\t$type\t$covd_bps\t$mut_cnt\t$bmr";
                }
            }

            $outMrFh->print( "$gene\t$type\t$covd_bps\t$mut_cnt\t$bmr" );
        }
        $outMrFh->close;
        $inMrFh->close;

        $this->{_GENE_MR_FILE} = $new_gene_mr_file;
    }
    #
    # If BMR modifiers were provided, then create another 
    # gene_mr_file with modified BMRs
    #
    if ( defined $this->{_BMR_MODIFIER_FILE} ) {
        my $inBmrModFh = IO::File->new( $this->{_BMR_MODIFIER_FILE} ) or die "Couldn't open $this->{_BMR_MODIFIER_FILE}. $!\n";
        my %bmr_modifier = ();
        while ( my $line = $inBmrModFh->getline ) {
            next if( $line =~ m/^#/ );
            chomp( $line );
            my ( $gene, $modifier ) = split( /\t/, $line );
            ( $modifier > 0 ) or die "$modifier is an invalid bmr-modifier. Please fix values in $this->{_BMR_MODIFIER_FILE}.\n";
            $bmr_modifier{$gene} = $modifier;
        }
        $inBmrModFh->close;
        #
        ## using File::Temp module instead of Genome::Sys
        # my ( $tmp_gene_mr_file ) = Genome::Sys->create_temp_file_path();
        #
        my ( undef, $tmp_gene_mr_file ) = tempfile();
        my $tmpFh = IO::File->new( $tmp_gene_mr_file, ">" ) or die "Temporary file could not be created. $!";
        my $inMrFh = IO::File->new( $this->{_GENE_MR_FILE} ) or die "Couldn't open $this->{_GENE_MR_FILE}. $!\n";
        while ( my $line = $inMrFh->getline ) {
            if ( $line =~ m/^#/ ) {
                $tmpFh->print( $line );
                next;
            }
            chomp( $line );
            my ( $gene, $type, $covd_bps, $mut_cnt, $bmr ) = split( /\t/, $line );
            $bmr = $bmr * $bmr_modifier{$gene} if( defined $bmr_modifier{$gene} );
            $tmpFh->print( "$gene\t$type\t$covd_bps\t$mut_cnt\t$bmr\n" );
        }
        $tmpFh->close;
        $inMrFh->close;
        $this->{_GENE_MR_FILE} = $tmp_gene_mr_file;
    }
    # rearrange gene mrs file for load balance 
    # of parallelation or SMG test
    my ( %gene_hash, %gene_mutations_hash, @temp_arr, $temp_index );
    my $new_gene_mrs_fh = IO::File->new( $this->{_GENE_MR_FILE} ) or die "Couldn't open $this->{_GENE_MR_FILE}. $!\n";
    map{ @_ = split /\t/; $gene_hash{ $_[0] } .= $_; if ( $_[1] eq "Overall" ) { $gene_mutations_hash{$_[0]} += $_[3]; } } $new_gene_mrs_fh->getlines;
    $new_gene_mrs_fh->close();
    # sorting by mutation number
    my @sorted_genes = sort{ $gene_mutations_hash{$b} <=> $gene_mutations_hash{$a} } keys %gene_mutations_hash;
    my ( $con, @tt, @ss, $k, ); $con = ""; $k = 1;
    map{ $con .= $_.","; if ( $_ % $this->{_PROCESSORS} == 0 ){ push(@tt, $con); $con = ""; } } (1..@sorted_genes);
    # pick up last one if need
    unless( @sorted_genes % $this->{_PROCESSORS} == 0 ){ push(@tt, $con); };
    map{ if ( $k % 2 == 0 ){ push( @ss, reverse(split /,/, $_) ); } else { push( @ss, (split /,/, $_) ); }; $k++; } @tt;
    my ( undef, $resorted_gene_mrs ) = tempfile();
    my $resorted_gene_mrs_fh = IO::File->new( $resorted_gene_mrs, ">" ) or die "Temporary file could not be created. $!";
    map{ $resorted_gene_mrs_fh->print( $gene_hash{ $sorted_genes[$_-1] } ); } @ss;
    $resorted_gene_mrs_fh->close();
    $this->{_GENE_MR_FILE} = $resorted_gene_mrs;
    # Create a temporary intermediate file to hold the p-values
    my ( undef, $pval_file ) = tempfile();
    my $pvalFh = IO::File->new( $pval_file, ">" ) or die "Temporary file could not be created. $!";
    # Call R for Fisher combined test, Likelihood ratio test, and convolution test on each gene
    my $smg_cmd = "R --slave --args < " . __FILE__ . ".R $this->{_GENE_MR_FILE} $pval_file smg_test $this->{_PROCESSORS} $this->{_SKIP_LOW_MR_GENES}";
    WIFEXITED( system $smg_cmd ) or croak "Couldn't run: $smg_cmd ($?)";
    # Call R for calculating FDR on the p-values calculated in the SMG test
    my $fdr_cmd = "R --slave --args < " . __FILE__ . ".R $pval_file $output_file_detailed calc_fdr $this->{_PROCESSORS} $this->{_SKIP_LOW_MR_GENES}";
    WIFEXITED( system $fdr_cmd ) or croak "Couldn't run: $fdr_cmd ($?)";
    # Parse the R output to identify the SMGs (significant 
    # by at least 2 of 3 tests)
    my $smgFh = IO::File->new( $output_file_detailed ) or die "Couldn't open $output_file_detailed. $!\n";
    my ( @newLines, @smgLines );
    my $header = "#Gene\t" . join( "\t", @mut_categs );
    $header .= "\tTot Muts\tCovd Bps\tMuts pMbp\tP-value FCPT\tP-value LRT\tP-value CT\tFDR FCPT\tFDR LRT\tFDR CT\n";
    while ( my $line = $smgFh->getline ) {
        chomp( $line );
        if ( $line =~ m/^Gene\tp.fisher\tp.lr\tp.convol\tfdr.fisher\tfdr.lr\tfdr.convol$/ ) {
            push( @newLines, $header );
            push( @smgLines, $header );
        } else {
            my ( $gene, @pq_vals ) = split( /\t/, $line );
            my ( $p_fcpt, $p_lrt, $p_ct, $q_fcpt, $q_lrt, $q_ct ) = @pq_vals;
            my @mut_cnts;
            foreach ( @mut_categs ) {
                # If a mutation count is a fraction, round down the digits after the decimal point
                push( @mut_cnts, (( $gene_muts{$gene}{$_} =~ m/\./ ) ? sprintf( "%.2f", $gene_muts{$gene}{$_} ) : $gene_muts{$gene}{$_} ));
            }
            my $mut_per_mbp = ( $gene_bps{$gene} ? sprintf( "%.2f", ( $gene_muts{$gene}{Overall} / $gene_bps{$gene} * 1000000 )) : 0 );
            push( @newLines, join( "\t", $gene, @mut_cnts, $gene_muts{$gene}{Overall}, $gene_bps{$gene}, $mut_per_mbp, @pq_vals ) . "\n" );
            # If the FDR of at least two of these tests is less than the maximum 
            # allowed, we consider it an SMG
            if ( ( $q_fcpt <= $this->{_MAX_FDR} && $q_lrt <= $this->{_MAX_FDR} ) || ( $q_fcpt <= $this->{_MAX_FDR} && $q_ct <= $this->{_MAX_FDR} ) ||
               ( $q_lrt <= $this->{_MAX_FDR} && $q_ct <= $this->{_MAX_FDR} )) {
                push( @smgLines, join( "\t", $gene, @mut_cnts, $gene_muts{$gene}{Overall}, $gene_bps{$gene}, $mut_per_mbp, @pq_vals ) . "\n" );
            }
        }
    }
    $smgFh->close;
    #### do expressed filtering 
    # load expression black gene list
    #my $black_gene_list_ref = undef;
    #$black_gene_list_ref = $this->{_SKIP_NON_EXPRESSED_GENES} ? $this->black_genes( $this->{_BLACK_GENE_LIST_FILE} ) : undef;
    my $black_gene_list_ref = $this->black_genes( $this->{_BLACK_GENE_LIST_FILE} );
    # Add per-gene SNV and Indel counts to the detailed R output, and 
    # make the header friendlier
    my $outDetFh = IO::File->new( $output_file_detailed, ">" ) or die "Couldn't open $output_file_detailed. $!\n";
    foreach ( @newLines ) {
        if ( /^\#Gene/ ) { 
            s/\n/\tExpression\n/; 
            $outDetFh->print( $_ );
        } else {
            my ($gene) = /^(.*?)\t/;
            ( defined $black_gene_list_ref->{$gene} ) ? s/\n/\tnon-expressed\n/ : s/\n/\texpressed\n/;
            $outDetFh->print( $_ ) unless( defined $black_gene_list_ref->{$gene} and $this->{_SKIP_NON_EXPRESSED_GENES} );
        }
    }
    $outDetFh->close;
    #
    # Do the same for only the genes that we consider SMGs
    my $outFh = IO::File->new( $this->{_OUTPUT_FILE}, ">" ) or die "Couldn't open $this->{_OUTPUT_FILE} . $!\n";
    #$outFh->print( @smgLines );
    foreach ( @smgLines ) {
        if ( /^\#Gene/ ) { 
            s/\n/\tExpression\n/; 
            $outFh->print( $_ );
        } else {
            my ($gene) = /^(.*?)\t/;
            ( defined $black_gene_list_ref->{$gene} ) ? s/\n/\tnon-expressed\n/ : s/\n/\texpressed\n/;
            $outFh->print( $_ ) unless( defined $black_gene_list_ref->{$gene} and $this->{_SKIP_NON_EXPRESSED_GENES} );
        }
    }
    $outFh->close;
    # Call R for qqplot 
    my $qq_cmd = "R --slave --args < " . __FILE__ . ".qqplot.R $output_file_detailed $this->{_QQ_PLOT_FILE}";
    WIFEXITED( system $qq_cmd ) or croak "Couldn't run: $qq_cmd ($?)";

    return 1;
}

## black gene list
sub black_genes {
    my ( $this, $black_genes_matrix_file, ) = @_;
    return undef unless ( defined $black_genes_matrix_file );
    my $fh_black_genes = IO::File->new( $black_genes_matrix_file ) or warn "Couldn't open black gene expression matrix file: $black_genes_matrix_file,. $!\n";
    my %genes_hash = map{ /^(.*?)\t/; ($1, 1) } <$fh_black_genes>;
    return \%genes_hash;
}

## usage
sub help_text {
    my $this = shift;
        return <<HELP

USAGE
 music2 smg --gene-mr-file=? --output-file=? [--max-fdr=?] [--skip-low-mr-genes]
    [--downsample-large-genes] [--bmr-modifier-file=?] [--processors=?]

SYNOPSIS
 music2 smg \
      --gene-mr-file output_dir/gene_mrs \
      --output-file output_dir/smgs

(A "gene-mr-file" can be generated using the tool "music2 bmr calc-bmr".)

REQUIRED INPUTS
  gene-mr-file
    File with per-gene mutation rates (Created using "music2 bmr calc-bmr") 
  output-file
    Output file that will list significantly mutated genes and their p-values 

OPTIONAL INPUTS
  max-fdr
    The maximum allowed false discovery rate for a gene to be considered an SMG 
    Default value '0.2' if not specified
  skip-low-mr-genes
    Skip testing genes with MRs lower than the background MR 
    Default value 'true' if not specified
  noskip-low-mr-genes
    Make skip-low-mr-genes 'false'
  skip-non-expressed-genes
    Skip non expressed genes based on TCGA data
    Default value 'true' if not specified
  nonskip-non-expressed-genes
    Make skip-non-expressed-genes 'false'
  downsample-large-genes
    Downscale #bps in large genes, and the #muts proportionally 
    Default value 'false' if not specified
  nodownsample-large-genes
    Make downsample-large-genes 'false'
  bmr-modifier-file
    Tab delimited multipliers per gene that modify BMR before testing [gene_name bmr_modifier] 
  black-gene-list-file 
    Expression matrix file for non-expressed genes based on TCGA expression data
  processors
    Number of cores to use (requires 'foreach' and 'doMC' R packages) 
    Default value '1' if not specified
  qq-plot-file
    output qq plot file for SMG test result.
    Default value 'smg_test_qq_plot.pdf'

DESCRIPTION
    This script runs R-based statistical tools to identify Significantly Mutated Genes (SMGs), when
    given per-gene mutation rates categorized by mutation type, and the overall background mutation
    rates (BMRs) for each of those categories (gene_mr_file, created using "music bmr calc-bmr").

    P-values and false discovery rates (FDRs) for each gene in gene_mr_file is calculated using
    three tests: Fisher's Combined P-value test (FCPT), Likelihood Ratio test (LRT), and the
    Convolution test (CT). For a gene, if its FDR for at least 2 of these tests is <= max_fdr, it
    will be output as an SMG. Another output file with prefix "_detailed" will have p-values and
    FDRs for all genes.


ARGUMENTS

    --skip-low-mr-genes

      Genes with consistently lower MRs than the BMRs across mutation categories, may show up in
      the results as an SMG (by CT or LRT). If such genes are not of interest, they may be assigned
      a p-value of 1. This should also speed things up. Genes with higher Indel or Truncation rates
      than the background will not be skipped even if the gene's overall MR is lower than the BMR.
      If bmr-modifiers are applied, this step uses the modified BMRs instead.

    --skip-non-expressed-genes

      We provided users one black gene list in which genes are not expressed based on TCGA expression
      data. Those genes can be exclused using this filter. Users have to provide black-gene-list-file
      when using this option. The black gene list file for each corresponding Ensemble version can 
      be downloaded from MuSiC2 website. 

    --downsample-large-genes

      Large genes have a higher power for detecting mutational significance than smaller genes. So
      a workaround is to downsample their covd bps and mutation count, keeing MR unchanged. A gene
      with more covd bps than Q3+(Q3-Q1)*4.5 is considered an outlier, where Q1 & Q3 are the lower
      and upper quartiles respectively, of #covd_bps in the gene-mr-file. For each outlier, the
      #covd bps in each mutation category is brought closer to Q3+(Q3-Q1)*4.5, and the #muts
      reduced proportionally, to keep the gene's MR mostly unchanged.

    --bmr-modifier-file

      The user can provide a BMR modifier for each gene in the ROI file, which is a multiplier for
      the categorized background mutation rates, before testing them against the gene's categorized
      mutation rates. Such a file can be used to correct for regional or systematic bias in
      mutation rates across the genome that may be correlated to CpG deamination or DNA repair
      processes like transcription-coupled repair or mismatch repair. Mutation rates have also been
      associated with DNA replication timing, where higher mutation rates are seen in late
      replicating regions. Note that the same per-gene multiplier is used on each mutation category
      of BMR. Any genes from the ROI file that are not in the BMR modifier file will be tested
      against unmodified overall BMRs per mutation category. BMR modifiers of <=0 are not
      permitted, because that's just silly.

HELP

}

1;

