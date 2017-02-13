package TGI::MuSiC2::LongGeneFilter;

#__STANDARD PERL PACKAGES
use strict;
use warnings;

#__LIBRARIES
use Statistics::Distributions;

use Carp;  # Use carp/croak/confess for improved error handling http://perldoc.perl.org/Carp.html
use IO::File;
use Getopt::Long;

# TODO: 
# * get rid of pos_genes, this will be done in a post-processing step
# * standardize input, output format.  Needs to be reentrant
# * Allow user to define whether output will be genes which remain or which are discarded
# * Allow values for y_thresh_l, etc. (6, 0.5) to be modified by user



# Evaluate command line arguments to obtain list of data files.
# if data_files defined, that is a comma-separated list of filenames
# if list_file defined, that file containing a list of filenames
# If both are defined that is an error
# If neither is defined that is an error
sub get_data_filenames {
    my ($this, $data_files, $list_file) = @_;

    die("Please define either --data-files or --list-file\n") if !($data_files xor $list_file);
    my @files;
    if ($list_file) {
        my $DF = IO::File->new( $list_file ) or die "Couldn't open $list_file. $!\n";
        while( my $line = $DF->getline ) {
            chomp $line;
            push(@files, $line);
        }
        $DF->close;
    } else {
        @files = split(/,/, $data_files);
    }
    return (\@files); 
}

#__GET LIST OF POSITIVE GENES (FOR CALLER ASSESSMENT)
sub get_pos_genes {
    my ($this, $pos_genes_file) = @_;
    my $pos_genes = {};
    my $DF = IO::File->new( $pos_genes_file ) or die "Couldn't open $pos_genes_file. $!\n";
    while( my $line = $DF->getline ) {
        chomp $line;
        $pos_genes->{$line} = 1;
    }
    $DF->close;
    return ($pos_genes);
}

sub print_header{
    my $this = shift;
    my $s = <<HEADER;
# CANCER-SPECIFIC ADJUSTMENT OF P-VAL THRESHOLD FOR LARGE GENES\n#
# music2 long-gene-filter \n#
# ************************************************************\n#
# PARAMETERS\n#\n
# default delineation between 'typical' and 'large' genes: $this->{X_THRESH_GENE_SIZE}
#    (note that this boundary may be modified for some cancers)
# highest allowable Y-threshold: $this->{Y_MAX}\n#
# P-value target above which mutational significance status
#      and gene size are no longer statistically related: $this->{PVALUE_THRESHOLD}\n#
# ************************************************************\n
HEADER
    print($s);
}


# Returns: number of genes tested and hash of filtered genes
sub process_file {
    my ($this, $file, $pos_genes) = @_;

#    my $eval_done = 0;  # this can go away
    my $genes_tested = 0;  # this is returned
    my $filtered_long_genes = {}; # this is returned

#__READ DATA IN THIS FILE
    warn "# processing data in file: $file\n";
#     print "# -------------------------------------------------------------\n";
    print "# processing data in file: $file\n";
    print "# -------------------------------------------------------------\n";
    my ($pts, $gene_data, $max_y_large_gene) = $this->_read_data_file_ ($file);

#__Y-THRESHOLD IS CURRENTLY FIXED FOR LEFT-PLANE GENES
    my $iters = 0;

    for (my $y_thresh_l = $this->{Y_THRESH_FIXED}; $y_thresh_l >= 6; $y_thresh_l -= 0.5) {
        warn "#   trying cutoff of Y = $y_thresh_l\n";

#__X-THRESHOLD THAT DELINEATES TYPICAL VS LARGE GENES
        for (my $x_thresh = $this->{X_THRESH_GENE_SIZE}; $x_thresh > $this->{X_INCR};
                $x_thresh -= $this->{X_INCR}) {

#__MOVE Y-THRESHOLD UPWARD FOR LARGE GENES UNTIL SIGNIFICANCE VANISHES
            for (my $y_thresh_r = $y_thresh_l; $y_thresh_r <= $this->{Y_MAX};
                    $y_thresh_r += $this->{Y_INCR}) {

#__GET 2x2 TABLE CELL COUNTS
                my ($top_left, $bot_left, $top_right, $bot_right) = _counts_ (
                        $pts, $x_thresh, $y_thresh_l, $y_thresh_r
                        );

#__PERFORM 2x2 TABLE SIGNIFICANCE TEST
                my ($pval, $grand_total) = _test_ ($top_left, $bot_left, $top_right, $bot_right);
                if ($genes_tested == 0) {
                    $genes_tested = $grand_total;
                }

#__TERMINATE ONCE SIGNIFICANCE VANISHES: CURRENT Y-THRESHOLD IS...
                if ($pval > $this->{PVALUE_THRESHOLD}) {
# Success! print out diagnostics
#__OUTPUT
                    my $p_val_typical_gene = 10**(-$y_thresh_l); #Log10
                    my $p_val_large_gene = 10**(-$y_thresh_r);   #Log10
                    print "   total gene evalluations for this cancer: $grand_total\n";
                    print "d  GENE SIZE: TYPICAL GENES <= $x_thresh <= LARGE GENES";
                    print "   (unchanged from default)" if
                        $x_thresh eq $this->{X_THRESH_GENE_SIZE};
                    print "\n";
                    print "   RECOMMENDED -LN(P-VALUE) CUT-OFFS\n";
                    print "d     typical genes: -LN(P) = $y_thresh_l " .
                        "(P-val = $p_val_typical_gene)";
                    print "   (unchanged from default)" if
                        $y_thresh_l eq $this->{Y_THRESH_FIXED};
                    print "\n";
                    print "d     large genes:   -LN(P) = $y_thresh_r " .
                        "(P-val = $p_val_large_gene)";
                    print "   (unchanged from default)" if
                        $y_thresh_r eq $this->{Y_THRESH_FIXED};
                    print "\n";
                    print "   CALCULATION META-DATA\n";
                    print "      total iterations = $iters\n";
                    print "      final 2x2 table:   $top_left,  $top_right\n";
                    print "                         $bot_left,  $bot_right\n";
                    print "      final P-val $pval > threshold $this->{PVALUE_THRESHOLD} " .
                        "(mutational significance status not statistically " .
                        "related to gene size)\n";

                    my ($old_tps, $old_fps, $old_tns, $old_fns) = (0, 0, 0, 0);
                    my ($new_tps, $new_fps, $new_tns, $new_fns) = (0, 0, 0, 0);
                    my ($original_fp, $current_fp, $true_pos) = ({}, {}, {});
                    foreach my $y (sort _numeric_ keys %{$gene_data}) {
                        foreach my $x (keys %{$gene_data->{$y}}) {
                            foreach my $gene_name (keys %{$gene_data->{$y}->{$x}}) {

#__ORIGINAL (FIXED-Y) FILTERING (NO LONG-GENE POST-HOC FILTER)
                                if ($y <= $y_thresh_l) {
                                    if (exists $pos_genes->{$gene_name}) {
                                        $old_fns++;
                                    } else {
                                        $old_tns++;
                                    }
                                } else {
                                    if (exists $pos_genes->{$gene_name}) {
                                        $old_tps++;
                                    } else {
                                        $old_fps++;
                                    }
                                }

#__NEW FILTERING RESULTS WITH LONG-GENE POST-HOC FILTER
                                if ($x > $x_thresh && $y > $y_thresh_l && $y < $y_thresh_r) {
                                    $filtered_long_genes->{$gene_name}++;
                                    $current_fp->{$gene_name} = "($x, $y)";
                                    if (exists $pos_genes->{$gene_name}) {
                                        $new_fns++;
                                    } else {
                                        $new_tns++;
                                    }
                                } elsif ($y <= $y_thresh_l) {
                                    $original_fp->{$gene_name} = "($x, $y)";
                                    if (exists $pos_genes->{$gene_name}) {
                                        $new_fns++;
                                    } else {
                                        $new_tns++;
                                    }
                                } else {
                                    $true_pos->{$gene_name} = "($x, $y)";
                                    if (exists $pos_genes->{$gene_name}) {
                                        $new_tps++;
                                    } else {
                                        $new_fps++;
                                    }
                                }
                            }
                        }
                    }
                    print "      OLD FILTERING: fp = $old_fps fn = $old_fns\n";
                    print "      NEW FILTERING: fp = $new_fps fn = $new_fns\n";
                    print "      additional genes now filtered as FPs under 'large gene' test (* = NOT part of 366 gene set)\n";
                    foreach my $gene_name (sort keys %{$current_fp}) {
                        print "d      ";
                        if (exists $pos_genes->{$gene_name}) {
                            print "  ";
                        } else {
                            print "* ";
                        }
                        print "$gene_name   $current_fp->{$gene_name}\n";
                    }
                    print "      genes kept as significant (* = NOT part of 366 gene set)\n";
                    foreach my $gene_name (sort keys %{$true_pos}) {
                        print "d      ";
                        if (exists $pos_genes->{$gene_name}) {
                            print "  ";
                        } else {
                            print "* ";
                        }
                        print "$gene_name   $true_pos->{$gene_name}\n";
                    }
                    print "\n";
                    return($genes_tested, $filtered_long_genes);
                } else {
                    $iters++;
                }
            }
        }
    }
    print "\n     NO CONVERGENCE USING CURRENT PARAMETERS\n\n";
    return($genes_tested, $filtered_long_genes);
}

sub new {
    my $class = shift;
    my $this = {};

    $this->{POS_GENE_FILE} = "LongGeneFilter.dat";
    $this->{DATA_FILES} = undef;
    $this->{LIST_FILE} = undef;
    $this->{X_THRESH_GENE_SIZE} = 5000;
    $this->{X_INCR} = 500;
    $this->{Y_THRESH_FIXED} = 8;
    $this->{Y_MAX} = 18;
    $this->{Y_INCR} = 0.1;
    $this->{PVALUE_THRESHOLD} = 0.005;

    bless $this, $class;
    $this->process();

    return $this;
}


sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (  # http://search.cpan.org/~chips/perl5.004_05/lib/Getopt/Long.pm
            'pos-gene-file=s'       => \$this->{POS_GENE_FILE},
            'data-files=s'           => \$this->{DATA_FILES},
            'list-file=s'           => \$this->{LIST_FILE},
            'x-thresh-gene-size=i'  => \$this->{X_THRESH_GENE_SIZE},
            'x-incr=i'              => \$this->{X_INCR},
            'y-thresh-fixed=f'      => \$this->{Y_THRESH_FIXED},
            'y-max=f'               => \$this->{Y_MAX},
            'y-incr=f'              => \$this->{Y_INCR},
            'pvalue-threshold=f'    => \$this->{PVALUE_THRESHOLD},
            'help' => \$help,
            );

    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->usage_text(); }

    my $pos_genes = $this->get_pos_genes($this->{POS_GENE_FILE});

    my @files = @{$this->get_data_filenames($this->{DATA_FILES}, $this->{LIST_FILE})};
    my $num_cancers = scalar @files;

    $this->print_header();

#####################
#  MAIN PROCESSING  #
#####################

#__PROCESS EACH FILE
    my $all_filtered_long_genes = {};
    my $total_gene_evals = 0;

    foreach my $file (@files) {
        my ($gene_evals, $filtered_long_genes) = $this->process_file($file, $pos_genes);
        $total_gene_evals += $gene_evals;
        foreach my $gene (keys %{$filtered_long_genes}) {  
            $all_filtered_long_genes->{$gene} += $filtered_long_genes->{$gene};
        }
    }

#__DIAGNOSTIC GRAND TALLY
    print "\n";
    print "TOTAL GENE EVALUATIONS OVER ALL $num_cancers CANCERS = $total_gene_evals\n\n";

#__UNION OF ALL GENES NETTED BY THE LONG GENE FILTER
    print "\nunionized list of all genes filtered as FPs under 'large gene' test over all cancers (* = part of 366 gene set)\n";
    foreach my $gene_name (sort keys %{$all_filtered_long_genes}) {
        if (exists $pos_genes->{$gene_name}) {
            print "* ";
        } else {
            print "  ";
        }
        print "$gene_name   ($all_filtered_long_genes->{$gene_name} cancers)\n";
    }
}

sub _numeric_ {$a <=> $b}

#   THE LINE "if (defined $y && $y) {" **IS** READING OCCURENCES OF 0.0
#   IN THE INPUT FILES BECAUSE, ALTHOUGH PERL PROCESSES NUMBERS AND
#   TEXT-THAT-RESEMBLES-A-NUMBER DIFFERENTLY (SEE BELOW EXAMPLE), IT *READS*
#   A FILE OF NUMBERS INITIALLY AS TEXT
#
#   bash-3.2$ perl             bash-3.2$ perl
#   $s = 0.0;                  $s = "0.0";   <---quotation marks on this one
#   if ($s) {                  if ($s) {
#     print "yes\n";              print "yes\n";
#   } else {                   } else {
#     print "no\n";               print "no\n";
#   }                          }
#   no                         yes


# 3 column format, tab separated
# gene_name in all upper-case
# x = gene size (integer) 
# y = -log10(MuSiC P-value)
sub _read_data_file_ {
    my ($this, $file) = @_;
    open (F, $file) || die "cant open $file";
    my ($pts, $gene_data, $max_y_large_gene) = ([], {}, 0);
    while (<F>) {

#__PARSE LINE
        next if /^#/;
        chomp;
        my ($gene_name, $x, $y) = split;

#__PROCESS A LEGITIMATE LINE
        if (defined $y && $y) {

#__STORE POINTS AND TRACK MAX Y-VALUE FOR LARGE GENES AND GENE NAMES
            push @{$pts}, [$x, $y];
            if ($x > $this->{X_THRESH_GENE_SIZE}) {
                $max_y_large_gene = $y if $y > $max_y_large_gene;
            }
            $gene_data->{$y}->{$x}->{$gene_name}++;
        }
    }
    close (F);
    return ($pts, $gene_data, $max_y_large_gene);
}


#__GET 2x2 TABLE CELL COUNTS
sub _counts_ {
    my ($pts, $x_thresh, $y_thresh_l, $y_thresh_r, $x_min) = @_;

#__DEFAULTS
    $y_thresh_r = $y_thresh_l unless $y_thresh_r;
    $x_min = 0 unless $x_min;

#__CREATE OBSERVED 2X2 TABLE THAT RESULTS FROM THESE THRESHOLD SETTINGS
    my ($top_left, $bot_left, $top_right, $bot_right) = (0, 0, 0, 0);
    foreach my $pt (@{$pts}) {
        my ($x, $y) = @{$pt};

#__LEFT-HALF-PLANE FOR TYPICAL GENES COMPARES TO FIXED Y-THRESHOLD
        if ($x <= $x_thresh) {
            if ($y <= $y_thresh_l) {
                $bot_left++;
            } else {
                $top_left++;
            }

#__AND RIGHT-HALF-PLANE FOR LARGE GENES TO PROVISIONAL Y-THRESHOLD
        } else {
            if ($y <= $y_thresh_r) {
                $bot_right++;
            } else {
                $top_right++;
            }
        }
    }
    return ($top_left, $bot_left, $top_right, $bot_right);
}


#__PERFORM 2x2 TABLE SIGNIFICANCE TEST
sub _test_ {
    my ($top_left, $bot_left, $top_right, $bot_right) = @_;

#__MARGINAL TOTALS FOR THESE THRESHOLDS
    my $row_top = $top_left + $top_right;
    my $row_bot = $bot_left + $bot_right;
    my $col_left = $top_left + $bot_left;
    my $col_righ = $top_right + $bot_right;

#__GRAND TOTAL
    my $grand_total = $bot_left + $top_left + $bot_right + $top_right;

#__SANITY CHECKING
#  die "tallying problem" unless $row_top + $row_bot == $grand_total;
#  die "tallying problem" unless $col_left + $col_righ == $grand_total;

#__THEORETICAL 2X2 TABLE OF EXPECTED VALUES FROM MARGINAL TOTALS
    my $top_left_expec = $col_left * $row_top / $grand_total;
    my $bot_left_expec = $col_left * $row_bot / $grand_total;
    my $top_right_expec = $col_righ * $row_top / $grand_total;
    my $bot_right_expec = $col_righ * $row_bot / $grand_total;

#__COMPUTE CHI-SQUARE STATISTIC
    my $chisq = ($top_left - $top_left_expec)**2 / $top_left_expec +
        ($top_right - $top_right_expec)**2 / $top_right_expec +
        ($bot_left - $bot_left_expec)**2 / $bot_left_expec +
        ($bot_right - $bot_right_expec)**2 / $bot_right_expec;

#__2X2 TABLE TEST HAS DOF OF 1
    my $dof = 1;

#__COMPUTE GOODNESS-OF-FIT P-VALUE OF OBSERVED VS EXPECTED
    my $pval = Statistics::Distributions::chisqrprob ($dof, $chisq);
    return ($pval, $grand_total);
}

sub usage_text {
    my $this = shift;
    return <<HELP

        MuSiC2 Long Gene Filter Module

        Find conditions for which significance status is no longer related to gene
        size --- this is done by iteratively raising the cut-off for longer genes
        as compared to shorter ones

        USAGE 
        music2 long-gene-filter --data-file=?|--list-file=? [--pos-gene-file=?]
        [x-thresh-gene-size=?] [x-incr=?] [y-thresh-fixed=?] [y-max=?]
        [y-incr=?] [pvalue-threshold=?]

        SEE ALSO

            music2 long-gene-filter --help 

HELP
}

sub help_text {
    my $this = shift;
    my $usage = usage_text();

    return $usage . <<HELP

        ARGUMENTS

        data-files, list-file

            Define list of data files to process in one of two ways:
            * data-files argument is a comma-separated list of data filenames
            * list-file argument is a file containing data filenames, one per line

        Either data-files or list-file is mandatory

            The data files are in TSV format with the following columns:
            * gene_name in all upper-case
            * gene size (integer) 
            * -log10(MuSiC P-value)

        pos-gene-file

            Filename of list of genes, one per line.
            These are the genes which should be kept in the analysis --- it is an error
            to filter them --- all other genes should be filtered.
            Allows for sensitivity & specificity estimation.

        x-thresh-gene-size

            Set the threshold between "typical" genes and "large" genes.  Default 5000

        x-incr

            Down-increment of threshold delimiting typical and large genes.  Default 500

        y-thresh-fixed

            Highest p-value for typical genes is also lower-bound for large genes (the
            actual p-value is exp(- y-thresh-fixed).  Default value 8 (corresponding to
            p-value = 0.00000001)

        y-max

            Set absolute lowest allowable p-value (highest y-threshold) to use. Default 18

        y-incr

            Increment of ln(p-value).  Default 0.1

        pvalue-threshold

            p-value threshold above which we say that mutational significance status
            is no longer related to gene size.  Default 0.005

        ===========
        DESCRIPTION
        ===========

        (1) Goodness of fit 2x2 table test for each cancer in the form of
        of the numbers of genes in each of 2 categories: (long or not long) and
        (significant p-value or not significant p-value)

                                | NOT LONG     LONG |
             -------------------+-------------------+---------------
             SIGNIF P-VALUE     |   a          b    |   a + b
             NOT SIGNIF P-VALUE |   c          d    |   c + d
             -------------------+-------------------+---------------
                                | a + c      b + d  | a + b + c + d

        NOTES:

        * note the protocal for using Statistics::Distributions as written in
        research/in_progress/integrated_analysis/scripts/Statistics/Normality.pm

        * from Statistics::Distributions internal documentation

        \$chisprob=Statistics::Distributions::chisqrprob (3,6.25);
        print "upper probability of the chi-square distribution (3 degrees "
            ."of freedom, chi-squared = 6.25): Q = 1-G = \$chisprob\n";

HELP
}

1;
