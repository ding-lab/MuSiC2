package TGI::MuSiC2::LongGeneFilter;

# standard perl packages
use strict;
use warnings;

# libraries
use Statistics::Distributions;

use Carp;  # Use carp/croak/confess for improved error handling http://perldoc.perl.org/Carp.html
use IO::File;
use Getopt::Long;

# TODO: 
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
    my ($this, $file) = @_;

# read data in this file
    print "# processing data in file: $file\n";
    print "# -------------------------------------------------------------\n";

    my ($pts, $gene_data, $max_y_large_gene) = $this->read_data_file($file);
    my $genes_tested = @{$pts};

    my $iters = 0;

# y-threshold is variable for left-plane genes
    #for (my $y_thresh_l = $this->{Y_THRESH_FIXED}; $y_thresh_l >= 6; $y_thresh_l -= 0.5) {
    # defaults: Y_THRESH_L_MIN = 6, Y_THRESH_L_DEC = 0.5
    for (my $y_thresh_l = $this->{Y_THRESH_FIXED}; $y_thresh_l >= $this->{Y_THRESH_L_MIN}; $y_thresh_l -= $this->{Y_THRESH_L_DEC}) {
        print"#   trying cutoff of Y = $y_thresh_l\n";

# x-threshold that delineates typical vs large genes
        for (my $x_thresh = $this->{X_THRESH_GENE_SIZE}; $x_thresh > $this->{X_INCR};
                $x_thresh -= $this->{X_INCR}) {

# move y-threshold upward for large genes until significance vanishes
            for (my $y_thresh_r = $y_thresh_l; $y_thresh_r <= $this->{Y_MAX};
                    $y_thresh_r += $this->{Y_INCR}) {

                if ($this->evaluate_long_filter($pts, $x_thresh, $y_thresh_l, $y_thresh_r)) {
                    my ($filtered_long_genes) = $this->get_filtered_genes($gene_data, $x_thresh, $y_thresh_l, $y_thresh_r);
                    print "      total iterations = $iters\n";
                    return ($genes_tested, $filtered_long_genes);
                } else {
                    $iters++;
                }
            }
        }
    }

    # TODO: how to deal with this error condition correctly?
    print "\n     NO CONVERGENCE USING CURRENT PARAMETERS\n\n";
    return($genes_tested, {});
}

sub evaluate_long_filter {
    my ($this, $pts, $x_thresh, $y_thresh_l, $y_thresh_r) = @_;

# get 2X2 table cell counts
    my ($top_left, $bot_left, $top_right, $bot_right) = _counts_( $pts, $x_thresh, $y_thresh_l, $y_thresh_r);

# perform 2X2 table significance test
    my ($pval) = _test_ ($top_left, $bot_left, $top_right, $bot_right);

# terminate once significance vanishes: current y-threshold is...
    if ($pval > $this->{PVALUE_THRESHOLD}) {

        # Success! print out diagnostics
        $this->print_diagnostics($x_thresh, $y_thresh_l, $y_thresh_r, $top_left, $bot_left, $top_right, $bot_right, $pval);
        return 1;
    } else {
        return 0;
    }
}

sub print_diagnostics {
    my ($this, $x_thresh, $y_thresh_l, $y_thresh_r, $top_left, $bot_left, $top_right, $bot_right, $pval) = @_;

    my $grand_total = $top_left + $bot_left + $top_right + $bot_right;
    my $p_val_typical_gene = 10**(-$y_thresh_l); #Log10
    my $p_val_large_gene = 10**(-$y_thresh_r);   #Log10

    print "   total gene evaluations for this cancer: $grand_total\n";
    print "d  GENE SIZE: TYPICAL GENES <= $x_thresh <= LARGE GENES";
    print "   (unchanged from default)" if $x_thresh eq $this->{X_THRESH_GENE_SIZE};
    print "\n";
    print "   RECOMMENDED -LN(P-VALUE) CUT-OFFS\n";
    print "d     typical genes: -LN(P) = $y_thresh_l (P-val = $p_val_typical_gene)";
    print "   (unchanged from default)" if $y_thresh_l eq $this->{Y_THRESH_FIXED};
    print "\n";
    print "d     large genes:   -LN(P) = $y_thresh_r (P-val = $p_val_large_gene)";
    print "   (unchanged from default)" if $y_thresh_r eq $this->{Y_THRESH_FIXED};
    print "\n";
    print "   CALCULATION META-DATA\n";
    print "      final 2x2 table:   $top_left,  $top_right\n";
    print "                         $bot_left,  $bot_right\n";
    print "      final P-val $pval > threshold $this->{PVALUE_THRESHOLD} " .
        "(mutational significance status not statistically " .
        "related to gene size)\n";
}

sub get_filtered_genes {    
    my ($this, $gene_data, $x_thresh, $y_thresh_l, $y_thresh_r) = @_;

    my $filtered_long_genes = {}; # this is returned

    my ($new_tps, $new_tns) = (0, 0);
    my ($original_fp, $current_fp, $true_pos) = ({}, {}, {});
    foreach my $y (sort _numeric_ keys %{$gene_data}) {
        foreach my $x (keys %{$gene_data->{$y}}) {
            foreach my $gene_name (keys %{$gene_data->{$y}->{$x}}) {

# new filtering results with long-gene post-hoc filter
                if ($x > $x_thresh && $y > $y_thresh_l && $y < $y_thresh_r) {
                    $filtered_long_genes->{$gene_name}++;
                    $current_fp->{$gene_name} = "($x, $y)";
# with pos_genes gone, all are new_fns/new_tps
                    $new_tns++;
                } elsif ($y <= $y_thresh_l) {
                    $original_fp->{$gene_name} = "($x, $y)";
                    $new_tns++;
                } else {
                    $true_pos->{$gene_name} = "($x, $y)";
                    $new_tps++;
                }
            }
        }
    }
    print "      p = $new_tps n = $new_tns\n";
    print "      additional genes now filtered as FPs under 'large gene' test\n";
    foreach my $gene_name (sort keys %{$current_fp}) {
        print "$gene_name   $current_fp->{$gene_name}\n";
    }
    print "      genes kept as significant\n";
    foreach my $gene_name (sort keys %{$true_pos}) {
        print "$gene_name   $true_pos->{$gene_name}\n";
    }
    print "\n";

    return($filtered_long_genes);
}


sub new {
    my $class = shift;
    my $this = {};

    $this->{DATA_FILES} = undef;
    $this->{LIST_FILE} = undef;
    $this->{X_THRESH_GENE_SIZE} = 5000;
    $this->{X_INCR} = 500;
    $this->{Y_THRESH_FIXED} = 8;
    $this->{Y_THRESH_L_MIN} = 6;
    $this->{Y_THRESH_L_DEC} = 0.5;
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
            'data-files=s'          => \$this->{DATA_FILES},
            'list-file=s'           => \$this->{LIST_FILE},
            'x-thresh-gene-size=i'  => \$this->{X_THRESH_GENE_SIZE},
            'x-incr=i'              => \$this->{X_INCR},
            'y-thresh-fixed=f'      => \$this->{Y_THRESH_FIXED},
            'y-thresh-l-min=f'      => \$this->{Y_THRESH_L_MIN},
            'y-thresh-l-dec=f'      => \$this->{Y_THRESH_L_DEC},
            'y-max=f'               => \$this->{Y_MAX},
            'y-incr=f'              => \$this->{Y_INCR},
            'pvalue-threshold=f'    => \$this->{PVALUE_THRESHOLD},
            'help' => \$help,
            );
    # defaults: Y_THRESH_L_MIN = 6, Y_THRESH_L_DEC = 0.5

    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->usage_text(); }

    my @files = @{$this->get_data_filenames($this->{DATA_FILES}, $this->{LIST_FILE})};
    my $num_cancers = scalar @files;

    $this->print_header();

#####################
#  MAIN PROCESSING  #
#####################

# PROCESS EACH FILE
    my $all_filtered_long_genes = {};
    my $total_gene_evals = 0;

    foreach my $file (@files) {
        my ($gene_evals, $filtered_long_genes) = $this->process_file($file);

        $total_gene_evals += $gene_evals;
        foreach my $gene (keys %{$filtered_long_genes}) {  
            $all_filtered_long_genes->{$gene} += $filtered_long_genes->{$gene};
        }
    }

# diagnostic grand tally
    print "\n";
    print "TOTAL GENE EVALUATIONS OVER ALL $num_cancers CANCERS = $total_gene_evals\n\n";

# union of all genes netted by the long gene filter
    print "\nunionized list of all genes filtered as FPs under 'large gene' test over all cancers\n";
    foreach my $gene_name (sort keys %{$all_filtered_long_genes}) {
        print "$gene_name   ($all_filtered_long_genes->{$gene_name} cancers)\n";
    }
}

sub _numeric_ {$a <=> $b}

#   the line "IF (DEFINED $Y && $Y) {" **is** reading occurences of 0.0
#   in the input files because, although perl processes numbers and
#   text-that-resembles-a-number differently (see below example), it *reads*
#   a file of numbers initially as text
#
#   bash-3.2$ perl             bash-3.2$ perl
#   $s = 0.0;                  $s = "0.0";   <---quotation marks on this one
#   if ($s) {                  if ($s) {
#     print "yes\n";              print "yes\n";
#   } else {                   } else {
#     print "no\n";               print "no\n";
#   }                          }
#   no                         yes


# Input data is 3 column format, tab separated:
# * gene_name in all upper-case
# * x = gene size (integer) 
# * y = -log10(MuSiC P-value)
# Returns ($pts, $gene_data, $max_y_large_gene):
# * pts is array of (X,Y) coordinates, with X gene size and Y=-log10(P-value)
# * gene_data is a hash of counts of X,Y,Gene_name
# * max_y_large_gene is the largest Y (lnP) value of any gene larger than X-threshold

sub read_data_file {
    my ($this, $file) = @_;
    open (F, $file) || die "cant open $file";
    my ($pts, $gene_data, $max_y_large_gene) = ([], {}, 0);
    while (<F>) {

# parse line
        next if /^#/;
        chomp;
        my ($gene_name, $x, $y) = split;

# process a legitimate line
        if (defined $y && $y) {

# store points and track max y-value for large genes and gene names
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


# get 2X2 table cell counts
sub _counts_ {
    my ($pts, $x_thresh, $y_thresh_l, $y_thresh_r, $x_min) = @_;

# defaults
    $y_thresh_r = $y_thresh_l unless $y_thresh_r;
    $x_min = 0 unless $x_min;

# create observed 2x2 table that results from these threshold settings
    my ($top_left, $bot_left, $top_right, $bot_right) = (0, 0, 0, 0);
    foreach my $pt (@{$pts}) {
        my ($x, $y) = @{$pt};

# left-half-plane for typical genes compares to fixed y-threshold
        if ($x <= $x_thresh) {
            if ($y <= $y_thresh_l) {
                $bot_left++;
            } else {
                $top_left++;
            }

# and right-half-plane for large genes to provisional y-threshold
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


# perform 2X2 table significance test
sub _test_ {
    my ($top_left, $bot_left, $top_right, $bot_right) = @_;

# marginal totals for these thresholds
    my $row_top = $top_left + $top_right;
    my $row_bot = $bot_left + $bot_right;
    my $col_left = $top_left + $bot_left;
    my $col_righ = $top_right + $bot_right;

# grand total is simply the sum of the arguments.  
    my $grand_total = $bot_left + $top_left + $bot_right + $top_right;

# sanity checking
#  die "tallying problem" unless $row_top + $row_bot == $grand_total;
#  die "tallying problem" unless $col_left + $col_righ == $grand_total;

# theoretical 2x2 table of expected values from marginal totals
    my $top_left_expec = $col_left * $row_top / $grand_total;
    my $bot_left_expec = $col_left * $row_bot / $grand_total;
    my $top_right_expec = $col_righ * $row_top / $grand_total;
    my $bot_right_expec = $col_righ * $row_bot / $grand_total;

# compute chi-square statistic
    my $chisq = ($top_left - $top_left_expec)**2 / $top_left_expec +
        ($top_right - $top_right_expec)**2 / $top_right_expec +
        ($bot_left - $bot_left_expec)**2 / $bot_left_expec +
        ($bot_right - $bot_right_expec)**2 / $bot_right_expec;

# 2x2 table test has dof of 1
    my $dof = 1;

# compute goodness-of-fit p-value of observed vs expected
    my $pval = Statistics::Distributions::chisqrprob ($dof, $chisq);
    return ($pval);
}

sub usage_text {
    my $this = shift;
    return <<HELP

        MuSiC2 Long Gene Filter Module

        Find conditions for which significance status is no longer related to gene
        size --- this is done by iteratively raising the cut-off for longer genes
        as compared to shorter ones

        USAGE 
        music2 long-gene-filter --data-file=?|--list-file=? 
        [x-thresh-gene-size=?] [x-incr=?] [y-thresh-fixed=?] [y-max=?]
        [y-incr=?] [pvalue-threshold=?]
        [y-thresh-l-min=?] [y-thresh-l-dec=?]

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

        x-thresh-gene-size

            Set the threshold between "typical" genes and "large" genes.  Default 5000

        x-incr

            Down-increment of threshold delimiting typical and large genes.  Default 500

        y-thresh-fixed, y-thresh-l-min=f, y-thresh-l-dec

            Highest p-value for typical genes is also lower-bound for large genes (the
            actual p-value is exp(- y-thresh-fixed).  Default value 8 (corresponding to
            p-value = 0.00000001).  TODO: describe thresh-l-min, thresh-l-dec

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
