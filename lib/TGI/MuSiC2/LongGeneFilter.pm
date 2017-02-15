package TGI::MuSiC2::LongGeneFilter;

# standard perl packages
use strict;
use warnings;

# libraries
use Statistics::Distributions;
use Scalar::Util qw(looks_like_number);


use Carp;  # Use carp/croak/confess for improved error handling http://perldoc.perl.org/Carp.html
use IO::File;
use Getopt::Long;

sub write_genes {
    my ($this, $fn, $genes) = @_;
    my $DF = IO::File->new( $fn, 'w' ) or die "Couldn't open $fn. $!\n";
    foreach my $gene (@{$genes}) {
        print $DF "$gene\n";
    }
    $DF->close;
}


# run filter algorithm by cycling through parameters and evaluate convergence.
# Write filtered genes to files if convergence occurs and return 1; return 0 if
# convergence fails.
sub apply_filter {
    my ($this, $pts, $gene_data, $pass_fn, $fail_fn, $lgf_fail_fn) = @_;
    my $iters = 0;

# y-threshold is variable for left-plane genes
    for (my $y_thresh_l = $this->{Y_THRESH_FIXED}; $y_thresh_l >= $this->{Y_THRESH_L_MIN}; $y_thresh_l -= $this->{Y_THRESH_L_DEC}) {
        print"#   trying y_thresh_l cutoff $y_thresh_l\n";

# x-threshold that delineates typical vs large genes
        for (my $x_thresh = $this->{X_THRESH_GENE_SIZE}; $x_thresh > $this->{X_INCR};
                $x_thresh -= $this->{X_INCR}) {

# move y-threshold upward for large genes until significance vanishes
            for (my $y_thresh_r = $y_thresh_l; $y_thresh_r <= $this->{Y_MAX};
                    $y_thresh_r += $this->{Y_INCR}) {

                if ($this->evaluate_long_filter($pts, $x_thresh, $y_thresh_l, $y_thresh_r)) {
                    my ($pass_genes, $fail_genes, $lgf_fail_genes) = $this->get_filtered_genes($gene_data, $x_thresh, $y_thresh_l, $y_thresh_r);

                    if ($pass_fn) { $this->write_genes($pass_fn, $pass_genes); }
                    if ($fail_fn) { $this->write_genes($fail_fn, $fail_genes); }
                    if ($lgf_fail_fn) { $this->write_genes($lgf_fail_fn, $lgf_fail_genes); }

                    print "   Total iterations = $iters\n";
                    return;
                } else {
                    $iters++;
                }
            }
        }
    }
    die("No convergence using current parameters\n");
}

# Evaluate whether filter converged with given parameters. If yes, print user diagnostics and return 1; otherwise return 0.
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

    print "\n   Total genes : $grand_total\n";
    print "   Large gene threshold: $x_thresh ";
    print "   (unchanged from default)" if $x_thresh eq $this->{X_THRESH_GENE_SIZE};
    print "\n";
    print "   Recommended -ln(P-value) cut-offs\n";
    print "      typical genes: -ln(P-val) = $y_thresh_l (P-val = $p_val_typical_gene)";
    print "   (unchanged from default)" if $y_thresh_l eq $this->{Y_THRESH_FIXED};
    print "\n";
    print "      large genes:   -ln(P-val) = $y_thresh_r (P-val = $p_val_large_gene)";
    print "   (unchanged from default)" if $y_thresh_r eq $this->{Y_THRESH_FIXED};
    print "\n";
    print "   Calculation meta-data\n";
    print "      final 2x2 table:   $top_left,  $top_right\n";
    print "                         $bot_left,  $bot_right\n";
    print "      final P-val $pval > threshold $this->{PVALUE_THRESHOLD} " .
        "(mutational significance status not statistically related to gene size)\n";
}

sub get_filtered_genes {    
    my ($this, $gene_data, $x_thresh, $y_thresh_l, $y_thresh_r) = @_;

    # list of genes which are unfiltered, filtered, and newly filtered by this
    # ("lgf") filter (because of an increased y_thresh_r).  Note that "newly
    # lgf filtered" is a subset of "filtered"
    my ($pass_genes, $fail_genes, $lgf_fail_genes) = ([], [], []); 

    foreach my $y (sort _numeric_ keys %{$gene_data}) {
        foreach my $x (keys %{$gene_data->{$y}}) {
            foreach my $gene_name (keys %{$gene_data->{$y}->{$x}}) {
                if ($x > $x_thresh && $y > $y_thresh_l && $y < $y_thresh_r) {
                    push(@{$lgf_fail_genes}, $gene_name);
                    push(@{$fail_genes}, $gene_name);
                } elsif ($y <= $y_thresh_l) {
                    push(@{$fail_genes}, $gene_name);
                } else {
                    push(@{$pass_genes}, $gene_name);
                }
            }
        }
    }
    return($pass_genes, $fail_genes, $lgf_fail_genes);
}


sub new {
    my $class = shift;
    my $this = {};

    $this->{DATA_FILE} = undef;
    $this->{SMG_FILE} = undef;
    $this->{SMG_COLUMN} = 9;
    $this->{GENE_SIZE_FILE} = undef;
    $this->{PASS_FILE} = undef;
    $this->{FAIL_FILE} = undef;
    $this->{LGF_FAIL_FILE} = undef;
    $this->{X_THRESH_GENE_SIZE} = 5000;
    $this->{X_INCR} = 500;
    $this->{Y_THRESH_FIXED} = 8;
    $this->{Y_THRESH_L_MIN} = 6;
    $this->{Y_THRESH_L_DEC} = 0.5;
    $this->{Y_MAX} = 18;
    $this->{Y_INCR} = 0.1;
    $this->{PVALUE_THRESHOLD} = 0.005;
    $this->{LOG_ZERO_P} = 25;

    bless $this, $class;
    $this->process();

    return $this;
}


sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (  # http://search.cpan.org/~chips/perl5.004_05/lib/Getopt/Long.pm
            'data-file=s'           => \$this->{DATA_FILE},
            'smg-file=s'            => \$this->{SMG_FILE},
            'smg-column=i'          => \$this->{SMG_COLUMN},
            'log-zero-p=f'          => \$this->{LOG_ZERO_P},
            'gene-size-file=s'      => \$this->{GENE_SIZE_FILE},
            'pass-file=s'           => \$this->{PASS_FILE},
            'fail-file=s'           => \$this->{FAIL_FILE},
            'lgf-fail-file=s'       => \$this->{LGF_FAIL_FILE},
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

    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->usage_text(); }

    my ($pts, $gene_data);
    if ($this->{DATA_FILE}) {
        ($pts, $gene_data) = $this->read_data_file($this->{DATA_FILE});
    } else {
        my $gene_sizes = $this->get_gene_sizes($this->{GENE_SIZE_FILE});
        my $smg_pval = $this->get_smg_pval($this->{SMG_FILE}, $this->{SMG_COLUMN});
        ($pts, $gene_data) = $this->merge_pval_size($gene_sizes, $smg_pval);
    }

    $this->print_header();
    $this->apply_filter($pts, $gene_data, $this->{PASS_FILE}, $this->{FAIL_FILE}, $this->{LGF_FAIL_FILE});

}

sub _numeric_ {$a <=> $b}  # used for string/number comparisons

# Read a TSV file with gene name, gene size in columns 1, 2.  Return hash with gene name
# as key and size as value.
sub get_gene_sizes {
    my ($this, $gene_size_fn) = @_;

    my $sizes = {};

    my $DF = IO::File->new( $gene_size_fn ) or die "Couldn't open $gene_size_fn. $!\n";
    while( my $line = $DF->getline ) {
        next if ($line =~ /^#/);
        chomp $line;
        my ($gene, $size) = split ' ', $line;
        die ("Size not a number in $gene_size_fn") if not looks_like_number($size);
        $sizes->{$gene} = $size;
    }
    $DF->close;

    return($sizes);
}

# Read a TSV file with gene name in first column and arbitrary number of additional columns.
# smg_column indicates column containing p-value data.
# Return hash with gene name as key, -log10(p-value) as data
# -log10(p-value) capped at --log-zero-p value
sub get_smg_pval {
    my ($this, $smg_file, $smg_column) = @_;

    # column defs are 1-index, perl is 0
    my $c = $smg_column - 1;

    my $pvals = {};

    my $DF = IO::File->new( $smg_file ) or die "Couldn't open $smg_file. $!\n";
    while(my $line = $DF->getline ) {
        next if ($line =~ /^#/);
        chomp $line;
        my @t = split ' ', $line;
        my $tlen = @t;
        die("Requested column ($smg_column) which does not exist in $smg_file\n") if ($c >= $tlen);
        my $pval = $t[$c];

        die ("P-value not a number in $smg_file") if not looks_like_number($pval);
        if ($pval == 0) {
            $pvals->{$t[0]} = $this->{LOG_ZERO_P};
        } else {
            $pvals->{$t[0]} = -log($pval)/log(10);# http://perldoc.perl.org/functions/log.html
        }
    }
    $DF->close;

    return($pvals);
}
    
# Add gene size info to pvalue data.  If size information for a given gene does not exist
# that is an error
#
# Returns ($pts, $gene_data):
# * pts is array of (X,Y) coordinates, with X gene size and Y=-log10(P-value)
# * gene_data is a hash of counts of X,Y,Gene_name
sub merge_pval_size {
    my ($this, $gene_sizes, $smg_pval) = @_;
    my ($pts, $gene_data) = ([], {});

    # iterate over all genes in smg_file
    foreach my $gene_name (keys %{$smg_pval}) {
        die("Size of gene $gene_name not defined in $this->{GENE_SIZE_FILE}\n") if not exists $gene_sizes->{$gene_name};
        my $x = $gene_sizes->{$gene_name};
        my $y = $smg_pval->{$gene_name};
        push @{$pts}, [$x, $y];
        $gene_data->{$y}->{$x}->{$gene_name}++;
    }
    return ($pts, $gene_data);
}


# Input data is 3 column format, tab separated:
# * gene_name in all upper-case
# * x = gene size (integer) 
# * y = -log10(MuSiC P-value)
# Returns ($pts, $gene_data):
# * pts is array of (X,Y) coordinates, with X gene size and Y=-log10(P-value)
# * gene_data is a hash of counts of X,Y,Gene_name

sub read_data_file {
    my ($this, $file) = @_;
    open (F, $file) || die "cant open $file";
    my ($pts, $gene_data) = ([], {});
    while (<F>) {

# parse line
        next if /^#/;
        chomp;
        my ($gene_name, $x, $y) = split;
        die ("Size not a number in $file") if not looks_like_number($x);
        die ("P-value not a number in $file") if not looks_like_number($y);
# store points 
        push @{$pts}, [$x, $y];
        $gene_data->{$y}->{$x}->{$gene_name}++;
    }
    close (F);
    return ($pts, $gene_data);
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

sub print_header{
    my $this = shift;
    my $s = <<HEADER;
# Adjustment of p-val threshold for large genes\n#
# music2 long-gene-filter \n#
# ************************************************************\n#
# Initial Parameters\n#
# default delineation between 'typical' and 'large' genes: $this->{X_THRESH_GENE_SIZE}
#    (note that this boundary may be modified for some cancers)
# highest allowable Y-threshold: $this->{Y_MAX}\n#
# P-value target above which mutational significance status
#      and gene size are no longer statistically related: $this->{PVALUE_THRESHOLD}\n#
# ************************************************************\n
HEADER
    print($s);
}

sub usage_text {
    my $this = shift;
    return <<HELP

        MuSiC2 Long Gene Filter Module

            Find conditions for which significance status is no longer related to gene
            size --- this is done by iteratively raising the cut-off for longer genes
            as compared to shorter ones

        USAGE 

            music2 long-gene-filter [--data-file=?] [--smg-file=?] [--smg-column=?]
            [--gene-size-file] [--pass-file=?] [--fail-file=?] [--lgf-fail-file=?]
            [--x-thresh-gene-size=?] [--x-incr=?] [--y-thresh-fixed=?] [--y-max=?]
            [--y-incr=?] [--pvalue-threshold=?] [log-zero-P=?]
            [--y-thresh-l-min=?] [--y-thresh-l-dec=?]

        SEE ALSO

            music2 long-gene-filter --help 

HELP
}

sub help_text {
    my $this = shift;
    my $usage = usage_text();

    return $usage . <<HELP

        ARGUMENTS

        data-file

            The data files is in TSV format with the following columns:
            * gene_name in all upper-case
            * gene size in base pairs (integer) 
            * -log10(MuSiC P-value)

            Must specify either data-file, or smg-file and gene-size file

        smg-file

            Data file in TSV format specifying MuSiC P-values, typically as
            generated by smg module.  First column specifies genes, Nth column specifies
            the P-values, with N defined by smg-column.  Genes must be unique.  Argument 
            required if data-file not defined.

        smg-column
            
            If reading smg-file, define column which contains P-values used for calculations.
            First column (containing genes) is 1.  Default 9.

        log-zero-P

            If smg-file has a P-value of 0, define value -log10(P-value) (since -infinity is 
            inconvenient).  Default 25.

        gene-size-file
    
            File in TSV format defining gene sizes in base pairs (col 2) for given gene (col 1).
            Genes must be unique, and all genes in smg-file must exist in gene-size-file. Argument
            requred if data-file not defined.

        pass-file, fail-file, lgf-fail-file

            Output filenames.  Genes which pass, fail, and fail because of Long Gene Filter
            (i.e., increased y_thresh_r) get written to these files respectively.  See below 
            for details.  Optional.

        x-thresh-gene-size

            Set the threshold between "typical" genes and "large" genes.  Default 5000

        x-incr

            Down-increment of threshold delimiting typical and large genes.  Default 500

        y-thresh-fixed

            Highest p-value for typical genes is also lower-bound for large genes (the
            actual p-value is exp(- y-thresh-fixed).  Default value 8 (corresponding to
            p-value = 0.00000001).  

        y-thresh-l-min=f, y-thresh-l-dec

            TODO: describe thresh-l-min, thresh-l-dec
            Defaults: 6, 0.5, resp.

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
        
        Visual description of algorithm:
         
                      |. . . . . . . . . . . . . . .
                      | . . . . . . . . . . . . . . 
                      |. . . . . . .================    y_thresh_r
               -ln(P) | . . . . . . " " " " " " " " 
                      |============= " " " " " " " "    y_thresh_l
                      | ' ' ' ' ' ' ' ' ' ' ' ' ' ' 
                      |' ' ' ' ' ' ' ' ' ' ' ' ' ' '
                      --------------]----------------- gene size
                                x_thresh
             
              Each gene exists as a point on a (gene size, -ln(P)) plane, and its significance
              is defined by the domain in which it falls on the plane.
              Absent the long gene filter, y_thresh_l = y_thresh_r, and all genes with -lnP < y_thresh_l
              are considered not significant.
             
              The long gene filter adds an x_thresh size cutoff, beyond which 
              y_thresh_r > y_thresh_l determines significance:  
             
              * Significant genes (indicated by . above, aka Pass_Genes) have 
                -lnP > y_thresh_l if size < x_thresh, and -lnP > y_thresh_r if size > x_thresh.  
             
              * Insignificant genes (indicated by ' and " above, aka Fail_Genes) are the complement 
                of this set.  
              
              Genes which are significant prior to this filter but insignificant after this
              filter is applied (indicated by " above) are a subset of Fail_Genes and are
              listed as New_Fail_Genes

              Algorithm proceeds by adjusting (in order), y_thresh_r, x_thresh, y_thresh_l,
              until the 2x2 table test no longer indicates a correlation between p-value and
              gene size. (TODO - clarify/confirm)
             
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
