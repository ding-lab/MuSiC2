package TGI::MuSiC2::Smg;
use Moose;
use Moose::Util::TypeConstraints;
use MooseX::StrictConstructor;
use Statistics::Descriptive;
use Carp;
use POSIX qw( WIFEXITED );
use IO::File;
use File::Temp qw/ tempfile /;
use Getopt::Long;
use Try::Tiny;
use namespace::autoclean;

subtype 'ExistingFileName',
        as 'Str',
        where   { -s $_[0] },
        message { "Gene mutation rate file not found or is empty: $_[0]" };

has gene_mr_file              => (is => 'ro', isa => 'ExistingFileName', required => 1);
has output_file               => (is => 'ro', isa => 'Str', required => 1);
has max_fdr                   => (is => 'ro', isa => 'Num', default => 0.2);
has processors                => (is => 'ro', isa => 'Int',  default => 1);
has skip_low_mr_genes         => (is => 'ro', isa => 'Bool', default => 1);
has skip_non_expressed_genes  => (is => 'ro', isa => 'Bool', default => 1);
has downsample_large_genes    => (is => 'ro', isa => 'Bool', default => undef);
has bmr_modifier_file         => (is => 'ro', isa => 'ExistingFileName');
has black_gene_list_file      => (is => 'ro', isa => 'Str');

sub option_flag {
    my ($options, $flag) = @_;

    my $option_name = $flag;
    $option_name =~ s/-/_/g;

    return (
        $flag     => sub { $options->{$option_name} = $_[1] },
        "no$flag" => sub { $options->{$option_name} = undef },
    );
}

sub option_flags {
    my ($options, @flags) = @_;
    return map { option_flag($options, $_) } @flags;
}

sub option_input {
    my ($options, $input) = @_;

    my ($option_name, $option_type) = split(/=/, $input);
    $option_name =~ s/-/_/g;

    return (
        $input => sub { $options->{$option_name} = $_[1] },
    );
}

sub option_inputs {
    my ($options, @inputs) = @_;
    return map { option_input($options, $_) } @inputs;
}

sub execute {
    my $class = shift;

    my $has_args = @ARGV ? 1 : undef;
    my %options = ();

    my $had_good_options = GetOptions (
        option_inputs(\%options, qw(
            gene-mr-file=s
            output-file=s
            max-fdr=f
            bmr-modifier-file=s
            black-gene-list-file=s
            processors=i
            help
        )),
        option_flags(\%options, qw(
            skip-low-mr-genes
            skip-non-expressed-genes
            downsample-large-genes
        )),
    );

    if (!$has_args || !$had_good_options || $options{help}) {
        print help_text();
        exit($had_good_options ? 0 : 2);
    }

    try {
        $class->new(%options)->execute_smg();
        exit(0);
    }
    catch {
        my $error_message = shift;
        $error_message =~ s{ at /.*$}{}s;
        print STDERR "$error_message\n";
        exit(2);
    };
}

sub execute_smg {
    my $self = shift;

    print "Beifang, the number of processors is: " . $self->processors . "\n";
    use Data::Dumper;
    print Dumper($self);
}

sub help_text {
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

__PACKAGE__->meta->make_immutable;
1;

