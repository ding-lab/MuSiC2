package TGI::MuSiC2::CalcBmrModifier;
##
# given genes mr file, calculate window based bmr modifier
##
use strict;
use warnings;
#
#
use IO::File;
use Getopt::Long;
use File::Temp qw/ tempfile /;

sub new {
    my $class = shift;
    my $this = {};

    $this->{_ROI_FILE} = undef;
    $this->{_GENE_MRS_FILE} = undef;
    $this->{_OUTPUT_MODIFIER_FILE} = "window_based_gene_bmr_modifier"; 
    $this->{_WINDOW_SIZE} = 1000000;

    bless $this, $class;
    $this->process();

    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (

        'roi-file=s'               => \$this->{_ROI_FILE},
        'gene-mrs-file=s'          => \$this->{_GENE_MRS_FILE},
        'output-modifier-file=s'   => \$this->{_OUTPUT_MODIFIER_FILE},

        'window-size=i'            => \$this->{_WINDOW_SIZE},

        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    #### processing ####
    #
    # Check on all the input data
    print STDERR "ROI file not found or is empty: $this->{_ROI_FILE}\n" unless( -s $this->{_ROI_FILE} );
    return undef unless( -s $this->{_ROI_FILE} );
    print STDERR "window based mutation rate file not found or is empty: $this->{_GENE_MRS_FILE}\n" unless( -s $this->{_GENE_MRS_FILE} );
    return undef unless( -s $this->{_GENE_MRS_FILE} );
    unless ( $this->{_WINDOW_SIZE} =~ /\d+/ ) {  
        print STDERR "Window size format is not valid !\n"; 
        return undef; 
    };
    #
    ## generate 0-based window based bmrs file 
    #
    my @cont = (); my %gene_mr_hash;
    my $gene_mrs_fh = IO::File->new( $this->{_GENE_MRS_FILE} ) or die "Can not open gene_mrs file. $!";
    while ( my $line = $gene_mrs_fh->getline ) {
        next if ($line =~ /^#Gene/);
        chomp($line); my @t = split /\t/, $line;
        $gene_mr_hash{'cov'}{$t[0]} += $t[2];
        $gene_mr_hash{'mut'}{$t[0]} += $t[3];
    }
    $gene_mrs_fh->close();

    foreach my $w (sort keys %{$gene_mr_hash{'mut'}}) {
        my $cont_line = "";
        my ($chr, $start) = split /:/, $w;
        $chr=~s/^CHR//;
        $start--;
        if ($gene_mr_hash{'mut'}{$w} >= 1) {
            $cont_line .= "$chr\t$start\t".($start+$this->{_WINDOW_SIZE})."\t$gene_mr_hash{'mut'}{$w}\t".int($gene_mr_hash{'cov'}{$w}/3)."\t".$gene_mr_hash{'mut'}{$w}/($gene_mr_hash{'cov'}{$w}/3)."\n"; ;
            push( @cont, $cont_line);
        }
    }
    my ( undef, $temp_0_based_win_bmrs ) = tempfile();
    my $win_bmrs_fh = IO::File->new( $temp_0_based_win_bmrs, ">" ) or die "Temporary file could not be created. $!";
    $win_bmrs_fh->print( @cont );
    $win_bmrs_fh->close();
    @cont = ();
    my $roi_fh = IO::File->new( $this->{_ROI_FILE} ) or die "Can not open ROI file. $!";
    while ( my $line = $roi_fh->getline ) {
        chomp( $line );
        my @t = split /\t/, $line;
        $t[1]--;
        push ( @cont, join("\t", @t)."\n" ); 
    }
    $roi_fh->close();
    my ( undef, $temp_0_based_roi_file ) = tempfile();
    my $new_rois_fh = IO::File->new( $temp_0_based_roi_file, ">" ) or die "Temporary file could not be created. $!";
    $new_rois_fh->print( @cont );
    $new_rois_fh->close();
    my ( undef, $temp_0_based_win_bmrs_sorted ) = tempfile();
    my ( undef, $temp_0_based_roi_file_sorted ) = tempfile();
    system( "joinx1.8 sort -i $temp_0_based_roi_file -o $temp_0_based_roi_file_sorted" );
    system( "joinx1.8 sort -i $temp_0_based_win_bmrs -o $temp_0_based_win_bmrs_sorted" );
    my ( undef, $temp_intersected_file ) = tempfile();
    ## do intersecting
    system( "joinx1.8 intersect --output-both -a $temp_0_based_roi_file_sorted -b $temp_0_based_win_bmrs_sorted -o $temp_intersected_file" );
    @cont = ();
    my $intersect_fh = IO::File->new( $temp_intersected_file ) or die "Temporary file could not be opened. $!";
    while ( my $line = $intersect_fh->getline ) { chomp( $line ); push( @cont, $line ); };
    $intersect_fh->close();
    my %contu = map{ @_ = split /\t/; ( join( "\t", @_[4..8] ), 1) } @cont;
    my ( $total_mr, $m, $c, $f2 ) = (0.0, 0, 0, 0);
    map{ @_ = split /\t/; $m += $_[3]; $c +=$_[4]; } keys %contu;
    $total_mr = $m/$c;
    my %mh = (); my %ch = ();
    map{ @_ = split /\t/; $_[9] = $_[9]/$total_mr; $f2 = abs($_[2] - $_[1]); $ch{$_[3]}+=$f2; $mh{$_[3]}+=($f2*$_[9]); } @cont;
    ## pour out modifier file
    #
    my $output_fh = IO::File->new( $this->{_OUTPUT_MODIFIER_FILE}, ">" ) or die "Output file could not be created. $!";
    map{ $output_fh->print( "$_\t", $mh{$_}/$ch{$_}, "\n"); }  sort keys %mh;
    $output_fh->close();
    #
    return 1;

}

## usage
sub help_text {
    my $this = shift;
        return <<HELP
USAGE
 music2 bmr calc-bmr-modifier --roi-file=? --gene-mrs-file=? [--output-modifier-file=?] [--window-size=?]

SYNOPSIS
 music2 bmr calc-bmr-modifier \
    --roi-file input_dir/roi_file \
    --gene-mrs-file input_dir/gene_mrs \

 music2 bmr calc-bmr-modifier \
    --roi-file input_dir/roi_file \
    --gene-mrs-file input_dir/gene_mrs \
    --window-size 1000000

REQUIRED INPUTS
  roi-file
    Tab delimited list of ROIs [chr start stop gene_name] (See DESCRIPTION) 
  gene-mrs-file
    Window based mutation rate file

OPTIONAL INPUTS
  output-modifier-file
    Output window based Backgroud Mutation Rate (BMR) modifier file  
    Default file name 'window_based_gene_bmr_modifier' if not specified
  window-size
    The window size for local gene/regional specific BMR calculation 
    Default value '1000000' if not specified

REQUIRED OUTPUTS
  No

DESCRIPTION
    Given window based gene BMRs file and roi file, this script calculates the BMR modifier for given ROI 
    regions.

ARGUMENTS

    --gene-mrs-file
       
      The window based mutation rate file which can be calculated by window based bmr calculation.

    --window-size

      The window size for local gene/regional specific BMR calculation. The typical size is 1M and
      user can choose appropriate window size based on input gene/region ROIs.

HELP

}

1;


