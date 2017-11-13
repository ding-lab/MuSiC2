package TGI::MuSiC2::Bmr;

use warnings;
use strict;

use IO::File;

use TGI::MuSiC2::CalcBmr; 
use TGI::MuSiC2::CalcCovg; 
use TGI::MuSiC2::CalcWigCovg; 
use TGI::MuSiC2::CalcCovgHelper;
use TGI::MuSiC2::CalcWindowRoi;
use TGI::MuSiC2::CalcWindowMaf;
use TGI::MuSiC2::CalcBmrModifier;

## process subcommands
sub new {
    my $class = shift;
    my $this = {};
    $this->{SUBCOMMAND} = shift;
    bless $this, $class;
    $this->process();

    return $this;
}

sub process {
    my $this = shift;
    my %cmds = map{
        ($_, 1)
    } qw( calc-bmr 
          calc-covg 
          calc-covg-helper 
          calc-wig-covg 
          calc-window-roi 
          calc-window-maf 
          calc-bmr-modifier );
    unless (defined $this->{SUBCOMMAND}) { die help_text(); };
    unless (exists $cmds{ $this->{SUBCOMMAND}}) {
        warn ' Please give valid sub command ! ', "\n";
        die help_text();
    }
    SWITCH:{
        $this->{SUBCOMMAND} eq 'calc-bmr'          && do { TGI::MuSiC2::CalcBmr->new();         last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-bmr-modifier' && do { TGI::MuSiC2::CalcBmrModifier->new(); last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-covg'         && do { TGI::MuSiC2::CalcCovg->new();        last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-covg-helper'  && do { TGI::MuSiC2::CalcCovgHelper->new();  last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-wig-covg'     && do { TGI::MuSiC2::CalcWigCovg->new();     last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-window-roi'   && do { TGI::MuSiC2::CalcWindowRoi->new();   last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-window-maf'   && do { TGI::MuSiC2::CalcWindowMaf->new();   last SWITCH; };

        $this->{SUBCOMMAND} eq 'help'       && do { die help_text(); last SWITCH; };
    }

}

## Bmr subcmds
sub help_text {
    my $this = shift;
        return <<HELP        
Sub-commands for music2 bmr:
calc-bmr             Calculates mutation rates given per-gene coverage (from     
                      "music2 bmr calc-covg"), and a mutation list                
calc-bmr-modifier    Calculates Backgroud Mutation Rate (BMR) modifier for gene/region specific
                       BMR estimation
calc-covg            Uses calcRoiCovg.c to count covered bases per-gene for each 
                      given tumor-normal pair of BAMs.                           
calc-covg-helper     Uses calcRoiCovg.c to count covered bases per-gene for a    
                      tumor-normal pair of BAMs.                                 
calc-wig-covg        Count covered bases per-gene for each given wiggle track    
                      format file.
calc-window-roi      Given gene/region ROIs, generates window based ROIs.
calc-window-maf      Given mutation list (MAF), generates window based MAF file.

HELP
}

1;

__END__
 
=head1 NAME

TGI::MuSiC2::Bmr - Sub-commands for music2 bmr.

=head1 SYNOPSIS

TGI::MuSiC2::Bmr - Sub-commands for music2 bmr.

=head1 DESCRIPTION

MuSiC2 - Mutational Significance in Cancer (Cancer Mutation Analysis) version 2. 

=head1 AUTHOR

Beifang Niu E<lt>beifang.cn@gmail.comE<gt>

=head1 SEE ALSO

https://github.com/ding-lab/MuSiC2

=head1 LICENSE

This library is free software with MIT licence; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

