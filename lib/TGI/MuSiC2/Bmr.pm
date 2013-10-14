package TGI::MuSiC2::Bmr;

use warnings;
use strict;

use IO::File;

use TGI::MuSiC2::CalcBmr; 
use TGI::MuSiC2::CalcCovg; 
use TGI::MuSiC2::CalcWigCovg; 
use TGI::MuSiC2::CalcCovgHelper;

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
    my %cmds = map{ ($_, 1) } qw( calc-bmr calc-covg calc-covg-helper calc-wig-covg );
    unless ( defined $this->{SUBCOMMAND} ) { die help_text(); };
    unless ( exists $cmds{ $this->{SUBCOMMAND} } ) {
        warn ' Please give valid sub command ! ', "\n";
        die help_text();
    }

    SWITCH:{
        $this->{SUBCOMMAND} eq 'calc-bmr'         && do { TGI::MuSiC2::CalcBmr->new();        last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-covg'        && do { TGI::MuSiC2::CalcCovg->new();       last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-covg-helper' && do { TGI::MuSiC2::CalcCovgHelper->new(); last SWITCH; };
        $this->{SUBCOMMAND} eq 'calc-wig-covg'    && do { TGI::MuSiC2::CalcWigCovg->new();    last SWITCH; };

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
calc-covg            Uses calcRoiCovg.c to count covered bases per-gene for each 
                      given tumor-normal pair of BAMs.                           
calc-covg-helper     Uses calcRoiCovg.c to count covered bases per-gene for a    
                      tumor-normal pair of BAMs.                                 
calc-wig-covg        Count covered bases per-gene for each given wiggle track    
                      format file.                                               

HELP
}

1;
