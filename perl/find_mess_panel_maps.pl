use strict;
use warnings;

use 5.022;

use Path::Tiny qw /path/;

my $base_dir = 'output/maxent/HIA_TEST_INV_BS';
my @dirs;
#say "Globbing $base_dir/*";
#my @dirs = map {path($_)->basename} glob ("$base_dir/*");

open my $fh, '<', 'output/maxent/species_dirs.txt' or die $!;
while (my $line = <$fh>) {
    chomp $line;
    push @dirs, $line;
}

#@dirs = qw /
#Abies_firma
#Acer_cappadocicum
#Acer_crataegifolium
#Acer_circinatum
#Acer_pectinatum
#Acer_oliverianum
#Acer_campestre
#Acer_japonicum
#Acer_griseum
#Abies_pinsapo
#/;

my @exists;
foreach my $dir (sort @dirs) {
    
    my $exists = -e "$base_dir/$dir/full/${dir}_mess_panel.png";
    if ($exists) {
        push @exists, $dir;
        say $dir;
    }
}


