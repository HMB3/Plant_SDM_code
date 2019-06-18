use strict;
use warnings;

use 5.022;
$|=1;
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



foreach my $dir (sort @dirs) {

    #say $species_name;
    my @files;
    #push @files, glob ("$base_dir/$dir/full/${dir}_*20?0*.tif");
    #push @files, glob ("$base_dir/$dir/full/${dir}_*20?0*_SUA_cell_count.csv");
    #push @files, glob ("$base_dir/$dir/full/${dir}_*20?0*_SUA_gain_loss_table.csv");
    push @files, glob ("$base_dir/$dir/full/${dir}_*20?0**.tif.aux.xml");
    
    #say $_ foreach @files;
    unlink $_ foreach @files;
}

