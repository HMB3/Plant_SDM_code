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


my @exists;
my %exists_hash;

if ($ARGV[0]) {
    my $existing_file = $ARGV[0];
    open my $ifh, '<', $existing_file or die $!;
    
    while (my $line = <$ifh>) {
        chomp $line;
        push @exists, $line;
        my $species = get_species_name($line);
        $exists_hash{$species}++;
    }
    $ifh->close;
}

foreach my $dir (sort @dirs) {
    my $species_name = get_species_name($dir);
    next if exists $exists_hash{$species_name};
    #say $species_name;
    #my @files = glob ("$base_dir/$dir/full/*2030*gain_loss*.tif");
    my @files = glob ("$base_dir/$dir/full/*_future_not_novel_*.tif");
    push @exists, @files;
    say $_ foreach @files;
}


sub get_species_name {
    my ($line) = @_;
    my $basename = path($line)->basename;
    my $species = $basename;
    if ($basename =~ m/^(.+)_gain_loss/) {
        $species = $1 // die "did not get species name from $line";
    }
    return $species;
}
