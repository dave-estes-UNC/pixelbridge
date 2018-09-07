#!/usr/bin/perl

use strict;
use warnings;

my ($pb, $rec_dir, $out_dir) = @ARGV;

if ((not defined $pb) or (not defined $rec_dir) or (not defined $out_dir)) {
  die "Usage: runcachesim.sh <path-to-pixelbridge> <recordings-dir> <output-dir>"
}

open my $VIDEOS, '<', "$rec_dir/videos";
my @videos = <$VIDEOS>;
close $VIDEOS;

open my $FILES, '<', "$rec_dir/files";
my @files = <$FILES>;
close $FILES;

my %tests = (
    "fb" => "--mode fb",
    "flat" => "--mode flat",
    "cache8" => "--mode cache --tc 10000 --bits 8",
    "cache4" => "--mode cache --tc 10000 --bits 4",
    "cache2" => "--mode cache --tc 10000 --bits 2",
    "dct1" => "--mode dct --quality 1",
    "dct4" => "--mode dct --quality 5",
    "dct100" => "--mode dct --quality 100"
    );

mkdir $out_dir;

my $logging_opt = "--logcosts fv --headless --csv";
#my $file_opt = "--start 200";
#my $video_opt = "--start 300 --frames 800";
#my $rewind_opt = "--start 300 --frames 1200 --rewind 800 200";
my $file_opt = "--start 200 -- frames 240";
my $video_opt = "--start 300 --frames 240";
my $rewind_opt = "--start 300 --frames 240 --rewind 120 60";

for my $test (keys %tests) {
    print "\nRunning: '$test' with arguments: $tests{$test}.\n\n";

    my $out_subdir = "$out_dir/$test";
    mkdir $out_subdir;

    foreach (@files) {
	chomp;
	print "\nRunning $pb $file_opt $logging_opt $rec_dir/$_ >> $out_subdir/$test.log\n\n";
	system "$pb $file_opt $logging_opt $rec_dir/$_ >> $out_subdir/$test.log";
    }
    system "for x in `ls costlog-*.json`; do gzip \$x; mv \$x.gz $out_subdir/; done";

    foreach (@videos) {
	chomp;
	print "\nRunning $pb $video_opt $logging_opt $rec_dir/$_ >> $out_subdir/$test.log\n\n";
	system "$pb $video_opt $logging_opt $rec_dir/$_ >> $out_subdir/$test.log";
    }
    system "for x in `ls costlog-*.json`; do gzip \$x; mv \$x.gz $out_subdir/; done";

    foreach (@videos) {
	chomp;
	print "\nRunning $pb $rewind_opt $logging_opt $rec_dir/$_ >> $out_subdir/$test.log\n\n";
	system "$pb $rewind_opt $logging_opt $rec_dir/$_ >> $out_subdir/$test.log";
    }
    system "for x in `ls costlog-*.json`; do gzip \$x; mv \$x.gz $out_subdir/; done";
}
