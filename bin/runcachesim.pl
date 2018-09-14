#!/usr/bin/perl

use strict;
use warnings;

my ($pb, $cs, $rec_dir, $out_dir) = @ARGV;

if ((not defined $pb) or (not defined $rec_dir) or (not defined $out_dir)) {
  die "Usage: runcachesim.sh <path-to-pixelbridge> <path-to-nddiCacheSim.py> <recordings-dir> <output-dir>"
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
    "dct4" => "--mode dct --quality 4",
    "dct100" => "--mode dct --quality 100"
    );

my $total = 8 * (33 + 7 + 7) * 2;
my $current = 1;

mkdir $out_dir;

my $logging_opt = "--logcosts fv --headless --csv";

# These represent the original options passed for the ACMMM2011 paper
#my $file_opt = "--start 200";
#my $video_opt = "--start 300 --frames 800";
#my $rewind_opt = "--start 300 --frames 1200 --rewind 800 200";

# For the cache analysis, we're only going to use 10 seconds of playback
my $file_opt = "--start 200 --frames 240";
my $video_opt = "--start 300 --frames 240";
my $rewind_opt = "--start 300 --frames 240 --rewind 120 60";

# For testing, just use a few frames
#DEBUG my $file_opt = "--start 200 --frames 1";
#DEBUG my $video_opt = "--start 300 --frames 1";
#DEBUG my $rewind_opt = "--start 300 --frames 10 --rewind 8 2";

sub run {
    my ($test, $opt, $out, $group, @items) = @_;

    # Run pixelbridge to collect the cost model statistics and individual memory accesses
    foreach (@items) {
	chomp;
	print "\n($current/$total) Running... $pb $opt $logging_opt $rec_dir/$_ >> $out/$test-$group.log\n\n";
	$current = $current + 1;
	system("$pb $opt $logging_opt $rec_dir/$_ >> $out/$test-$group.log") == 0
	    or die "FATAL: pixelbridge crashed.";
    }

    # Find the .json files containing the individual memory accesses. Archive them as gzips
    # and output into the CSV along with the cost model statistics.
    my $item = "";
    my $json = "";
    open my $LOG, '<', "$out/$test-$group.log";
    open my $CSV, '>', "$out/$test-$group.csv";
    print $CSV "Frames,Commands Sent,Bytes Transmitted,IV Num Reads,IV Bytes Read,IV Num Writes,IV Bytes Written,CP Num Reads,CP Bytes Read,CP Num Writes,CP Bytes Written,FV Num Reads,FV Bytes Read,FV Num Writes,FV Bytes Written,Pixels Mapped,Pixels Blended,PSNR,File Name,L1-Hit,L1-Miss,L1-Load,L1-Store,L1-Evict,L2-Hit,L2-Miss,L2-Load,L2-Store,L2-Evict,L3-Hit,L3-Miss,L3-Load,L3-Store,L3-Evict,Mem-Hit,Mem-Miss,Mem-Load,Mem-Store,Mem-Evict\n";
    while (<$LOG>) {
	if (/,$rec_dir\/(.*)/) {
	    chomp;
	    print $CSV $_;
	    $item = $1;
	} elsif (/(costlog-.*\.json)/) {
	    $json = $1;
	    print "\n($current/$total) Running... $cs $json\n\n";
	    $current = $current + 1;
	    my @cache_results = `$cs $json`;
	    #DEBUG my @cache_results = `$cs --limit 2 $json`;
	    foreach my $result (@cache_results) {
		print $result;
		if ($result =~ /\( *(\d*)B\).*\( *(\d*)B\).*\( *(\d*)B\).*\( *(\d*)B\).*\( *(\d*)B\)/) {
		    print $CSV ",$1,$2,$3,$4,$5";
		}
	    }
	    print $CSV "\n";
	    system "gzip $json";
	    system "mv ${json}.gz $out";
	}
    }
    close $CSV;
    close $LOG;
}

for my $test (keys %tests) {
    print "\nRunning: '$test' with arguments: $tests{$test}.\n\n";

    my $out_subdir = "$out_dir/$test";
    mkdir $out_subdir;

    run($test, $tests{$test}." ".$file_opt, $out_subdir, "files", @files);
    run($test, $tests{$test}." ".$video_opt, $out_subdir, "video", @videos);
    run($test, $tests{$test}." ".$rewind_opt, $out_subdir, "rewind", @videos);
}
