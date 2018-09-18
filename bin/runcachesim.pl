#!/usr/bin/perl

use strict;
use warnings;
use threads;
use threads::shared;
use Thread::Queue;

my $usage =
    "\n".
    "Usage: runcachesim.sh [--pb <path-to-pixelbridge>] [--cs <path-to-nddiCacheSim.py>]\n".
    "                      [--tests <comma-separated-list>] [--sets <comma-separated-list>]\n".
    "                      <recordings-dir> <output-dir>\n".
    "  - Must provide either --pb, --cs, or both. If just --cs, then it assumes that there are .json.gz files in the out directory already.\n".
    "  - If --tests is not used, then all tests are run. Similar for --sets.\n".
    "  - --tests takes a comma-separated list from: fb, flat, cache8, cache4, cache2, dct1, dct4, dct100.\n".
    "  - --sets takes a comma-separated list from: files (desktop recordings), video (video clips), rewind (video clips with section rewinding).\n".
    "\n";

my ($pb, $cs, %opt_tests, %opt_sets, $rec_dir, $out_dir);

while (my $arg = shift @ARGV) {
    if ($arg eq "--help") {
        die $usage;
    } elsif ($arg eq "--pb") {
        $pb = shift @ARGV;
    } elsif ($arg eq "--cs") {
        $cs = shift @ARGV;
    } elsif ($arg eq "--tests") {
        %opt_tests = map { $_ => 1 } split "," , shift @ARGV;
    } elsif ($arg eq "--sets") {
        %opt_sets = map { $_ => 1} split "," , shift @ARGV;
    } elsif (not defined $rec_dir) {
        $rec_dir = $arg;
    } elsif (not defined $out_dir) {
        $out_dir = $arg;
    } else {
        die $usage;
    }

}

if (((not defined $pb) and (not defined $cs)) or (not defined $rec_dir) or (not defined $out_dir)) {
    die $usage;
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
my $count = 0;

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

my $workQ = Thread::Queue->new();
my $printQ = Thread::Queue->new();

my @workerThreads;
my %workerThreadsBusy :shared;

# worker thread routine
sub worker {
    # Thread will loop until no more work
    while (defined(my $work = $workQ->dequeue())) {
        $workerThreadsBusy{threads->tid()} = 1;
        my ($workCount,$workJson, $workCsv, $workOutDir, $workOutFile) = split('\|\|', $work);

        if (not -e "./$workJson") {
            system "gunzip -k $workOutDir/${workJson}.gz";
            system "mv $workOutDir/${workJson} .";
        }

        print "\n($workCount/$total) Running... $cs --ijson $workJson\n\n";
        my @cache_results = `$cs --ijson $workJson`;
        #DEBUG my @cache_results = `$cs --ijson --limit 2 $workJson`;
        foreach my $result (@cache_results) {
            print $result;
            if ($result =~ /\( *(\d*)B\).*\( *(\d*)B\).*\( *(\d*)B\).*\( *(\d*)B\).*\( *(\d*)B\)/) {
                $workCsv = "$workCsv,$1,$2,$3,$4,$5";
            }
        }
        $printQ->enqueue("$workCsv||$workOutFile");
        if (not -e "$workOutDir/${workJson}.gz") {
            system "gzip $workJson";
            system "mv ${workJson}.gz $workOutDir";
        } else {
            system "rm $workJson";
        }
        $workerThreadsBusy{threads->tid()} = 0;
    }
};

# start worker threads and print thread
for (1..4) {
    my $thr = threads->create(\&worker);
    push @workerThreads, $thr;
    $workerThreadsBusy{$thr->tid()} = 0;
}

my $printThread = threads->create(
    sub {
        while (defined(my $work = $printQ->dequeue())) {
            my ($workCsv, $workOutFile) = split('\|\|', $work);
            open my $OUT, '>>', $workOutFile;
            print $OUT "$workCsv\n";
            close $OUT;
        }
    }
    );

sub waitWorkerThreadsIdle {
    my $busy = 0;
    while (($workQ->pending() > 0) || $busy) {
        $busy = 0;
        for my $tid (keys %workerThreadsBusy) {
            $busy = $busy || $workerThreadsBusy{$tid};
        }
        sleep(1);
    }
};

sub run {
    my ($test, $opt, $out, $group, @items) = @_;
    my $logfile = "$out/$test-$group.log";

    #
    # First Pass: PixelBridge
    #
    if (defined $pb) {
        # Run pixelbridge to collect the cost model statistics and individual memory accesses
        foreach (@items) {
            chomp;
            $count = $count + 1;
            print "\n($count/$total) Running... $pb $opt $logging_opt $rec_dir/$_ >> $logfile\n\n";
            system("$pb $opt $logging_opt $rec_dir/$_ >> $logfile") == 0
                or die "FATAL: pixelbridge crashed.";
            if (not defined $cs) {
                system("gzip *.json");
                system("mv *.gz $out");
            }
        }
    }

    #
    # Second Pass: Cache Sim
    #
    if (defined $cs) {
        # Find the .json files containing the individual memory accesses. Archive them as gzips.
        # and output into the CSV along with the cost model statistics.
        my $json = "";
        my $result = "";
        open my $LOG, '<', $logfile;

        # Print the CSV headers
        my $csvfile = "$out/$test-$group.csv";
        open my $CSV, '>', $csvfile;
        print $CSV "Frames,Commands Sent,Bytes Transmitted,IV Num Reads,IV Bytes Read,IV Num Writes,IV Bytes Written,CP Num Reads,CP Bytes Read,CP Num Writes,CP Bytes Written,FV Num Reads,FV Bytes Read,FV Num Writes,FV Bytes Written,Pixels Mapped,Pixels Blended,PSNR,File Name,L1-Hit,L1-Miss,L1-Load,L1-Store,L1-Evict,L2-Hit,L2-Miss,L2-Load,L2-Store,L2-Evict,L3-Hit,L3-Miss,L3-Load,L3-Store,L3-Evict,Mem-Hit,Mem-Miss,Mem-Load,Mem-Store,Mem-Evict\n";
        close $CSV;

        # Process the log and enqueue cache sim work for the worker threads
        while (<$LOG>) {
            if (/,$rec_dir\/(.*)/) {
                chomp;
                $result = $_;
            } elsif (/(costlog-.*\.json)/) {
                $json = $1;
                $count = $count + 1;
                $workQ->enqueue("$count||$json||$result||$out||$csvfile");
            }
        }

        close $LOG;

        # Wait for the worker/print queues to drain before running the next set of pixelbridge experiments
        waitWorkerThreadsIdle();
    }
}

for my $test (keys %tests) {
    if (not defined %opt_tests or exists $opt_tests{$test}) {
        print "\nRunning: '$test' with arguments: $tests{$test}.\n\n";

        my $out_subdir = "$out_dir/$test";
        mkdir $out_subdir;

        run($test, $tests{$test}." ".$file_opt, $out_subdir, "files", @files) if (not defined %opt_sets or exists $opt_sets{"files"});
        run($test, $tests{$test}." ".$video_opt, $out_subdir, "video", @videos) if (not defined %opt_sets or exists $opt_sets{"video"});
        run($test, $tests{$test}." ".$rewind_opt, $out_subdir, "rewind", @videos) if (not defined %opt_sets or exists $opt_sets{"rewind"});
    }
}

# Shutdown the queues and threads
$workQ->end();
foreach my $thr (@workerThreads) {
    $thr->join();
}
$printQ->end();
$printThread->join();
