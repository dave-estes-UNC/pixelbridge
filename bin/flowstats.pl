#!/usr/bin/perl

my $SPAN_LENGTH = 240;

my @span, @spanTotals;
my $count = 0, $spanTotal = 0;
my $low, $high;

while (<STDIN>) {
    if (!/^FlowCSV,\d+,(\d+),(.*)$/) { next; }

    $frame = $1;
    $flow = $2;

    if ($flow =~ /inf/) {
	$flow = 0;
    }

    if ($count < $SPAN_LENGTH) {
	$spanTotal += $flow;
	push @span, $flow;
	$count++;
	if ($count == $SPAN_LENGTH) {
	    push @spanTotals, $spanTotal;
	    $low = $high = $spanTotal;
	}
    } else {
	$spanTotal -= shift @span;
	$spanTotal += $flow;
	push @span, $flow;
	push @spanTotals, $spanTotal;
	if ($spanTotal > $high) {
	    $high = $spanTotal;
	    #print "New High: " . ($frame - $SPAN_LENGTH) . " $flow\n";
	} elsif ($spanTotal < $low) {
	    $low = $spanTotal;
	    #print "New Low: " . ($frame - $SPAN_LENGTH) . " $flow\n";
	}
    }
    #print "$frame $flow $spanTotal\n";
}

my $lFrame, $mFrame = 2, $hFrame;
my $median = $low + ($high - $low) / 2;
foreach my $i (0 .. $#spanTotals) {
    if ($low == $spanTotals[$i]) {
	$lFrame = $i + 2;
    }
    if ($high == $spanTotals[$i]) {
	$hFrame = $i + 2;
    }
    if (abs ($median - $spanTotals[$i]) < abs ($median - $spanTotals[$mFrame-2])) {
	$mFrame = $i + 2;
	#print "New Median: $mFrame " . ($median - $spanTotals[$i]) . "\n";
    }
}

print "Low: $lFrame, $low\n";
print "High: $hFrame, $high\n";
print "Median: $mFrame, $median\n";
