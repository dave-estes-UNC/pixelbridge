#!/bin/sh

./build/objdir-release/pixelbridge --mode flow --csv ../recordings/acmmm2011/video/bourne-10-action.m4v > output/bourne-10-action_flow.csv
./build/objdir-release/pixelbridge --mode flow --csv ../recordings/acmmm2011/video/bourne-2-moderate.m4v > output/bourne-2-moderate_flow.csv
./build/objdir-release/pixelbridge --mode flow --csv ../recordings/acmmm2011/video/bourne-13-dialog.m4v > output/bourne-13-dialog_flow.csv
./build/objdir-release/pixelbridge --mode flow --csv ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_flow.csv
./build/objdir-release/pixelbridge --mode flow --csv ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_flow.csv
./build/objdir-release/pixelbridge --mode flow --csv ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_flow.csv
./build/objdir-release/pixelbridge --mode flow --csv ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_flow.csv
./build/objdir-release/pixelbridge --mode flow --csv ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_flow.csv

./flowstats.pl < output/bourne-10-action_flow.csv
./flowstats.pl < output/bourne-2-moderate_flow.csv
./flowstats.pl < output/bourne-13-dialog_flow.csv
./flowstats.pl < output/captain-720_flow.csv
./flowstats.pl < output/captain-1080_flow.csv
./flowstats.pl < output/limitless-720_flow.csv
./flowstats.pl < output/limitless-1080_flow.csv
./flowstats.pl < output/Elysium-4k_flow.csv
