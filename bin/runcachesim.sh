#!/bin/sh

# Runs all of the experiments while simulating cachec. This includes the recordings, the videos
# with normal playback, and the videos with rewind playbayback

if [ $# -ne 3 ]
then
	echo "Usage: runcachesim.sh <path-to-pixelbridge> <recordings-dir> <output-dirt>"
	exit 1
fi

PB=$1
REC_DIR=$2
OUT_DIR=$3

# For testing and exiting early
for x in `cat $REC_DIR/files`; do $PB --mode fb --start 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/fb.log
#for x in `cat $REC_DIR/videos`; do $PB --mode fb --start 300 --frames 800 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/video-fb.log
#for x in `cat $REC_DIR/videos`; do $PB --mode fb --start 300 --frames 1200 --rewind 800 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/rewind-fb.log
for x in `ls costlog-*.json`; do gzip $x; mv $x.gz $OUT_DIR; done
exit;

# Framebuffer Mode
for x in `cat $REC_DIR/files`; do $PB --mode fb --start 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/fb.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --start 300 --frames 800 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/video-fb.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --start 300 --frames 1200 --rewind 800 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/rewind-fb.log

# Flat Tiled Mode
for x in `cat $REC_DIR/files`; do $PB --mode flat --start 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/flat.log
for x in `cat $REC_DIR/videos`; do $PB --mode flat --start 300 --frames 800 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/video-flat.log
for x in `cat $REC_DIR/videos`; do $PB --mode flat --start 300 --frames 1200 --rewind 800 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/rewind-flat.log

# Cached Tiled Mode (8)
for x in `cat $REC_DIR/files`; do $PB --mode cache --tc 10000 --bits 8 --start 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/cache8.log
for x in `cat $REC_DIR/videos`; do $PB --mode cache --tc 10000 --bits 8 --start 300 --frames 800 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/video-cache8.log
for x in `cat $REC_DIR/videos`; do $PB --mode cache --tc 10000 --bits 8 --start 300 --frames 1200 --rewind 800 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/rewind-cache8.log

# Cached Tiled Mode (4)
for x in `cat $REC_DIR/files`; do $PB --mode cache --tc 10000 --bits 4 --psnr --start 200 $REC_DIR/$x; done > $OUT_DIR/cache4.log
for x in `cat $REC_DIR/videos`; do $PB --mode cache --tc 10000 --bits 4 --psnr --start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-cache4.log
for x in `cat $REC_DIR/videos`; do $PB --mode cache --tc 10000 --bits 4 --psnr --start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-cache4.log

# Cached Tiled Mode (2)
for x in `cat $REC_DIR/files`; do $PB --mode cache --tc 10000 --bits 2 --psnr --start 200 $REC_DIR/$x; done > $OUT_DIR/cache2.log
for x in `cat $REC_DIR/videos`; do $PB --mode cache --tc 10000 --bits 2 --psnr --start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-cache2.log
for x in `cat $REC_DIR/videos`; do $PB --mode cache --tc 10000 --bits 2 --psnr --start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-cache2.log

# DCT Tiled Mode (1)
for x in `cat $REC_DIR/files`; do $PB --mode dct --quality 1 --psnr --start 200 $REC_DIR/$x; done > $OUT_DIR/dct1.log
for x in `cat $REC_DIR/videos`; do $PB --mode dct --quality 1 --psnr --start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-dct1.log
for x in `cat $REC_DIR/videos`; do $PB --mode dct --quality 1 --psnr --start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-dct1.log

# DCT Tiled Mode (4)
for x in `cat $REC_DIR/files`; do $PB --mode dct --quality 4 --psnr --start 200 $REC_DIR/$x; done > $OUT_DIR/dct4.log
for x in `cat $REC_DIR/videos`; do $PB --mode dct --quality 4 --psnr --start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-dct4.log
for x in `cat $REC_DIR/videos`; do $PB --mode dct --quality 4 --psnr --start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-dct4.log

# DCT Tiled Mode (100)
for x in `cat $REC_DIR/files`; do $PB --mode dct --quality 100 --psnr --start 200 $REC_DIR/$x; done > $OUT_DIR/dct100.log
for x in `cat $REC_DIR/videos`; do $PB --mode dct --quality 100 --psnr --start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-dct100.log
for x in `cat $REC_DIR/videos`; do $PB --mode dct --quality 100 --psnr --start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-dct100.log

# Count Mode (Perfect Pixel Latching)
for x in `cat $REC_DIR/files`; do $PB --mode count --start 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/count.log
for x in `cat $REC_DIR/videos`; do $PB --mode count for--start 300 --frames 800 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/video-count.log
for x in `cat $REC_DIR/videos`; do $PB --mode count –start 300 --frames 1200 --rewind 800 200 --logcosts fv $REC_DIR/$x; done > $OUT_DIR/rewind-count.log

# Comment this out to run the blending experiments
exit

# Blending Temporal
for x in `cat $REC_DIR/files`; do $PB --mode fb --blend t --start 200 $REC_DIR/$x; done > $OUT_DIR/blend-temporal.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --blend t--start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-blend-temporal.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --blend t--start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-blend-temporal.log

# Blending Coefficient Plane
for x in `cat $REC_DIR/files`; do $PB --mode fb --blend cp --start 200 $REC_DIR/$x; done > $OUT_DIR/blend-coefficient-plane.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --blend cp--start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-blend-coefficient-plane.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --blend cp--start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-blend-coefficient-plane.log

# Blending Frame Volume
for x in `cat $REC_DIR/files`; do $PB --mode fb --blend fv --start 200 $REC_DIR/$x; done > $OUT_DIR/blend-frame-volume.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --blend fv--start 300 --frames 800 $REC_DIR/$x; done > $OUT_DIR/video-blend-frame-volume.log
for x in `cat $REC_DIR/videos`; do $PB --mode fb --blend fv--start 300 --frames 1200 --rewind 800 200 $REC_DIR/$x; done > $OUT_DIR/rewind-frame-volume.log
