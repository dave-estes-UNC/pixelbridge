#!/bin/sh

./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v > $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 2:1,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 4:1,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 8:1,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 16:1,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 2:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 4:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 8:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 16:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 4:4,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 8:4,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 16:4,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 4:1,2:1,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 8:1,4:1,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 16:1,8:1,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 4:1,2:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 8:1,4:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 16:1,8:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 4:2,2:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 8:2,4:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 16:2,8:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 4:4,2:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 8:4,4:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
./build/objdir-release/pixelbridge --frames 100 --mode dct --dctscales 16:4,8:2,1:8 --dctdelta $1 --psnr ../recordings/acmmm2011/video/bourne-10-action.m4v >> $2
