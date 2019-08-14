#!/bin/sh

if [ $# -ne 4 ]
then
	echo "Usage: runauto.sh <start> <frames> <budget> <infile>"
	exit 1
fi

echo "GroupCSV,$4,$1,$2,1:8,$3,Snap,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 1:8 --dctsnap --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,1:8,$3,Trim,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 1:8 --dcttrim --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,1:8,$3,Snap,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 1:8 --dctsnap --dctplanes 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,1:8,$3,Trim,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 1:8 --dcttrim --dctplanes 0 --psnr --csv $4

echo "GroupCSV,$4,$1,$2,16:1 1:8,$3,Snap,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 16:1,1:8 --dctsnap --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,16:1 1:8,$3,Trim,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 16:1,1:8 --dcttrim --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,16:1 1:8,$3,Snap,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 16:1,1:8 --dctsnap --dctplanes 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,16:1 1:8,$3,Trim,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 16:1,1:8 --dcttrim --dctplanes 0 --psnr --csv $4

echo "GroupCSV,$4,$1,$2,2:2 1:8,$3,Snap,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 2:2,1:8 --dctsnap --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,2:2 1:8,$3,Trim,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 2:2,1:8 --dcttrim --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,2:2 1:8,$3,Snap,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 2:2,1:8 --dctsnap --dctplanes 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,2:2 1:8,$3,Trim,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 2:2,1:8 --dcttrim --dctplanes 0 --psnr --csv $4

echo "GroupCSV,$4,$1,$2,4:3 1:8,$3,Snap,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 4:3,1:8 --dctsnap --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,4:3 1:8,$3,Trim,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 4:3,1:8 --dcttrim --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,4:3 1:8,$3,Snap,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 4:3,1:8 --dctsnap --dctplanes 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,4:3 1:8,$3,Trim,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 4:3,1:8 --dcttrim --dctplanes 0 --psnr --csv $4

echo "GroupCSV,$4,$1,$2,8:2 4:2 1:8,$3,Snap,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 8:2,4:2,1:8 --dctsnap --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,8:2 4:2 1:8,$3,Trim,Delta"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 8:2,4:2,1:8 --dcttrim --dctdelta 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,8:2 4:2 1:8,$3,Snap,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 8:2,4:2,1:8 --dctsnap --dctplanes 0 --psnr --csv $4
echo "GroupCSV,$4,$1,$2,8:2 4:2 1:8,$3,Trim,Planes"
./build/objdir-release/pixelbridge --start $1 --frames $2 --mode dct --dctbudget $3 --dctscales 8:2,4:2,1:8 --dcttrim --dctplanes 0 --psnr --csv $4
