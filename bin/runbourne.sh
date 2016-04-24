#!/bin/sh

##
## Movie DVD
##

# Low Flow
if [ ! -e output/bourne-13-dialog_LowFlow_LowCompression.csv]
then
    ./runauto.sh 5131 240 260640 ../recordings/acmmm2011/video/bourne-13-dialog.m4v > output/bourne-13-dialog_LowFlow_LowCompression.csv
fi
if [ ! -e output/bourne-13-dialog_LowFlow_MediumCompression.csv ]
then
    ./runauto.sh 5131 240 130320 ../recordings/acmmm2011/video/bourne-13-dialog.m4v >output/bourne-13-dialog_LowFlow_MediumCompression.csv
fi
if [ ! -e output/bourne-13-dialog_LowFlow_HighCompression.csv ]
then
    ./runauto.sh 5131 240 81450 ../recordings/acmmm2011/video/bourne-13-dialog.m4v > output/bourne-13-dialog_LowFlow_HighCompression.csv
fi

# Median Flow
if [ ! -e output/bourne-10-action_MedianFlow_LowCompression.csv ]
then
    ./runauto.sh 5258 240 260640 ../recordings/acmmm2011/video/bourne-10-action.m4v > output/bourne-10-action_MedianFlow_LowCompression.csv
fi
if [ ! -e output/bourne-10-action_MedianFlow_MediumCompression.csv ]
then
    ./runauto.sh 5258 240 130320 ../recordings/acmmm2011/video/bourne-10-action.m4v > output/bourne-10-action_MedianFlow_MediumCompression.csv
fi
if [ ! -e output/bourne-10-action_MedianFlow_HighCompression.csv ]
then
    ./runauto.sh 5258 240 81450 ../recordings/acmmm2011/video/bourne-10-action.m4v > output/bourne-10-action_MedianFlow_HighCompression.csv
fi

# High Flow
if [ ! -e output/bourne-10-action_HighFlow_LowCompression.csv ]
then
    ./runauto.sh 785 240 260640 ../recordings/acmmm2011/video/bourne-10-action.m4v > output/bourne-10-action_HighFlow_LowCompression.csv
fi
if [ ! -e output/bourne-10-action_HighFlow_MediumCompression.csv ]
then
    ./runauto.sh 785 240 130320 ../recordings/acmmm2011/video/bourne-10-action.m4v > output/bourne-10-action_HighFlow_MediumCompression.csv
fi
if [ ! -e output/bourne-10-action_HighFlow_HighCompression.csv ]
then
    ./runauto.sh 785 240 81450 ../recordings/acmmm2011/video/bourne-10-action.m4v > output/bourne-10-action_HighFlow_HighCompression.csv
fi
