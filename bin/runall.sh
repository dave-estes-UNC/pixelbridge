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


##
## Trailer 720
##

# Low Flow
if [ ! -e output/captain-720_LowFlow_LowCompression.csv ]
then
    ./runauto.sh 2 240 696320 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_LowFlow_LowCompression.csv
fi
if [ ! -e output/captain-720_LowFlow_MediumCompression.csv ]
then
    ./runauto.sh 2 240 348160 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_LowFlow_MediumCompression.csv
fi
if [ ! -e output/captain-720_LowFlow_HighCompression.csv ]
then
    ./runauto.sh 2 240 217600 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_LowFlow_HighCompression.csv
fi

# Median Flow
if [ ! -e output/captain-720_MedianFlow_LowCompression.csv ]
then
    ./runauto.sh 1178 240 696320 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_MedianFlow_LowCompression.csv
fi
if [ ! -e output/captain-720_MedianFlow_MediumCompression.csv ]
then
    ./runauto.sh 1178 240 348160 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_MedianFlow_MediumCompression.csv
fi
if [ ! -e output/captain-720_MedianFlow_HighCompression.csv ]
then
    ./runauto.sh 1178 240 217600 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_MedianFlow_HighCompression.csv
fi

# High Flow
if [ ! -e output/captain-720_HighFlow_LowCompression.csv ]
then
    ./runauto.sh 3075 240 696320 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_HighFlow_LowCompression.csv
fi
if [ ! -e output/captain-720_HighFlow_MediumCompression.csv ]
then
    ./runauto.sh 3075 240 348160 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_HighFlow_MediumCompression.csv
fi
if [ ! -e output/captain-720_HighFlow_HighCompression.csv ]
then
    ./runauto.sh 3075 240 217600 ../recordings/acmmm2011/video/captain-720.mov > output/captain-720_HighFlow_HighCompression.csv
fi


##
## Movie 720p
##

# Low Flow
if [ ! -e output/limitless-720_LowFlow_LowCompression.csv ]
then
    ./runauto.sh 1109 240 696320 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_LowFlow_LowCompression.csv
fi
if [ ! -e output/limitless-720_LowFlow_MediumCompression.csv ]
then
    ./runauto.sh 1109 240 348160 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_LowFlow_MediumCompression.csv
fi
if [ ! -e output/limitless-720_LowFlow_HighCompression.csv ]
then
    ./runauto.sh 1109 240 217600 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_LowFlow_HighCompression.csv
fi

# Median Flow
if [ ! -e output/limitless-720_MedianFlow_LowCompression.csv ]
then
    ./runauto.sh 358 240 696320 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_MedianFlow_LowCompression.csv
fi
if [ ! -e output/limitless-720_MedianFlow_MediumCompression.csv ]
then
    ./runauto.sh 358 240 348160 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_MedianFlow_MediumCompression.csv
fi
if [ ! -e output/limitless-720_MedianFlow_HighCompression.csv ]
then
    ./runauto.sh 358 240 217600 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_MedianFlow_HighCompression.csv
fi

# High Flow
if [ ! -e output/limitless-720_HighFlow_LowCompression.csv ]
then
    ./runauto.sh 25 240 696320 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_HighFlow_LowCompression.csv
fi
if [ ! -e output/limitless-720_HighFlow_MediumCompression.csv ]
then
    ./runauto.sh 25 240 348160 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_HighFlow_MediumCompression.csv
fi
if [ ! -e output/limitless-720_HighFlow_HighCompression.csv ]
then
    ./runauto.sh 25 240 217600 ../recordings/acmmm2011/video/limitless-720.mov > output/limitless-720_HighFlow_HighCompression.csv
fi


##
## Trailer 1080
##

# Low Flow
if [ ! -e output/captain-1080_LowFlow_LowCompression.csv ]
then
    ./runauto.sh 3 240 1566720 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_LowFlow_LowCompression.csv
fi
if [ ! -e output/captain-1080_LowFlow_MediumCompression.csv ]
then
    ./runauto.sh 3 240 783360 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_LowFlow_MediumCompression.csv
fi
if [ ! -e output/captain-1080_LowFlow_HighCompression.csv ]
then
    ./runauto.sh 3 240 489600 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_LowFlow_HighCompression.csv
fi

# Median Flow
if [ ! -e output/captain-1080_MedianFlow_LowCompression.csv ]
then
    ./runauto.sh 2749 240 1566720 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_MedianFlow_LowCompression.csv
fi
if [ ! -e output/captain-1080_MedianFlow_MediumCompression.csv ]
then
    ./runauto.sh 2749 240 783360 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_MedianFlow_MediumCompression.csv
fi
if [ ! -e output/captain-1080_MedianFlow_HighCompression.csv ]
then
    ./runauto.sh 2749 240 489600 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_MedianFlow_HighCompression.csv
fi

# High Flow
if [ ! -e output/captain-1080_HighFlow_LowCompression.csv ]
then
    ./runauto.sh 2936 240 1566720 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_HighFlow_LowCompression.csv
fi
if [ ! -e output/captain-1080_HighFlow_MediumCompression.csv ]
then
    ./runauto.sh 2936 240 783360 ../recordings/acmmm2011/video/captain-1080.mov >output/captain-1080_HighFlow_MediumCompression.csv
fi
if [ ! -e output/captain-1080_HighFlow_HighCompression.csv ]
then
    ./runauto.sh 2936 240 489600 ../recordings/acmmm2011/video/captain-1080.mov > output/captain-1080_HighFlow_HighCompression.csv
fi


##
## Movie 1080
##

# Low Flow
if [ ! -e output/limitless-1080_LowFlow_LowCompression.csv ]
then
    ./runauto.sh 1013 240 1566720 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_LowFlow_LowCompression.csv
fi
if [ ! -e output/limitless-1080_LowFlow_MediumCompression.csv ]
then
    ./runauto.sh 1013 240 783360 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_LowFlow_MediumCompression.csv
fi
if [ ! -e output/limitless-1080_LowFlow_HighCompression.csv ]
then
    ./runauto.sh 1013 240 489600 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_LowFlow_HighCompression.csv
fi

# Median Flow
if [ ! -e output/limitless-1080_MedianFlow_LowCompression.csv ]
then
    ./runauto.sh 352 240 1566720 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_MedianFlow_LowCompression.csv
fi
if [ ! -e output/limitless-1080_MedianFlow_MediumCompression.csv ]
then
    ./runauto.sh 352 240 783360 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_MedianFlow_MediumCompression.csv
fi
if [ ! -e output/limitless-1080_MedianFlow_HighCompression.csv ]
then
    ./runauto.sh 352 240 489600 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_MedianFlow_HighCompression.csv
fi

# High Flow
if [ ! -e output/limitless-1080_HighFlow_LowCompression.csv ]
then
    ./runauto.sh 1400 240 1566720 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_HighFlow_LowCompression.csv
fi
if [ ! -e output/limitless-1080_HighFlow_MediumCompression.csv ]
then
    ./runauto.sh 1400 240 783360 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_HighFlow_MediumCompression.csv
fi
if [ ! -e output/limitless-1080_HighFlow_HighCompression.csv ]
then
    ./runauto.sh 1400 240 489600 ../recordings/acmmm2011/video/limitless-1080.mov > output/limitless-1080_HighFlow_HighCompression.csv
fi


##
## Trailer 4K
##

# Low Flow
if [ ! -e output/Elysium-4k_LowFlow_LowCompression.csv ]
then
    ./runauto.sh 2930 240 4147200 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_LowFlow_LowCompression.csv
fi
if [ ! -e output/Elysium-4k_LowFlow_MediumCompression.csv ]
then
    ./runauto.sh 2930 240 3110400 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_LowFlow_MediumCompression.csv
fi
if [ ! -e output/Elysium-4k_LowFlow_HighCompression.csv ]
then
    ./runauto.sh 2930 240 2073600 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_LowFlow_HighCompression.csv
fi

# Median Flow
if [ ! -e output/Elysium-4k_MedianFlow_LowCompression.csv ]
then
    ./runauto.sh 2495 240 4147200 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_MedianFlow_LowCompression.csv
fi
if [ ! -e output/Elysium-4k_MedianFlow_MediumCompression.csv ]
then
    ./runauto.sh 2495 240 3110400 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_MedianFlow_MediumCompression.csv
fi
if [ ! -e output/Elysium-4k_MedianFlow_HighCompression.csv ]
then
    ./runauto.sh 2495 240 2073600 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_MedianFlow_HighCompression.csv
fi

# High Flow
if [ ! -e output/Elysium-4k_HighFlow_LowCompression.csv ]
then
    ./runauto.sh 2647 240 4147200 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_HighFlow_LowCompression.csv
fi
if [ ! -e output/Elysium-4k_HighFlow_HediumCompression.csv ]
then
    ./runauto.sh 2647 240 3110400 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_HighFlow_HediumCompression.csv
fi
if [ ! -e output/Elysium-4k_HighFlow_HighCompression.csv ]
then
    ./runauto.sh 2647 240 2073600 ../recordings/2015/Elysium_trailer_3840x2160_3176f_24fps_RMN_QP24.mkv > output/Elysium-4k_HighFlow_HighCompression.csv
fi
