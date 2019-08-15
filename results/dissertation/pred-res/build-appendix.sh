#!/bin/bash

csvpath=$1
appendix=$2
images=$3

echo "\chapter{Prediction Residual Decoder Results}" > $appendix
echo "\label{chap:pred_res}" >> $appendix

python pred-res-grapher.py "Bourne DVD Low Flow and Low Compression" $csvpath/bourne-13-dialog_LowFlow_LowCompression.csv $csvpath/bourne-13-dialog_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD Low Flow and Medium Compression" $csvpath/bourne-13-dialog_LowFlow_MediumCompression.csv $csvpath/bourne-13-dialog_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD Low Flow and High Compression" $csvpath/bourne-13-dialog_LowFlow_HighCompression.csv $csvpath/bourne-13-dialog_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD Median Flow and Low Compression" $csvpath/bourne-10-action_MedianFlow_LowCompression.csv $csvpath/bourne-10-action_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD Median Flow and Medium Compression" $csvpath/bourne-10-action_MedianFlow_MediumCompression.csv $csvpath/bourne-10-action_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD Median Flow and High Compression" $csvpath/bourne-10-action_MedianFlow_HighCompression.csv $csvpath/bourne-10-action_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD High Flow and Low Compression" $csvpath/bourne-10-action_HighFlow_LowCompression.csv $csvpath/bourne-10-action_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD High Flow and Medium Compression" $csvpath/bourne-10-action_HighFlow_MediumCompression.csv $csvpath/bourne-10-action_flow.csv $3 >> $appendix
python pred-res-grapher.py "Bourne DVD High Flow and High Compression" $csvpath/bourne-10-action_HighFlow_HighCompression.csv $csvpath/bourne-10-action_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 Low Flow and Low Compression" $csvpath/captain-720_LowFlow_LowCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 Low Flow and Medium Compression" $csvpath/captain-720_LowFlow_MediumCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 Low Flow and High Compression" $csvpath/captain-720_LowFlow_HighCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 Median Flow and Low Compression" $csvpath/captain-720_MedianFlow_LowCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 Median Flow and Medium Compression" $csvpath/captain-720_MedianFlow_MediumCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 Median Flow and High Compression" $csvpath/captain-720_MedianFlow_HighCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 High Flow and Low Compression" $csvpath/captain-720_HighFlow_LowCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 High Flow and Medium Compression" $csvpath/captain-720_HighFlow_MediumCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 720 High Flow and High Compression" $csvpath/captain-720_HighFlow_HighCompression.csv $csvpath/captain-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 Low Flow and Low Compression" $csvpath/limitless-720_LowFlow_LowCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 Low Flow and Medium Compression" $csvpath/limitless-720_LowFlow_MediumCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 Low Flow and High Compression" $csvpath/limitless-720_LowFlow_HighCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 Median Flow and Low Compression" $csvpath/limitless-720_MedianFlow_LowCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 Median Flow and Medium Compression" $csvpath/limitless-720_MedianFlow_MediumCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 Median Flow and High Compression" $csvpath/limitless-720_MedianFlow_HighCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 High Flow and Low Compression" $csvpath/limitless-720_HighFlow_LowCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 High Flow and Medium Compression" $csvpath/limitless-720_HighFlow_MediumCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 720 High Flow and High Compression" $csvpath/limitless-720_HighFlow_HighCompression.csv $csvpath/limitless-720_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 Low Flow and Low Compression" $csvpath/captain-1080_LowFlow_LowCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 Low Flow and Medium Compression" $csvpath/captain-1080_LowFlow_MediumCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 Low Flow and High Compression" $csvpath/captain-1080_LowFlow_HighCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 Median Flow and Low Compression" $csvpath/captain-1080_MedianFlow_LowCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 Median Flow and Medium Compression" $csvpath/captain-1080_MedianFlow_MediumCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 Median Flow and High Compression" $csvpath/captain-1080_MedianFlow_HighCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 High Flow and Low Compression" $csvpath/captain-1080_HighFlow_LowCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 High Flow and Medium Compression" $csvpath/captain-1080_HighFlow_MediumCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Captain America 1080 High Flow and High Compression" $csvpath/captain-1080_HighFlow_HighCompression.csv $csvpath/captain-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 Low Flow and Low Compression" $csvpath/limitless-1080_LowFlow_LowCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 Low Flow and Medium Compression" $csvpath/limitless-1080_LowFlow_MediumCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 Low Flow and High Compression" $csvpath/limitless-1080_LowFlow_HighCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 Median Flow and Low Compression" $csvpath/limitless-1080_MedianFlow_LowCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 Median Flow and Medium Compression" $csvpath/limitless-1080_MedianFlow_MediumCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 Median Flow and High Compression" $csvpath/limitless-1080_MedianFlow_HighCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 High Flow and Low Compression" $csvpath/limitless-1080_HighFlow_LowCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 High Flow and Medium Compression" $csvpath/limitless-1080_HighFlow_MediumCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Limitless 1080 High Flow and High Compression" $csvpath/limitless-1080_HighFlow_HighCompression.csv $csvpath/limitless-1080_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k Low Flow and Low Compression" $csvpath/Elysium-4k_LowFlow_LowCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k Low Flow and Medium Compression" $csvpath/Elysium-4k_LowFlow_MediumCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k Low Flow and High Compression" $csvpath/Elysium-4k_LowFlow_HighCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k Median Flow and Low Compression" $csvpath/Elysium-4k_MedianFlow_LowCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k Median Flow and Medium Compression" $csvpath/Elysium-4k_MedianFlow_MediumCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k Median Flow and High Compression" $csvpath/Elysium-4k_MedianFlow_HighCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k High Flow and Low Compression" $csvpath/Elysium-4k_HighFlow_LowCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k High Flow and Medium Compression" $csvpath/Elysium-4k_HighFlow_HediumCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
python pred-res-grapher.py "Elysium 4k High Flow and High Compression" $csvpath/Elysium-4k_HighFlow_HighCompression.csv $csvpath/Elysium-4k_flow.csv $3 >> $appendix
