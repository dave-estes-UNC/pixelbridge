import sys
import shutil

#.m4v .mov .mkv
def gettir(filename, start, count):
    start = str(start)
    count = str(count)
    infile = open(filename)
    outfile = open(filename + '.new', 'w')
    for line in infile:
        line = line.replace('.mov,', '.mov,'+start+','+count+',') 
        line = line.replace('.m4v,', '.m4v,'+start+','+count+',') 
        line = line.replace('.mkv,', '.mkv,'+start+','+count+',') 
        outfile.write(line)
    outfile.close()
    infile.close()
    shutil.move(filename + '.new', filename)

gettir('Elysium-4k_HighFlow_HediumCompression.csv', 2647, 240)
gettir('Elysium-4k_HighFlow_HighCompression.csv', 2647, 240)
gettir('Elysium-4k_HighFlow_LowCompression.csv', 2647, 240)

gettir('Elysium-4k_LowFlow_HighCompression.csv', 2930, 240)
gettir('Elysium-4k_LowFlow_LowCompression.csv', 2930, 240)
gettir('Elysium-4k_LowFlow_MediumCompression.csv', 2930, 240)

gettir('Elysium-4k_MedianFlow_HighCompression.csv', 2495, 240)
gettir('Elysium-4k_MedianFlow_LowCompression.csv', 2495, 240)
gettir('Elysium-4k_MedianFlow_MediumCompression.csv', 2495, 240)

gettir('bourne-10-action_HighFlow_HighCompression.csv', 785, 240)
gettir('bourne-10-action_HighFlow_LowCompression.csv', 785, 240)
gettir('bourne-10-action_HighFlow_MediumCompression.csv', 785, 240)

gettir('bourne-10-action_MedianFlow_HighCompression.csv', 5258, 240)
gettir('bourne-10-action_MedianFlow_LowCompression.csv', 5258, 240)
gettir('bourne-10-action_MedianFlow_MediumCompression.csv', 5258, 240)

gettir('bourne-13-dialog_LowFlow_HighCompression.csv', 5131, 240)
gettir('bourne-13-dialog_LowFlow_LowCompression.csv', 5131, 240)
gettir('bourne-13-dialog_LowFlow_MediumCompression.csv', 5131, 240)

gettir('captain-1080_HighFlow_HighCompression.csv', 2936, 240)
gettir('captain-1080_HighFlow_LowCompression.csv', 2936, 240)
gettir('captain-1080_HighFlow_MediumCompression.csv', 2936, 240)

gettir('captain-1080_LowFlow_HighCompression.csv', 3, 240)
gettir('captain-1080_LowFlow_LowCompression.csv', 3, 240)
gettir('captain-1080_LowFlow_MediumCompression.csv', 3, 240)

gettir('captain-1080_MedianFlow_HighCompression.csv', 2749, 240)
gettir('captain-1080_MedianFlow_LowCompression.csv', 2749, 240)
gettir('captain-1080_MedianFlow_MediumCompression.csv', 2749, 240)

gettir('captain-720_HighFlow_HighCompression.csv', 3075, 240)
gettir('captain-720_HighFlow_LowCompression.csv', 3075, 240)
gettir('captain-720_HighFlow_MediumCompression.csv', 3075, 240)

gettir('captain-720_LowFlow_HighCompression.csv', 2, 240)
gettir('captain-720_LowFlow_LowCompression.csv', 2, 240)
gettir('captain-720_LowFlow_MediumCompression.csv', 2, 240)

gettir('captain-720_MedianFlow_HighCompression.csv', 1178, 240)
gettir('captain-720_MedianFlow_LowCompression.csv', 1178, 240)
gettir('captain-720_MedianFlow_MediumCompression.csv', 1178, 240)

gettir('limitless-1080_HighFlow_HighCompression.csv', 1400, 240)
gettir('limitless-1080_HighFlow_LowCompression.csv', 1400, 240)
gettir('limitless-1080_HighFlow_MediumCompression.csv', 1400, 240)

gettir('limitless-1080_LowFlow_HighCompression.csv', 1013, 240)
gettir('limitless-1080_LowFlow_LowCompression.csv', 1013, 240)
gettir('limitless-1080_LowFlow_MediumCompression.csv', 1013, 240)

gettir('limitless-1080_MedianFlow_HighCompression.csv', 352, 240)
gettir('limitless-1080_MedianFlow_LowCompression.csv', 352, 240)
gettir('limitless-1080_MedianFlow_MediumCompression.csv', 352, 240)

gettir('limitless-720_HighFlow_HighCompression.csv', 25, 240)
gettir('limitless-720_HighFlow_LowCompression.csv', 25, 240)
gettir('limitless-720_HighFlow_MediumCompression.csv', 25, 240)

gettir('limitless-720_LowFlow_HighCompression.csv', 1109, 240)
gettir('limitless-720_LowFlow_LowCompression.csv', 1109, 240)
gettir('limitless-720_LowFlow_MediumCompression.csv', 1109, 240)

gettir('limitless-720_MedianFlow_HighCompression.csv', 358, 240)
gettir('limitless-720_MedianFlow_LowCompression.csv', 358, 240)
gettir('limitless-720_MedianFlow_MediumCompression.csv', 358, 240)
