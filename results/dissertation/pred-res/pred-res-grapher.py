#!/usr/bin/python

# Usage: python pred-res-grapher.py <clip-name> <result.csv> <flow.csv>
#
# Will output a separate graph and latex snippet for every result in the file.
# The graph will plot psnr for all 240 frames as a line and then the flow
# for each frame x (flow from x-1 to x) will be plotted as a bar underneath.

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv) != 4:
    sys.exit("Wrong number of arguments. See usage in header comment.")

clipname = sys.argv[1]


#
# plotdata
#
def plotdata(title, config, flow, psnr, filename):
    frames = np.arange(len(flow))
    if len(flow) != len(psnr):
        sys.exit("flow and psnr data not equal length.")
    
    fig, ax1 = plt.subplots()
    #ax1.set_title(title)
    
    ax1.set_xlabel('frames (f)')
    ax1.set_ylabel('flow', color='grey')
    ax1.set_yscale('log')
    ax1.bar(frames, flow, width=1.0, color='lightgrey')
    ax1.tick_params(axis='y', labelcolor='grey')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('PSNR', color='blue')
    ax2.plot(frames, psnr, color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig('images/' + filename)
    #plt.show()

#
# printlatex
#
def printlatex(title, config, filename, newfig, lastfig, lf):
    if newfig:
        print("\\begin{figure}")
        print("\\centering")

    if filename != '':
        print("\\begin{minipage}{.5\\textwidth}")
        print("\\centering")
        print("\\epsfig{file=images/" + filename + ", width=3.0in}")
        print("\\caption{PSNR vs. optical flow of the " + title + "clip using the configuration " + config + ".}")
        print("\\label{fig:chart" + filename.replace('.eps', '') + "}")
        print("\\end{minipage}" + ('' if lf else '%'))
    
    if lastfig:
        print("\\end{figure}")



#
# Parse the flow
#
flowdata = []
flowfile = open(sys.argv[3])
for line in flowfile:
    cells = line.rstrip().split(',')
    if cells[0] == 'FlowCSV':
        flowdata.append(float(cells[3]))
flowfile.close()


#
# Parse the results
#
psnrdata = []
filename = '';
config = '';
start = 0
last = 0
fignum = 0;
resultfile = open(sys.argv[2])
for line in resultfile:
    cells = line.rstrip().split(',')
    if cells[0] == 'GroupCSV':
        # If we have psnrdata, then plot it.
        if len(psnrdata) > 0:
            subflowdata = flowdata[start:last]
            filename = clipname + '_' + config + '.eps';
            filename = filename.replace(' ', '')
            filename = filename.replace(':', '-')
            filename = filename.replace(',', '_')
            plotdata(clipname, config, subflowdata, psnrdata, filename)
            printlatex(clipname, config, filename,
                           True if fignum % 6 == 0 else False,
                           True if fignum % 6 == 5 else False,
                           True if fignum % 2 == 1 else False)
            fignum = fignum + 1

        psnrdata = []
        filename = cells[1]
        config = cells[4] + ',' + cells[6] + ',' + cells[7]
        start = int(cells[2]) - 2
        last = start + int(cells[3])
    elif cells[0] == 'RenderCSV':
        psnrdata.append(float(cells[3]))

if fignum % 6 != 0:
    printlatex('', '', '', False, True, False);
resultfile.close()

