#!/usr/bin/python

# Usage: python pred-res-grapher.py <result.csv> <flow.csv>
#
# Will output a separate graph and latex snippet for every result in the file.
# The graph will plot psnr for all 240 frames as a line and then the flow
# for each frame x (flow from x-1 to x) will be plotted as a bar underneath.

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv) != 3:
    sys.exit("Wrong number of arguments. See usage in header comment.")

#
# plotdata
#
def plotdata(name, config, flow, psnr):
    frames = np.arange(len(flow))
    if len(flow) != len(psnr):
        sys.exit("flow and psnr data not equal length.")
    
    fig, ax1 = plt.subplots()
    
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
    plt.savefig("out.eps")
    plt.show()


flowdata = []
#
# Parse the flow
#
flowfile = open(sys.argv[2])
for line in flowfile:
    cells = line.rstrip().split(',')
    if cells[0] == 'FlowCSV':
        flowdata.append(float(cells[3]))
flowfile.close()
#print(flowdata)
#print("Flow Parsing Complete. " + str(len(flowdata)) + " frames.")

#
# Parse the results
#
psnrdata = []
name = '';
config = '';
start = 0
last = 0
resultfile = open(sys.argv[1])
for line in resultfile:
    cells = line.rstrip().split(',')
    if cells[0] == 'GroupCSV':
        # If we have psnrdata, then plot it.
        if len(psnrdata) > 0:
            subflowdata = flowdata[start:last]
            plotdata(name, config, subflowdata, psnrdata)
            exit(1)
        psnrdata = []
        name = cells[1]
        config = cells[4] + ',' + cells[6] + ',' + cells[7]
        start = int(cells[2]) - 2
        last = start + int(cells[3])
    elif cells[0] == 'RenderCSV':
        psnrdata.append(float(cells[3]))

resultfile.close()

