For the dissertation, we added a separate chart for every results: 1080 in total.
pred-res-grapher.ph parsed the results and flow and produced graphs and latex snippets.
The result and csv files are in the results for ACMMM2015.

pred-res-grapher.ph is driven by build-appendix.sh, which will create
all of the graphs and generate the latex for the prediction residual
appendix.

From this directory, you can supply the path to the csv file, the
appendix filename and the path to the dissertation images. E.g.

    $ ./build-appendix.sh ../../ACMMM2015/pred-res ../../../../disseration/appendices/pred_res.tex ../../../../disseration/images
