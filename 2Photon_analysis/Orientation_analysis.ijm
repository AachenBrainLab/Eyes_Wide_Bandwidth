dir1 = getDirectory("C:/Users/balla/Desktop/MC_bandwidth_dist_analysis");
dir2 = getDirectory("C:/Users/balla/Desktop/MC_bandwidth_dist_analysis");
list = getFileList(dir1);
setBatchMode(true);
for (i=0; i<list.length; i++) {
 showProgress(i+1, list.length);
 open(dir1+list[i]);
run("Directionality", "method=[Local gradient orientation] nbins=180 histogram_start=0 histogram_end=180 display_table");
 saveAs("Results", dir2 + list[i] + ".csv");
 close();
}


