Sub-directory contains MATLAB scripts to process and display Franks data from 
"Bolding KA,Franks KM 2018. Recurrent cortical circuits implement concentration-invariant odor coding. Science. Vol: 361(6407)"

1) Must first download data from CRCNS.org at: https://dx.doi.org/10.6080/K00C4SZB , click on THY1 link and download: processedTHY1.tar.gz AND ExperimentCatalog_THY1.txt
2) Run the script: getPlotTHY1.m, assumes ExperimentCatalog_THY1.txt is in same directory, and all data from processedTHY1.tar.gz is in sub-folder named 'processed'

Two scripts to plot the awake head-fixed freely breathing data from same paper.  
Assumes have data from all 8 mice in 170608.mat, -09.mat, -13.mat, -14.mat, -18.mat, -19.mat, -21.mat, -22.mat 
* Time-varying spiking statistics (PSTH, Var, Cov) in 100ms half overlapping windows, for both Ethyl Butyrate (black) and Hexanoal (green) run: plot_timeStats_DtoF.m 
* Shows increases from spontaneous to evoked (time-averaged) for same stats (PSTH, Var, Cov) and odors, run: plot_ssIncr_AtoC.m . 
Also displays p-values reported in our paper.
