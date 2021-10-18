# sst-drift
This repository contains the core routine to fit a model of SST temporal evolution to SST observations from surface drifters along their trajectories as described in Elipot et al. 2021 (to be submitted).

The test code is `main.m` which load some test data in `level0_sst_data_drifter_id55366.mat`, run the estimation with the scripts `lowessatx.m` and `ssteval.m`, and plot some results. Some statistics are calculated thanks to `weightedMedian.n` (Haase 2021).

Sven Haase (2021). weighted median (https://www.mathworks.com/matlabcentral/fileexchange/23077-weighted-median), MATLAB Central File Exchange. Retrieved October 18, 2021.
