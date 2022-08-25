# sst-drift
This repository contains the core routine to fit a model of SST temporal evolution to SST observations from surface drifters along their trajectories as described in Elipot et al. 2022 (revised for Scientific Data). Preprint available on arXiv with doi: [10.48550/arXiv.2201.08289](
https://doi.org/10.48550/arXiv.2201.08289)

The test code is `main.m` which loads some test data in `level0_sst_data_drifter_id55366.mat`, runs the estimation with the scripts `lowessatx.m` and `ssteval.m`, and plots some results. Some statistics are calculated thanks to `weightedMedian.m` (Haase 2021).

Sven Haase (2021). weighted median (https://www.mathworks.com/matlabcentral/fileexchange/23077-weighted-median), MATLAB Central File Exchange. Retrieved October 18, 2021.
