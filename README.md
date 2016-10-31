Estimation of effective connectivity (EC) and input statistics (variances in matrix Sigma) to interpret fMRI/BOLD time series

related to paper: Gilson M, Moreno-Bote R, Ponce-Alvarez A, Ritter P, Deco G. Estimation of Directed Effective Connectivity from fMRI Functional Connectivity Hints at Asymmetries of Cortical Connectome. PLoS Comput Biol 2016, 12: e1004762; dx.doi.org/10.1371/journal.pcbi.1004762

The scripts are written in python 2.7 and use numpy, scipy and matplotlib.

Here we have 22 subjects with 2 conditions: viewing a black screen or a movie (Ponce-Alvarez A, Deco G, Hagmann P, Romani GL, Mantini D, Corbetta M (2015) Resting-State Temporal Synchronization Networks Emerge from Connectivity Topology and Heterogeneity. PLoS Comput Biol 11(2): e1004100. doi:10.1371/journal.pcbi.1004100). The fMRI covariance matrices are stored in FC_emp.npy (indices: subject; condition = black screen / movie; time shift = 1 / 2; ROI; ROI). The corresponding time constants for the BOLD autocovariances are in tau_x.npy (indices: subject; ROI).

File optimization_movie.py: recovers the effective connectivity as well as the input covariance matrix for the FC matrices and structural connectivity (SC) matrix given by DTI corresponding to Hagmann et al. (PLoS Biol 2008), which is a skeleton for EC; saves results for individual EC and Sigma in 'res_movie' directory (in EC_mod.npy and Sigma_mod.npy files, together with the covariance matrices FC0 and FCtau with no time shift and a time shift equal to 1 or 2, cf i_tau; also saves the mean over all subjects); the algorithm to tune Sigma has been improved compared to the paper Gilson et al.

File plot_results_movie.py: check the consistency over individual EC and relationship with FC0
