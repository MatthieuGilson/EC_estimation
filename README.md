Estimation of effective connectivity (EC) and input statistics (variances in matrix Sigma) to interpret fMRI/BOLD time series

related to paper: Gilson M, Moreno-Bote R, Ponce-Alvarez A, Ritter P, Deco G. Estimation of Directed Effective Connectivity from fMRI Functional Connectivity Hints at Asymmetries of Cortical Connectome. PLoS Comput Biol 2016, 12: e1004762; dx.doi.org/10.1371/journal.pcbi.1004762

File optimization_movie.py: recovers the connectivity as well as the input covariance matrix for the fMRI covariance matrices in FC_emp.npy (indices: subject, condition = black screen / movie, time shift, ROI, ROI) and structural connectivity (SC) matrix given by DTI corresponding to Hagmann et al. (PLoS Biol 2008), which is a skeleton for EC; saves results for individual EC and Sigma in 'res_movie' directory (in EC_mod.npy and Sigma_mod.npy files, together with the covariance matrices FC0 and FCtau with no time shift and a time shift equal to 1 or 2, cf. i_tau; also saves the mean over all subjects); the algorithm to tune Sigma has been improved compared to the paper Gilson et al.

File plot_results_movie.py: check the consistency over individual EC and relationship with FC0
