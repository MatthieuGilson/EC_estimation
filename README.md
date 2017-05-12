Estimation of effective connectivity (EC) and input statistics (variances in matrix Sigma) to interpret fMRI/BOLD time series

related to draft:
Gilson M, Deco G, Friston K, Hagmann P, Mantini D, Betti V, Romani GL, Corbetta M. Effective connectivity inferred from fMRI transition dynamics during movie viewing points to a balanced reconfiguration of cortical interactions. http://biorxiv.org/content/early/2017/02/20/110015

model in: Gilson M, Moreno-Bote R, Ponce-Alvarez A, Ritter P, Deco G. Estimation of Directed Effective Connectivity from fMRI Functional Connectivity Hints at Asymmetries of Cortical Connectome. PLoS Comput Biol 2016, 12: e1004762; dx.doi.org/10.1371/journal.pcbi.1004762

The scripts are written in python 2.7 and use numpy, scipy and matplotlib.
Run scripts in order:
1) optimization_movie.py
2) plot_results_movie.py

Here we have 22 subjects with 2 conditions: viewing a black screen or a movie (Ponce-Alvarez A, Deco G, Hagmann P, Romani GL, Mantini D, Corbetta M (2015) Resting-State Temporal Synchronization Networks Emerge from Connectivity Topology and Heterogeneity. PLoS Comput Biol 11(2): e1004100. doi:10.1371/journal.pcbi.1004100). The fMRI time series are stored in rest_movie_ts.npy (indices: subject; condition = black screen / movie with 2 sessions each; ROI index; time in TRs).

File optimization_movie.py: recovers the effective connectivity as well as the input covariance matrix (Sigma) for the BOLD time series and structural connectivity (SC) matrix given by DTI corresponding to Hagmann et al. (PLoS Biol 2008), which is a skeleton for EC; saves results for individual EC and Sigma in 'res_movie' directory (in EC_mod.npy and Sigma_mod.npy files, together with the model covariance matrices FC0 and FCtau that aim to reproduce the empirical counterparts in FC_emp, the corresponding time shift equal to 1 or 2 TR is determined by i_tau); the algorithm to tune Sigma has been improved compared to the paper Gilson et al. 2016.

File plot_results_movie.py: check the consistency over individual EC and the relationship with FC0; check significant changes in EC and Sigma (input variances of network model) between rest and movie.
