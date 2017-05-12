import os, sys
import numpy as np
import scipy.linalg as spl
import scipy.stats as stt
import matplotlib.pyplot as pp
import pickle as pkl

os.system('clear')

obj_dir = './'
res_dir = 'res_movie/'
if not os.path.exists(res_dir):
	print 'create directory:', res_dir
	os.makedirs(res_dir)

print obj_dir, ';', res_dir

graph_format = 'png'

##################
# fMRI time series
ts_emp = np.load('rest_movie_ts.npy')

n_sub = 22 # number of subjects
n_cond = 2 # rest and movie conditions
n_session = 2 # 2 sessions per condition
N = 66 # number of ROIs
T = 300 # number of TRs of the recording

# time shifts for FC: 0, 1 and 2 TR
v_tau = np.arange(3,dtype=float)
n_tau = v_tau.size

# FC = spatiotemporal covariance of BOLD signals (average of the 2 sessions)
FC_emp = np.zeros([n_sub,n_cond,n_tau,N,N])
for i_sub in range(n_sub):
	for i_run in range(n_cond*n_session):
		ts_emp[i_sub,i_run,:,:] -= np.outer(ts_emp[i_sub,i_run,:,:].mean(1),np.ones(T)) # center the time series
for i_sub in range(n_sub):
	for i_cond in range(n_cond):
		for i_session in range(n_session):
			for i_tau in range(n_tau):
				FC_emp[i_sub,i_cond,i_tau,:,:] += np.tensordot(ts_emp[i_sub,i_cond*2+i_session,:,0:T-n_tau+1],ts_emp[i_sub,i_cond*2+i_session,:,i_tau:T-n_tau+1+i_tau],axes=(1,1)) / float((T-n_tau)*n_session)

FC_emp *= 0.5/FC_emp[:,0,0,:,:].mean()
print 'max FC value (most of the distribution should be between 0 and 1):', FC_emp.mean()

# time constant for BOLD autocovariances
slopes = np.zeros([n_sub,n_cond,N])
for i_sub in range(n_sub):
	for i_cond in range(n_cond):
		for i in range(N):
			ac_tmp = np.maximum(FC_emp[i_sub,i_cond,:,i,i],1e-10) # autocovariance for time shifts in v_tau; with lower bound to avoid negative values (cf. log)
			slopes[i_sub,i_cond,i] = np.polyfit(v_tau,np.log(ac_tmp),1)[0] # slope of autocovariance for ROI i

tau_x = -1./slopes.mean(2) # inverse of negative slope of autocovariance


#################
# structural data
SC_anat = np.load(obj_dir+'SC_anat.npy')

lim_SC = 0. # limit DTI value to determine SC (only connections with larger values are tuned)

# mask for existing connections for EC
mask_EC = np.zeros([N,N],dtype=bool) # EC weights to tune
mask_EC[SC_anat>lim_SC] = True
for i in range(N):
	mask_EC[i,i] = False # no self connection
	mask_EC[i,N-1-i] = False # additional interhemispheric connections
print 'EC density:', mask_EC.sum()/float(N*(N-1))

# diagonal mask for input noise matrix (here, no input cross-correlation)
mask_Sigma = np.eye(N,dtype=bool)


##############
# optimization
w_C = 1. # weight reference to set max in optimization (to increase if too many saturated estimated weights)

# optimzation rates (to avoid explosion of activity, Sigma is tuned quicker)
epsilon_EC = 0.0005
epsilon_Sigma = 0.05

min_val_EC = 0. # minimal value for tuned EC elements
max_val_EC = 1. # maximal value for tuned EC elements
min_val_Sigma = 0. # minimal value for tuned Sigma elements


i_tau = 1 # time shift for optimization (in TR; can be 1 or 2)
tau = v_tau[i_tau]
print 'opt with time shift', tau, 'TR'

EC_mod = np.zeros([n_sub,n_cond,N,N])
Sigma_mod = np.zeros([n_sub,n_cond,N,N])
FC0_mod = np.zeros([n_sub,n_cond,N,N])
FCtau_mod = np.zeros([n_sub,n_cond,N,N])

for i_sub in range(n_sub):
	for i_cond in range(n_cond):
		print
		print 'sub', i_sub, '; cond', i_cond

		# initial EC
		EC = np.zeros([N,N]) # initial connectivity
		Sigma = np.eye(N)  # initial noise

		# record best fit (matrix distance between model and empirical FC)
		best_dist = 1e10

		# objective FC matrices (empirical)
		FC0_obj = FC_emp[i_sub,i_cond,0,:,:]
		FCtau_obj = FC_emp[i_sub,i_cond,i_tau,:,:]

		stop_opt = False
		i_opt = 0
		while not stop_opt:

			# calculate Jacobian of dynamical system
			J = -np.eye(N)/tau_x[i_sub,:].mean() + EC

			# calculate FC0 and FCtau for model
			FC0 = spl.solve_lyapunov(J,-Sigma)
			FCtau = np.dot(FC0,spl.expm(J.T*tau))

			# matrices of model error
			Delta_FC0 = FC0_obj-FC0
			Delta_FCtau = FCtau_obj-FCtau

			# calculate error between model and empirical data for FC0 and FC_tau (matrix distance)
			dist_FC_tmp = 0.5*(np.sqrt((Delta_FC0**2).sum()/(FC0_obj**2).sum())+np.sqrt((Delta_FCtau**2).sum()/(FCtau_obj**2).sum()))

			# calculate Pearson correlation between model and empirical data for FC0 and FC_tau
			Pearson_FC_tmp = 0.5*(stt.pearsonr(FC0.reshape(-1),FC0_obj.reshape(-1))[0]+stt.pearsonr(FCtau.reshape(-1),FCtau_obj.reshape(-1))[0])

			# record best model parameters
			if dist_FC_tmp<best_dist:
				best_dist = dist_FC_tmp
				best_Pearson = Pearson_FC_tmp
				i_best = i_opt
				EC_mod_tmp = np.array(EC)
				Sigma_mod_tmp = np.array(Sigma)
				FC0_mod_tmp = np.array(FC0)
				FCtau_mod_tmp = np.array(FCtau)
			else:
				stop_opt = i_opt>100

			# Jacobian update
			Delta_J = np.dot(np.linalg.pinv(FC0),Delta_FC0+np.dot(Delta_FCtau,spl.expm(-J.T*tau))).T/tau

			# update EC (recurrent connectivity)
			EC[mask_EC] += epsilon_EC * Delta_J[mask_EC]
			EC[mask_EC] = np.clip(EC[mask_EC],min_val_EC,max_val_EC)

			# update Sigma (input variances)
			Delta_Sigma = -np.dot(J,Delta_FC0)-np.dot(Delta_FC0,J.T)
			Sigma[mask_Sigma] += epsilon_Sigma * Delta_Sigma[mask_Sigma]
			Sigma[mask_Sigma] = np.maximum(Sigma[mask_Sigma],min_val_Sigma)

			# check for stop
			if not stop_opt:
				if (i_opt)%50==0:
					print 'opt step:', i_opt
					print 'dist FC:', dist_FC_tmp, '; Pearson FC:', Pearson_FC_tmp
				i_opt += 1
			else:
				print 'stop at step', i_opt, 'with best FC dist:', best_dist, '; best FC Pearson:', best_Pearson

		EC_mod[i_sub,i_cond,:,:] = EC_mod_tmp
		Sigma_mod[i_sub,i_cond,:,:] = Sigma_mod_tmp
		FC0_mod[i_sub,i_cond,:,:] = FC0_mod_tmp
		FCtau_mod[i_sub,i_cond,:,:] = FCtau_mod_tmp


# save results
np.save(res_dir+'FC_emp.npy',FC_emp) # empirical spatiotemporal FC
np.save(res_dir+'mask_EC.npy',mask_EC) # mask of optimized connections
np.save(res_dir+'mask_Sigma.npy',mask_Sigma) # mask of optimized Sigma elements

np.save(res_dir+'EC_mod.npy',EC_mod) # estimated EC matrices
np.save(res_dir+'Sigma_mod.npy',Sigma_mod) # estimated Sigma matrices
np.save(res_dir+'FC0_mod.npy',FC0_mod) # model FC0
np.save(res_dir+'FCtau_mod.npy',FCtau_mod) # model FCtau (tau = 1 or 2 TR)

# various compiled results (if needed)
if False:
	np.save(res_dir+'mean_FC0_rest.npy',FC_emp[:,0,0,:,:].mean(0))
	np.save(res_dir+'mean_FC1_rest.npy',FC_emp[:,0,1,:,:].mean(0))
	np.save(res_dir+'mean_EC_rest.npy',EC_mod[:,0,:,:].mean(0))
	np.save(res_dir+'mean_Sigma_rest.npy',Sigma_mod[:,0,:,:].mean(0))
	np.save(res_dir+'mean_FC0_movie.npy',FC_emp[:,1,0,:,:].mean(0))
	np.save(res_dir+'mean_FC1_movie.npy',FC_emp[:,1,1,:,:].mean(0))
	np.save(res_dir+'mean_EC_movie.npy',EC_mod[:,1,:,:].mean(0))
	np.save(res_dir+'mean_Sigma_movie.npy',Sigma_mod[:,1,:,:].mean(0))
