import os, sys
import numpy as np
import scipy.linalg as spl
import scipy.stats as stt
import matplotlib.pyplot as pp
import pickle as pkl

os.system('clear')

obj_dir = 'obj_movie/'
res_dir = 'res_movie/'
if not os.path.exists(res_dir):
	print 'create directory:', res_dir
	os.makedirs(res_dir)

print obj_dir, '/', res_dir

graph_format = 'png'

# network parameters
n_sub = 22
n_cond = 2

N = 66 # size of network

v_tau = np.arange(2,dtype=float)
n_tau = v_tau.size

w_C = 1. # weight reference to set max in optimization (to increase if too many saturated estimated weights)
noise_level = 0.5

# optimzation steps and rate (to avoid explosion of activity, diagonal tuned quickly)
n_opt = 2000
epsilon_EC = 0.0005
epsilon_Sigma = 0.05

i_tau = 1 # time shift for optimization (in TR; here can be 1 or 2)

min_val_EC = 0. # maximal value for EC
max_val_EC = 1. # maximal value for EC
min_val_Sigma = 0. # minimal value for Sigma

lim_SC = 0. # limit DTI value to determine SC (only connections with larger values are tuned)


# load empirical data
SC_anat = np.load(obj_dir+'SC_anat.npy')
FC_emp = np.load(obj_dir+'FC_emp.npy')
tau_x = np.load(obj_dir+'tau_x.npy')


# mask for existing connections for EC
mask_EC = np.zeros([N,N],dtype=bool) # EC weights to tune
mask_EC[SC_anat>lim_SC] = True
for i in range(N):
	mask_EC[i,i] = False # no self connection
	mask_EC[i,N-1-i] = False # additional interhemispherical connections
print 'EC density:', mask_EC.sum()/float(N*(N-1))

# diagonal mask for input noise matrix (so far, tune noise matrix for diagonal elements only)
mask_Sigma = np.eye(N,dtype=bool)


np.save(res_dir+'mask_EC.npy',mask_EC)
np.save(res_dir+'mask_Sigma.npy',mask_Sigma)


# optimization
print '*opt*'
tau = v_tau[i_tau]
print 'tau:', tau

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
		EC[np.logical_not(mask_EC)] = 0
		Sigma = np.eye(N)  # initial noise * (0.5+0.5*np.random.rand())

		# record best fit
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

			# update EC
			EC[mask_EC] += epsilon_EC * Delta_J[mask_EC]
			EC[mask_EC] = np.clip(EC[mask_EC],min_val_EC,max_val_EC)

			# update noise
			Delta_Sigma = -np.dot(J,Delta_FC0)-np.dot(Delta_FC0,J.T)
			Sigma[mask_Sigma] += epsilon_Sigma * Delta_Sigma[mask_Sigma]
			Sigma[mask_Sigma] = np.maximum(Sigma[mask_Sigma],min_val_Sigma)

			# check for stop
			if i_opt<n_opt-1 and not stop_opt:
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
np.save(res_dir+'EC_mod.npy',EC_mod)
np.save(res_dir+'Sigma_mod.npy',Sigma_mod)
np.save(res_dir+'FC0_mod.npy',FC0_mod)
np.save(res_dir+'FCtau_mod.npy',FCtau_mod)

np.save(res_dir+'mean_FC0_blackscreen.npy',FC_emp[:,0,0,:,:].mean(0))
np.save(res_dir+'mean_FC1_blackscreen.npy',FC_emp[:,0,1,:,:].mean(0))
np.save(res_dir+'mean_EC_blackscreen.npy',EC_mod[:,0,:,:].mean(0))
np.save(res_dir+'mean_Sigma_blackscreen.npy',Sigma_mod[:,0,:,:].mean(0))
np.save(res_dir+'mean_FC0_movie.npy',FC_emp[:,1,0,:,:].mean(0))
np.save(res_dir+'mean_FC1_movie.npy',FC_emp[:,1,1,:,:].mean(0))
np.save(res_dir+'mean_EC_movie.npy',EC_mod[:,1,:,:].mean(0))
np.save(res_dir+'mean_Sigma_movie.npy',Sigma_mod[:,1,:,:].mean(0))
