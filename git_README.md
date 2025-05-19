# README 

## correlation functions used in Markov state models 
To calculate the state-to-state correlation functions used in Markov state models (MSMs) use ```MSM_correlation``` library. The file contains a collection to calculate the correlation matrix element for a range of lagtimes and different molecular dynamics (MD) simulations. 
You need to install a python environment with the packages ```MDTraj```, ```NUMPY```, ```matplptlib``` and ```deeptime```. As input to the functions provide the MD trajectories form different thermostats. You will get files and visualisations for mean correlations vs lagtime for all thermostats. 

### You need to:
- define MD data and analysis output directories 
```python
# MD data 
input=''
# analysis output
output=''
```
- parametrise the model for lag times and write out frequency of the MD simulation 
```python
# lagtime range 
lags=np.array(range(1,600, 10))
# MD output frequence 
nstxout=100
```
- initialise functions 
- load MD data (use ```mdtraj```)
- project MD data into the torsion angle space (use  ```torions```)
- discretise the MD trajectories via ```cluster_UnitGrid_C```
- create the MSM and write out the state-to-state correlation functions (use ```get_correlation```)
- evaluate the mean over a number of runs
```python
trajectory = md.load(input+traj_file + '.trr', top=input+gro_file)
    ITS=[]
for i in range(5):
    residue_indices = trajectory.topology.select(
        'resid '+str(i*28)+' to '+str(27+i*28)
        )
    subset_traj = trajectory.atom_slice(residue_indices)
    phi, psi = torions(subset_traj, number=28)
    CV_traj = np.array([phi, psi]).swapaxes(1,0)
    dtrajs = cluster_UnitGrid_C(phi, psi) 
    correlations=get_correlation(
            dtraj=dtrajs.astype(int), lagtimes=lags,
            countmode='sliding', reversible=True, 
            stationary_distribution_constraint=None)
    ITS.append(correlations)
ITS=np.array(ITS)
ITS = ITS/ITS.max()
correlations_mean= np.mean(ITS, axis=0)
correlations_std = np.std(ITS, axis=0)
```
- plot the correlation values as a function of the lag time (use ```plot_correlation```)
```python
ITS = [corr_vr, corr_nh, corr_be, corr_sd_001, corr_sd_01, corr_sd_1, corr_sd_10, corr_sd_100]
STD = [std_vr, std_nh, std_be, std_sd_001, std_sd_01, std_sd_1, std_sd_10, std_sd_100]
Lables = [r'VR', r'NH', r'BE', r'$\tau_T=0.01$ ps', r'$\tau_T=0.1$ ps', r'$\tau_T=1$ ps', r'$\tau_T=10$ ps', r'$\tau_T=100$ ps']
Colors = [vr_color, nh_color, be_color, sd_001_color, sd_01_color, sd_1_color, sd_10_color, sd_100_color]
plot_correlation(lags, ITS, STD, nstxout, Colors, Lables)
```

You can find detailed infos about the input and output of all functions in the documentation string of each method. README 

## building a Markov state models for pentane 
To create a Markov state models (MSMs) use ```MSM``` library. The file contains a collection to create and visualize a deeptime MSM for pentane. You need to install a python environment with the packages ```MDTraj```, ```NUMPY```, ```matplptlib``` and ```deeptime```. As input to the functions provide the MD trajectories for different thermostats. You will get files and visualisations for implied timescales and eigenvectors. 

### You need to:
- define MD data and analysis output directories 
```python
# MD data 
input='' # input directory
# analysis output
output=''
```
- initialise functions 
- load MD data (use ```mdtraj```)
- parametrise the model for lag time, number of eigenvector and write out frequency of the MD simulation 
- project MD data into the torsion angle space (use  ```torions```)
- discretise the MD trajectories via ```cluster_UnitGrid_C```
- create the MSM and get the eigenvectors for a number of processes, e.g. ```k=0``` via ```get_eigenvectors```
- visualize them with ```plot_eigenvectors2D```
```python
input='' # input directory
traj_file = filename_sd_1 # change for integrator
gro_file = traj_file + '.gro'
trajectory = md.load(input+traj_file + '.trr', top=input+gro_file)
phi, psi = torions(trajectory, number=140)
lag = 50
nstxout=100
CV_traj = np.array([psi, phi]).swapaxes(1,0)
dtrajs = cluster_UnitGrid_C(phi, psi)
plot_assignment_uniform(CV_traj, dtrajs, scip=10)

evecs,lcs=get_eigenvectors(
    dtraj=dtrajs.astype(int), lagtime=lag, reversible=True,
    stationary_distribution_constraint=None, 
    countmode='sliding', c)

plot_eigenvectors2D(evecs, CV_traj[:,0],CV_traj[:,1], dtrajs, 0, np.absolute(evecs[0]).max())
plot_eigenvectors2D(evecs, CV_traj[:,0],CV_traj[:,1], dtrajs, 1, np.absolute(evecs[1:]).max())
plot_eigenvectors2D(evecs, CV_traj[:,0],CV_traj[:,1], dtrajs, 2, np.absolute(evecs[1:]).max())
```
- evaluate the implied timescales for a list of lag times and thermostats (ue ```get_implied_timescales```)
- plot the implied timescales for different processes ```k=0``` (use ````plot_ITS_bar_broken````)
```python
lags=np.array(range(1,600, 10))
ITS = [its_nh, its_be, its_vr,  np.zeros_like(its_be),its_sd_001, its_sd_01,its_sd_1, its_sd_10, its_sd_100 ]
STD = [std_nh, std_be, std_vr,  np.zeros_like(its_be),std_sd_001, std_sd_01, std_sd_1, std_sd_10, std_sd_100 ]
Lables = [r'NH',  r'BE', r'VR',  ' ', r'$\tau_T=0.01$ ps',r'$\tau_T=0.1$ ps',r'$\tau_T=1$ ps', r'$\tau_T=10$ ps', r'$\tau_T=100$ ps']
Colors = [nh_color, be_color, vr_color, 'white', sd_001_color, sd_01_color, sd_1_color, sd_10_color, sd_100_color]
ITS_001 = [its_sd_001]
STD_001 = [std_sd_001]
Lables_001 = [r'$\tau_T=0.01$ ps']
Colors_001 = [sd_001_color]
plot_ITS_bar_broken(100,
            ITS, STD,
            Colors, Lables,
            ITS_001, STD_001,
            Colors_001, Lables_001, 
            20,k=0)
```
You can find detailed infos about the input and output of all functions in the documentation string of each method. 