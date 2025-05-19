'''
@author: J.-L. Schaefer<joana-lysiane DOT schaefer AT fu-berlin DOT de>

Collection to calculate the correlation matrix element 
for a range of lagtimes and different MD simulations.

Dependencies:
    MDTraj
    NUMPY
    matplptlib
    deeptime

Input:
    MD trajectories for different thermostats 
Output:
    files and visualisations for 
    mean correlation vs lagtime for all thermostats 
'''
# Imports 
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.markov import TransitionCountEstimator

# MD data 
input=''
# analysis output
output=''

# MD data name
filename_vr = 'vr_1' 
filename_nh = 'nh' 
filename_be = 'be' 
filename_sd_001 = 'sd_001'
filename_sd_01 = 'sd_01' 
filename_sd_1 = 'sd_1'
filename_sd_10 = 'sd_10'
filename_sd_100 = 'sd_100' 

Lables = [filename_vr, 
          filename_nh, 
          filename_be, 
          filename_sd_001, 
          filename_sd_01, 
          filename_sd_1, 
          filename_sd_10, 
          filename_sd_100]

# lagtime range 
lags=np.array(range(1,600, 10))
# MD output frequence 
nstxout=100

###############################################
# functions 
###############################################
def get_correlation(dtraj,
                    lagtimes,
                    countmode,
                    reversible,
                    stationary_distribution_constraint,
                    ):
    """ To get correlation matrix element C_{5,5} of a 
    deeptime MSM.
    """
    # instantiate the MaximumLikelihoodMSM estimator for MSM 
    estimator = MaximumLikelihoodMSM(reversible=reversible,
                                     stationary_distribution_constraint=stationary_distribution_constraint)
    # create MSM models to analyse for timescales
    models = []
    for lagtime in lagtimes:
        # collect statistics from discret trajectory using deeptime’s transition count estimator
        count_estimator = TransitionCountEstimator(lagtime=lagtime, 
                                                   count_mode=countmode)
        # use transition count estimator to fit a count model
        counts = count_estimator.fit(dtraj).fetch_model() 
        correlation = counts.count_matrix[4,4]
        models.append(correlation)
    correlations = np.array(models)
    return correlations   

def torions(traj, number=1):
    ''' to give torsion angle along plane right (4, 5, 7, 8)
    or left (7, 8, 10, 11) of C3 in pentane
    Attention atom index is one lower then in gro or top file
    for pentane in water number = 1
    for pentane in pentane number = 140
    '''
    ind_psi = np.array([np.array([7,8,10,11])+i*17 for i in range(number)])
    ind_phi = np.array([np.array([4, 5, 7, 8])+i*17 for i in range(number)])
    psi=[]
    phi=[]
    for i in range(number):
        psi_degrees = md.compute_dihedrals(traj, [ind_psi[i]])
        phi_degrees = md.compute_dihedrals(traj, [ind_phi[i]])
        phi_degrees = np.degrees(phi_degrees)
        psi_degrees = np.degrees(psi_degrees)
        phi_degrees = np.mod(phi_degrees + 180, 360) - 180
        psi_degrees = np.mod(psi_degrees + 180, 360) - 180
        phi.append(phi_degrees)
        psi.append(psi_degrees)
    phi = np.array(phi).flatten()
    psi = np.array(psi).flatten()
    return phi, psi 

def cluster_UnitGrid_C(x, y):
    '''
    Discretizes the (x, y) coordinates into grid regions based on specific intervals for both dimensions.

    Args:
        x (array-like): 1D or 2D array of x-coordinates, where each value lies in the range [-180, 180].
        y (array-like): 1D or 2D array of y-coordinates, where each value lies in the range [-180, 180].
    
    Returns:
        numpy.ndarray: A 1D array of cluster IDs corresponding to each (x, y) point in the input arrays.
                       The cluster IDs are integers in the range [0, 8] based on the intervals.
    
    Example:
        >>> x = np.array([-100, 50, 170])
        >>> y = np.array([100, -50, -170])
        >>> cluster_UnitGrid(x, y)
        array([4, 3, 8])
    '''
   
    intervals = [-180, -70, 70, 180]
    x_bins = np.digitize(x, bins=intervals) - 1
    y_bins = np.digitize(y, bins=intervals) - 1
    cluster_ids = x_bins * 3 + y_bins  
    return cluster_ids

def plot_correlation(MSM_lagtimes, MSM_ITS_mean, MSM_ITS_std,MD_nout,Colors,Lables):

    fig, ax = plt.subplots(1, 1, figsize=(5.25,4.))
    # MD traj in fs convert to ps
    lagtimes = MSM_lagtimes * MD_nout * 0.001 
    its_data = MSM_ITS_mean
    its_std = MSM_ITS_std 

    for idx, integrator in enumerate(its_data):
        if Lables[idx] == 'BE':
            plt.plot(lagtimes, integrator[:], linewidth=1.5, ls='--',
                        color=Colors[idx], label=Lables[idx])
        else:
            plt.plot(lagtimes, integrator[:], linewidth=1.5, 
                        color=Colors[idx], label=Lables[idx])

        ax.fill_between(lagtimes,integrator[:]-its_std[idx],
                        integrator[:]+its_std[idx], 

                        color=Colors[idx],alpha=.15)

    ax.set_xlim(lagtimes.min(), lagtimes.max())
    
    ax.set_xlabel('Lag time in ps', fontsize=19)
    ax.set_ylabel('Crorrelation', fontsize=19)
    
    ax.tick_params(axis="both",which="major",labelsize=19)
    fig.legend(loc='lower left', bbox_to_anchor=(0.01, 0.00001), borderaxespad=5,fontsize=14, ncol=1)

    plt.tight_layout()

###############################################
# run over all MD simulations 
###############################################
for traj_file in Lables:
    # Load the .xtc file using MDTraj
    gro_file = traj_file + '.gro'
    trajectory = md.load(input+traj_file + '.trr', top=input+gro_file)
    ITS=[]
    for i in range(5):
        residue_indices = trajectory.topology.select(
            'resid '+str(i*28)+' to '+str(27+i*28)
            )
        subset_traj = trajectory.atom_slice(residue_indices)
        print(subset_traj)
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
    np.save(output+traj_file+'_'+'correlations_mean', correlations_mean)
    np.save(output+traj_file+'_'+'correlations_std', correlations_std)

###############################################
# visualisation 
###############################################
corr_vr= np.load(output+filename_vr+'_'+'correlations_mean.npy')
corr_nh= np.load(output+filename_nh+'_'+'correlations_mean.npy')
corr_be= np.load(output+filename_be+'_'+'correlations_mean.npy')
corr_sd_001= np.load(output+filename_sd_001+'_'+'correlations_mean.npy')
corr_sd_01= np.load(output+filename_sd_01+'_'+'correlations_mean.npy')
corr_sd_1= np.load(output+filename_sd_1+'_'+'correlations_mean.npy')
corr_sd_10= np.load(output+filename_sd_10+'_'+'correlations_mean.npy')
corr_sd_100= np.load(output+filename_sd_100+'_'+'correlations_mean.npy')

std_vr= np.load(output+filename_vr+'_'+'correlations_std.npy')
std_nh= np.load(output+filename_nh+'_'+'correlations_std.npy')
std_be= np.load(output+filename_be+'_'+'correlations_std.npy')
std_sd_001= np.load(output+filename_sd_001+'_'+'correlations_std.npy')
std_sd_01= np.load(output+filename_sd_01+'_'+'correlations_std.npy')
std_sd_1= np.load(output+filename_sd_1+'_'+'correlations_std.npy')
std_sd_10= np.load(output+filename_sd_10+'_'+'correlations_std.npy')
std_sd_100= np.load(output+filename_sd_100+'_'+'correlations_std.npy')

vr_color="green"
vr_color='#00509E' 
vr_color_="dashed"
be_color= '#008B8B' 
be_linestyle="dotted"
nh_color= '#003366' 
nh_linestyle="solid"
sd_001_color='#D32F2F' 
sd_01_color='#800080' 
sd_1_color='#00509E' 
sd_10_color='#006400' 
sd_100_color='#FFCC00' 

ITS = [corr_vr, corr_nh, corr_be, corr_sd_001, corr_sd_01, corr_sd_1, corr_sd_10, corr_sd_100]
STD = [std_vr, std_nh, std_be, std_sd_001, std_sd_01, std_sd_1, std_sd_10, std_sd_100]
Lables = [r'VR', r'NH', r'BE', r'$\tau_T=0.01$ ps', r'$\tau_T=0.1$ ps', r'$\tau_T=1$ ps', r'$\tau_T=10$ ps', r'$\tau_T=100$ ps']
Colors = [vr_color, nh_color, be_color, sd_001_color, sd_01_color, sd_1_color, sd_10_color, sd_100_color]

plot_correlation(lags, ITS, STD, nstxout, Colors, Lables)

