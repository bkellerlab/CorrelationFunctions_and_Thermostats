'''
@author: J.-L. Schaefer <joana-lysiane DOT schaefer AT fu-berlin DOT de>

Collection to create and visualize a deeptime MSM 
of pentane. 

Dependencies:
    MDTraj
    NUMPY
    matplptlib
    deeptime

Input:
    MD trajectories for different thermostats 
Output:
    files and visualisations for 
    implied timescales and eigenvectors
'''
# Imports 
import mdtraj as md
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm._markov_state_model import MarkovStateModel
from deeptime.markov.tools.estimation import largest_connected_set

# MD data 
input='' # input directory
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
# number of evaluated eigenvalues
nvecs = 3
# MD output frequence 
nstxout=100

###############################################
# functions 
###############################################
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
    '''Discretizes the (x, y) coordinates into grid regions based on specific intervals for both dimensions.

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

def plot_assignment_uniform(CV_traj, dtrajs, scip, ngrid=3, x_boundaries = [-180, -70, 70, 180], y_boundaries=[-180, -70, 70, 180]):
    """Plot a 2D scatter plot of the torsion angle trajectory with clustering assignments.

    Args:
        CV_traj (ndarray): A 2D array of shape (n_samples, 2) representing the phi,psi 
        trajectory of over time. 
        dtrajs (array-like): A 1D array of discretized trajectory data, where each entry corresponds
                              to the cluster ID assigned.
        scip (int): The step size (downsampling factor) used to reduce the number of points 
                    visualized for better readability.
        ngrid (int, optional): The grid size, specifying the number of clusters along each axis. 
                               Default is 3 (3x3 grid, resulting in 9 clusters) we can also chose a
                               (2x2 grid, 4 clusters) .

    Returns:
        None: Displays a 2D scatter plot with the trajectory, grid, and cluster assignments. 
    """
    fig, ax = plt.subplots(figsize=(8, 8)) 
    if ngrid==3:
        for x in x_boundaries:
            ax.axvline(x, color='k', linestyle='--', alpha=0.5)
        for y in y_boundaries:
            ax.axhline(y, color='k', linestyle='--', alpha=0.5)
    
    colors = plt.cm.jet(np.linspace(0, 1, ngrid*ngrid))
    
    scatter = ax.scatter(CV_traj[:,0], CV_traj[:,1], c='k', marker='o')
    scatter = ax.scatter(CV_traj[:,0][::scip], CV_traj[:,1][::scip], c=dtrajs[::scip], cmap='jet', marker='o')
    
    ax.set_xlim([-180, 180])
    ax.set_ylim([-180, 180])
    ax.set_xlabel(r'$\phi$',fontsize=16)
    ax.set_ylabel(r'$\psi$',fontsize=16)
    ax.tick_params(axis="both",which="major",labelsize=16)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.ax.tick_params(labelsize=16) 
    cbar.set_label('Cluster ID',fontsize=16)
    plt.grid(True)
    plt.show()

def get_implied_timescales(dtraj,
                            lagtimes,
                            k,
                            countmode,
                            reversible,
                            stationary_distribution_constraint,
                            ):
    """Calculate implied timescales from a discrete trajectory using deeptime.

    Args:
        dtraj (int list): The discrete trajectory, where each element represents a state.
        lagtimes (list or array): Lagtimes at which to compute the implied timescales.
        k (int): The number of MSM eigenvalues.
        countmode (str): The mode to be used by the transition count estimator 
                         (default., 'sliding').
        reversible (bool): If True, the estimator will assume a reversible MSM.
        stationary_distribution_constraint (bool): If True, apply a constraint 
                                                   on the stationary distribution.

    Returns:
        np.ndarray: Implied timescales for each lagtime in `lagtimes`. 

    Example:
        >>> dtraj = np.array([0, 1, 2, 1, 0, 2, 1])  # example discrete trajectory
        >>> lagtimes = [10, 20, 30]  # example lagtimes to compute implied timescales
        >>> k = 5  # number of eigenvalues
        >>> countmode = 'sliding'
        >>> reversible = True
        >>> stationary_distribution_constraint = None
        >>> timescales = get_implied_timescales(dtraj, lagtimes, k, countmode, 
                                                reversible, stationary_distribution_constraint)
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
        # re-estimate a MSM based on the count model
        msm = estimator.fit(counts).fetch_model()
        # create MSM model based on transition matrix
        msm_model = MarkovStateModel(msm.transition_matrix, 
                                     lagtime=lagtime,
                                     n_eigenvalues=k)
        # collect timescales from the MSM model
        models.append(msm_model.timescales())
    its_data = np.array(models)
    return its_data   

def get_eigenvectors(dtraj,
                     lagtime,
                     reversible,
                     stationary_distribution_constraint,
                     countmode,
                     k,
                     right=False
                     ):
    """Compute the eigenvectors of the transition matrix from a discrete trajectory 
    using deeptime.

    Args:
        dtraj (int list): The discrete trajectory, where each element represents a state.
        lagtime (int): Lagtime at which to compute the MSM.
        k (int): The number of MSM eigenvalues.
        countmode (str): The mode to be used by the transition count estimator 
                         (default., 'sliding').
        reversible (bool): If True, the estimator will assume a reversible MSM.
        stationary_distribution_constraint (bool): If True, apply a constraint 
                                                   on the stationary distribution.
        right (bool, optional): If True, returns the **right** eigenvectors. If False, 
                                 returns the **left** eigenvectors. Default is False.

    Returns:
        tuple: A tuple containing:
            - eigenvecs (ndarray): The eigenvectors (left or right) of the MSM's 
                                    transition matrix. Shape is (n_states, k), where 
                                    `n_states` is the number of states and `k` is the 
                                    number of eigenvectors requested.
            - lcs (list): A list of the largest connected sets (LCS) based on the 
                          transition count matrix.

    Example:
        >>> dtraj = np.array([0, 1, 2, 1, 0, 2, 1])  # example discrete trajectory
        >>> lagtime = 10  # example lagtime for eigenvector calculation
        >>> k = 5  # number of eigenvectors to compute
        >>> right = True  # compute right eigenvectors
        >>> eigenvectors, lcs = get_eigenvectors(dtraj, lagtime, reversible=True, 
                                                 stationary_distribution_constraint=False, 
                                                 countmode='symmetric', k=k, right=right)
    """
    # instantiate the MaximumLikelihoodMSM estimator for MSM 
    estimator = MaximumLikelihoodMSM(
            reversible=reversible,
            stationary_distribution_constraint=stationary_distribution_constraint)

    # collect statistics from discret trajectory using deeptime’s transition count estimator
    count_estimator = TransitionCountEstimator(lagtime=lagtime, 
                                               count_mode=countmode)
    # use transition count estimator to fit a count model
    counts = count_estimator.fit(dtraj).fetch_model() 

    # re-estimate a MSM based on the count model
    msm = estimator.fit(counts).fetch_model()
    # create MSM model based on transition matrix
    msm_model = MarkovStateModel(msm.transition_matrix, 
                                 lagtime=lagtime,
                                 n_eigenvalues=k)
    
    if right:
        # collect right eigenvectors from the MSM model
        eigenvecs = msm_model.eigenvectors_right(k)
    else:
        # collect left eigenvectors from the MSM model
        eigenvecs = msm_model.eigenvectors_left(k)
    # store the largest connected sets
    lcs = largest_connected_set(counts.count_matrix)
    return eigenvecs, lcs

def plot_eigenvectors2D(MSM_eigenvector, psi,phi, MSM_traj, n_vec, limit_cb, scip=None):
    '''Plot a 2D scatter plot of eigenvector values over the two torsion coordinates.

    Args:
        MSM_eigenvector (ndarray): The eigenvectors from the MSM.
        psi (ndarray): The CV 1 values.
        phi (ndarray): The CV 2 values.
        MSM_traj (ndarray): The MSM trajectory data, used to build the eigenvector.
        n_vec (int): The index of the eigenvector to plot.
        limit_cb (float): The limit for the colorbar range.
        scip (int, optional): Step size for plotting data points (default is None).

    Returns:
        None: Displays a 2D scatter plot with a colorbar representing the population.
    '''
    if abs(np.sum(np.sign(MSM_eigenvector[n_vec])))==MSM_eigenvector[n_vec].shape[0]:
        cmin = 0
        cmap=mpl.colormaps['Reds']
        eigenvector=np.concatenate([abs(MSM_eigenvector[n_vec])])
    else:
        cmin = -limit_cb
        cmap=mpl.colormaps['RdBu_r'] 
        eigenvector=np.concatenate([MSM_eigenvector[n_vec]])

    f, ax2 = plt.subplots(1, 1, figsize=(7, 5))
    im = ax2.scatter(
        x=psi[::scip],
        y=phi[::scip],
        c=eigenvector[MSM_traj],
        s=10,
        cmap=cmap,
        vmin=cmin, vmax=limit_cb
    )
    ax2.set_facecolor("lightgrey")
    ax2.grid()
    ax2.set_xlabel(r'$\phi$', fontsize=20)
    ax2.set_ylabel(r'$\psi$', fontsize=20)
    ax2.set_xticks([-100,0,100])
    ax2.set_yticks([-100,0,100])
    ax2.tick_params(axis='both', labelsize=18)
    for spine in ax2.spines.values():
        spine.set_visible(False)
    # Create the colorbar
    cbar = plt.colorbar(im, ax=ax2) 
    cbar.ax.tick_params(labelsize=18)   
    plt.show() 

def plot_ITS_bar_broken(MD_nout, MSM_ITS_mean_1, MSM_ITS_std_1, Colors_1, Lables_1,
                        MSM_ITS_mean_2, MSM_ITS_std_2, Colors_2, Lables_2, cut, k=2):
    # Convert MD traj in fs to ps
    scale_factor = MD_nout * 0.001
    its_data_1 = np.array(MSM_ITS_mean_1) * scale_factor
    its_std_1 = np.array(MSM_ITS_std_1) * scale_factor
    its_data_2 = np.array(MSM_ITS_mean_2) * scale_factor
    its_std_2 = np.array(MSM_ITS_std_2) * scale_factor

    # Prepare data
    data = [
        (its_data_1, its_std_1, Colors_1, Lables_1),
        (its_data_2, its_std_2, Colors_2, Lables_2)
    ]
    
    # Setup figure and axes
    fig, (ax, ax1) = plt.subplots(2, 1, figsize=(4, 3.), sharex=True)
    fig.subplots_adjust(hspace=0.0001)
    
    axes = [ax1, ax]
    LAB = [' ', '', '', '', ' ', '', ' ', '', '']

    # Initialize a list for the bar patches to use in the legend
    legend_patches = []

    # Process each dataset (1st and 2nd)
    for axis, (its_data, its_std, Colors, Lables) in zip(axes, data):
        means = [np.mean(its[cut:, k]) for its in its_data]
        std_devs = [np.mean(std[cut:, k]) for std in its_std]
        
        # Plot error bars and bars
        for mean, std_dev, label, color in zip(means, std_devs, Lables, Colors):
            if std_dev==0:
                pass
            else:
                axis.errorbar(label, mean, yerr=std_dev, fmt='-', capsize=5, capthick=2, color='k')
            if label=='VR':
                bar = axis.bar(label, mean, color=color, label=label,hatch='/', width=0.99)
            else:
                bar = axis.bar(label, mean, color=color, label=label, width=0.99)
            if axis==ax1:
                legend_patches.append(Patch(color=color, label=label))


    # Adjust x-axis labels and ticks
    ax1.tick_params(axis='x', rotation=45)
    ax1.set_xticklabels(LAB, ha="right")
    ax1.tick_params(axis="both", which="major", labelsize=19)
    ax.tick_params(axis="both", which="major", labelsize=19)

    # Add vertical axis label
    fig.text(0.0009, 0.5, 'Implied timescale in ps', ha='center', va='center', rotation='vertical', fontsize=19)

    # Remove spines between axes
    ax.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax.tick_params(labelbottom=False)
    ax.tick_params(axis='x', which='both', bottom=False)  # Hide tick labels at the top

    # Add slanted lines to separate the subplots
    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([1, 0], [0, 0], transform=ax.transAxes, **kwargs)
    ax.plot([1, 0], [1, 1], transform=ax1.transAxes, **kwargs)

    # Set y-limits for the axes
    ax.set_ylim(242, 258) #(165, 178)
    ax1.set_ylim(0, 63)

    # Add legend for bars only
    #fig.legend(legend_patches, np.concatenate([Lables_1]), loc="upper right",  fontsize=15,
    #           bbox_to_anchor=(.7, .97), borderaxespad=0.1)
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

def plot_ITS_bar_broken_process2(MD_nout, MSM_ITS_mean_1, MSM_ITS_std_1, Colors_1, Lables_1,
                        MSM_ITS_mean_2, MSM_ITS_std_2, Colors_2, Lables_2, cut, k=2):
    # Convert MD traj in fs to ps
    scale_factor = MD_nout * 0.001
    its_data_1 = np.array(MSM_ITS_mean_1) * scale_factor
    its_std_1 = np.array(MSM_ITS_std_1) * scale_factor
    its_data_2 = np.array(MSM_ITS_mean_2) * scale_factor
    its_std_2 = np.array(MSM_ITS_std_2) * scale_factor

    # Prepare data
    data = [
        (its_data_1, its_std_1, Colors_1, Lables_1),
        (its_data_2, its_std_2, Colors_2, Lables_2)
    ]
    
    # Setup figure and axes
    fig, (ax, ax1) = plt.subplots(2, 1, figsize=(4.25, 3.25), sharex=True)

    fig.subplots_adjust(hspace=0.0001)
    
    axes = [ax1, ax]
    LAB = [r'NH',  r'BE', r'VR',  ' ', r'$0.01$',r'$0.1$',r'$1$', r'$10$', r'$100$']

    # Initialize a list for the bar patches to use in the legend
    legend_patches = []

    # Process each dataset (1st and 2nd)
    for axis, (its_data, its_std, Colors, Lables) in zip(axes, data):
        means = [np.mean(its[cut:, k]) for its in its_data]
        std_devs = [np.mean(std[cut:, k]) for std in its_std]
        
        # Plot error bars and bars
        for mean, std_dev, label, color in zip(means, std_devs, Lables, Colors):
            if std_dev==0:
                pass
            else:
                axis.errorbar(label, mean, yerr=std_dev, fmt='-', capsize=5, capthick=2, color='k')
            if label=='VR':
                bar = axis.bar(label, mean, color=color, label=label,hatch='/', width=0.99)
            else:
                bar = axis.bar(label, mean, color=color, label=label, width=0.99)
            if axis==ax1:
                legend_patches.append(Patch(color=color, label=label))


    # Adjust x-axis labels and ticks
    ax1.tick_params(axis='x', rotation=90)
    ax1.set_xticklabels(LAB, ha="center")
    ax1.tick_params(axis="both", which="major", labelsize=19)
    ax.tick_params(axis="both", which="major", labelsize=19)

    # Add vertical axis label
    fig.text(0.0009, 0.5, 'Implied timescale in ps', ha='center', va='center', rotation='vertical', fontsize=19)

    # Remove spines between axes
    ax.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    #ax1.xaxis.tick_top()
    #ax.xaxis.tick_top()
    ax.tick_params(labelbottom=False)
    ax.tick_params(axis='x', which='both', bottom=False)  # Hide tick labels at the top

    # Add slanted lines to separate the subplots
    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([1, 0], [0, 0], transform=ax.transAxes, **kwargs)
    ax.plot([1, 0], [1, 1], transform=ax1.transAxes, **kwargs)

    # Set y-limits for the axes
    ax.set_ylim(165, 178) #change if you want to plot different its range
    ax1.set_ylim(0, 63)

    # Add legend for bars only
    #fig.legend(legend_patches, np.concatenate([Lables_1]), loc="upper right",  fontsize=15,
    #           bbox_to_anchor=(.7, .97), borderaxespad=0.1)
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

###############################################
# MSM eigenvectors for a number of pentane 
# molecules of a specific MD setup 
###############################################
input='' # input directory
traj_file = filename_sd_1 # change for integrator
gro_file = traj_file + '.gro'
trajectory = md.load(input+traj_file + '.trr', top=input+gro_file)
phi, psi = torions(trajectory, number=140)
lag = 50
lags=np.array(range(1,600, 10))
nstxout=100
CV_traj = np.array([psi, phi]).swapaxes(1,0)
dtrajs = cluster_UnitGrid_C(phi, psi)
plot_assignment_uniform(CV_traj, dtrajs, scip=10)

evecs,lcs=get_eigenvectors(
    dtraj=dtrajs.astype(int), lagtime=lag, reversible=True,
    stationary_distribution_constraint=None, 
    countmode='sliding', k=3)

plot_eigenvectors2D(evecs, CV_traj[:,0],CV_traj[:,1], dtrajs, 0, np.absolute(evecs[0]).max())
plot_eigenvectors2D(evecs, CV_traj[:,0],CV_traj[:,1], dtrajs, 1, np.absolute(evecs[1:]).max())
plot_eigenvectors2D(evecs, CV_traj[:,0],CV_traj[:,1], dtrajs, 2, np.absolute(evecs[1:]).max())

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
        its=get_implied_timescales(
            dtraj=dtrajs.astype(int), lagtimes=lags, k=nvecs, 
            countmode='sliding', reversible=True, 
            stationary_distribution_constraint=None)
        ITS.append(its)
    ITS=np.array(ITS)
    its_mean= np.mean(ITS, axis=0)
    its_std = np.std(ITS, axis=0)
    np.save(output+traj_file+'_'+'its_mean', its_mean)
    np.save(output+traj_file+'_'+'its_std', its_std)

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

its_vr= np.load(output+filename_vr+'_'+'its_mean.npy')
its_nh= np.load(output+filename_nh+'_'+'its_mean.npy')
its_be= np.load(output+filename_be+'_'+'its_mean.npy')
its_sd_001= np.load(output+filename_sd_001+'_'+'its_mean.npy')
its_sd_01= np.load(output+filename_sd_01+'_'+'its_mean.npy')
its_sd_1= np.load(output+filename_sd_1+'_'+'its_mean.npy')[:,0:2]
its_sd_10= np.load(output+filename_sd_10+'_'+'its_mean.npy')
its_sd_100= np.load(output+filename_sd_100+'_'+'its_mean.npy')

std_vr= np.load(output+filename_vr+'_'+'its_std.npy')
std_nh= np.load(output+filename_nh+'_'+'its_std.npy')
std_be= np.load(output+filename_be+'_'+'its_std.npy')
std_sd_001= np.load(output+filename_sd_001+'_'+'its_std.npy')
std_sd_01= np.load(output+filename_sd_01+'_'+'its_std.npy')
std_sd_1= np.load(output+filename_sd_1+'_'+'its_std.npy')[:,0:2]
std_sd_10= np.load(output+filename_sd_10+'_'+'its_std.npy')
std_sd_100= np.load(output+filename_sd_100+'_'+'its_std.npy')

# eigenvalues process 1
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

#  eigenvalues process 2 
plot_ITS_bar_broken_process2(100,
            ITS, STD,
            Colors, Lables,
            ITS_001, STD_001,
            Colors_001, Lables_001, 
            20,k=1)
