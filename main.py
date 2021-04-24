import numpy as np
import matplotlib.pyplot as plt
import source
import os
from netCDF4 import Dataset
from os.path import join
from source import plot_precip_distribution, prepare_data_for_precip_distribution 
    # these are in source/precip_distribution.py
from source import sizes_of_cluster, size_struct_to_array, power_law, power_law_plot, plot_dict_update, plot_dict_init
    # these are in source/classes.py
from source import load_obj, save_obj
    # these are in source/pickle_io.py


def CESM_precip_clusters(data_path, **kwargs):
        
    lat_start = -15; lat_end = 15
    lon_start = 160; lon_end = 222.2; # central Pacific, consistent with the longitude span of 62.2 degrees in CRM

    if 'p_threshold' in kwargs:
        p_threshold = kwargs.get('p_threshold')    
    else:
        p_threshold = 0.7 * 24
    thr_str = '_{0:.1f}'.format(p_threshold/24)
    lat = Dataset(join(data_path, 'precip_data_CESM/precip_h_001.nc')).variables['latitude'][:]
    lon = Dataset(join(data_path, 'precip_data_CESM/precip_h_001.nc')).variables['longitude'][:]
    y_ind = np.where(np.logical_and(lat >= lat_start, lat <= lat_end).reshape(-1, 1))[0]
    x_ind = np.where(np.logical_and(lon >= lon_start, lon <= lon_end).reshape(-1, 1))[0]
    area = Dataset(join(data_path, 'precip_data_CESM/area.nc')).variables['area'][0, y_ind, x_ind]

    n_ensembles = [1] # for Paul's figure
    Sizes_h, Powers_h = [], []

    for n in n_ensembles:
        precip_h = Dataset(join(data_path, 'precip_data_CESM/precip_h_{0:03d}.nc'.format(n))).variables['precip'][:, y_ind, x_ind]
        p_ind_h = precip_h > p_threshold
        temp_Sizes_h, temp_Powers_h = sizes_of_cluster(p_ind_h, precip_h, area, PERIODIC_X = False)
        Sizes_h = Sizes_h + temp_Sizes_h
        Powers_h = Powers_h + temp_Powers_h
        del p_ind_h, precip_h
    
    save_obj(Sizes_h, data_path, 'precip_data_CESM/Sizes_h' + thr_str)
    save_obj(Powers_h, data_path, 'precip_data_CESM/Powers_h' + thr_str)

    # calculate and save the total area of the domain in meters^2
    save_obj(np.sum(area), data_path, 'precip_data_CESM/total_area')

def CRM_precip_clusters(data_path, **kwargs):

    R = 6371000.0
    lat_start = -15; lat_end = 15
    lon_start = 0; lon_end = 62.2;
    if 'p_threshold' in kwargs:
        p_threshold = kwargs.get('p_threshold')                    
    else:                            
        p_threshold = 0.7 * 24
    thr_str = '_{0:.1f}'.format(p_threshold/24)
    lat = Dataset(join(data_path, 'precip_data_CRM/precip_0K.nc')).variables['latitude'][:] / R / np.pi * 180
    lon = Dataset(join(data_path, 'precip_data_CRM/precip_0K.nc')).variables['longitude'][:] / R / np.pi * 180
    y_ind = np.where(np.logical_and(lat >= lat_start, lat <= lat_end))[0]
    x_ind = np.where(np.logical_and(lon >= lon_start, lon <= lon_end))[0]
    area = Dataset(join(data_path, 'precip_data_CRM/area.nc')).variables['area'][0, y_ind, x_ind]

    precip_0K_1 = Dataset(join(data_path, 'precip_data_CRM/precip_0K.nc')).variables['precip'][:, y_ind, x_ind]
    precip_0K_2 = Dataset(join(data_path, 'precip_data_CRM/precip_0K_Yani_run.nc')).variables['precip'][:, y_ind, x_ind]
    precip_0K = np.concatenate((precip_0K_1, precip_0K_2))
    p_ind_0K = precip_0K > p_threshold
    Sizes_0K, Powers_0K = sizes_of_cluster(p_ind_0K, precip_0K, area, PERIODIC_X = True)
    save_obj(Sizes_0K, data_path, 'precip_data_CRM/Sizes_0K' + thr_str)
    save_obj(Powers_0K, data_path, 'precip_data_CRM/Powers_0K' + thr_str)
    del precip_0K_1, precip_0K_2, p_ind_0K, precip_0K

    # calculate and save the total area of the domain in meters^2
    save_obj(np.sum(area), data_path, 'precip_data_CRM/total_area')

def TRMM_precip_clusters(data_path, **kwargs):
    
    lat_start = -15; lat_end = 15
    lon_start = 160; lon_end = 222.2; # central Pacific, consistent with the longitude span of 62.2 degrees in CRM

    if 'p_threshold' in kwargs:
        p_threshold = kwargs.get('p_threshold')
    else:
        p_threshold = 0.7 * 24
    thr_str = '_{0:.1f}'.format(p_threshold/24)
    lat = Dataset(join(data_path, 'precip_data_TRMM/precip_TRMM_30S-30N.nc')).variables['latitude'][:]
    lon = Dataset(join(data_path, 'precip_data_TRMM/precip_TRMM_30S-30N.nc')).variables['longitude'][:]
    y_ind = np.where(np.logical_and(lat >= lat_start, lat <= lat_end))[0]
    x_ind = np.where(np.logical_and(lon >= lon_start, lon <= lon_end))[0]
    area = Dataset(join(data_path, 'precip_data_TRMM/area.nc')).variables['area'][0, y_ind, x_ind]
    
    precip = Dataset(join(data_path, 'precip_data_TRMM/precip_TRMM_30S-30N.nc')).variables['precip'][:, y_ind, x_ind]
    p_ind = precip > p_threshold
    
    Sizes, Powers = sizes_of_cluster(p_ind, precip, area, PERIODIC_X = False)
    
    del p_ind, precip
    save_obj(Sizes, data_path, 'precip_data_TRMM/Sizes' + thr_str)
    save_obj(Powers, data_path, 'precip_data_TRMM/Powers' + thr_str)

    # calculate and save the total area of the domain in meters^2
    save_obj(np.sum(area), data_path, 'precip_data_TRMM/total_area')

def precip_clusters():

    # a wrapper function to select the precipitation data, and calculate precipitation clusters
    
    # define the lower threshold of precip to include points in precip clusters
    p_threshold = 0.7 * 24;
    thr_str = '_{0:.1f}'.format(p_threshold/24)
    data_path = './data/'

    # CRM precipitation
    CRM_precip_clusters(data_path, p_threshold = p_threshold)

    # CESM precipitation
    CESM_precip_clusters(data_path, p_threshold = p_threshold)

    # TRMM precipitation
    TRMM_precip_clusters(data_path, p_threshold = p_threshold)

if __name__ == '__main__':

    # select precipitation data and calculate precipitation clusters
    precip_clusters()

    data_path = './data/'
    # collect precipitation
    prepare_data_for_precip_distribution(data_path)

    ## Plots for the paper in Phylosophical Transactions of the Royal Society B ##
    plot_path = './figures/'
    if not os.path.isdir(plot_path):
        os.mkdir(plot_path)

    ## Plot 1: compare CRM and TRMM at CESM resolution
    
    nc_filename_CRM  = join(data_path, 'precip_data_CRM/precip_CRM_at_CESM.nc')
    nc_filename_TRMM = join(data_path, 'precip_data_TRMM/precip_TRMM_at_CESM.nc')
    nc_filename_CESM = join(data_path, 'precip_data_CESM/precip_CESM.nc')
    precip_CRM_at_CESM  = Dataset(nc_filename_CRM ).variables['precip'][:]
    precip_TRMM_at_CESM = Dataset(nc_filename_TRMM).variables['precip'][:]
    precip_CESM         = Dataset(nc_filename_CESM).variables['precip'][:]
    area_CESM           = Dataset(join(data_path, 'precip_data_CESM/area.nc')).variables['area'][0, 16:48, 0:50]
    area_CRM            = np.ones(area_CESM.shape)
    datasets = [precip_CRM_at_CESM, precip_TRMM_at_CESM, precip_CESM]
    areas    = [area_CRM, area_CESM, area_CESM]
    
    # set colors
    colors = ['tab:blue', 'tab:orange', 'tab:red', 'tab:green', 'tab:purple', 'tab:brown']
    color_ind = [3, 2, 0]
    colors_plot = [colors[i] for i in color_ind]
    
    labels = [r'Control simulation', r'TRMM (observations)', r'CESM-LE (GCM)']
    filename = plot_path + 'precip_distributions.pdf'
    figsize = (10, 4.2)
    fig = plt.figure(figsize = figsize)
    
    plot_precip_distribution(fig, datasets, areas, colors_plot, labels, 
            min_size = 1, bin_size = 30, title = r'(a) Precipitation rate', LEGEND = True, LATEX = False, ylim = (1e-7, 0.5))
    
    ## Plot 2: distribution of cluster sizes
    
    class model_params:
        def __init__(self, obj_name, size_name, power_name, model_name, color, 
                fit_start_s, fit_end_s, fit_start_p, fit_end_p, marker, area_name):
            self.obj_name       = obj_name
            self.size_name      = size_name
            self.power_name     = power_name
            self.model_name     = model_name
            self.color          = color
            self.fit_start_s    = fit_start_s
            self.fit_end_s      = fit_end_s
            self.fit_start_p    = fit_start_p
            self.fit_end_p      = fit_end_p
            self.marker         = marker
            self.area_name      = area_name

    model_params_dict = {}
    thr_str = '_0.7'
    
    model_params_dict.update({'SAM_0K': model_params('./data/precip_data_CRM' , 'Sizes_0K'+thr_str, 'Powers_0K'+thr_str, 
        r'SAM' , colors_plot[0], 0, 17, 0, 17, ',', 'total_area')})
    model_params_dict.update({'TRMM'  : model_params('./data/precip_data_TRMM', 'Sizes'   +thr_str, 'Powers'   +thr_str, 
        r'TRMM', colors_plot[1], 0, 17, 0, 17, ',', 'total_area')})
    model_params_dict.update({'CESM'  : model_params('./data/precip_data_CESM', 'Sizes_h' +thr_str, 'Powers_h' +thr_str, 
        r'CESM', colors_plot[2], 0, 17, 0, 17, ',', 'total_area')})

    fig_name_size  = plot_path + 'multi_cluster_size'  + thr_str + '.pdf'
    fig_name_power = plot_path + 'multi_cluster_power' + thr_str + '.pdf'

    model_names = ['SAM_0K', 'TRMM', 'CESM']
    counter = 4;
    plot_dict_size  = plot_dict_init(fig_name_size , r'Cluster size (km$^2$)', r'Frequency distribution (km$^{-4}$)')
    plot_dict_power = plot_dict_init(fig_name_power, r'Cluster power (GW)'   , r'Frequency distribution (GW$^{-1}$ km$^{-2}$)')
    Rs = np.zeros((len(model_names))) 
        # log likelihood ratio between power law and lognormal, positive means power law is more likely
    p_values = np.zeros((len(model_names))) # p value of the sign of R
    sigmas = np.zeros((len(model_names))) # standard error of the exponent
    MLE = False

    for i in range(len(model_names)):
        counter = counter - 1
        model = model_params_dict[model_names[i]]
        
        obj_name = model.obj_name
        total_area = load_obj(obj_name, model.area_name)
        Sizes = load_obj(obj_name, model.size_name)
        Powers = load_obj(obj_name, model.power_name)
        sizes2  = size_struct_to_array(Sizes) / 1e6 # change from m^2 to km^2
        powers2 = size_struct_to_array(Powers)
        N = len(Sizes) * total_area / 1e6 # change from m^2 to km^2
            # normalize by the total area and time steps to get the number density of observing a cluster 
            # of certain size per km^2 at any instance. 
            # This would integrate to the frequency of findning a cluster per km^2
        # cluster size
        fit_start_s, fit_end_s = model.fit_start_s, model.fit_end_s
        reg, hist, bins_x, bins, sigmas[i], Rs[i], p_values[i] = power_law(sizes2, fit_start_s, fit_end_s, 
                absolute_prob = True, N = N, MLE = MLE, LIMITS = True)
        
        # cluster power of precip events
        fit_start_p, fit_end_p = model.fit_start_p, model.fit_end_p
        reg_p, hist_p, bins_p_x, bins_p, sigma, R, p = power_law(powers2, fit_start_p, fit_end_p, 
                absolute_prob = True, N = N, MLE = MLE)
        
        plot_dict_size = plot_dict_update(plot_dict_size,
                    reg, hist.reshape(-1, 1), bins_x.reshape(-1, 1), model.color, model.marker, model.model_name, fit_start_s, fit_end_s)
        plot_dict_power = plot_dict_update(plot_dict_power, 
                    reg_p, hist_p.reshape(-1, 1), bins_p_x.reshape(-1, 1), model.color, model.marker, model.model_name, fit_start_p, fit_end_p)

    power_law_plot(fig, plot_dict_size, LINE = True, title = r'(b) Cluster size', LATEX = False, ylim = (3e-15, 2e-7))
 
    fig.savefig(filename, bbox_inches = 'tight')
    fig.show()
    

