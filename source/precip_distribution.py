import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import join
import scipy as sp

def CESM_precip(data_path, **kwargs):

    obj_name = join(data_path, 'precip_data_CESM')
    lat_start = -15; lat_end = 15
    lon_start = 160; lon_end = 222.2; # central Pacific, consistent with the longitude span of 62.2 degrees in CRM

    lat = Dataset(join(obj_name, 'precip_h_001.nc')).variables['latitude'][:]
    lon = Dataset(join(obj_name, 'precip_h_001.nc')).variables['longitude'][:]
    time = Dataset(join(obj_name, 'precip_h_001.nc')).variables['time'][:]
    y_ind = np.where(np.logical_and(lat >= lat_start, lat <= lat_end).reshape(-1, 1))[0]
    x_ind = np.where(np.logical_and(lon >= lon_start, lon <= lon_end).reshape(-1, 1))[0]
    area = Dataset(join(obj_name, 'area.nc')).variables['area'][0, y_ind, x_ind]

    #n_ensembles = [1, 2, 3, 4, 5, 35]
    n_ensembles = [1] # for Paul's figure

    precip = np.zeros((len(n_ensembles)*len(time), len(y_ind), len(x_ind)))
    for i in range(len(n_ensembles)):
        n = n_ensembles[i]
        precip[i*len(time):(i+1)*len(time), :, :] = \
            Dataset(join(obj_name, 'precip_h_{0:03d}.nc'.format(n))).variables['precip'][:, y_ind, x_ind]
    
    return precip, area, lon[x_ind], lat[y_ind], time

def CRM_precip(data_path, **kwargs):

    obj_name = join(data_path, 'precip_data_CRM')
    R = 6371000.0
    lat_start = -15; lat_end = 15
    lon_start = 0; lon_end = 62.2;

    lat = Dataset(join(obj_name, 'precip_0K.nc')).variables['latitude'][:] / R / np.pi * 180
    lon = Dataset(join(obj_name, 'precip_0K.nc')).variables['longitude'][:] / R / np.pi * 180
    y_ind = np.where(np.logical_and(lat >= lat_start, lat <= lat_end))[0]
    x_ind = np.where(np.logical_and(lon >= lon_start, lon <= lon_end))[0]
    area = Dataset(join(obj_name, 'area.nc')).variables['area'][0, y_ind, x_ind]

    precip_0K_1 = Dataset(join(obj_name, 'precip_0K.nc')).variables['precip'][:, y_ind, x_ind]
    precip_0K_2 = Dataset(join(obj_name, 'precip_0K_Yani_run.nc')).variables['precip'][:, y_ind, x_ind]
    precip_0K = np.concatenate((precip_0K_1, precip_0K_2))

    time = np.arange(precip_0K.shape[0])

    return precip_0K, area, lon[x_ind] + 160, lat[y_ind], time

def TRMM_precip(data_path, **kwargs):

    obj_name = join(data_path, 'precip_data_TRMM')
    lat_start = -15; lat_end = 15
    if 'lon_start' in kwargs:
        lon_start = kwargs.get('lon_start')
    else:
        lon_start = 160 # central Pacific, consistent with the span of 62.2 degrees in CRM
    if 'lon_end' in kwargs:
        lon_end = kwargs.get('lon_end')
    else:
        lon_end = 222.2 # central Pacific, consistent with the span of 62.2 degrees in CRM

    lat = Dataset(join(obj_name, 'precip_TRMM_30S-30N.nc')).variables['latitude'][:]
    lon = Dataset(join(obj_name, 'precip_TRMM_30S-30N.nc')).variables['longitude'][:]
    time = Dataset(join(obj_name, 'precip_TRMM_30S-30N.nc')).variables['time'][:]
    y_ind = np.where(np.logical_and(lat >= lat_start, lat <= lat_end))[0]
    x_ind = np.where(np.logical_and(lon >= lon_start, lon <= lon_end))[0]
    area = Dataset(join(obj_name, 'area.nc')).variables['area'][0, y_ind, x_ind]

    if 'Nt' in kwargs:
        Nt = kwargs.get('Nt')
    else:
        Nt = len(Dataset(join(obj_name, 'precip_TRMM_30S-30N.nc')).variables['time'])
    t_ind = np.arange(0, Nt, 1)
        
    precip = Dataset(join(obj_name, 'precip_TRMM_30S-30N.nc')).variables['precip'][t_ind, y_ind, x_ind]

    return precip, area, lon[x_ind], lat[y_ind], time

def area_interpolation(precip_source, lon_source, lat_source, lon_target, lat_target, **kwargs):

    if len(precip_source.shape) == 2:
        precip_source = precip_source.reshape(np.concatenate(([1], precip_source.shape)))

    dlon_source = lon_source[1] - lon_source[0]
    dlat_source = lat_source[1] - lat_source[0]
    lon_source_bnds = np.concatenate((lon_source - dlon_source/2, [lon_source[-1] + dlon_source/2]))
    lat_source_bnds = np.concatenate((lat_source - dlat_source/2, [lat_source[-1] + dlat_source/2]))

    Nt = precip_source.shape[0]
    Nlon_t = len(lon_target)
    Nlat_t = len(lat_target)
    Nlon_s = len(lon_source)
    Nlat_s = len(lat_source)
    
    # get weight matrix
    #weights = np.zeros((Nlat_t, Nlon_t, Nlat_s, Nlon_s))
    precip_target = np.zeros((Nt, Nlat_t, Nlon_t))
    temp_weights = np.zeros((Nlat_s, Nlon_s))
    dlon_target = lon_target[1] - lon_target[0]
    dlat_target = lat_target[1] - lat_target[0]
    
    for i in range(Nlon_t):
        print('{0:d} out of {1:d}'.format(i, Nlon_t))
        lon_center = lon_target[i]
        lon_l = lon_target[i] - dlon_target/2 # left bound
        lon_r = lon_target[i] + dlon_target/2 # right bound
        dis_l = lon_source_bnds[1:] - lon_l
        dis_r = lon_r - lon_source_bnds[0:-1]
        lon_ind = np.logical_and(dis_l > 0, dis_r > 0)
        dis_lon = np.minimum(dis_l.clip(0, dlon_source), dis_r.clip(0, dlon_source))
        for j in range(Nlat_t):
            lat_center = lat_target[j]
            lat_b = lat_target[j] - dlat_target/2 # bottom bound
            lat_u = lat_target[j] + dlat_target/2 # upper bound
            # figure out which boxes from the source grid lies within the target grid
            dis_b = lat_source_bnds[1:] - lat_b
            dis_u = lat_u - lat_source_bnds[0:-1] 
            lat_ind = np.logical_and(dis_b > 0, dis_u > 0)
            dis_lat = np.minimum(dis_b.clip(0, dlat_source), dis_u.clip(0, dlat_source))
            if 'AREA' in kwargs and kwargs.get('AREA'):
                lat_overlap_u = (lat_b + np.cumsum(dis_lat)) * (dis_lat != 0).astype('float')
                lat_overlap_b = lat_overlap_u - dis_lat
                temp_weights[:, :] = np.outer(np.sin(lat_overlap_u/180*np.pi) - np.sin(lat_overlap_b/180*np.pi), dis_lon)
            else:
                temp_weights[:, :] = np.outer(dis_lat, dis_lon)
            temp_weights = temp_weights / np.sum(temp_weights)
            #weights[j, i, :, :] = temp_weights 

            for t in range(Nt):
                precip_target[t, j, i] = np.sum(precip_source[t, :, :] * temp_weights)

    '''
    # get interpolated precip
    for t in range(Nt):
        if np.mod(t, 10) == 0:
            print('{0:d} out of {1:d}'.format(t, Nt))
        for i in range(Nlon_t):
            for j in range(Nlat_t):
                precip_target[t, j, i] = np.sum(precip_source[t, :, :] * weights[j, i, :, :])
    '''

    '''
    # test case
    dlat_source = 2.; dlon_source = 3.
    dlat_target = 9.; dlon_target = 8.
    lon_source = np.arange(dlon_source/2, 360 + dlon_source/2, dlon_source)
    lat_source = np.arange(-90 + dlat_source/2, 90 + dlat_source/2, dlat_source)
    lon_target = np.arange(dlon_target/2, 360 + dlon_target/2, dlon_target)
    lat_target = np.arange(-90 + dlat_target/2, 90 + dlat_target/2, dlat_target)
    Lon_source, Lat_source = np.meshgrid(lon_source, lat_source)
    precip_source = np.cos(Lat_source / 180 * np.pi) * np.sin(Lon_source / 180 * np.pi / 2)
    precip_target = area_interpolation(precip_source, lon_source, lat_source, lon_target, lat_target)
    precip_target = precip_target[0, :, :]
    
    plt.subplot(1, 2, 1)
    plt.imshow(precip_source)
    plt.colorbar()
    plt.subplot(1, 2, 2)
    plt.imshow(precip_target)
    plt.colorbar()
    plt.show()

    np.sum(precip_target * dlat_target * dlon_target) - np.sum(precip_source * dlat_source * dlon_source)
    '''

    return precip_target

def plot_precip_distribution(fig, datasets, areas, colors_plot, labels, **kwargs):

    Nd = len(datasets)

    if 'bin_size' in kwargs:
        bin_size = kwargs.get('bin_size')
    else:
        bin_size = 30
    if 'max_size' in kwargs:
        max_size = kwargs.get('max_size')
    else:
        max_size = np.max([np.max(dataset) for dataset in datasets])
    if 'min_size' in kwargs:
        min_size = kwargs.get('min_size')
    else:
        min_size = np.max([np.min(dataset[dataset != 0]) for dataset in datasets])
        min_size = max((min_size, 1))

    # treat 0 points separately
    normalize_factor = np.zeros((Nd, 1))
    for c in range(Nd):
        normalize_factor[c] = 1 - np.sum(datasets[c] < min_size) / float(np.prod(datasets[c].shape))

    log10_range = np.ceil(np.log10(max_size/min_size)*100)/100
    bins = np.power(10, np.linspace(np.log10(min_size), log10_range + np.log10(min_size), bin_size))
    bins_x = np.sqrt(bins[1:] / bins[0:-1]) * bins[0:-1]
    
    hist_count = np.zeros((bin_size - 1, Nd))
    for c in range(Nd):
        precip = datasets[c]
        weights = np.tile(areas[c], (precip.shape[0], 1, 1))
        hist_count_0, _ = np.histogram(precip, bins = bins, weights = weights, density = True)
        hist_count[:, c] = hist_count_0 * normalize_factor[c]
        #hist_count_0, _ = np.histogram(precip, bins = bins, density = True)
        #hist_count[:, c] = hist_count_0 / (bins[1:] - bins[0:-1]) \
        #        / float(sum(hist_count_0)) * normalize_factor[c]

    # get 99.9 and 99.99 percentile
    percentile_999  = np.zeros(Nd)
    percentile_9999 = np.zeros(Nd)
    from scipy import interpolate
    for c in range(Nd):
        CDF = 1 - normalize_factor[c] + np.cumsum(hist_count[:, c] * (bins[1:] - bins[0:-1]))
        f = interpolate.interp1d(CDF, bins_x)
        percentile_999[c] = f(0.999)
        percentile_9999[c] = f(0.9999)

    if 'LATEX' in kwargs and kwargs.get('LATEX'):
        plt.rc('text', usetex=True)
        plt.rcParams["font.family"] = "Times New Roman"

    # begin plotting
    ax = fig.add_subplot(121)
    
    p = []
    hist_count[hist_count == 0] = np.nan
    for c in range(Nd):
        temp_p = plt.plot(bins_x, hist_count[:, c], colors_plot[c], linewidth = 1)
        p.append(temp_p)
        # plot dots for percentiles
        f = interpolate.interp1d(np.log(bins_x), np.log(hist_count[:, c]))
        plt.scatter(percentile_999 [c], np.exp(f(np.log(percentile_999 [c]))), s=10, marker='o', c=colors_plot[c])
        #plt.scatter(percentile_9999[c], np.exp(f(np.log(percentile_9999[c]))), s=10, marker='s', c=colors_plot[c])

    if 'LEGEND' in kwargs:
        LEGEND = kwargs.get('LEGEND')
    else:
        LEGEND = true
    if LEGEND:
        ax.legend(handles = [a[0] for a in p], labels = labels, frameon = False, loc = 'lower left')

    plt.xscale('log')
    plt.yscale('log')
    if 'ylim' in kwargs:
        plt.ylim(kwargs.get('ylim'))
        
    # reduce number of minor ticks in y
    import matplotlib.ticker as mtick
    locmaj = mtick.LogLocator(base=10,numticks=10)
    ax.yaxis.set_major_locator(locmaj)
    locmin = mtick.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
    ax.yaxis.set_minor_locator(locmin)    
    ax.yaxis.set_minor_formatter(mtick.NullFormatter())

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xlabel(r'Precipitation rate (mm day$^{-1}$)')
    plt.ylabel(r'Frequency distribution (mm$^{-1}$ day)')
    if 'title' in kwargs:
        plt.title(kwargs.get('title'))
    
def prepare_data_for_precip_distribution(data_path):
    
    # get precipitation data
    precip_CESM, area_CESM, lon_CESM, lat_CESM, time_CESM = CESM_precip(data_path)
    precip_CRM , area_CRM , lon_CRM , lat_CRM , time_CRM  = CRM_precip (data_path)
    precip_TRMM, area_TRMM, lon_TRMM, lat_TRMM, time_TRMM = TRMM_precip(data_path)

    nc_filename_CRM  = join(data_path, 'precip_data_CRM/precip_CRM_at_CESM.nc' )
    nc_filename_TRMM = join(data_path, 'precip_data_TRMM/precip_TRMM_at_CESM.nc')
    nc_filename_CESM = join(data_path, 'precip_data_CESM/precip_CESM.nc'        )
    
    precip_CRM_at_CESM  = area_interpolation(precip_CRM, lon_CRM, lat_CRM, lon_CESM, lat_CESM, AREA = False)
    writeNetCDF_2D(precip_CRM_at_CESM, lat_CESM, lon_CESM, time_CRM, 
            nc_filename_CRM, 'precip', 'SAM 0K precip at CESM resolution')
    precip_TRMM_at_CESM = area_interpolation(precip_TRMM, lon_TRMM, lat_TRMM, lon_CESM, lat_CESM, AREA = True)
    writeNetCDF_2D(precip_TRMM_at_CESM, lat_CESM, lon_CESM, time_TRMM, 
            nc_filename_TRMM, 'precip', 'TRMM precip at CESM resolution')
    writeNetCDF_2D(precip_CESM, lat_CESM, lon_CESM, time_CESM,
            nc_filename_CESM, 'precip', 'CESM precip that matches TRMM and SAM')

    '''
    # compare area weighting using cos(lat) and non-area weighting
    precip_TRMM_at_CESM = Dataset('precip_data_TRMM/precip_TRMM_at_CESM.nc').variables['precip']
    precip_TRMM_at_CESM_original = area_interpolation(precip_TRMM, lon_TRMM, lat_TRMM, lon_CESM, lat_CESM, AREA = False)
    np.max(np.abs(precip_TRMM_at_CESM - precip_TRMM_at_CESM_original)) # this number is 0.408, and the mean is 0.001
    '''

