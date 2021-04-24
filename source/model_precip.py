import numpy as np
import os
from os.path import join
import fnmatch
from netCDF4 import Dataset
#from pickle_io import load_obj, save_obj
from matplotlib import colors
'''
from model_precip import writeNetCDF_2D, block_mean
'''

def writeNetCDF_2D(var, lat, lon, time, nc_filename, varname, file_description, **kwargs):
    
    # write NetCDF file
    if os.path.isfile(nc_filename):
        os.remove(nc_filename)
    nc_file = Dataset(nc_filename, 'w', format='NETCDF4')
    nc_file.description = file_description

    Ntime = len(time)
    Nlat = len(lat)
    Nlon = len(lon)
    
    # dimensions
    nc_file.createDimension('longitude', Nlon)
    nc_file.createDimension('latitude', Nlat)
    nc_file.createDimension('time', Ntime)

    # variables
    nc_time = nc_file.createVariable('time', 'f4', ('time',))
    nc_time.units = 'days since 1998-01-01 00:00:00'
    if 'calendar' in kwargs:
        nc_time.calendar = kwargs.get('calendar')
    else:
        nc_time.calendar = 'standard'
    nc_lon = nc_file.createVariable('longitude', 'f4', ('longitude',))
    nc_lon.units = 'degrees_east'
    nc_lon.long_name = 'longitude'
    nc_lat = nc_file.createVariable('latitude', 'f4', ('latitude',))
    nc_lat.units = 'degrees_north'
    nc_lat.long_name = 'latitude'
    nc_var = nc_file.createVariable(varname, 'f8', ('time', 'latitude', 'longitude'), fill_value = float('nan'))

    # assign variables
    nc_time[:] = time
    nc_lon[:] = lon
    nc_lat[:] = lat
    nc_var[:] = var

    nc_file.close()

def block_mean(var, n, axis):

    shapes = np.array(var.shape)
    if np.mod(shapes[axis], n) != 0:
        raise Exception('the block size has to be divisible by the size of that axis!')
    new_shapes = shapes[:]
    new_shapes[axis] = new_shapes[axis] / n
    new_shapes = np.insert(new_shapes, axis + 1, n)
    var = var.reshape(new_shapes)
    var = np.nanmean(var, axis = axis + 1)

    return var

def cal_precip_CRM_Yani(**kwargs):

    if 'obj_name' in kwargs:
        obj_name = kwargs.get('obj_name')
    else:
        obj_name = 'precip_data_CRM'
    if not os.path.isdir(obj_name):
        os.mkdir(obj_name)

    path = '/net/aimsir/archive1/janniy/Paul_runs/'

    # get 2D file names
    filename_2D_0K_1 = join(path, 'ctl_12km_576x1440x48_ctl_288.2Dcom_1.nc')
    filename_2D_0K_2 = join(path, 'ctl_12km_576x1440x48_ctl_288.2Dcom_2.nc')

    # select region between 30S to 30N
    R = 6371000.0 # in m
    dx = 12e3; dy = 12e3; dt = 1./8. # time is in days
    Omega = 7.2921e-5 # in rad/s
    x = np.arange(0, 576) * dx
    y = np.arange(0, 1440) * dy
    time = np.arange(1, 4801) * dt
    y_center = (np.arange(-720, 720) + 0.5) * dy
    lat = y_center / R / np.pi * 180.
    y_start = np.min(np.where(lat > -30)) # should be 442
    y_end   = np.max(np.where(lat <  30)) # should be 997
    y_ind   = np.arange(y_start, y_end + 1, 1)
    n       = 2 # accumulation indices, from 3-hourly to 6-hourly

    print('selecting CRM precipitation, Yani run')
    dataset_0K_1 = Dataset(filename_2D_0K_1)
    dataset_0K_2 = Dataset(filename_2D_0K_2)
    precip_0K_1 = dataset_0K_1.variables['Prec'][:, y_ind, :]
            # 100 days spin-in time according to Bill in his description of data
    precip_0K_2 = dataset_0K_2.variables['Prec'][:, y_ind, :]
    time_0K = time
    precip_0K = np.concatenate((precip_0K_1, precip_0K_2))
    del precip_0K_1, precip_0K_2
    # change from 3-hourly rates to 6-hourly rates
    precip_0K = block_mean(precip_0K, n, 0)
    time_0K = time_0K[1::2]
    # write NetCDF files
    varname = 'precip'
    nc_filename = join(obj_name, 'precip_0K_Yani_run.nc')
    file_description = 'CRM precipitation data'
    writeNetCDF_2D(precip_0K, y_center[y_ind], x, time_0K, nc_filename, varname, file_description)
    del precip_0K


def cal_precip_CRM(**kwargs):

    if 'obj_name' in kwargs:
        obj_name = kwargs.get('obj_name')
    else:
        obj_name = 'precip_data_CRM'
    if not os.path.isdir(obj_name):
        os.mkdir(obj_name)

    path_0K = '/net/aimsir/archive1/pog/bill_crm_data/qobskm12x576/'
    path_4K = '/net/aimsir/archive1/pog/bill_crm_data/qobs4Kkm12x576/'

    # get 2D file names
    filename_2D_0K_1 = join(path_0K, 'qobskm12x576_576x1440x48_ctl_288.2Dcom_1.nc')
    filename_2D_0K_2 = join(path_0K, 'qobskm12x576_576x1440x48_ctl_288.2Dcom_2.nc')
    filename_2D_4K_1 = join(path_4K, 'qobs4Kkm12x576_576x1440x48_ctl_288.2Dcom_1.nc')
    filename_2D_4K_2 = join(path_4K, 'qobs4Kkm12x576_576x1440x48_ctl_288.2Dcom_2.nc')

    # select region between 30S to 30N
    R = 6371000.0 # in m
    dx = 12e3; dy = 12e3; dt = 1./8. # time is in days
    Omega = 7.2921e-5 # in rad/s
    x = np.arange(0, 576) * dx
    y = np.arange(0, 1440) * dy
    time = np.arange(1, 5601) * dt
    y_center = (np.arange(-720, 720) + 0.5) * dy
    lat = y_center / R / np.pi * 180.
    y_start = np.min(np.where(lat > -30)) # should be 442
    y_end   = np.max(np.where(lat <  30)) # should be 997
    y_ind   = np.arange(y_start, y_end + 1, 1)
    n       = 2 # accumulation indices, from 3-hourly to 6-hourly

    print('selecting CRM 0K precipitation')
    dataset_0K_1 = Dataset(filename_2D_0K_1)
    dataset_0K_2 = Dataset(filename_2D_0K_2)
    precip_0K_1 = dataset_0K_1.variables['Prec'][800:, y_ind, :] 
            # 100 days spin-in time according to Bill in his description of data
    precip_0K_2 = dataset_0K_2.variables['Prec'][:, y_ind, :]
    time_0K = time[800:]
    precip_0K = np.concatenate((precip_0K_1, precip_0K_2))
    del precip_0K_1, precip_0K_2
    # change from 3-hourly rates to 6-hourly rates
    precip_0K = block_mean(precip_0K, n, 0)
    time_0K = time_0K[1::2]
    # write NetCDF files
    varname = 'precip'
    nc_filename = join(obj_name, 'precip_0K.nc')
    file_description = 'CRM precipitation data'
    writeNetCDF_2D(precip_0K, y_center[y_ind], x, time_0K, nc_filename, varname, file_description)
    del precip_0K
    
    print('selecting CRM 4K precipitation')
    dataset_4K_1 = Dataset(filename_2D_4K_1)
    dataset_4K_2 = Dataset(filename_2D_4K_2)
    precip_4K_1 = dataset_4K_1.variables['Prec'][800:, y_ind, :]
    precip_4K_2 = dataset_4K_2.variables['Prec'][:, y_ind, :]
    time_4K = time[800:]
    precip_4K = np.concatenate((precip_4K_1, precip_4K_2))
    del precip_4K_1, precip_4K_2
    # change from 3-hourly rates to 6-hourly rates
    precip_4K = block_mean(precip_4K, n, 0)
    time_4K = time_4K[1::2]
    # write NetCDF files
    nc_filename = join(obj_name, 'precip_4K.nc')
    file_description = 'CRM precipitation data'
    writeNetCDF_2D(precip_4K, y[y_ind], x, time_4K, nc_filename, varname, file_description)
    del precip_4K
    
    # calculate areas of each grid point
    area = np.ones((len(y_ind), len(x)))
    area = area * 12000**2 # in m^2
    nc_filename = join(obj_name, 'area.nc')
    file_description = 'CRM area data'
    writeNetCDF_2D(area, y_center[y_ind], x, [0], nc_filename, 'area', file_description)

def cal_precip_CESM(**kwargs):

    if 'obj_name' in kwargs:
        obj_name = kwargs.get('obj_name')
    else:
        obj_name = 'precip_data_CESM'
    if not os.path.isdir(obj_name):
        os.mkdir(obj_name)

    lat_start, lat_end = -30, 30
    data_path = '/net/aimsir/archive1/ziweili/CESM_LENS/data/'
    n_ensembles = [1, 2, 3, 4, 5, 35]
    file_description = 'CESM precipitation data'
    time_ind = np.arange(365*4*(2002 - 1990), 365*(2005 - 1990 + 1)*4) # years 2002-2005 for Paul's figure

    # get lat, lon arrays
    filename_0 = data_path + 'b.e11.B20TRC5CNBDRD.f09_g16.' + \
            '001.cam.h2.PRECT.1990010100Z-2005123118Z.nc'
    dataset_0 = Dataset(filename_0)
    lat = np.array(dataset_0.variables['lat'])
    lon = np.array(dataset_0.variables['lon'])
    y_start = np.min(np.where(lat > lat_start))
    y_end   = np.max(np.where(lat < lat_end))
    y_ind   = np.arange(y_start, y_end + 1, 1)
    varname = 'precip'
    
    # historical data
    for n_ensemble in n_ensembles:
        n_str = '{0:03d}'.format(n_ensemble)
        print('CESM historical: ' + n_str)
        filename_h = data_path + 'b.e11.B20TRC5CNBDRD.f09_g16.' + \
                n_str + '.cam.h2.PRECT.1990010100Z-2005123118Z.nc'
        dataset_h = Dataset(filename_h)
        time = np.array(dataset_h.variables['time'])
        time = time - (1998 - 1850) * 365 # move from days since 1850 to 1998
        precip_h = dataset_h.variables['PRECT'][time_ind, y_ind, :] * 86400.0 * 1000 
                # get years 1991-2000, and transform from m/s to mm/day
        nc_filename = join(obj_name, 'precip_h_' + n_str + '.nc')
        writeNetCDF_2D(precip_h, lat[y_ind], lon, time[time_ind], 
                nc_filename, varname, file_description, calendar = 'no_leap')
    
    # calculate areas of each grid point
    R = 6371000.0 # in m
    dlon = lat[1] - lat[0]
    dlat = lon[1] - lon[0]
    Lon, Lat = np.meshgrid(lon, lat)
    area = R**2 * np.cos(Lat/180.*np.pi) * dlat/180.*np.pi * dlon/180.*np.pi # in m^2
    # save as a NetCDF file
    nc_filename = join(obj_name, 'area.nc')
    file_description = 'CESM area data'
    writeNetCDF_2D(area[y_ind, :], lat[y_ind], lon, [0], nc_filename, 'area', file_description)


def cal_precip_TRMM(**kwargs):
    
    if 'obj_name' in kwargs:
        obj_name = kwargs.get('obj_name')
    else:
        obj_name = 'precip_data_TRMM'
    if not os.path.isdir(obj_name):
        os.mkdir(obj_name)

    from pyhdf.SD import SD, SDC # documentation: http://fhs.github.io/pyhdf/
    path = '/net/aimsir/archive1/ziweili/cluster_pdf/TRMM_3B42/'
    
    # get years
    filenames_1 = []
    years = np.arange(2002, 2005 + 1, 1) # 2002-2005 for figure requested by Paul
    for year in years:
        filenames_1 = filenames_1 + [name for name in os.listdir(path) if 
                fnmatch.fnmatch(name, '3B42.' + str(year) + '????.??.*.HDF')]
    # get months
    filenames_2 = []
    months = np.arange(1, 13, 1) # Jan. - Dec.
    for month in months:
        filenames_2 = filenames_2 + [name for name in filenames_1 if 
                fnmatch.fnmatch(name, '3B42.????' + '{0:02d}'.format(month) + '??.??.*.HDF')]
    
    filenames_3 = []
    hours = ['00', '03', '06', '09', '12', '15', '18', '21'] # get 3-hourly, and average to be 6-hourly
    hour_interval = 3.;
    for hour in hours:
        filenames_3 = filenames_3 + [name for name in filenames_2 if
                fnmatch.fnmatch(name, '3B42.????????.' + hour + '.*.HDF')]
    
    filenames = [path + name for name in sorted(filenames_3)]
    del filenames_1, filenames_2, filenames_3

    filename = filenames[0]
    yx_shape = SD(filename, SDC.READ).select('precipitation').get().T.shape
    
    # save as NetCDF files
    precip = np.zeros((len(filenames), yx_shape[0], yx_shape[1]), dtype = 'float16')
    for counter in range(len(filenames)):
        filename = filenames[counter]
        if (counter + 1) % 100 == 0:
            print('counter = {0:d}/{1:d}'.format(counter+1, len(filenames)))
        file = SD(filename, SDC.READ)
        precip_obj = file.select('precipitation')
        precip[counter, :, :] = precip_obj.get().T * 24.
            # from mm/h to mm/day
        file.end()
    
    dlat = 0.25
    dlon = 0.25
    lat = np.arange(-49.875, 49.875 + dlat, dlat)
    lon = np.arange(-180, 180, dlon)
    
    # select latitudes between 30S and 30N
    #lat_start = -30; lat_end = 30
    lat_start = -50; lat_end = 50
    y_start = np.min(np.where(lat >= lat_start))
    y_end = np.max(np.where(lat <= lat_end))
    y_ind   = np.arange(y_start, y_end + 1, 1)
    lat = lat[y_ind]

    nc_filename = join(obj_name, 'precip_TRMM_{0:d}S-{0:d}N_3hourly.nc'.format(np.abs(lat_start), np.abs(lat_end)))
    file_description = 'TRMM precipitation data'
    varname = 'precip'
    time = (np.arange(precip.shape[0])+1) * hour_interval / 24. + (years[0] - 1998) * 365 + np.ceil((years[0] - 2000.) / 4.)
            # convert to days since 1998-01-01, will be recognized by ncview
            # the third term accounts for leap years of 2000, 2004 and 2008 (2000/400 is an integer)
    precip = precip[:, y_ind, :]
    precip[precip == -float('inf')] = float('nan')

    # change from 3-hourly rates to 6-hourly rates
    #n = 2
    n = 1
    if n == 1:
        precip_in = precip
    else:
        precip_in = block_mean(precip, n, 0)
        time = time[1::n]
        

    # need to roll the data to get longitudinal degrees larger than zero
    var = np.zeros(precip_in.shape)
    for t in range(precip_in.shape[0]):
        if (t + 1) % 100 == 0:
            print('t = {0:d}'.format(t+1))
        var[t, :, :] = np.roll(precip_in[t, :, :], -720, axis = 1)
    lon2 = np.mod(np.roll(lon, -720), 360)
   
    # write NetCDF file
    writeNetCDF_2D(var, lat, lon2, time, nc_filename, varname, file_description)

    # calculate areas of each grid point
    R = 6371000.0 # in m
    Lon, Lat = np.meshgrid(lon2, lat)
    area = R**2 * np.cos(Lat/180.*np.pi) * dlat/180.*np.pi * dlon/180.*np.pi # in m^2
    nc_filename = join(obj_name, 'area.nc')
    file_description = 'TRMM area data'
    writeNetCDF_2D(area, lat, lon2, [0], nc_filename, 'area', file_description)


