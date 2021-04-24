import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import powerlaw

class grid_point:
    def __init__(self, x, y, ind, bonds_ind, precip, area):
        self.sub = np.array([y, x])
        self.ind = ind
            # index for this grid point
        self.bonds_ind = bonds_ind
            # 1-by-4 vector for the indices of points that it has connection to        
        self.precip = precip
        self.area = area

class cluster:
    
    def __init__(self, n_x, n_y):
        self.n_x = n_x
        self.n_y = n_y
        self.subs = np.array([], dtype = 'int')
        self.inds = np.array([], dtype = 'int')
        self.center_sub = []
        self.sizes = np.array([], dtype = 'float')
        self.precips = np.array([], dtype = 'float')
        self.toappend = np.array([], dtype = 'int')
        self.cluster_power = []
    
    def get_toappend(self, point, ind2point):
        inds = point.bonds_ind
        # find the nearest neighbor bonds and get rid of duplicate points
        self.toappend = np.unique(np.append(self.toappend, inds[inds != -1]))
        # remove points already assigned to the cluster
        self.toappend = self.toappend[ind2point[self.toappend] != -1]
    
    def append_point(self, point, ind2point):
        if len(self.subs) == 0: # if self.subs is empty
            self.subs = point.sub
        elif len(self.subs.shape) == 1:
            self.subs = np.stack((self.subs, point.sub), axis = 0)
        else:
            self.subs = np.concatenate((self.subs, [point.sub]), axis = 0)
        self.inds = np.append(self.inds, point.ind)
        self.precips = np.append(self.precips, point.precip)
        self.sizes = np.append(self.sizes, point.area)
        #self.size += 1
        ind2point[point.ind] = -1
        self.get_toappend(point, ind2point)
        return ind2point
    
    def get_center_sub(self):
        self.center_sub = np.mean(self.subs, axis = 0)
    
    def get_cluster_power(self):
        self.cluster_power = sum(self.precips * self.sizes) / 1.e3 * 1.e3 / 24. / 1.4e6  # from mm/day*m^2 to GW (1.4e6 kg/h)
        return self.cluster_power

def find_cluster(p_ind_t, precip_t, area, **kwargs):
    
    # a nearest-neighbor cluster algorithm
    # area is the area of each grid point, the same shape as precip_t and p_ind_t
    
    from classes import cluster, grid_point
    n_x = p_ind_t.shape[1]
    n_y = p_ind_t.shape[0]

    if 'PERIODIC_X' in kwargs:
        PERIODIC_X = kwargs.get('PERIODIC_X')
    else:
        PERIODIC_X = True
    
    # define connectivity arrays: right, left, top, and bottom
    if PERIODIC_X:
        #print('periodic in x')
        ind_r = np.logical_and(p_ind_t, p_ind_t[:, np.append(np.array(range(n_x-1))+1, 0)]) # right
        ind_l = np.logical_and(p_ind_t, p_ind_t[:, np.append(n_x-1, np.array(range(n_x-1)))]) # left
    else:
        ind_r = np.append(np.logical_and(p_ind_t[:, range(n_x-1)], p_ind_t[:, np.array(range(n_x-1))+1]),
                          np.zeros((n_y, 1), dtype=bool), axis=1) # right
        ind_l = np.append(np.zeros((n_y, 1), dtype=bool), 
                          np.logical_and(p_ind_t[:, np.array(range(n_x-1))+1], p_ind_t[:, range(n_x-1)]), axis=1) # left

    ind_t = np.append(np.logical_and(p_ind_t[range(n_y-1), :], p_ind_t[np.array(range(n_y-1))+1, :]),
                      np.zeros((1, n_x), dtype=bool), axis=0) # top
    ind_b = np.append(np.zeros((1, n_x), dtype=bool),
                      np.logical_and(p_ind_t[np.array(range(n_y-1))+1, :], p_ind_t[range(n_y-1), :]), axis=0) # bottom
    
    # get whether there's bond or not for a grid point
    bonds = np.concatenate((ind_t.reshape(n_x*n_y, 1),
                            ind_b.reshape(n_x*n_y, 1),
                            ind_l.reshape(n_x*n_y, 1),
                            ind_r.reshape(n_x*n_y, 1)), axis = 1)
    Y, X = np.unravel_index(range(n_x*n_y), (n_y, n_x))

    points = []
    ind2point = np.zeros((n_x*n_y, ), dtype = 'int') - 1
    # initiate precipitating grid points indicated by p_ind_t
    for i in np.where(p_ind_t.reshape(-1,))[0]:
        temp_subs = np.array([[Y[i]+1, Y[i]-1, Y[i], Y[i]], [X[i], X[i], X[i]-1, X[i]+1]])
        # get the index for the points that this grid point is connected to
        temp_bonds_ind = (bonds[i, :]*2 - 1) * np.ravel_multi_index(temp_subs, (n_y, n_x), mode='wrap')
            # the y-direction boundary is dealt with properly because bonds[i, :] accounts for y boundaries
        temp_bonds_ind[temp_bonds_ind < 0] = -1
        points.append(grid_point(x = X[i], y = Y[i], ind = i,
                    bonds_ind = temp_bonds_ind, precip = precip_t[Y[i], X[i]], area = area[Y[i], X[i]]))
        ind2point[np.ravel_multi_index([Y[i], X[i]], (n_y, n_x))] = len(points) - 1

    clusters = []
    while any(ind2point != -1):
        # if there are still grid points not assigned to a cluster, start a new cluster
        new_cluster = cluster(n_x, n_y)
        temp_point = points[min(ind2point[np.where(ind2point != -1)[0]])]
        ind2point = new_cluster.append_point(temp_point, ind2point)
        c_finished = False
        while not c_finished:
            # if there are no further points to append, finish the current cluster
            if len(new_cluster.toappend) == 0:
                c_finished = True
            else:
                # get indices of points that need to be appended
                temp_points_ind = ind2point[new_cluster.toappend]
                for i in temp_points_ind:
                    ind2point = new_cluster.append_point(points[i], ind2point)

        clusters.append(new_cluster)
    return clusters

def sizes_of_cluster(p_ind, precip, area, **kwargs):

    Sizes = []
    Powers = []
    Nt = p_ind.shape[0]

    count = 0
    for t in range(Nt):
        count += 1
        if count % 10 == 0:
            print('t = {0:d} out of {1:d}'.format(count, Nt))
        clusters = find_cluster(p_ind[t, :, :], precip[t, :, :], area, **kwargs)
        sizes  = np.array([sum(cluster.sizes) for cluster in clusters])
        powers = np.array([cluster.get_cluster_power() for cluster in clusters])
        Sizes.append(sizes)
        Powers.append(powers)

    return Sizes, Powers

def size_struct_to_array(Sizes):
    sizes2 = np.zeros(np.int(1e9), dtype = 'float')
    ind_s1 = 0
    for t in range(len(Sizes)):
        len_1 = len(Sizes[t])
        sizes2[ind_s1 : ind_s1 + len_1] = Sizes[t]
        ind_s1 = ind_s1 + len_1
    sizes2 = sizes2[0:ind_s1]
    return sizes2[0:ind_s1]

class Reg:
    def __init__(self, coef_, intercept_):
        self.coef_ = np.array([[coef_]])
        self.intercept_ = intercept_

def power_law(sizes2, fit_start, fit_end, **kwargs):
    
    min_size, max_size = min(sizes2), max(sizes2)
    log10_range = np.ceil(np.log10(max_size/min_size) * 100) / 100
    temp_bins = np.power(10, np.linspace(np.log10(min_size), log10_range + np.log10(min_size), 25))
    temp_increment = np.ceil((temp_bins[1:] - temp_bins[0:-1]) / min_size)
    if 'bins' in kwargs:
        bins = kwargs.get('bins')
    else:
        bins = np.zeros((25,))
        bins[0] = min_size
        for i in range(len(bins) - 1):
            bins[i + 1] = bins[i] + min_size * temp_increment[i]
    bins_x = np.sqrt(bins[1:] / bins[0:-1]) * bins[0:-1]
    hist_count_0, _ = np.histogram(sizes2, bins = bins)

    if 'absolute_prob' in kwargs and kwargs.get('absolute_prob'):
        hist_count = hist_count_0 / (bins[1:] - bins[0:-1]) / kwargs.get('N')
    else:
        hist_count = hist_count_0 / (bins[1:] - bins[0:-1]) / float(sum(hist_count_0))

    if fit_end > len(bins) - 2:
        fit_end = len(bins) - 2
    
    # get fitting indices according to fit_start, fit_end, and the requirement that hist_count is none-zero
    fit_ind = np.zeros(bins_x.shape, dtype=bool)
    fit_ind[fit_start:fit_end] = True; fit_ind[hist_count == 0] = False
    
    # check if maximum likelihood estimation is specified, if not true, use linear regression
    if 'MLE' in kwargs and kwargs.get('MLE') == False:
        # do linear regression
        reg = LinearRegression().fit(np.log(bins_x[fit_ind]).reshape(-1, 1),
                                     np.log(hist_count[fit_ind]).reshape(-1, 1))
        R = 0; p = 0; sigma = 0
    else:
        if 'LIMITS' in kwargs and kwargs.get('LIMITS') == False:
            fit = powerlaw.Fit(sizes2)
        else:
            fit = powerlaw.Fit(sizes2, xmin = bins[fit_start], xmax = bins[fit_end+1])
        C = fit.power_law._pdf_continuous_normalizer * \
                float(sum(hist_count_0[fit_ind])) / float(sum(hist_count_0))
        R, p = fit.distribution_compare('power_law', 'lognormal')
        if 'absolute_prob' in kwargs and kwargs.get('absolute_prob'):
            C = C * float(sum(hist_count_0)) / kwargs.get('N')
        reg = Reg(-fit.alpha, np.log(C))
        sigma = fit.sigma

    return reg, hist_count, bins_x, bins, sigma, R, p

def power_law_plot(fig, plot_dict, **kwargs):
    
    markersize = 1
    if 'LATEX' in kwargs and kwargs.get('LATEX'):
        plt.rc('text', usetex=True)
        plt.rcParams["font.family"] = "Times New Roman"
    ax = fig.add_subplot(122)
    handles = []
    
    if 'LINE' in kwargs :
        LINE = kwargs.get('LINE')
    else:
        LINE = True
    
    for i in range(len(plot_dict.get('reg'))):
        
        reg     = plot_dict.get('reg')[i]
        color   = plot_dict.get('colors')[i]
        bins_x  = plot_dict.get('bins_x')[i]
        marker  = plot_dict.get('markers')[i]
        hist    = plot_dict.get('hist')[i]
        fit_start = plot_dict.get('fit_start')[i]
        fit_end = plot_dict.get('fit_end')[i]
        # plot binned data PDF
        p1 = ax.plot(bins_x, hist, '-' + marker,
                color = color, markersize = markersize,
                markeredgewidth = 0.9, markeredgecolor = color, markerfacecolor = 'none',
                label = plot_dict.get('run_legend')[i], linewidth = 1)
        
        if LINE:
            # plot regression line
            p2 = ax.plot((bins_x[fit_start], bins_x[fit_end]),
                    np.exp(reg.coef_[0, 0]*np.log((bins_x[fit_start], bins_x[fit_end])) + reg.intercept_),
                    '--', color = p1[0].get_color(), linewidth = 1)
            if 'LATEX' in kwargs and kwargs.get('LATEX'):
                p2[0].set_label(r'$n = {0:2.2f}$'.format(reg.coef_[0, 0]))
            else:
                p2[0].set_label(r'$n = {0:2.2f}$'.format(reg.coef_[0, 0]))
            '''
            ax.plot(bins_x[fit_start:fit_end], hist[fit_start:fit_end], marker,
                    color = p1[0].get_color(), markersize = markersize,
                    markeredgewidth=0.0)
            '''
            handles.append(p2[0])

    labels = [h.get_label() for h in handles]
    plt.xlabel(plot_dict.get('xlabel'))
    plt.ylabel(plot_dict.get('ylabel'))
    if 'ylim' in kwargs:
        plt.ylim(kwargs.get('ylim'))
    ax.legend(handles = handles, labels = labels, frameon = False, loc = 'lower left')
    import matplotlib.ticker as mtick
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
   
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.yscale('log')
    if 'xscale' in kwargs:
        plt.xscale(kwargs.get('xscale'))
    else:
        plt.xscale('log')
 
    # solution: https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large
    locmaj = mtick.LogLocator(base=10,numticks=15)
    ax.yaxis.set_major_locator(locmaj)
    locmin = mtick.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=15)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(mtick.NullFormatter())
    
    if 'title' in kwargs:
        plt.title(kwargs.get('title'))

def plot_dict_init(fig_name, xlabel, ylabel):
    
    plot_dict = {'reg': [], 'hist': [], 'bins_x': [], 'colors': [], 'markers': [], 'run_legend': [],
                'fit_start': [], 'fit_end': [],
                'xlabel': xlabel, 'ylabel': ylabel, 'filename': fig_name}
    
    return plot_dict

def plot_dict_update(plot_dict, reg, hist, bins_x, color, marker, legend, fit_start, fit_end):
    
    plot_dict['reg']       = plot_dict['reg']         + [reg]
    plot_dict['hist']      = plot_dict['hist']        + [hist]
    plot_dict['bins_x']    = plot_dict['bins_x']      + [bins_x]
    plot_dict['colors']    = plot_dict['colors']      + [color]
    plot_dict['markers']   = plot_dict['markers']     + [marker]
    plot_dict['run_legend']= plot_dict['run_legend']  + [legend]
    plot_dict['fit_start'] = plot_dict['fit_start']   + [fit_start]
    plot_dict['fit_end']   = plot_dict['fit_end']     + [fit_end]
    
    return plot_dict



