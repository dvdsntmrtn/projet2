import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy import stats

from varexpe import Expe, Variable

def plot_EBM(my_axis, my_EBM_name, datas_EBM, ybeg, my_color):
    y = datas_EBM[my_EBM_name, 'normalized_T']
    x = np.arange(ybeg, ybeg+len(y))
    my_axis.plot(x, y, linestyle='-', color=my_color)
    y = datas_EBM[my_EBM_name, 'normalized_T0']
    my_axis.plot(x, y, linestyle='--', color=my_color)

def plot_figure_0(figname, datas_EBM, figtype):
    fig, ax = plt.subplots()
    ax.set_xlim([1, 1000])

    ybeg = 1
    # plot_EBM(ax, 'abrupt-4xCO2', datas_EBM, ybeg, 'black')
    plot_EBM(ax, 'abrupt-2xCO2', datas_EBM, ybeg, 'black')
    plot_EBM(ax, 'stab-2xCO2-dab', datas_EBM, ybeg, 'red')
    plot_EBM(ax, 'expo-2xCO2', datas_EBM, ybeg, 'blue')

    fig.savefig(figname+'.'+figtype)


def plot_expe_var(my_axis, my_expe, datasets, varname, vartype, xmax, line=False, final_mean=False, start_year=0):
    x = np.arange(my_expe.ybeg, my_expe.yend+1)
    if len(datasets[(my_expe.expid(), varname, vartype)].shape) > 1:
        y = datasets[(my_expe.expid(), varname, vartype)][:, 0, 0].values
    else:
        y = datasets[(my_expe.expid(), varname, vartype)]
    if line:
        my_axis.plot(x, y, linestyle='-', color = my_expe.color) #, label=my_expe.name)
    else:
        my_axis.plot(x, y, marker = my_expe.marker, color = my_expe.color, label=my_expe.name, linewidth=1.0, markersize=2)
    if final_mean:
        x = np.arange(my_expe.ybeg, xmax)
        my_axis.plot(x, np.full_like(x, np.mean(y[-50:]), dtype=np.double), color = my_expe.color, linestyle='--', linewidth=0.7)
        # my_axis.axhline(y = np.mean(y[-50:]), xmin=1999, xmax=2800, color = my_expe.color, linestyle='--', linewidth=0.7)
    if start_year > 0:
        xx = np.array([my_expe.ybeg])
        yy = np.array([start_year])
        my_axis.plot(xx, yy, marker = '*', color = my_expe.color, markersize=12)
        # my_axis.axvline(x = my_expe.ybeg, ymin=0, ymax=start_year, color = my_expe.color, linestyle=':', linewidth=0.7)

def plot_figure_1(dictexpes, datasets, varname, vartype, figtype='png', plot_EBM = False):
    fig, ax = plt.subplots()
    
    ax.set_xlim([1850, 2800])

    xmin, xmax = plt.xlim()

    plot_expe_var(ax, dictexpes['CNRM-CM6-1_piControl_r1'], datasets, varname, vartype, xmax)
    plot_expe_var(ax, dictexpes['CNRM-CM6-1_abrupt-2xCO2_r1'], datasets, varname, vartype, xmax)    
    plot_expe_var(ax, dictexpes['CNRM-CM6-1_abrupt-4xCO2_r1'], datasets, varname, vartype, xmax)

    e_4x = dictexpes['CNRM-CM6-1_abrupt-4xCO2_r1']
    z_4x = datasets[(e_4x.expid(), varname, vartype)]
    exp = dictexpes['CNRM-CM6-1_stab-1p4xCO2-dab_r1']
    z_ = z_4x[exp.ybeg-1850, 0, 0]
    plot_expe_var(ax, dictexpes['CNRM-CM6-1_stab-1p4xCO2-dab_r1'], datasets, varname, vartype, xmax, final_mean=True, start_year=z_)
    exp = dictexpes['CNRM-CM6-1_stab-2xCO2-dab_r1']
    z_ = z_4x[exp.ybeg-1850, 0, 0]
    plot_expe_var(ax, dictexpes['CNRM-CM6-1_stab-2xCO2-dab_r1'], datasets, varname, vartype, xmax, final_mean=True, start_year=z_)
    exp = dictexpes['CNRM-CM6-1_stab-2p8xCO2-dab_r1']
    z_ = z_4x[exp.ybeg-1850, 0, 0]
    plot_expe_var(ax, dictexpes['CNRM-CM6-1_stab-2p8xCO2-dab_r1'], datasets, varname, vartype, xmax, final_mean=True, start_year=z_)

    #tmp plot_expe_var(ax, dictexpes['CNRM-CM6-1_stab-2xCO2-tab8x_r1'], datasets, varname, vartype, xmax)
    plot_expe_var(ax, dictexpes['CNRM-CM6-1_expo-2xCO2_r1'], datasets, varname, vartype, xmax)
    
    plot_expe_var(ax, dictexpes['EBM_abrupt-4xCO2_r1'], datasets, varname, vartype, xmax, line = True)

    if plot_EBM:
        plot_expe_var(ax, dictexpes['EBM_abrupt-4xCO2_r1'], datasets, varname, vartype, xmax, line = True)
        plot_expe_var(ax, dictexpes['EBM_abrupt-2xCO2_r1'], datasets, varname, vartype, xmax, line = True)
        plot_expe_var(ax, dictexpes['EBM_stab-2xCO2-dab_r1'], datasets, varname, vartype, xmax, line = True)
        plot_expe_var(ax, dictexpes['EBM_stab-1p4xCO2-dab_r1'], datasets, varname, vartype, xmax, line = True)
        plot_expe_var(ax, dictexpes['EBM_stab-2p8xCO2-dab_r1'], datasets, varname, vartype, xmax, line = True)
        plot_expe_var(ax, dictexpes['EBM_expo-2xCO2_r1'], datasets, varname, vartype, xmax, line = True)
        plot_expe_var(ax, dictexpes['EBM_stab-2xCO2-tab8x_r1'], datasets, varname, vartype, xmax, line = True)
    
    # Add 0_axis
    ax.axhline(y=0, color='k', linewidth=0.5)
    #
    plt.xlabel('years')
    if varname == 'tas_Amon':
        plt.ylabel('delta tas (K)')
        ax.legend(loc='lower right', fontsize=8)
    elif varname == 'rnet':
        plt.ylabel('delta net TOA (W/m2)')
        ax.legend(loc='upper right', fontsize=8)
    # ax.grid(True)
    fig.savefig('evol_'+varname+'.'+figtype)


def plot_figure_1b(dictexpes, datasets, figtype='png'):
    
    fig, ax = plt.subplots()
    
    # plot_with_linregress(ax, dictexpes['CNRM-CM6-1_abrupt-2xCO2_r1'], datasets, 'k', 0.4, 10)
    # plot_with_linregress(ax, dictexpes['CNRM-CM6-1_abrupt-4xCO2_r1'], datasets, 'k', 0.8, 10)
    
    plot_scatters(ax, dictexpes['CNRM-CM6-1_abrupt-2xCO2_r1'], datasets)
    plot_scatters(ax, dictexpes['CNRM-CM6-1_abrupt-4xCO2_r1'], datasets)
    plot_scatters(ax, dictexpes['CNRM-CM6-1_stab-2p8xCO2-dab_r1'], datasets)
    plot_scatters(ax, dictexpes['CNRM-CM6-1_stab-1p4xCO2-dab_r1'], datasets)
    plot_scatters(ax, dictexpes['CNRM-CM6-1_stab-2xCO2-dab_r1'], datasets)
    plot_scatters(ax, dictexpes['CNRM-CM6-1_expo-2xCO2_r1'], datasets)
    #tmp plot_scatters(ax, dictexpes['CNRM-CM6-1_stab-2xCO2-tab8x_r1'], datasets)

    plot_line(ax, dictexpes['EBM_abrupt-2xCO2_r1'], datasets)
    plot_line(ax, dictexpes['EBM_abrupt-4xCO2_r1'], datasets)
    plot_line(ax, dictexpes['EBM_stab-2xCO2-dab_r1'], datasets)
    plot_line(ax, dictexpes['EBM_stab-2p8xCO2-dab_r1'], datasets)
    plot_line(ax, dictexpes['EBM_stab-1p4xCO2-dab_r1'], datasets)
    plot_line(ax, dictexpes['EBM_expo-2xCO2_r1'], datasets)
    #tmp plot_line(ax, dictexpes['EBM_stab-2xCO2-tab8x_r1'], datasets)

    ymin = -4
    ymax = 8
    ax.set_ylim([ymin,ymax])
    
    # Add 0_axis
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)
    # ax.grid(True)
    plt.xlabel('delta tas (K)')
    plt.ylabel('delta net TOA (W/m2)')
    ax.legend(loc='upper right', fontsize=8)
    fig.savefig('NT_plot'+'.'+figtype)

def plot_zonmean(my_axis, my_expe, datasets, varname, vartype, last_years, dashed=False):  
    ilats = np.arange(0,128)    
    d = datasets[(my_expe.expid(), varname, vartype)][:,ilats,:]
    x = d[:,ilats,:]['lat']
    # t = np.mean(d[-20:,ilats,0].values, axis=0)
    t = np.mean(d[-last_years:,ilats,0].values, axis=0)
    if dashed:
        my_axis.plot(x, t / np.mean(t), marker = my_expe.marker, color = my_expe.color, label=my_expe.name+' : '+'last-'+str(last_years)+'y', linestyle='--')
    else:
        my_axis.plot(x, t / np.mean(t), marker = my_expe.marker, color = my_expe.color, label=my_expe.name+' : '+'last-'+str(last_years)+'y')

def plot_zonmean_between(my_axis, my_expe, datasets, varname, vartype, year0, year1, dashed=False):  
    ilats = np.arange(0,128)    
    d = datasets[(my_expe.expid(), varname, vartype)][:,ilats,:]
    x = d[:,ilats,:]['lat']
    t = np.mean(d[year0:year1,ilats,0].values, axis=0)
    if dashed:
        # my_axis.semilogy(x, t / np.mean(t), marker = my_expe.marker, color = my_expe.color, label=str(year0)+'-'+str(year1), linestyle='--')
        my_axis.plot(x, t / np.mean(t), marker = my_expe.marker, color = my_expe.color, label=my_expe.name+' : '+str(year0)+'-'+str(year1), linestyle='--')
    else:
        # my_axis.semilogy(x, t / np.mean(t), marker = my_expe.marker, color = my_expe.color, label='y'+str(year0)+'-y'+str(year1))
        my_axis.plot(x, t / np.mean(t), marker = my_expe.marker, color = my_expe.color, label=my_expe.name+' : y'+str(year0)+'-y'+str(year1))


def plot_figure_3(dictexpes, datasets, figtype='png'):
    fig, ax = plt.subplots()
    varname = 'tas_Amon'
    vartype = 'anom'
    
    # plot_zonmean(ax, dictexpes['CNRM-CM6-1_piControl_r1'], datasets, varname, vartype)
    plot_zonmean_between(ax, dictexpes['CNRM-CM6-1_abrupt-4xCO2_r1'], datasets, varname, vartype, 50, 150)
    plot_zonmean(ax, dictexpes['CNRM-CM6-1_abrupt-4xCO2_r1'], datasets, varname, vartype, 50, dashed=True)
    plot_zonmean_between(ax, dictexpes['CNRM-CM6-1_abrupt-2xCO2_r1'], datasets, varname, vartype, 50, 150)
    #tmp plot_zonmean(ax, dictexpes['CNRM-CM6-1_abrupt-2xCO2_r1'], datasets, varname, vartype, 50, dashed=True)
    plot_zonmean(ax, dictexpes['CNRM-CM6-1_stab-2xCO2-dab_r1'], datasets, varname, vartype, 50)
    plot_zonmean(ax, dictexpes['CNRM-CM6-1_stab-1p4xCO2-dab_r1'], datasets, varname, vartype, 50)
    plot_zonmean(ax, dictexpes['CNRM-CM6-1_stab-2p8xCO2-dab_r1'], datasets, varname, vartype, 50)
    #tmp plot_zonmean(ax, dictexpes['CNRM-CM6-1_stab-2xCO2-tab8x_r1'], datasets, varname, vartype, 50)
    #tmp plot_zonmean(ax, dictexpes['CNRM-CM6-1_expo-2xCO2_r1'], datasets, varname, vartype, 50)
    
    plt.ylabel('delta tas / ECS')
    plt.xlabel('latitude')
    ax.legend(loc='upper left', fontsize=8)
    
    fig.savefig('figure3'+'.'+figtype)



def plot_with_linregress(my_ax, datasets, expid, my_color, my_alpha, my_xmax):
    x = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
    y = datasets[(expid, 'rnet', 'anom')][:, 0, 0].values
    my_ax.plot(x, y, '.', color=my_color, alpha=my_alpha)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    xx = np.arange(0, my_xmax)
    yy = slope*xx + intercept
    my_ax.plot(xx, yy, linestyle='-', color=my_color, alpha=my_alpha)

def scatters(my_ax, datasets, expid, my_color, my_alpha, my_marker_size):
    if len(datasets[(expid, 'tas_Amon', 'anom')].shape) > 1:
        x = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
        y = datasets[(expid, 'rnet', 'anom')][:, 0, 0].values
    else:
        x = datasets[(expid, 'tas_Amon', 'anom')]
        y = datasets[(expid, 'rnet', 'anom')]
    my_ax.plot(x, y, 'x', color=my_color, alpha=my_alpha, markersize=my_marker_size)

def scatters_by_decade(my_ax, datasets, expid, my_color, my_alpha, my_label):
    if len(datasets[(expid, 'tas_Amon', 'anom')].shape) > 1:
        x = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
        y = datasets[(expid, 'rnet', 'anom')][:, 0, 0].values        
    else:
        x = datasets[(expid, 'tas_Amon', 'anom')]
        y = datasets[(expid, 'rnet', 'anom')]
    xx = np.mean(np.reshape(x[:10*np.int(np.size(x)/10)], [np.int(np.size(x)/10), 10]), axis=1)
    yy = np.mean(np.reshape(y[:10*np.int(np.size(y)/10)], [np.int(np.size(y)/10), 10]), axis=1)
    my_ax.plot(xx, yy, 'o', color=my_color, alpha=my_alpha, label=my_label)

def F_and_lbda_from_piClim(datasets, expid_ctl, expid_pert, ny):
    x0 = np.mean(datasets[(expid_ctl, 'tas_Amon', 'anom')][:ny,0,0].values, axis=0)
    y0 = np.mean(datasets[(expid_ctl, 'rnet', 'anom')][:ny,0,0].values, axis=0)
    x1 = np.mean(datasets[(expid_pert, 'tas_Amon', 'anom')][:ny,0,0].values, axis=0)
    y1 = np.mean(datasets[(expid_pert, 'rnet', 'anom')][:ny,0,0].values, axis=0)
    lbda = - (y1 - y0) / (x1 - x0)
    F = y0 + lbda * x0
    return F, lbda

def plot_scatters(my_axis, my_expe, my_datasets):
    scatters_by_decade(my_axis, my_datasets, my_expe.expid(), my_expe.color, 1.0, my_expe.name)    
    scatters(my_axis, my_datasets, my_expe.expid(), my_expe.color, 0.3, 2)

def plot_line(my_axis, my_expe, my_datasets):
    expid = my_expe.expid()
    if len(my_datasets[(expid, 'tas_Amon', 'anom')].shape) > 1:
        x = my_datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
        y = my_datasets[(expid, 'rnet', 'anom')][:, 0, 0].values        
    else:
        x = my_datasets[(expid, 'tas_Amon', 'anom')]
        y = my_datasets[(expid, 'rnet', 'anom')]
    my_axis.plot(x, y, linestyle='-', color=my_expe.color)
