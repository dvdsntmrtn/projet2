# import os
# os.environ['CLIMAF_CACHE'] = 'chemin_temporaire'

from climaf.api import *
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import xarray as xr
import numpy as np
from scipy import stats

from varexpe import Expe, Variable
from ebm import *
from auxilliary import *
from auxilliary_plot import *

#-- Define specific CLIMAF projects
define_my_CLIMAF_projects()

# -- Define dict of expes and of vars
my_dictexps_1 = dict_expes_stab_article('tmpCMIP6')
my_dictexps_2 = dict_expes_stab_article('zonCMIP6')
my_dictvars = dict_vars_NT()

dict_vars_tas = {}
v = Variable(name='tas', table='Amon')
dict_vars_tas[v.name] = v

# -- Load data
datasets_1 = load_datas(my_dictexps_1, my_dictvars, verbose=True)
datasets_2 = load_datas(my_dictexps_2, dict_vars_tas, verbose=True)

# -- EBM
# T = {}
# N = {}
expid = 'CNRM-CM6-1_abrupt-4xCO2_r1'
n_y = 149
TT = datasets_1[(expid, 'tas_Amon', 'anom')][:n_y, 0, 0].values
NN = datasets_1[(expid, 'rnet', 'anom')][:n_y, 0, 0].values
EBM_4x = derive_EBM_I(TT, NN, 4.0, nyr_s=50)
print(EBM_4x)

y_deb = 1850
n_years = 1000

datas_EBM = {}
EBM_academic = EBM()
print('-----------------------------------')
print(EBM_academic)

def add_EBM_datasets(my_dictexps, my_datasets,my_name, my_EBM, my_forcage, y_deb, n_years, my_color):
    output_EBM = analytical_EBM(my_EBM, my_forcage, n_years, lprint = False)
    e = Expe(project='', model='EBM', name=my_name, ybeg=y_deb, yend=y_deb+n_years, color=my_color)
    my_dictexps[e.expid()] = e
    my_datasets[e.expid(), 'tas_Amon', 'anom'] = output_EBM['T']
    my_datasets[e.expid(), 'rnet', 'anom'] = output_EBM['N']

    

forcage_ab4x = FORCING(typ='abrupt', xCO2_infty = 4.0)
add_EBM_datasets(my_dictexps_1, datasets_1, 'abrupt-4xCO2', EBM_4x, forcage_ab4x, y_deb, n_years, my_dictexps_1['CNRM-CM6-1_abrupt-4xCO2_r1'].color)

forcage_ab2x = FORCING(typ='abrupt', xCO2_infty = 2.0)
add_EBM_datasets(my_dictexps_1, datasets_1, 'abrupt-2xCO2', EBM_4x, forcage_ab2x, y_deb, n_years, my_dictexps_1['CNRM-CM6-1_abrupt-2xCO2_r1'].color)

forcage_dab2x = FORCING(typ='double_abrupt', xCO2_0 = 4.0, t_0 = 287, xCO2_infty = 2.0)
add_EBM_datasets(my_dictexps_1, datasets_1, 'stab-2xCO2-dab', EBM_4x, forcage_dab2x, y_deb, n_years, my_dictexps_1['CNRM-CM6-1_stab-2xCO2-dab_r1'].color)

forcage_dab1p4x = FORCING(typ='double_abrupt', xCO2_0 = 4.0, t_0 = 119, xCO2_infty = 1.4)
add_EBM_datasets(my_dictexps_1, datasets_1, 'stab-1p4xCO2-dab', EBM_4x, forcage_dab1p4x, y_deb, n_years, my_dictexps_1['CNRM-CM6-1_stab-1p4xCO2-dab_r1'].color)

forcage_dab2p8x = FORCING(typ='double_abrupt', xCO2_0 = 4.0, t_0 = 574, xCO2_infty = 2.8)
add_EBM_datasets(my_dictexps_1, datasets_1, 'stab-2p8xCO2-dab', EBM_4x, forcage_dab2p8x, y_deb, n_years, my_dictexps_1['CNRM-CM6-1_stab-2p8xCO2-dab_r1'].color)

def expo_EBM(my_EBM, my_tau_e, x_infty, n_y):
    optim_x0 = optimal_x0(my_tau_e, my_EBM.parameters()['tau_s'], x_infty)
    forcage = FORCING(typ='expo', xCO2_0 = optim_x0, xCO2_infty = x_infty, tau_e = my_tau_e)
    return analytical_EBM(my_EBM, forcage, n_y)

def optimal_expo_EBM(my_EBM, x_infty, n_y):
    optim_x0 = optimal_x0(my_EBM.parameters()['tau_e'], my_EBM.parameters()['tau_s'], x_infty)
    forcage = FORCING(typ='expo', xCO2_0 = optim_x0, xCO2_infty = x_infty, tau_e = my_EBM.parameters()['tau_e'])
    return analytical_EBM(my_EBM, forcage, n_y)

x_infty = 2.0
my_EBM = EBM_4x
optim_x0 = optimal_x0(my_EBM.parameters()['tau_e'], my_EBM.parameters()['tau_s'], x_infty)
print(optim_x0)
forcage_expo2x = FORCING(typ='expo', xCO2_0 = optim_x0, xCO2_infty = x_infty, tau_e = my_EBM.parameters()['tau_e'])
add_EBM_datasets(my_dictexps_1, datasets_1, 'expo-2xCO2', EBM_4x, forcage_expo2x, y_deb, n_years, my_dictexps_1['CNRM-CM6-1_expo-2xCO2_r1'].color)

forcage_2x_tab8x = FORCING(typ='triple_abrupt', xCO2_0 = 4.0, t_0 = 150, xCO2_m = 8.0, t_m = 150+74, xCO2_infty = 2.0)
add_EBM_datasets(my_dictexps_1, datasets_1, 'stab-2xCO2-tab8x', EBM_4x, forcage_2x_tab8x, y_deb, n_years, my_dictexps_1['CNRM-CM6-1_stab-2xCO2-tab8x_r1'].color)
 
# T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
# N[expid] = datasets[(expid, 'rnet', 'anom')][:, 0, 0].values

# exit()


def optimal_dab_EBM(my_EBM, x0, x_infty, n_y):
    optim_t_0 = optimal_t0(my_EBM.parameters()['tau_s'], x0, x_infty)
    print('optimal_t_0 (%6.2f, %6.2f) = %4.0f'%(x0,x_infty,optim_t_0))
    forcage = FORCING(typ='double_abrupt', xCO2_0=x0, t_0=optim_t_0, xCO2_infty=x_infty)
    return analytical_EBM(my_EBM, forcage, n_y)

def add_datas_EBM(my_datas_EBM, my_name, my_analytical_EBM):
    my_datas_EBM[my_name, 'normalized_T'] = my_analytical_EBM['normalized_T']
    my_datas_EBM[my_name, 'normalized_T0'] = my_analytical_EBM['normalized_T0']

add_datas_EBM(datas_EBM, 'abrupt-2xCO2', analytical_EBM(EBM_academic, forcage_ab2x, n_years))
add_datas_EBM(datas_EBM, 'abrupt-4xCO2', analytical_EBM(EBM_academic, forcage_ab4x, n_years))
add_datas_EBM(datas_EBM, 'stab-2xCO2-dab', optimal_dab_EBM(EBM_academic, 4.0, 2.0, n_years))
add_datas_EBM(datas_EBM, 'expo-2xCO2', optimal_expo_EBM(EBM_academic, 2.0, n_years))

plot_figure_0('toto', datas_EBM, 'png')

exit()

plot_figure_1(my_dictexps_1, datasets_1, 'tas_Amon', 'anom', 'jpg')
plot_figure_1(my_dictexps_1, datasets_1, 'rnet', 'anom', 'jpg')

plot_figure_1b(my_dictexps_1, datasets_1, 'jpg')

plot_figure_3(my_dictexps_2, datasets_2, 'jpg')

exit()

plot_figure_1(my_dictexps, datasets, 'hc2000_HOMOImon', 'anom')
plot_figure_1(my_project, datasets, 'hc700_HOMOImon', 'anom')
plot_figure_1(my_project, datasets, 'heatc_HOMOImon', 'anom')
plot_figure_1(my_project, datasets, 'hcont3000_Emon', 'anom')

exit()



exit()

ECS = {}
T = {}
N = {}

expid = 'CNRM-CM6-1_abrupt-4xCO2_r1'
TT = datasets[(expid, 'tas_Amon', 'anom')][:149, 0, 0].values
NN = datasets[(expid, 'rnet', 'anom')][:149, 0, 0].values
EBM_4x_I = derive_EBM_I(TT, NN, 4.0, nyr_s=50)
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
N[expid] = datasets[(expid, 'rnet', 'anom')][:, 0, 0].values

print(EBM_4x_I)

t_0 = optimal_t0(EBM_4x_I.parameters()['tau_s'], 4.0, 1.5)
print('t_0_1p5x = %4.0f'%t_0)

t_0 = optimal_t0(EBM_4x_I.parameters()['tau_s'], 8.0, 2.0)
print('t_0_2x = %4.0f'%t_0)

t_0_4x = optimal_t0(EBM_4x_I.parameters()['tau_s'], 8.0, 4.0)
print('t_0_4x = %4.0f'%t_0_4x)

def optimal_dab_EBM(my_EBM, x0, x_infty, n_y):
    optim_t_0 = optimal_t0(my_EBM.parameters()['tau_s'], x0, x_infty)
    print('optimal_t_0 (%6.2f, %6.2f) = %4.0f'%(x0,x_infty,optim_t_0))
    forcage = FORCING(typ='double_abrupt', xCO2_0=x0, t_0=optim_t_0, xCO2_infty=x_infty)
    return analytical_EBM(my_EBM, forcage, n_y)

def optimal_dab2_EBM(my_EBM, x0, x_infty, n_y):
    optim_t_0 = 50
    print('___optimal_t_0 (%6.2f, %6.2f) = %4.0f'%(x0,x_infty,optim_t_0))
    forcage = FORCING(typ='double_abrupt', xCO2_0=x0, t_0=optim_t_0, xCO2_infty=x_infty)
    return analytical_EBM(my_EBM, forcage, n_y)

def optimal_expo_EBM(my_EBM, x_infty, n_y):
    optim_x0 = optimal_x0(my_EBM.parameters()['tau_e'], my_EBM.parameters()['tau_s'], x_infty)
    forcage = FORCING(typ='expo', xCO2_0 = optim_x0, xCO2_infty = x_infty, tau_e = my_EBM.parameters()['tau_e'])
    return analytical_EBM(my_EBM, forcage, n_y)

def expo_EBM(my_EBM, my_tau_e, x_infty, n_y):
    optim_x0 = optimal_x0(my_tau_e, my_EBM.parameters()['tau_s'], x_infty)
    forcage = FORCING(typ='expo', xCO2_0 = optim_x0, xCO2_infty = x_infty, tau_e = my_tau_e)
    return analytical_EBM(my_EBM, forcage, n_y)


exit()


expid = 'CNRM-CM6-1_stab-1p4xCO2-dab_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[1969-1850:1969-1850+len(y)]
ax.plot(x, y, 'o', color='m', alpha=0.8)


nyy = 50
ttot = 300

fig, ax = plt.subplots()

years = np.arange(1850, 5000)

expid = 'CNRM-CM6-1_abrupt-4xCO2_r1'
y = T[expid]
x = years[0:len(y)]
ax.plot(x, y, 'o', color='k', alpha=0.8)

T149 = y[149]
xx = 4.0

my_forcage_ab4x = FORCING(typ='abrupt', xCO2_infty=4.0)
output_ab4x_I = analytical_EBM(EBM_4x_I, my_forcage_ab4x, ttot, lprint=False)

expid = 'CNRM-CM6-1_stab-1p4xCO2-dab_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[1969-1850:1969-1850+len(y)]
ax.plot(x, y, 'o', color='m', alpha=0.8)

expid = 'CNRM-CM6-1_stab-2p8xCO2-dab_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[2424-1850:2424-1850+len(y)]
ax.plot(x, y, 'o', color='m', alpha=0.8)

expid = 'CNRM-CM6-1_stab-2xCO2-dab_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[2137-1850:2137-1850+len(y)]
ax.plot(x, y, 'o', color='m', alpha=0.8)

expid = 'CNRM-CM6-1_expo-4xCO2_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[0:len(y)]
##### 
ax.plot(x, y, 'o', color='m', alpha=0.8)

expid = 'CNRM-CM6-1_expo-2xCO2_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[0:len(y)]
##### ax.plot(x, y, 'o', color='m', alpha=0.8)

expid = 'CNRM-CM6-1_stab-2xCO2-tab8x_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[2000-1850:2000-1850+len(y)]
ax.plot(x, y, '.', color='k', alpha=0.8)

expid = 'CNRM-CM6-1_stab-4xCO2-tab16x_r1'
T[expid] = datasets[(expid, 'tas_Amon', 'anom')][:, 0, 0].values
y = T[expid]
x = years[2000-1850:2000-1850+len(y)]
ax.plot(x, y, '+', color='k', alpha=0.8)

my_forcage_ab4x = FORCING(typ='abrupt', xCO2_infty=4.0)
output_ab4x_I = analytical_EBM(EBM_4x_I, my_forcage_ab4x, ttot, lprint=False)
ECS = EBM_4x_I.F / EBM_4x_I.lbda
print('ECS = %5.2f'%ECS)
y = output_ab4x_I['T']
print(1+np.argmax(y/ECS > 0.99))
print(1+np.argmax(y/ECS > 0.95))
x = years[0:len(y)]
#_ plt.plot(x, y, linewidth=2, linestyle='-', color='r', alpha=0.4)


my_forcage_ab8x = FORCING(typ='abrupt', xCO2_infty=2.0)
output_ab8x_I = analytical_EBM(EBM_4x_I, my_forcage_ab8x, ttot, lprint=False)
y = output_ab8x_I['T']
x = years[0:len(y)]
#_ 
#_ plt.plot(x, y, linestyle='-', color='g', alpha=0.4)

_x_0 = 4.0
_t_0 = 150
_t_m = 150 + 51
_x_infty = 2.0
_tau_s = EBM_4x_I.parameters()['tau_s']
_x_m = np.exp((np.log(_x_0)*(1-np.exp(_t_0/_tau_s))+np.log(_x_infty)*np.exp(_t_m/_tau_s)) / (np.exp(_t_m/_tau_s) - np.exp(_t_0/_tau_s)))
print('-----> x_m =%5.2f'%_x_m)

_x_infty = 2.0
_x_m = 16.0
__t_m = _tau_s * np.log((np.log(_x_0) - (np.log(_x_0)-np.log(_x_m))*np.exp(_t_0/_tau_s))/(np.log(_x_m) - np.log(_x_infty)))
print('<<<<----------------------------<<<<')
print('-----> t_m =%5.2f'%__t_m)

# exit()


my_forcage_test = FORCING(typ='triple_abrupt', xCO2_0=_x_0, t_0=_t_0, xCO2_m=_x_m, t_m=_t_m, xCO2_infty=_x_infty)
output_ab8x_I = analytical_EBM(EBM_4x_I, my_forcage_test, ttot, lprint=True)
y = output_ab8x_I['T']
x = years[0:len(y)]
##### 
plt.plot(x, y, linestyle='-', color='r', alpha=0.8)
# plt.plot(x, output_ab8x_I['T0'], linestyle='--', color='r', alpha=0.3)

_x_m = 8.0
_t_m = 150 + 74
_x_infty = 2.0
my_forcage_test = FORCING(typ='triple_abrupt', xCO2_0=_x_0, t_0=_t_0, xCO2_m=_x_m, t_m=_t_m, xCO2_infty=_x_infty)
output_ab8x_I = analytical_EBM(EBM_4x_I, my_forcage_test, ttot, lprint=True)
y = output_ab8x_I['T']
x = years[0:len(y)]
##### 
plt.plot(x, y, linestyle='-', color='b', alpha=0.8)
# plt.plot(x, output_ab8x_I['T0'], linestyle='--', color='r', alpha=0.8)


y = optimal_dab_EBM(EBM_4x_I, 4, 2, ttot)['T']
x = years[0:len(y)]
##### plt.plot(x, y, linestyle='-', color='b', alpha=0.7)
print('mean =%5.2f'%np.mean(y[-nyy:]))

forcage = FORCING(typ='double_abrupt', xCO2_0=4, t_0=239, xCO2_infty=2)
y = analytical_EBM(EBM_4x_I, forcage, ttot)['T']
y0 = analytical_EBM(EBM_4x_I, forcage, ttot)['T0']
x = years[0:len(y)]
# plt.plot(x, y, linestyle='-', color='m', alpha=0.7)
# plt.plot(x, y0, linestyle='--', color='m', alpha=0.7)
print('mean =%5.2f'%np.mean(y[-nyy:]))

y = optimal_dab_EBM(EBM_4x_I, 4, np.sqrt(2), ttot)['T']
x = years[0:len(y)]
##### plt.plot(x, y, linestyle='-', color='b', alpha=0.7)
print('mean =%5.2f'%np.mean(y[-nyy:]))

# y = optimal_dab_EBM(EBM_4x_I, 4, 2*np.sqrt(2), ttot)['T']
# x = years[0:len(y)]
##### plt.plot(x, y, linestyle='-', color='b', alpha=0.7)
# print('mean =%5.2f'%np.mean(y[-nyy:]))

y = optimal_dab_EBM(EBM_4x_I, 4, 2, ttot)['T']
x = years[0:len(y)]
##### plt.plot(x, y, linestyle='-', color='b', alpha=0.7)
print('mean =%5.2f'%np.mean(y[-nyy:]))

y = optimal_expo_EBM(EBM_4x_I, 4, ttot)['T']
x = years[0:len(y)]
##### plt.plot(x, y, linestyle='-', color='r', alpha=0.5)
print('mean =%5.2f'%np.mean(y[-nyy:]))

y = optimal_expo_EBM(EBM_4x_I, 2, ttot)['T']
x = years[0:len(y)]
##### plt.plot(x, y, linestyle='-', color='r', alpha=0.5)
print('mean =%5.2f'%np.mean(y[-nyy:]))

y = optimal_expo_EBM(EBM_4x_I, 2.8, ttot)['T']
x = years[0:len(y)]
# plt.plot(x, y, linestyle='-', color='r', alpha=0.5)
print('mean =%5.2f'%np.mean(y[-nyy:]))

y = expo_EBM(EBM_4x_I, 150, 2.0, ttot)['T']
x = years[0:len(y)]
# plt.plot(x, y, linestyle='--', color='b')
print('mean =%5.2f'%np.mean(y[-nyy:]))

y = expo_EBM(EBM_4x_I, 100, 2.0, ttot)['T']
x = years[0:len(y)]
# plt.plot(x, y, linestyle='--', color='b')
print('mean =%5.2f'%np.mean(y[-nyy:]))

y = expo_EBM(EBM_4x_I, 75, 2.0, ttot)['T']
x = years[0:len(y)]
# plt.plot(x, y, linestyle='--', color='b')
print('mean =%5.2f'%np.mean(y[-nyy:]))

# Add 0_axis
ax.axhline(y=0, color='k', linewidth=0.5)
# ax.axvline(x=0, color='k', linewidth=0.5)

ax.set_xlim([1850,2100])

ax.grid(True)
fig.savefig('evol_T.png')



exit()
