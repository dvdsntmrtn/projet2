from climaf.api import *
import xarray as xr
import numpy as np
from varexpe import Expe, Variable

# Useless ?
def avg_var(modvar) :
    modvar_spcavg = space_average(modvar)
    serie = ccdo(modvar_spcavg,operator='yearavg')
    return serie

def avg_zon(x) :
    u = ccdo(x, operator='zonmean')
    y = ccdo(u, operator='yearavg')
    return y


def define_my_CLIMAF_projects():
    '''
    Define new CLIMAF's project by setting :
      - attributes
      - data localization
    '''
    name = 'myCMIP6'
    cproject(name, 'model', 'experiment', 'member', 'variable', 'table')
    rootdirs = []
    rootdirs.append('/cnrm/amacs/USERS/stmartin/NO_SAVE/CMIP6/sorties/*/${model}_${experiment}_r${member}i1p1*')
    file_pattern = '${variable}_${table}_${model}_${experiment}_r${member}i1p1*_YYYYMM-YYYYMM.nc'
    dataloc(project = name, organization = 'generic', url = [rdir+'/'+file_pattern for rdir in rootdirs])

    name = 'CPL6214'
    cproject(name, 'model', 'experiment', 'member', 'variable', 'table', 'gridtype')
    rootdir = '/cnrm/amacs/USERS/stmartin/data3/sorties/PRE6'
    locadirs = ['${model}_${experiment}/X/monthly/${model}_${experiment}_arpsfx_monthly_${variable}_YYYY-YYYY.nc', '${model}_${experiment}/X/${model}_${experiment}_arpsfx_monthly_${variable}_YYYY-YYYY.nc']
    dataloc(project = name, organization = 'generic', url=[rootdir+'/'+ldir for ldir in locadirs])

    # temporary ?
    name = 'tmpCMIP6'
    cproject(name, 'model', 'experiment', 'member', 'variable', 'table')
    rootdirs = []
    rootdirs.append('/home/stmartin/work/lxamacs/data')
    rootdirs.append('/cnrm/amacs/USERS/stmartin/data3/sorties/STAB/CLIMAF')
    file_pattern = '${model}_${experiment}_r${member}_${variable}_${table}.gmean.annual.nc'    
    dataloc(project = name, organization = 'generic', url = [rdir+'/'+file_pattern for rdir in rootdirs])
    
    name = 'zonCMIP6'
    cproject(name, 'model', 'experiment', 'member', 'variable', 'table')
    rootdirs = []
    rootdirs.append('/home/stmartin/work/lxamacs/data')
    rootdirs.append('/cnrm/amacs/USERS/stmartin/data3/sorties/STAB/CLIMAF')
    file_pattern = '${model}_${experiment}_r${member}_${variable}_${table}.zmean.annual.nc'    
    dataloc(project = name, organization = 'generic', url = [rdir+'/'+file_pattern for rdir in rootdirs])


def dict_expes_stab_article(project_name):
    '''
    Define the specific dict of Expe-s for stabilization study
    '''
    dict_allexpes = {}
    eCTL = Expe(project=project_name, model='CNRM-CM6-1', name='piControl', member=1, ybeg=1850, yend=2349, marker=',', color='silver')
    eCTL.expe_control = eCTL  
    dict_allexpes[eCTL.expid()] = eCTL
    e = Expe(project=project_name, model='CNRM-CM6-1', name='abrupt-4xCO2', member=1, ybeg=1850, yend=2762, expe_control=eCTL, marker=',', color='dimgray')
    dict_allexpes[e.expid()] = e
    e = Expe(project=project_name, model='CNRM-CM6-1', name='abrupt-2xCO2', member=1, ybeg=1850, yend=2189, expe_control=eCTL, marker=',', color='black')
    dict_allexpes[e.expid()] = e
    e = Expe(project=project_name, model='CNRM-CM6-1', name='stab-1p4xCO2-dab', member=1, ybeg=1969, yend=2203, expe_control=eCTL, marker='.', color='purple')
    dict_allexpes[e.expid()] = e
    e = Expe(project=project_name, model='CNRM-CM6-1', name='stab-2xCO2-dab', member=1, ybeg=2137, yend=2277, expe_control=eCTL, marker='.', color='orangered')
    dict_allexpes[e.expid()] = e
    e = Expe(project=project_name, model='CNRM-CM6-1', name='stab-2p8xCO2-dab', member=1, ybeg=2424, yend=2576, expe_control=eCTL, marker='.', color='orange')
    dict_allexpes[e.expid()] = e
    e = Expe(project=project_name, model='CNRM-CM6-1', name='expo-2xCO2', member=1, ybeg=1850, yend=2017, expe_control=eCTL, marker='.', color='dodgerblue')
    dict_allexpes[e.expid()] = e
    e = Expe(project=project_name, model='CNRM-CM6-1', name='stab-2xCO2-tab8x', member=1, ybeg=2000, yend=2033, expe_control=eCTL, marker='.', color='green')
    dict_allexpes[e.expid()] = e
    return dict_allexpes


def dict_vars_T():
    '''
    Define the specific dict of Variable-s for N-T plot study
    '''
    dict_vars = {}
    v = Variable(name='tas', table='Amon')
    dict_vars[v.name] = v
    return dict_vars


def dict_vars_NT():
    '''
    Define the specific dict of Variable-s for N-T plot study
    '''
    dict_vars = {}
    v = Variable(name='tas', table='Amon')
    dict_vars[v.name] = v
    v = Variable(name='rsdt', table='Amon')
    dict_vars[v.name] = v
    v = Variable(name='rsut', table='Amon')
    dict_vars[v.name] = v
    v = Variable(name='rlut', table='Amon')
    dict_vars[v.name] = v
    return dict_vars


def dict_vars_heatc():
    '''
    Define the specific dict of Variable-s for heat contents
    '''
    dict_vars = {}
    v = Variable(name='hc2000', table='HOMOImon')
    dict_vars[v.name] = v
    v = Variable(name='hc700', table='HOMOImon')
    dict_vars[v.name] = v
    v = Variable(name='heatc', table='HOMOImon')
    dict_vars[v.name] = v
    v = Variable(name='hcont300', table='Emon')
    dict_vars[v.name] = v
    return dict_vars


def load_datas(dictexpes, dictvars, verbose=False):
    '''
    Get data from : 
    - a specific dict of Expe-s : dictexpes
    - a specific dict of Var-s : dictvars
    > Return a dictionary with a 3-element tuple (exp.expid(), varname, vartype='brut/mean/anom') keys
    '''
    datasets = {}
    # -- read datas
    for exp in dictexpes.values():
        if verbose:
            print exp.name, exp.member
        for var in dictvars.values():
            f = ds(project=exp.project, variable=var.name, table=var.table, gridtype=var.grid, model=exp.model, experiment=exp.name, member=exp.member, period=exp.period())
            if f.listfiles() is not None:
                xds = xr.open_dataset(cfile(f))
                datasets[(exp.expid(), var.varid(), 'brut')] = xds[var.name]
                if verbose:
                    print(f.listfiles())
            else:
                print var.name, ' not in ', exp.name
    # -- add 'temporal' mean value
    for exp in dictexpes.values():    
        for var in dictvars.values():
            if (datasets.has_key((exp.expid(), var.varid(), 'brut'))):
                datasets[(exp.expid(), var.varid(), 'mean')] = np.mean(datasets[(exp.expid(), var.varid(), 'brut')], axis=0)
    # -- add anomaly-from-the-control-temporal-mean value
    # -- add new var 'rnet' from rsdt, rsut and rlut values
    for exp in dictexpes.values():        
        for var in dictvars.values():
            if (datasets.has_key((exp.expid(), var.varid(), 'brut'))):
                datasets[(exp.expid(), var.varid(), 'anom')] = datasets[(exp.expid(), var.varid(), 'brut')] - datasets[(exp.expe_control.expid(), var.varid(), 'mean')]
        for x in ['brut', 'mean', 'anom']:
            if (datasets.has_key((exp.expid(), 'rsdt_Amon', 'brut'))):
                datasets[(exp.expid(), 'rnet', x)] = datasets[(exp.expid(), 'rsdt_Amon', x)] - datasets[(exp.expid(), 'rsut_Amon', x)] - datasets[(exp.expid(), 'rlut_Amon', x)]
    return datasets
