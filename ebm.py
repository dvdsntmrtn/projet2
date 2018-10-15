# Version 1.0.1

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn import linear_model

#--------------------------------------------------------------------------------------------------
class EBM(object):

    def __init__(self, F=7.65, lbda=1.25, c=9.844, c_0=133.325, gam=0.70, epsi=1.0, xCO2=4):
        self.F = F
        self.lbda = lbda
        self.c = c
        self.c_0 = c_0
        self.gam = gam
        self.epsi = epsi
        self.xCO2 = xCO2
        
    def __str__(self):
        param = self.parameters()
        xstr = 'F = %5.2f W/m2; lbda = %5.2f W/m2/K : ECS = %5.2f K'%(self.F, self.lbda, self.F/self.lbda)
        xstr += '\n'
        xstr += 'c = %5.1f; c_0 = %5.1f'%(self.c, self.c_0)
        xstr += '\n'
        xstr += 'gam = %4.2f; epsi = %4.2f'%(self.gam, self.epsi)
        xstr += '\n'
        xstr += 'a_f = %4.2f; a_s = %4.2f'%(param['a_f'], param['a_s'])
        xstr += '\n'
        xstr += 'tau_f = %4.2f; tau_s = %4.2f; tau_e = %4.2f'%(param['tau_f'], param['tau_s'], param['tau_e'])
        xstr += '\n'
        xstr += 'phi_f = %4.2f; phi_s = %4.2f'%(param['phi_f'], param['phi_s'])
        return xstr
    
    def parameters(self):
        #
        gammprim = self.epsi * self.gam
        c0prim   = self.epsi * self.c_0
        #
        b     = (self.lbda + gammprim)/self.c + gammprim/c0prim
        bstar = (self.lbda + gammprim)/self.c - gammprim/c0prim
        delta = b*b - 4*(self.lbda*gammprim)/(self.c*c0prim)
        # Constants
        phi_f = self.c/(2*gammprim)*(bstar - np.sqrt(delta))
        phi_s = self.c/(2*gammprim)*(bstar + np.sqrt(delta))
        # Relaxation times
        tau_f = (self.c*c0prim)/(2*self.lbda*gammprim)*(b - np.sqrt(delta))
        tau_s = (self.c*c0prim)/(2*self.lbda*gammprim)*(b + np.sqrt(delta))
        # ECS contributions
        a_f   =  (phi_s*tau_f*self.lbda)/(self.c*(phi_s - phi_f))
        a_s   = -(phi_f*tau_s*self.lbda)/(self.c*(phi_s - phi_f))
        # Outputs
        param = {}
        param['tau_f'] = tau_f
        param['tau_s'] = tau_s
        param['a_f']   = a_f
        param['a_s']   = a_s
        param['phi_f'] = phi_f
        param['phi_s'] = phi_s
        param['tau_e'] = a_f*tau_s + a_s*tau_f
        return param

    def plot(self, myForcing, n_years, lplotFNH = False):
        output = analytical_EBM(self, myForcing, n_years)
        x = np.arange(n_years+1)
        plt.figure()
        plt.grid(True)        
        plt.plot(x, output['T'], linestyle='-', color='k')
        plt.plot(x, output['T0'], linestyle='--', color='k')
        plt.plot(x, output['termI'], linestyle='-', color='r')
        plt.plot(x, output['termF'], linestyle='-', color='g')
        plt.plot(x, output['termS'], linestyle='-', color='b')
        plt.plot(x, output['termE'], linestyle='-', color='c')        
        plt.plot(x, output['term0F'], linestyle='--', color='g')
        plt.plot(x, output['term0S'], linestyle='--', color='b')
        plt.plot(x, output['term0E'], linestyle='--', color='c')
        if lplotFNH:
            plt.plot(x, output['N'], '.', color='r', alpha=0.8)
            plt.plot(x, output['F'], '.', color='g', alpha=0.8)
            plt.plot(x, output['H'], '.', color='b', alpha=0.8)

        

#--------------------------------------------------------------------------------------------------
class FORCING(object):
    def __init__(self, typ='abrupt', t_0=0, xCO2_0=0, tau_e=0, xCO2_infty=4, k_=1/140, t_m=0, xCO2_m=0):
        self.typ = typ
        self.t_0 = t_0
        self.xCO2_0 = xCO2_0
        self.tau_e = tau_e
        self.xCO2_infty = xCO2_infty
        self.t_m = t_m
        self.xCO2_m = xCO2_m
        self.k_ = k_


#--------------------------------------------------------------------------------------------------
def optimal_t0(tau_s, x_0, x_infty):
    return np.int(-tau_s*np.log(1.0 - np.log(x_infty)/np.log(x_0)))

def optimal_x0(tau_e, tau_s, x_infty):
    return np.exp(tau_s/tau_e*np.log(x_infty))

def analytical_EBM(myEBM, myFORCING, n_years=150, lprint=False):

    # temporal axis
    t_ = np.arange(n_years+1, dtype=np.float)
    
    # myEBM 
    forc = myEBM.F   
    lbda = myEBM.lbda
    c = myEBM.c
    c_0 = myEBM.c_0
    gam = myEBM.gam
    epsi = myEBM.epsi
    xCO2 = myEBM.xCO2

    paramEBM = myEBM.parameters()
    phi_f = paramEBM['phi_f']
    phi_s = paramEBM['phi_s']
    tau_f = paramEBM['tau_f']
    tau_s = paramEBM['tau_s']
    a_f = paramEBM['a_f']
    a_s = paramEBM['a_s']

    # myFORCING
    typforcing = myFORCING.typ
    t_0 = myFORCING.t_0
    xCO2_0 = myFORCING.xCO2_0
    tau_e = myFORCING.tau_e
    xCO2_infty = myFORCING.xCO2_infty
    t_m = myFORCING.t_m
    xCO2_m = myFORCING.xCO2_m
    k_ = myFORCING.k_

    # from xCO2 to forcing
    F_infty = forc / np.log(xCO2) * np.log(xCO2_infty)
    F_m     = forc / np.log(xCO2) * np.log(xCO2_m)
    F_0     = forc / np.log(xCO2) * np.log(xCO2_0)

    output = {}
    
    if typforcing == 'abrupt':
        if lprint:
            print('-----------------------------------------------')
            print('ABRUPT : F = F_infty')
            print('-----------------------------------------------')
            print('F_infty = %5.2f W/m2 (xCO2_infty = %4.1f)'%(F_infty, xCO2_infty))
            print('-----------------------------------------------')

        F = np.full_like(t_, F_infty, dtype=np.float)

        termI  = F / lbda       
        termF  = - c_step(F_infty, lbda, a_f)*np.exp(-t_/tau_f)
        termS  = - c_step(F_infty, lbda, a_s)*np.exp(-t_/tau_s)
        term0F = - c_step(F_infty, lbda, phi_f*a_f)*np.exp(-t_/tau_f)
        term0S = - c_step(F_infty, lbda, phi_s*a_s)*np.exp(-t_/tau_s)
        termE  = np.zeros(n_years+1)
        term0E = np.zeros(n_years+1)
    
    elif typforcing == 'linear':
        if lprint:
            print('-----------------------------------------------')
            print('LINEAR : F = k_ * t_ * F_infty')
            print('-----------------------------------------------')
            print('F_infty = %5.2f W/m2 (xCO2_infty = %4.1f)'%(F_infty, xCO2_infty))
            print('k_ = %5.2f'%k_)
            print('-----------------------------------------------')

        F = k_ * t_ * F_infty

        termI  = F / lbda
        termF  = c_linear(F_infty, k_, lbda, a_f, tau_f)*(np.exp(-t_/tau_f) - 1)
        termS  = c_linear(F_infty, k_, lbda, a_s, tau_s)*(np.exp(-t_/tau_s) - 1)
        term0F = c_linear(F_infty, k_, lbda, phi_f*a_f, tau_f)*(np.exp(-t_/tau_f) - 1)
        term0S = c_linear(F_infty, k_, lbda, phi_s*a_s, tau_s)*(np.exp(-t_/tau_s) - 1)
        termE  = np.zeros(n_years+1)
        term0E = np.zeros(n_years+1)
    
    elif typforcing == 'double_abrupt' :
        if lprint:
            print('-----------------------------------------------')
            print('DOUBLE ABRUPT : F[:t_0]=F_0 and F[t_0:]=F_infty')
            print('-----------------------------------------------')
            print('t_0 = %i'%t_0)
            print('F_0 = %5.2f W/m2 (xCO2_0 = %4.1f)'%(F_0, xCO2_0))
            print('F_infty = %5.2f W/m2 (xCO2_infty = %4.1f)'%(F_infty, xCO2_infty))
            print('-----------------------------------------------')

        tt_1 = np.arange(0, t_0+1)
        tt_2 = np.arange(t_0+1, n_years+1)

        F = np.zeros(n_years+1)
        F[tt_1] = F_0
        F[tt_2] = F_infty
        
        termI  = np.zeros(n_years+1)
        termF  = np.zeros(n_years+1)
        termS  = np.zeros(n_years+1)
        term0F = np.zeros(n_years+1)
        term0S = np.zeros(n_years+1)
        termE  = np.zeros(n_years+1)
        term0E = np.zeros(n_years+1)

        # Part 1 : tt_1 = [0:t_0]
        output_ab0   = analytical_EBM(myEBM, FORCING(typ='abrupt', xCO2_infty=xCO2_0), t_0)
        termI[tt_1]  = output_ab0['termI']
        termF[tt_1]  = output_ab0['termF']
        termS[tt_1]  = output_ab0['termS']
        term0F[tt_1] = output_ab0['term0F']
        term0S[tt_1] = output_ab0['term0S']

        # Part 2 : tt_2 = [t_0+1:n_years]
        termI[tt_2]  = F_infty/lbda
        termF[tt_2]  = c_2step(F_0, F_infty, t_0, lbda, a_f, tau_f)*np.exp(-t_[tt_2]/tau_f)
        termS[tt_2]  = c_2step(F_0, F_infty, t_0, lbda, a_s, tau_s)*np.exp(-t_[tt_2]/tau_s)
        term0F[tt_2] = c_2step(F_0, F_infty, t_0, lbda, phi_f*a_f, tau_f)*np.exp(-t_[tt_2]/tau_f)
        term0S[tt_2] = c_2step(F_0, F_infty, t_0, lbda, phi_s*a_s, tau_s)*np.exp(-t_[tt_2]/tau_s)

    elif typforcing == 'triple_abrupt' :
        if lprint:
            print('-----------------------------------------------')
            print('TRIPLE ABRUPT : F[:t_0]=F_0, F[t_0:t_m]=F_m') 
            print('                F[t_m:]=F_infty')
            print('-----------------------------------------------')
            print('t_0 = %i'%t_0)
            print('F_0 = %5.2f W/m2 (xCO2_0 = %4.1f)'%(F_0, xCO2_0))
            print('t_m = %i'%t_m)
            print('F_m = %5.2f W/m2 (xCO2_m = %4.1f)'%(F_m, xCO2_m))
            print('F_infty = %5.2f W/m2 (xCO2_infty = %4.1f)'%(F_infty, xCO2_infty))
            print('-----------------------------------------------')

        tt_1 = np.arange(0, t_0+1)
        tt_2 = np.arange(t_0+1, t_m+1)
        tt_3 = np.arange(t_m+1, n_years+1)

        F = np.zeros(n_years+1)
        F[tt_1] = F_0
        F[tt_2] = F_m
        F[tt_3] = F_infty
        
        termI  = np.zeros(n_years+1)
        termF  = np.zeros(n_years+1)
        termS  = np.zeros(n_years+1)
        term0F = np.zeros(n_years+1)
        term0S = np.zeros(n_years+1)
        termE  = np.zeros(n_years+1)
        term0E = np.zeros(n_years+1)

        # Part 1 : tt_1 = [0:t_0]
        output_ab0   = analytical_EBM(myEBM, FORCING(typ='abrupt', xCO2_infty=xCO2_0), t_0)
        termI[tt_1]  = output_ab0['termI']
        termF[tt_1]  = output_ab0['termF']
        termS[tt_1]  = output_ab0['termS']
        term0F[tt_1] = output_ab0['term0F']
        term0S[tt_1] = output_ab0['term0S']

        # Part 2 : tt_2 = [t_0+1:t_m]
        termI[tt_2]  = F_m/lbda
        termF[tt_2]  = c_2step(F_0, F_m, t_0, lbda, a_f, tau_f)*np.exp(-t_[tt_2]/tau_f)
        termS[tt_2]  = c_2step(F_0, F_m, t_0, lbda, a_s, tau_s)*np.exp(-t_[tt_2]/tau_s)
        term0F[tt_2] = c_2step(F_0, F_m, t_0, lbda, phi_f*a_f, tau_f)*np.exp(-t_[tt_2]/tau_f)
        term0S[tt_2] = c_2step(F_0, F_m, t_0, lbda, phi_s*a_s, tau_s)*np.exp(-t_[tt_2]/tau_s)

        # Part 3 : tt_3 = [t_m+1:n_years]
        termI[tt_3]  = F_infty/lbda
        termF[tt_3]  = c_3step(F_0, F_m, F_infty, t_0, t_m, lbda, a_f, tau_f)*np.exp(-t_[tt_3]/tau_f)
        termS[tt_3]  = c_3step(F_0, F_m, F_infty, t_0, t_m, lbda, a_s, tau_s)*np.exp(-t_[tt_3]/tau_s)
        term0F[tt_3] = c_3step(F_0, F_m, F_infty, t_0, t_m, lbda, phi_f*a_f, tau_f)*np.exp(-t_[tt_3]/tau_f)
        term0S[tt_3] = c_3step(F_0, F_m, F_infty, t_0, t_m, lbda, phi_s*a_s, tau_s)*np.exp(-t_[tt_3]/tau_s)

    elif typforcing == 'expo' :
        if lprint:
            print('-----------------------------------------------')
            print('EXPONENTIAL : F = F_infty+(F_0-F_infty)*exp(-t_/tau_e)')
            print('-----------------------------------------------')
            print('F_infty = %5.2f W/m2 (xCO2_infty = %4.1f)'%(F_infty, xCO2_infty))
            print('F_0 = %5.2f W/m2 (xCO2_0 = %4.1f)'%(F_0, xCO2_0))
            print('tau_e = %i'%tau_e)
            print('-----------------------------------------------')
       
        F = F_infty + (F_0 - F_infty)*np.exp(-t_/tau_e)
        
        if tau_e == tau_s or tau_e == tau_f :
            print('>> Error')
            return None

        termI  = np.full_like(t_, F_infty/lbda, dtype=np.float)
        termS  = -(c_step(F_infty, lbda, a_s) + c_expo(F_0, F_infty, tau_e, lbda, a_s, tau_s))*np.exp(-t_/tau_s)
        termF  = -(c_step(F_infty, lbda, a_f) + c_expo(F_0, F_infty, tau_e, lbda, a_f, tau_f))*np.exp(-t_/tau_f)
        termE  = (c_expo(F_0, F_infty, tau_e, lbda, a_s, tau_s) + c_expo(F_0, F_infty, tau_e, lbda, a_f, tau_f))*np.exp(-t_/tau_e)
        term0S = -(c_step(F_infty, lbda, phi_s*a_s) + c_expo(F_0, F_infty, tau_e, lbda, phi_s*a_s, tau_s))*np.exp(-t_/tau_s)
        term0F = -(c_step(F_infty, lbda, phi_f*a_f) + c_expo(F_0, F_infty, tau_e, lbda, phi_f*a_f, tau_f))*np.exp(-t_/tau_f)
        term0E = (c_expo(F_0, F_infty, tau_e, lbda, phi_s*a_s, tau_s) + c_expo(F_0, F_infty, tau_e, lbda, phi_f*a_f, tau_f))*np.exp(-t_/tau_e)
    
    T = termI + termS + termF + termE
    T0 = termI + term0S + term0F + term0E
    H = gam*(T - T0)
    N = F - lbda*T - (epsi-1)*H
    
    output['F'] = F
    output['T'] = T
    output['T0'] = T0
    output['N'] = N
    output['H'] = H
    
    output['termI'] = termI
    output['termS'] = termS
    output['termF'] = termF
    output['termE'] = termE
    output['term0S'] = term0S
    output['term0F'] = term0F
    output['term0E'] = term0E

    output['normalized_T'] = T / (F_infty/lbda)
    output['normalized_T0'] = T0 / (F_infty/lbda)
    
    return output


#--------------------------------------------------------------------------------------------------
# Auxiliary functions for analytical_EBM
def c_step(F_infty, lbda, x_n):
    return x_n * F_infty/lbda

def c_linear(F_infty, k_, lbda, x_n, tau_n):
    return k_ * F_infty/lbda * tau_n * x_n

def c_2step(F_0, F_infty, t_0, lbda, x_n, tau_n):
    return x_n/lbda*((F_0 - F_infty)*np.exp(t_0/tau_n) - F_0)

def c_3step(F_0, F_m, F_infty, t_0, t_m, lbda, x_n, tau_n):
    return x_n/lbda*((F_m - F_infty)*np.exp(t_m/tau_n) + (F_0 - F_m)*np.exp(t_0/tau_n) - F_0)

def c_expo(F_0, F_infty, tau_e, lbda, x_n, tau_n):
    pi_n = 1.0 / (1.0 - tau_n/tau_e)
    return x_n/lbda*pi_n*(F_0 - F_infty)
#--------------------------------------------------------------------------------------------------


def derive_EBM_I(T, N, xCO2, nyr_s = 30, nyr_f = 10):
    
    n_years = T.size
    
    #--------------------------------------------------------------------------------------------------
    # 1. Estimation of F and lambda from Gregory plot (linear fit of N)
    #--------------------------------------------------------------------------------------------------
    slope, intercept, r_value, p_value, std_err = stats.linregress(T, N)
    lbda = - slope
    forc = intercept
    T_eq = forc / lbda
    
    #--------------------------------------------------------------------------------------------------
    # 2. Estimation of a_s and tau_s from fit of log(1 - T/T_eq) : see Eq. (17) 
    #--------------------------------------------------------------------------------------------------
    t_i   = nyr_s
    x_    = np.arange(t_i, n_years)
    y_    = np.log(1 - T[t_i-1:n_years-1] / T_eq)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_, y_)
    tau_s = -1.0 / slope
    a_s   = np.exp(intercept)
    
    #--------------------------------------------------------------------------------------------------
    # 3. Estimation of tau_f by averaging over the first 10 years : see Eq. (18)
    #--------------------------------------------------------------------------------------------------
    a_f = 1 - a_s
    t_i = nyr_f
    t_  = np.arange(1, t_i)
    tau = t_ / (np.log(a_f) - np.log(1 - T[0:t_i-1]/T_eq - a_s*np.exp(-t_/tau_s)))
    tau_f = np.mean(tau)
    
    #--------------------------------------------------------------------------------------------------
    # 4. Derivation of c, c0 and gamma : see Eqs (19), (20) and (21)
    #--------------------------------------------------------------------------------------------------
    c   = lbda / (a_f / tau_f + a_s / tau_s)
    c_0 = lbda*(a_f*tau_f + a_s*tau_s) - c
    gam = c_0 / (a_s*tau_f + a_f*tau_s)

    output = EBM(F=forc, lbda=lbda, c=c, c_0=c_0, gam=gam, epsi=1.0, xCO2=xCO2)

    return output


def derive_EBM_II(T, N, xCO2):
    
    n_years = T.size
    
    #--------------------------------------------------------------------------------------------------
    # 0. set param to the EBM-1 values
    #--------------------------------------------------------------------------------------------------
    
    EBM_0 = derive_EBM_I(T, N, xCO2)
    
    forcage = FORCING(typ='abrupt', xCO2_infty=xCO2)

    datas_EBM = analytical_EBM(EBM_0, forcage, n_years)
    T0_EBM = datas_EBM['T0']
    T_EBM = datas_EBM['T']
    H_EBM = datas_EBM['H']
    
    n_iters = 10
    for iter in range(n_iters):
        
        X = np.c_[T, H_EBM[1:]]
        regr = linear_model.LinearRegression()
        regr.fit(X, N)
        forc = regr.intercept_
        lbda = - regr.coef_[0]
        epsi = 1 - regr.coef_[1]
        
        # print('=====>')       
        # print(forc)
        # print(lbda)
        # print(epsi)
        # print('=====>')       

        T_eq = forc / lbda
        
        t_i   = 80
        x_    = np.arange(t_i, n_years)
        y_    = np.log(1 - T[t_i-1:n_years-1] / T_eq)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_, y_)    
        tau_s = -1.0 / slope
        a_s   = np.exp(intercept)
        
        a_f = 1 - a_s
        t_i = 6
        t_  = np.arange(1, t_i)
        tau = t_ / (np.log(a_f) - np.log(1 - T[0:t_i-1]/T_eq - a_s*np.exp(-t_/tau_s)))
        tau_f = np.mean(tau)
        
        c   = lbda / (a_f / tau_f + a_s / tau_s)
        c_0 = lbda*(a_f*tau_f + a_s*tau_s) - c
        gam = c_0 / (a_s*tau_f + a_f*tau_s)

        # print(c)
        # print(c_0)
        
        c_0 = c_0 / epsi
        gam = gam / epsi
        
        myEBM = EBM(F=forc, lbda=lbda, c=c, c_0=c_0, gam=gam, epsi=epsi, xCO2=xCO2)
        
        datas_EBM = analytical_EBM(myEBM, forcage, n_years)
        T0_EBM = datas_EBM['T0']
        T_EBM = datas_EBM['T']
        H_EBM = datas_EBM['H']

    output = EBM(F=forc, lbda=lbda, c=c, c_0=c_0, gam=gam, epsi=epsi, xCO2=xCO2)

    return output

def derive_ECS_from_EBM(iEBM, myT, myN, xCO2):
    if iEBM == 1:
        ebm = derive_EBM_I(myT, myN, xCO2)
        return ebm.F / ebm.lbda
    elif iEBM == 2:
        ebm = derive_EBM_II(myT, myN, xCO2)
        return ebm.F / ebm.lbda
    else:
        print('%i is not known')%iEBM
        return None

def derive_ECS_from_stab(myT, n_years):
    return np.mean(myT[-n_years:])
