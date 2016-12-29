# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:07:01 2012
Check Pierce equation
@author: Wei
"""

import numpy as np
import matplotlib.pylab as plt
from scipy import special, optimize
import sys, math
sys.path.append(r'S:\projecten\QSIDE\simulations\Project_scripts\mainClass')  # used for check
import QSIDE  # used for check

sys.path.append(r"S:\software\intern\script_wgw\Multiple_diffraction")
from MDF_v2 import DoubleDiffraction2D, dist2D

def errSin(p, x, fresnelS):
    return fitSin(p, x) - fresnelS

def fitSin(p, x):
    if type(x)==list:
        x = np.array(x)
    fitS = p[0]/(x+0.47)*np.sin(p[1]*x**2+p[2]) + p[3]
    return fitS
    
def errCos(p, x, fresnelC):
    return fitCos(p, x) - fresnelC

def fitCos(p, x):
    if type(x)==list:
        x = np.array(x)
    fitC = p[0]/(x+0.47)*np.cos(p[1]*x**2+p[2]) + p[3]
    return fitC
    
def errCos_bySin(p, x, fresnelC):
    return fitCos_bySin(p, x) - fresnelC

def fitCos_bySin(p, x):
    if type(x)==list:
        x = np.array(x)
    fitC = 1./(np.pi*(x+p))*np.sin(np.pi/2*x**2)+0.5
    return fitC

def errSin_byCos(p, x, fresnelS):
    return fitSin_byCos(p, x) - fresnelS
    
def fitSin_byCos(p, x):
    if type(x)==list:
        x = np.array(x)
    fitS = -1./(np.pi*(x+p))*np.cos(np.pi/2*x**2)+0.5
    return fitS

class IntSinCos():
    def __init__(self, p=0.05):
        """ p is the extra correction for the equation, by adding p, we avoid 
            the singularity when x==0. C=0.5+1/(pi*x+p)*sin(0.5*pi*x**2); 
            S = 0.5-1/(pi*x)*cos(0.5*pi*x**2). 
            After fitting we get that for C p=2.04; for S p=1.46. 
            To vanish it in the f^2+g^2=C^2+S^2-S-C+0.5, p should equals each 
            other in both S and C. we use p=0.05 in default. When x=0 the Fresnel
            integral will be zero but f^2+g^2=0.5, from this we get p=0.05
        """
        self.p = p
        
    def _intC(self, x):
        if type(x)==list:
            x = np.array(x)
        fitC = 1./(np.pi*(x+self.p))*np.sin(np.pi/2*x**2)+0.5
        return fitC
    
    def _intS(self, x):
        if type(x)==list:
            x = np.array(x)
        fitS = -1./(np.pi*(x+self.p))*np.cos(np.pi/2*x**2)+0.5
        return fitS     
        
def fIntC(x):
    if type(x)==list:
        x = np.array(x)
    fitC = 0.37/(x+0.37)*np.sin(np.pi/2*x**2)+0.5
    return fitC          
    
def fIntS(x):
    if type(x)==list:
        x = np.array(x)
    fitS = -0.37/(x+0.37)*np.cos(np.pi/2*x**2)+0.5
    return fitS 
    
def fcos(x):
    p = [-0.4,  1.6,  1.2,  0.5]
    fitC = p[0]/(x+0.4)*np.cos(p[1]*x**2+p[2]) + p[3]
    return fitC

def fsin(x):
    p = [-0.4,  1.6,  1.2,  0.5]
    fitS = p[0]/(x+0.4)*np.sin(p[1]*x**2+p[2]) + p[3]
    return fitS
    
def fitf(x):
    x = np.abs(x)
    S = fsin(x)  # fresnel sin
    C = fcos(x)  # fresnel cosEquation~(\ref{Lpi}) would underestimate $A_{bar}$ a little for the interested noise-map height 4m. As a result, to calculate another time for the contribution of the image source is recommended.
    fx = (0.5-S)*np.cos(0.5*np.pi*x**2) - (0.5-C)*np.sin(0.5*np.pi*x**2)
    return fx    
    
def fitg(x):
    x = np.abs(x)
    S = fsin(x)  # fresnel sin
    C = fcos(x)  # fresnel cos
    gx = (0.5-C)*np.cos(0.5*np.pi*x**2) + (0.5-S)*np.sin(0.5*np.pi*x**2)
    return gx 
    
def fg(x):
    x = np.abs(x)
    SC = special.fresnel(x)
    S = SC[0]  # fresnel sin
    C = SC[1]  # fresnel cos
    fx = (0.5-S)*np.cos(0.5*np.pi*x**2) - (0.5-C)*np.sin(0.5*np.pi*x**2)
    gx = (0.5-C)*np.cos(0.5*np.pi*x**2) + (0.5-S)*np.sin(0.5*np.pi*x**2)
    return [fx, gx]
    
def f(x):
    x = np.abs(x)
    SC = special.fresnel(x)
    S = SC[0]  # fresnel sin
    C = SC[1]  # fresnel cos
    fx = (0.5-S)*np.cos(0.5*np.pi*x**2) - (0.5-C)*np.sin(0.5*np.pi*x**2)    
    return fx

def g(x):
    x = np.abs(x)
    SC = special.fresnel(x)
    S = SC[0]  # fresnel sin
    C = SC[1]  # fresnel cos
    gx = (0.5-C)*np.cos(0.5*np.pi*x**2) + (0.5-S)*np.sin(0.5*np.pi*x**2)
    return gx 
    
def S(x):
    x = np.abs(x)
    SC = special.fresnel(x)
    S = SC[0]  # fresnel sin
    return S
    
def C(x):
    x = np.abs(x)
    SC = special.fresnel(x)
    C = SC[1]  # fresnel cos
    return C

def x056S(x):
    p = 0.56
    fitS = -0.4/(x+p)*np.sin(1.6*x**2+1.5) + 0.5
    return fitS    
    
def x056C(x):
    p = 0.56
    fitC = -0.4/(x+p)*np.cos(1.6*x**2+1.5) + 0.5
    return fitC    
    
def x13piS(x):
    p = [1./np.pi, 1./np.pi]
    fitS = -p[0]/(x+p[1])*np.cos(0.5*np.pi*x**2) + 0.5
    return fitS

def x13piC(x):
    p = [1./np.pi, 1./np.pi]
    fitC = p[0]/(x+p[1])*np.sin(0.5*np.pi*x**2) + 0.5
    return fitC 
    
def series_S(x):
    return 0.5 - 1./(np.pi*x)*np.cos(0.5*np.pi*x**2.)

def series_C(x):
    return 0.5 + 1./(np.pi*x)*np.sin(0.5*np.pi*x**2.)

def x04(p, x):
    y = (p/(np.pi*(x+p)))**2.
    return y

def err_x04(p, x, trueValue):
    return x04(p,x)-trueValue
    
def fit_x04(x):
    trueValue = S(x)**2.+C(x)**2. - C(x) - S(x) + 0.5
    p0 = 0.38
    p1, success = optimize.leastsq(err_x04, p0, args=(x, trueValue))
    print "p1 for [0.4/(p1+x)]**2: ", p1
    plt.plot(x, trueValue, '-')
    plt.plot(x, x04(p1, x), '--')
    plt.grid()

def fit_fresnelIntSC(x):
    fresnelSC = special.fresnel(x)
    C = fresnelSC[1]
    p0 = 0.1
    p1, success = optimize.leastsq(errCos_bySin, p0, args=(x, C))
    print 'P1 for fresnelC: ', p1
    plt.plot(x, C, '-')
    plt.plot(x, fitCos_bySin(p1, x), '--')
    plt.legend(['FresnelC', 'Fitted cos(x)'], loc='best')
    plt.xlabel('x')
    plt.ylabel('C(x)')
    plt.grid()
    
    S = fresnelSC[0]
    p2, success = optimize.leastsq(errSin_byCos, p0, args=(x, S))
    print 'p2 for fresnelS: ', p2
    plt.figure()
    plt.plot(x, S, '-')
    plt.plot(x, fitSin_byCos(p2, x), '--')
    plt.legend(['FresnelS', 'Fitted sin(x)'], loc='best')
    plt.xlabel('x')
    plt.ylabel('S(x)')
    plt.grid()
        
def fit_fresenlSC(x):
    fresnelSC = special.fresnel(x)
    C = fresnelSC[1]
    p0 = [-0.36124879,  1.83658169,  1.05923556,  0.50227021]
    p1, success = optimize.leastsq(errCos, p0, args=(x, C))
    print 'P1 for fresnelC: ', p1
    plt.plot(x, C, '-')
    plt.plot(x, fitCos(p1, x), '--')
    plt.legend(['FresnelC', 'Fitted cos(x)'], loc='best')
    plt.xlabel('x')
    plt.ylabel('C(x)')
    plt.grid()
    
    S = fresnelSC[0]
    p2, success = optimize.leastsq(errSin, p0, args=(x, S))
    print 'p2 for fresnelS: ', p2
    plt.figure()
    plt.plot(x, S, '-')
    plt.plot(x, fitSin(p2, x), '--')
    plt.legend(['FresnelS', 'Fitted sin(x)'], loc='best')
    plt.xlabel('x')
    plt.ylabel('S(x)')
    plt.grid()


def plot_fresenlSC(x):
    fresnelSC = special.fresnel(x)
    C = fresnelSC[1]
    fitFcos = fcos(x)
    plt.plot(x, C, '-')
    plt.plot(x, fitFcos, '--')
    plt.legend(['FresnelC', 'Fitted cos(x)'], loc='best')
    plt.xlabel('x')
    plt.ylabel('C(x)')
    plt.grid()
    
    S = fresnelSC[0]
    fitFsin = fsin(x)
    plt.figure()
    plt.plot(x, S, '-')
    plt.plot(x, fitFsin, '--')
    plt.legend(['FresnelS', 'Fitted sin(x)'], loc='best')
    plt.xlabel('x')
    plt.ylabel('S(x)')
    plt.grid()

def fresnel_int_compare():
    x = np.linspace(0.01, 5., 500)
    C1 = x13piC(x)
    S1 = x13piS(x)
    C2 = x056C(x)
    S2 = x056S(x)
    C3 = fIntC(x)
    S3 = fIntS(x)
    C4 = series_C(x)
    S4 = series_S(x)
    Ct = C(x)
    St = S(x)
    plt.figure(1)
    plt.plot(x, C1, '-',  x, C3, '-+',  x, C4, '-.', x, Ct, '--')
    plt.legend(['1/[pi(1/pi+x)]*sin', '0.37/(0.37+x)*sin','1/(pi*x)*sin', 'theory'], loc='best')
    plt.xlabel('x')
    plt.ylabel('C(x)')
    plt.figure(2)
    plt.plot(x, S1, '-', x, S3, '-+',   x, S4, '-.', x, St, '--')
    plt.legend(['-1/[pi(1/pi+x)]*cos', '-0.37/(0.37+x)*cos','-1/(pi*x)*cos', 'theory'], loc='best')
    plt.ylim([-0.5, 1.])
    plt.xlabel('x')
    plt.ylabel('S(x)')

def fitfX056(x):
    x = np.abs(x)
    S = x056S(x)  # fresnel sin
    C = x056C(x)  # fresnel cos
    fx = (0.5-S)*np.cos(0.5*np.pi*x**2) - (0.5-C)*np.sin(0.5*np.pi*x**2)
    return fx    
    
def fitgX056(x):
    x = np.abs(x)
    S = x056S(x)  # fresnel sin
    C = x056C(x)  # fresnel cos
    gx = (0.5-C)*np.cos(0.5*np.pi*x**2) - (0.5-S)*np.sin(0.5*np.pi*x**2)
    return gx 

def fitfOfficial(x):
    x = np.abs(x)
    S = fIntS(x)  # fresnel sin
    C = fIntC(x)  # fresnel cos
    fx = (0.5-S)*np.cos(0.5*np.pi*x**2) - (0.5-C)*np.sin(0.5*np.pi*x**2)
    return fx   

def fitgOfficial(x):
    x = np.abs(x)
    S = fIntS(x)  # fresnel sin
    C = fIntC(x)  # fresnel cos
    gx = (0.5-C)*np.cos(0.5*np.pi*x**2) + (0.5-S)*np.sin(0.5*np.pi*x**2)
    return gx     

def FDTD_diffr_compare():
    # double diffraction
    dire2 = r'S:\projecten\QSIDE\simulations\Project_scripts\PieceEquation\noGround'
    simFile = 'Wi10_Hi11'
    obj2 = QSIDE.FDTD_2D(dire2, 'Wi10_Hi11', 'yverloop1')
    [rcvSplFF1, fr] = obj2._simMinusFF('oct')    
    spos = [5.2, 1.0]
    fr = [31.5, 63, 125, 250, 500, 1000]
    phis = math.atan(4.8/10.0)
    rpositions = np.loadtxt(dire2 + '\\' + simFile + '\\' + 'yverloop1.positions.txt')
    for n in range(rpositions.shape[0]):
        rpos = rpositions[n, :]
        phir = math.atan((rpos[0]-20.0)/abs((11-rpos[1])))
        mObj = DoubleDiffraction2D()
        [dLevels, distance] = mObj._levelrefFF2NoGround(spos, rpos, \
        [10., 11.], [20., 11.], 1.5*np.pi,1.5*np.pi, phis, phir, fr)    
        [dLevelsX037, distance] = mObj._levelrefFF2NoGroundX037(spos, rpos, \
        [10., 11.], [20., 11.], 1.5*np.pi,1.5*np.pi, phis, phir, fr)  
        [dLevelsX13pi, distance] = mObj._levelrefFF2NoGroundX13pi(spos, rpos, \
        [10., 11.], [20., 11.], 1.5*np.pi,1.5*np.pi, phis, phir, fr) 
        Lc = 10.0*np.log10(distance/dist2D(spos, rpos))
        plt.figure()
        print 'dLevels: ', dLevels
        print 'rcvSplFF1[:,n]+Lc: ', rcvSplFF1[:,n]+Lc
        plt.semilogx(fr, dLevels, '--', fr, dLevelsX037, '-+', fr, dLevelsX13pi, '-.', fr, rcvSplFF1[:,n]+Lc, '-')
        plt.xlim([25, 1300])
        plt.legend(['PI(0.4/(0.4+X))^2', 'PI(0.37/(0.37+X))^2', 'PI(1/[pi(1/pi+X)])^2', 'FDTD'], loc='best')
        
def FDTD_diffr_compare2():
    # double diffraction
    dire2 = r'S:\projecten\QSIDE\simulations\Project_scripts\PieceEquation\noGround'
    simFile = 'Wi4_Hi11'
    obj2 = QSIDE.FDTD_2D(dire2, 'Wi10_Hi11', 'yverloop1')
    [rcvSplFF1, fr] = obj2._simMinusFF('oct')    
    spos = [5.2, 1.0]
    fr = [31.5, 63, 125, 250, 500, 1000]
    phis = math.atan(4.8/10.0)
    rpositions = np.loadtxt(dire2 + '\\' + simFile + '\\' + 'yverloop1.positions.txt')
    for n in range(rpositions.shape[0]):
        rpos = rpositions[n, :]
        phir = math.atan((rpos[0]-14.0)/abs((11-rpos[1])))
        mObj = DoubleDiffraction2D()
        [dLevels, distance] = mObj._levelrefFF2NoGround(spos, rpos, \
        [10., 11.], [14., 11.], 1.5*np.pi,1.5*np.pi, phis, phir, fr)    
        [dLevelsX056, distance] = mObj._levelrefFF2NoGroundX037(spos, rpos, \
        [10., 11.], [14., 11.], 1.5*np.pi,1.5*np.pi, phis, phir, fr)  
        Lc = 10.0*np.log10(distance/dist2D(spos, rpos))
        plt.figure()
        print 'dLevels: ', dLevels
        print 'rcvSplFF1[:,n]+Lc: ', rcvSplFF1[:,n]+Lc
        plt.semilogx(fr, dLevels, '--', fr, dLevelsX056, '-+', fr, rcvSplFF1[:,n]+Lc, '-')
        plt.xlim([25, 1300])
        plt.legend(['p=0.4', 'p=0.56', 'FDTD'], loc='best')
        
        
def theory_diffr_compare():
    gama = np.sqrt(40.0/3.0)
    B = np.sqrt(3.0/4.0)
    MvS = np.sqrt(3.0)*np.cos(2.0/3.0*np.pi/4) - np.sqrt(3)/2.0
    
    angles = []
    pDiffr = []
    pDiffrFit =  []
    pDFO = []
    p0 = 0.37
    for thetaL in np.linspace(0, np.pi/2-0.01, 200):    
        MvL = np.sqrt(3)*np.cos(2.0/3.0*thetaL) - np.sqrt(3)/2.0
        Ys = gama*MvS
        YL = gama*MvL      
        if Ys>YL:   
            pDiffrFit.append(10*np.log10((p0/(Ys+p0))**2*(p0/(B*YL+p0))**2.))
            pDFO.append(10*np.log10((1./(np.pi*Ys))**2*(1./(np.pi*B*YL))**2))
            pDiffr.append(10*np.log10((f(Ys)**2 + g(Ys)**2)*(f(B*YL)**2 + g(B*YL)**2)))
            angles.append(thetaL)
        elif Ys<YL:
            pDiffrFit.append(10*np.log10((p0/(YL+p0))**2*(p0/(B*Ys+p0))**2.))
            pDFO.append(10*np.log10((1./(np.pi*YL))**2*(1./(np.pi*B*Ys))**2))
            pDiffr.append(10*np.log10((f(YL)**2 + g(YL)**2)*(f(B*Ys)**2 + g(B*Ys)**2)))
            angles.append(thetaL)
        else:
            pass
    
    plt.figure()
    plt.plot(angles, pDiffr, '-', angles, pDiffrFit, '--', angles, pDFO, '-.')
    plt.xlabel('Angle')
    plt.ylabel('10log(|Pdiffr/Pat,L|^2) [dB]')
    plt.xticks([0,np.pi/12., np.pi/6, np.pi/4, np.pi/3, np.pi/12+np.pi/3, np.pi/2], [0, 15, 30, 45, 60, 75,90])
    plt.legend(['theory', '0.37/(x+0.37)', '1/(pi*x)'], loc='best')
    plt.grid()    
    
    plt.figure(2)
    plt.plot(angles, np.array(pDiffr)-np.array(pDiffr), '-', angles, np.array(pDiffr)-np.array(pDiffrFit), '--',  \
    angles, np.array(pDiffr)-np.array(pDFO), '-.')
    np.savetxt("theo_theo_angles.txt", angles)
    np.savetxt("theo_theo_values.txt", np.array(pDiffr)-np.array(pDiffr))
    np.savetxt("theo_fit_difference.txt", np.array(pDiffr)-np.array(pDiffrFit))
    np.savetxt("theo_official_difference.txt", np.array(pDiffr)-np.array(pDFO))
    plt.xlabel('Angle')
    plt.ylabel('Difference [dB]')
    plt.xticks([0,np.pi/12., np.pi/6, np.pi/4, np.pi/3, np.pi/12+np.pi/3, np.pi/2], [0, 15, 30, 45, 60, 75,90])
    plt.legend(['theory', '0.37/(x+0.37)', '1/(pi*x)'], loc='best')
    plt.ylim([-8, 5])
    plt.grid()    
    plt.show()

        
def theory_diffr_compare2():
    
    pDiffr = []
    pDiffrFit =  []
    pDFO = []
    p0 = 0.37
    x = np.linspace(0.01, 4, 200) 
    pDiffrFit = (p0/(x+p0))**2
    pDFO = (1./(np.pi*x))**2
    pDiffr = C(x)**2. + S(x)**2. - C(x) - S(x) + 0.5        
    
    plt.figure()
    plt.plot(x, pDiffr, '-', x, pDiffrFit, '--', x, pDFO, '-.')
    plt.xlabel('x')
    plt.ylabel('10log(|Pdiffr/Pat,L|^2) [dB]')
    plt.legend(['theory', '0.37/(x+0.37)', '1/(pi*x)'], loc='best')
    plt.grid()    
    plt.ylim([0,1])
    
    plt.figure(2)
    plt.plot(x, pDiffr-pDiffr, '-', x, pDiffr-pDiffrFit, '--',  \
    x, pDiffr-pDFO, '-.')
    plt.xlabel('x')
    plt.ylabel('Difference [dB]')
    plt.legend(['theory', '0.37/(x+0.37)', '1/(pi*x)'], loc='best')
    plt.grid()    
    plt.show()   
    plt.ylim([-0.1, 0.1])

def replot_abar_ground():
    # double diffraction
    dire2 = r'S:\projecten\QSIDE\simulations\Project_scripts\PieceEquation\withGround'
    simFile = ['Wi10_Hi11', 'Wi20_Hi11', 'Wi40_Hi11', 'Wi80_Hi11']
    w = [10., 20., 40., 80.]
    for fnum, sf in enumerate(simFile):
        obj2 = QSIDE.FDTD_2D(dire2, sf,  'yverloop1')
        [rcvSplFF1, fr] = obj2._simMinusFF('oct')    
        spos = [5.2, 0.0]
        fr = [31.5, 63, 125, 250, 500, 1000]
        phis = math.atan(4.8/11.0)
        rpositions = np.loadtxt(dire2 + '\\' + sf + '\\' + 'yverloop1.positions.txt')
        n = 3 # indicate the receiver height 4.4m
        rpos = rpositions[n, :]
        phir = math.atan(4.5/6.6)
        mObj = DoubleDiffraction2D()
        [dLevels, distance] = mObj._levelrefFF2Ground(spos, rpos, \
        [10., 11.], [10.+w[fnum], 11.], 1.5*np.pi, 1.5*np.pi, phis, phir, fr) 
        Lc = 10.0*np.log10(distance/dist2D(spos, rpos))
        plt.figure()
        print 'dLevels: ', dLevels
        print 'rcvSplFF1[:,n]+Lc: ', rcvSplFF1[:,n]
        plt.semilogx(fr, dLevels, '--', fr, rcvSplFF1[:,n]+Lc, '-')
        plt.xlim([25, 1300])
        plt.legend(['Calculated', 'FDTD'], loc='best')
    
if __name__=='__main__':
#    fit_fresnelIntSC(np.linspace(1., 10., 500))
#    fit_x04(np.linspace(0, 5, 1000))
#    fresnel_int_compare()
    theory_diffr_compare()
#    theory_diffr_compare2()
#    FDTD_diffr_compare()
#    FDTD_diffr_compare2()
#    replot_abar_ground()
    plt.show()

