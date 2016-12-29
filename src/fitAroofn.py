# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 13:27:26 2012
quantify the effect of roof shape only. 
There are no canyons in this analysis
The roof shape effect with canyons are fitted in acan_aroof.py file

MODIFID FROM S:\projecten\QSIDE\simulations\Project_scripts\Roofshapes
fit_roofshape_05.py and fit_roofshape_06.py

@author: Wei
"""

import sys
import numpy as np
import math
from scipy import optimize
import matplotlib.pylab as plt
sys.path.append(r'S:\software\intern\script_wgw\Multiple_diffraction')
from MDF_v2 import AttenuationBar

def dist2D(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

def fitFun(p, Abar):
    ''' fitting function '''
    return p[0]*Abar+p[1]

def errFun(p, Abar, triDiff):
    ''' error function '''
    return fitFun(p, Abar) - triDiff

def fitFun2(p, Y):
    ''' fitting function '''
    return p[0]*10.0*np.log10((abs(p[1])/(abs(Y)+abs(p[1])))**2)

def errFun2(p, Y, triDiff):
    ''' error function '''
    return fitFun2(p,Y) - triDiff

def fitFun3(p, fr):
    ''' fitting function '''
    return p[0]*10.0*np.log10(fr)+p[1]

def errFun3(p, fr, triDiff):
    ''' error function '''
    return fitFun3(p,fr) - triDiff


class FitAroof():
    def __init__(self):        
        self.h = 12.0   # height of barrier
        self.ry = 4.   # height fo receiver
        self.fr = [31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000]
        self.waveLen = 340.0/np.array(self.fr)
        self.sybols = ['-*', '-+', '-^', '->', '-<', '-o']
        self.sybolsf = ['--*', '--+', '--^', '-->', '--<', '--o']
    def _fitSlopeRoofNoCanyon(self):
        widths = [10, 20, 40, 80]
        dist2barrier = [0.5, 1., 2., 3.]
        LdT, LdF, LdTminusLdF = [], [], []
        for n, w in enumerate(widths):
#            plt.subplot(2,2,n+1)
            for d in dist2barrier:
                rx = 10.0+w+d
                spos = [w-d, 0]
                rpos = [rx, self.ry]
                N1point = [10.0, self.h]
                N2point = [10+w, self.h]
                Nmpoint = [(N1point[0]+N2point[0])/2.0, N1point[1]+5.0]
                edgeVertex = [N1point, Nmpoint, N2point]
                phis = math.atan(d/self.h)
                phir = math.atan(d/(self.h-rpos[1]))                
                alfa = math.atan(5.0/w*2)
                alfa2 = np.pi - 2*alfa
                beta = [1.5*np.pi-alfa, 2*np.pi-alfa2, 1.5*np.pi-alfa]
                barObj = AttenuationBar()
                [Lt,Dt] = barObj._barLevelrefFFGround(spos, rpos, edgeVertex, beta, phis, phir, self.fr)
                LdT.append(Lt)
                [Lf, Df] = barObj._barLevelrefFFGround(spos, rpos, [N1point, N2point], [1.5*np.pi, 1.5*np.pi], phis, phir, self.fr)
                LdF.append(Lf)
        Abar = np.reshape(LdT, [-1, 1])              
        LdFarray = np.reshape(LdF, [-1, 1])
        p0 = [0.3, 2]
        assert (len(Abar)==len(LdFarray))
        p1, success = optimize.leastsq(errFun, p0,\
            args=(Abar[:,0],  Abar[:, 0]-LdFarray[:, 0]))
        print '\n\nP[0]*Abar+p[1], p: ', p1, "\n\n"
        
        # modify on 31/10/2012. to calculate mse
        trueValue = Abar[:, 0]-LdFarray[:, 0]
        np.savetxt('trueValue.txt', trueValue)
        np.savetxt('Abar.txt', Abar)
        
        #plot results        
#        widths = [10, 20, 40, 80]
#        dist2barrier = [0.5, 1, 2, 3]
        for n, w in enumerate(widths):
#            plt.subplot(2,2,n+1)
            plt.figure(n, (7, 5))
            lgs = []
            for m, d in enumerate(dist2barrier):
                lgs.append(str(d)+'m')
                lgs.append(str(d)+'m-fit')
                rx = 10.0+w+d
                spos = [w-d, 0]
                rpos = [rx, self.ry]
                N1point = [10.0, self.h]
                N2point = [10+w, self.h]
                Nmpoint = [(N1point[0]+N2point[0])/2.0, N1point[1]+5.0]
                edgeVertex = [N1point, Nmpoint, N2point]
                phis = math.atan(d/self.h)
                phir = math.atan(d/(self.h-rpos[1]))                
                alfa = math.atan(5.0/w*2)
                alfa2 = np.pi - 2*alfa
                beta = [1.5*np.pi-alfa, 2*np.pi-alfa2, 1.5*np.pi-alfa]
                barObj = AttenuationBar()
                [LdT, Dt] = barObj._barLevelrefFFGround(spos, rpos, edgeVertex, beta, phis, phir, self.fr)
                [LdF, Df] = barObj._barLevelrefFFGround(spos, rpos, [N1point, N2point], [1.5*np.pi, 1.5*np.pi], phis, phir, self.fr)                 
                plt.semilogx(self.fr, LdT-LdF, self.sybols[m], self.fr, fitFun(p1, LdT), self.sybolsf[m])
            plt.xlim([25, 10000])
            plt.grid()
            plt.ylabel("Lslope - Lflat dB")
            plt.xlabel("Frequency Hz")
            plt.xticks(self.fr, self.fr)
            plt.title("width"+str(w))
            plt.legend(lgs, loc="best")     
            plt.savefig(str(w)+'m.pdf')
    def _BKawai(self, rs, rr, w):
        """ rs -> distance from source to the diffracted edge
            rr -> distance from the diffracted edge to the receiver
            w -> width of the diffraction edge
        """
        return w*(w+rr+rs)/((w+rs)*(w+rr))
    def _AoverLambda(self, p1, p2, L, waveLen):
        """ A in Kawai's paper page 232"""
        return 2*p1*p2/(L*waveLen)
        
    def _fitSlopeRoofNoCanyon_DB(self):
        widths = [10, 20, 40, 80, 160, 320]
        dist2barrier = [0.5, 2, 5, 10, 20]
        LdT, LdF, LdTminusLdF, B, A= [], [], [], [], []
        Y = []
        for n, w in enumerate(widths):
            for d in dist2barrier:
                rx = 10.0+w+d
                spos = [w-d, 0]
                rpos = [rx, self.ry]
                N1point = [10.0, self.h]
                N2point = [10+w, self.h]
                Nmpoint = [(N1point[0]+N2point[0])/2.0, N1point[1]+5.0]
                edgeVertex = [N1point, Nmpoint, N2point]
                phis = math.atan(d/self.h)
                phir = math.atan(d/(self.h-rpos[1]))                
                alfa = math.atan(5.0/w*2)
                alfa2 = np.pi - 2*alfa
                beta = [1.5*np.pi-alfa, 2*np.pi-alfa2, 1.5*np.pi-alfa]
                barObj = AttenuationBar()
                [Lt,Dt] = barObj._barLevelrefFFGround(spos, rpos, edgeVertex, beta, phis, phir, self.fr)
                LdT.append(Lt)
                [Lf, Df] = barObj._barLevelrefFFGround(spos, rpos, [N1point, N2point], [1.5*np.pi, 1.5*np.pi], phis, phir, self.fr)
                LdF.append(Lf)
                
                # for fitfun2
                w23 = dist2D(Nmpoint, rpos)
                beta = 2.0*np.pi-alfa2
                nu = np.pi/beta
                Mv = (np.cos(nu*np.pi) - np.cos(nu*beta))/(nu*np.sin(nu*np.pi))
                rs = dist2D(spos, N1point)
                rr = dist2D(rpos, N2point)
                B = (self._BKawai(rs+w23, rr, w23)+self._BKawai(rs, rr+w23, w23))/2
                gama = self._AoverLambda(w23+rs, w23+rr, rs+rr+2.0*w23, self.waveLen)                
                Y.append(B*gama*Mv)                
        Y = np.reshape(Y, [-1, 1])
        Abar = np.reshape(LdT, [-1, 1])              
        LdFarray = np.reshape(LdF, [-1, 1])
        p0 = [0.3, 2]
        assert (len(Abar)==len(LdFarray))
        p1, success = optimize.leastsq(errFun2, p0,\
            args=(Y[:,0],  Abar[:, 0]-LdFarray[:, 0]))
        print '\n\nP[0]/sqrt(2)*(p[1]/(Y+p[1])), p: ', p1, "\n\n"
        
        # modify on 31/10/2012. to calculate mse
        np.savetxt('BY.txt', Y[:,0])
        
        #plot results
        for n, w in enumerate(widths):
            plt.figure()
            lgs = []
            for d in dist2barrier:
                lgs.append(str(d)+'m')
                lgs.append(str(d)+'m-fit')
                rx = 10.0+w+d
                spos = [w-d, 0]
                rpos = [rx, self.ry]
                N1point = [10.0, self.h]
                N2point = [10+w, self.h]
                Nmpoint = [(N1point[0]+N2point[0])/2.0, N1point[1]+5.0]
                edgeVertex = [N1point, Nmpoint, N2point]
                phis = math.atan(d/self.h)
                phir = math.atan(d/(self.h-rpos[1]))                
                alfa = math.atan(5.0/w*2)
                alfa2 = np.pi - 2*alfa
                beta = [1.5*np.pi-alfa, 2*np.pi-alfa2, 1.5*np.pi-alfa]
                barObj = AttenuationBar()
                [LdT,Dt] = barObj._barLevelrefFFGround(spos, rpos, edgeVertex, beta, phis, phir, self.fr)
                [LdF, Df] = barObj._barLevelrefFFGround(spos, rpos, [N1point, N2point], [1.5*np.pi, 1.5*np.pi], phis, phir, self.fr)
                
                # for fitfun2
                w23 = dist2D(Nmpoint, rpos)
                beta = 2.0*np.pi-alfa2
                nu = np.pi/beta
                Mv = (np.cos(nu*np.pi) - np.cos(nu*beta))/(nu*np.sin(nu*np.pi))
                rs = dist2D(spos, N1point)
                rr = dist2D(rpos, N2point)
                B = (self._BKawai(rs+w23, rr, w23)+self._BKawai(rs, rr+w23, w23))/2
                gama = self._AoverLambda(w23+rs, w23+rr, rs+rr+2.0*w23, self.waveLen)    
                Y = B*gama*Mv
                plt.semilogx(self.fr, LdT-LdF, '-', self.fr, fitFun2(p1, Y), '--')
                
#                print "LdT-LdF: ", LdT-LdF
#            plt.xlim([25, 10000])
            plt.grid()
            plt.ylabel("Lslope - Abar [dB]")
            plt.xlabel("Abar [dB]")
            plt.title("width"+str(w))
            plt.legend(["0.5m", "0.5m Fit", "1m", "1m Fit", \
            "2m", "2m Fit","3m","3m Fit",], loc="best")
        
    
        

if __name__=='__main__':   
    obj = FitAroof()
    obj._fitSlopeRoofNoCanyon()
    plt.show()
