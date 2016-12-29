# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 17:38:18 2013
calculate multiple reflections by using the united diffraction theory
@author: Wei
"""

import numpy as np
import math
import sys
sys.path.append(r'S:\software\intern\script_wgw\special')
from specialUnstandardLib import  trigamma      
sys.path.append(r'S:\software\intern\script_wgw\Multiple_diffraction')
import MDF_v2 as MDF
sys.path.append(r'S:\projecten\QSIDE\simulations\Project_scripts\mainClass')  # used for check
import QSIDE
import matplotlib.pylab as plt


class MultipleReflection(object):
    def __init__(self):
        ''' according to the analysis, suppose the reflections are a great number
            then all the values are the limitation when rs->oo. So the necessary
            functions are start with "_lim"
        '''
        pass
    def _limB(self, propaPath):
        ''' used to calculate the B values when the rs approaches to infinite
        
            propaPath = length of the propagation path. since the distance from
            source to the first diffraction point is infinite(rs->oo), the 
            propaPath MUST NOT include this. which is a little different with
            the propaPath in MDF_v2. 
        '''
        limBs = []
        for n in xrange(len(propaPath)-1):
            limBs.append(propaPath[n]/(sum(propaPath[n::])))
        limBs.append(1.)
        return limBs
        
    def _limA(self, propaPath):
        ''' used to calculate the A values when the rs approaches to infinite
        '''
        limAs = []
        for n in xrange(len(propaPath)):
            limAs.append(2.*sum(propaPath[n::]))
        return limAs
    
    def _limMvAlfa(self, nupi, infiniteAngle):
        ''' calculate the Mv(alfa) when rs->oo. When rs->oo, the angle between 
            the left edge of the diffraction piont and the line connecting the 
            diffraction point and the source will approach to alfai1, and the 
            angle between the left edge of the diffraction piont and the line 
            connecting the diffraction point and the receiver is alfai2. 
            
            nupi = diffraction angle. e.g. rectangle is 1.5*pi
            infiniteAngle = alfai2 - alfai1 
        '''
        assert(len(nupi)==len(infiniteAngle))
        limMvAlfa = []
        for n, a in enumerate(infiniteAngle):
            nu = np.pi/nupi[n]             
            limMvAlfa.append(abs((math.cos(nu*np.pi)-math.cos(nu*a))/(nu*math.sin(nu*np.pi))))
        return limMvAlfa 
    
    def _Xr(self, thetaBig, A, nuPi):
        ''' calculate X- in Kawai's paper
            thetaBig = the angle in equ(44)
            A = self._AKawai(p1, p2, L)
            nuPi = beta. the diffraction angle
            when compare X1-, X2- and X3-, wave number k will be a constant, so 
            here the we remove the k when calculate Xminus
        '''        
        if thetaBig<(np.pi-nuPi):
            Nminus = -1.0
        if thetaBig>=(np.pi-nuPi) and thetaBig<=(np.pi+nuPi):
            Nminus = 0.0
        if thetaBig>(np.pi+nuPi):
            Nminus = 1.0
        return 2.0*A*np.cos((2.0*Nminus*nuPi - thetaBig)/2.0)**2.0    
        
    def _limX(self, limB, limA, limMvAlfa, infiniteAngle, nupi, waveLen):        
        XMinus = []
        for n in range(len(limA)):
            XMinus.append(self._Xr(infiniteAngle[n], limA[n], nupi[n]))        
        Bk = np.array([1.]*len(limB))
        for xc in xrange(len(XMinus)-1):
            if XMinus[xc]<=XMinus[xc+1]:
                Bk[xc] = Bk[xc]*limB[xc]
                Bk[xc+1] = Bk[xc+1]*1.0
            else:
                Bk[xc] = Bk[xc]*1.0
                Bk[xc+1] = Bk[xc+1]*limB[xc]            
        limX = []
        for m in xrange(len(limA)):            
            limX.append(np.sqrt(Bk[m]*limA[m]/waveLen)*limMvAlfa[m])        
        return limX
    
    def _C1s(self, limX):
        N = len(limX)
        C1s = 1.
        for n in xrange(N):
            limXi = limX[n]
            C1s *= (0.4/(limXi+0.4))**2
        C1s *= 2.* (1./np.sqrt(2))**N
        return C1s        
        
    def _C2s(self, ws1, ws, propaPath):
        ''' calculate C2s. Since rs->oo, the propaPath does not include rs
            ws1 = distance from souce to the left facade(opposite direction to
            the source)
            ws = the width of the canyon
        '''
        C2s = ws1 - ws + sum(propaPath)
        return C2s
        
    def _sumPdiffSquare(self, ws1, ws, propaPath, nupi, infiniteAngle, waveLen):        
        limB = self._limB(propaPath)
        limA = self._limA(propaPath)
        limMvAlfa = self._limMvAlfa(nupi, infiniteAngle)
        limX = self._limX(limB, limA, limMvAlfa, infiniteAngle, nupi, waveLen)
        C1s = self._C1s(limX)
        C2s = self._C2s(ws1, ws, propaPath)        
        sumPdiffSquare = C1s/((4*np.pi)**2)*trigamma(1+C2s/(2*ws))/(4*ws**2)
        return sumPdiffSquare        
            
    def _realSourceOverHi(self, spos, rpos, propaPath, nupi, thetaSn, thetan, fr):
        objN = MDF.MultidiffractionN()
        [atteN, distN] = objN._calDiffrN(spos, rpos, propaPath, nupi, thetaSn, thetan, fr)
        return  [atteN, distN]    
        
    def _execute0(self):
        # case 6-12
        gds = ['k-', 'm-', 'c-', 'r-', 'b-', 'g-', 'y-']
        gds2 = ['k--', 'm--', 'c--', 'r--', 'b--', 'g--', 'y--']
        drct = r'S:\projecten\QSIDE\simulations\Project_scripts\PieceEquation\noGround'
        folder = 'rectangle_null6-12_64000'
        rID = 'yverloop1'
        fr = [31.5, 63, 125, 250, 500, 1000]
        obj = QSIDE.FDTD_2D(drct, folder, rID)
        [rcvSplFF1, fr2] = obj._simMinusFF('oct')
        
        spos = [10., 1.]
        rposen = np.loadtxt(drct+'\\'+folder+'\\'+rID+'.positions.txt')
        ws = 20.
        ws1 = 10.
        N0p = [0., 6.]
        N1p = [20., 12.]
        N2p = [30., 12.]    
        fout = open('6-12energy2.txt', 'wb')         
        for n in xrange(rposen.shape[0]):
            rpos = rposen[n, :]
            [atteN, distN] = self._realSourceOverHi(spos, rpos, [MDF.dist2D(spos, N1p), 10., MDF.dist2D(N2p, rpos)],\
                [1.5*np.pi, 1.5*np.pi], [math.atan(10./11.), 0], [1.5*np.pi, np.pi+math.atan((12.-rpos[1])/(rpos[0]-30.))], fr)
            pdiff1 = 10.**(0.1*atteN)*1./((4*np.pi)**2)*(1./distN**2)
            
            # 1st order imagesource <-
            [atteN2, distN2] = self._realSourceOverHi([-ws1, 1.], rpos, [MDF.dist2D([-ws1, 1], N1p), 10., MDF.dist2D(N2p, rpos)],\
                [1.5*np.pi, 1.5*np.pi], [math.atan(30./11.), 0], [1.5*np.pi, np.pi+math.atan((12.-rpos[1])/(rpos[0]-30.))], fr)
            pdiff2 = 10.**(0.1*atteN2)*1./((4*np.pi)**2)*(1./distN2**2)            
            
             # 1st order imagesource ->
            [atteN3, distN3] = self._realSourceOverHi([30., 1.], rpos, [MDF.dist2D([30., 1.], N1p), 10., MDF.dist2D(N2p, rpos)],\
                [2.*np.pi, 1.5*np.pi], [math.atan(11./10.), 0], [2.*np.pi, np.pi+math.atan((12.-rpos[1])/(rpos[0]-30.))], fr)
            pdiff3 = 10.**(0.1*atteN3)*1./((4*np.pi)**2)*(1./distN3**2)   
            
            r3r = MDF.dist2D(N2p, rpos)
            alfa2 = np.pi+math.atan((12.-rpos[1])/(rpos[0]-30.))  # angle difference in close to the receiver
            alfa12 = math.atan(6./20.)
            
            r3r = MDF.dist2D(N2p, rpos)
            alfa2 = np.pi+math.atan((12.-rpos[1])/(rpos[0]-30.))  # angle difference in close to the receiver
            alfa12 = math.atan(6./20.)
            sumPdiffSquare1 = self._sumPdiffSquare(ws1, ws, [MDF.dist2D(N0p, N1p), 10., r3r], \
                [2.*np.pi, 1.5*np.pi, 1.5*np.pi], [alfa12-2.*np.pi, np.pi, alfa2], 340./np.array(fr))
            
            r3r = MDF.dist2D(N2p, rpos)
            alfa3 = np.pi+math.atan((12.-rpos[1])/(rpos[0]-30.))  # angle difference in close to the receiver
            sumPdiffSquare2 = self._sumPdiffSquare(ws1, ws, [10., r3r], \
                [2.*np.pi, 1.5*np.pi], [2.*np.pi, alfa3], 340./np.array(fr))
            
            sumAllPdiff = pdiff1 + pdiff2 + pdiff3 + sumPdiffSquare1 + sumPdiffSquare2
            totalLevelAtte = 10.*np.log10(sumAllPdiff*(4*np.pi)**2 *distN**2)         
            
            plt.semilogx(fr, totalLevelAtte, gds2[n], \
            fr, rcvSplFF1[:, n]+10.0*np.log10(distN/MDF.dist2D(spos, rpos)),gds[n])
        fout.close()  


def calcualte_argin(rs, rr, w, ts, tr, waveLen):
    L = rs + rr + w
    Ys = np.sqrt(3)*(np.cos(2.0/3.0*ts)-0.5) * np.sqrt(2*rs*(w+rr)/(waveLen*L))
    Yr = np.sqrt(3)*(np.cos(2.0/3.0*tr)-0.5) * np.sqrt(2*rr*(w+rs)/(waveLen*L))
    Ygr = []
    Ysm = []
    for n, wlambda in enumerate(waveLen):
        if Ys[n] > Yr[n]:
            Ygr.append(np.float(Ys[n]))   # greater Y
            Ysm.append(np.float(np.sqrt(3.0)*(np.cos(2.0/3.0*tr)-0.5) * np.sqrt(2*rr*w/(wlambda*(w+rr)))))   # smaller Y
        else:
            Ygr.append(np.float(Yr[n]))   # greater Y
            Ysm.append(np.float(np.sqrt(3.0)*(np.cos(2.0/3.0*ts)-0.5) * np.sqrt(2*rs*w/(wlambda*(w+rs)))))   # smaller Y
    return [Ygr, Ysm]
    
def QSIDE_fit(Hs, Hi, Hr, N1point, N2point, Ws,  Wi,  Wr, spos, rpos, fr, rho1=0.97, rho2=0.97):
        """ N1point is building vertix in source canyon (x, y)
            N2point is building vertex in receiver canyon (x, y)
            spos is the source position
            rpos is the receiver position. 
            Be sure spos, rpos, N1point and N2point are under the same coordinate system
        """        
        # modified on Feb 1st 2012, June 25h
        betas = 1.5*np.pi
        betar = 1.5*np.pi
        sHi = spos[1]  # source height
        rHi = rpos[1]  # receiver height
        h1 = Hi-spos[1]
        h2 = Hi-rpos[1]
        xs = N1point[0]-spos[0]
        xr = rpos[0]-N2point[0]
        rs = np.sqrt(xs**2.+h1**2.)
        rr = np.sqrt(xr**2.+h2**2.)
        distSR2 = np.abs(rpos[1]-spos[1])**2.+(xs+Wi+xr)**2.
        phis = math.acos(h1/rs)
        phir = math.acos(h2/rr)        
        waveLength = 340./np.array(fr)
        if sHi==0:
            sHi = 0.01
        if rHi==0:
            rHi = 0.01
        if type(fr)==list:
            fr = np.array(fr) 
        if N1point[1]!=N2point[1]:
            N1point[1] = (N1point[1]+N2point[1])/2.0
            N2point[1] = (N1point[1]+N2point[1])/2.0
        
        # write the arguments out. modified on April 3rd, 2013
        [Ygr, Ysm] = calcualte_argin(rs, rr, abs(N2point[0]-N1point[0]), phis, phir, waveLength)
        
        THKbar = MDF.DoubleDiffraction2D()
        [minusAbarFlat, d] = THKbar._levelrefFF2Ground(spos, rpos, N1point, N2point, betas, betar, phis, phir, fr)
    
    #============ this comment is a simpified immplement of Dick's method
        Wi = np.abs(N1point[0]-N2point[0])
    
        # reflection coefficient alpha=0.97, beta=0.97
        C = 1.5*Ws+Wi+1.5*Wr
        C1s = 1./(3.31*math.cos(phir)*np.sqrt(rr/waveLength)+1.)**2.
        C1r = 1./(3.31*math.cos(phis)*np.sqrt(rs/waveLength)+1.)**2.
        C3s = 3.31*h1*np.sqrt(Wi/waveLength)+0.5*Ws+rr+Wi
        C3r = 3.31*h2*np.sqrt(Wi/waveLength)+0.5*Wr+rs+Wi
            
        if (Hs-sHi)/h1<=1./3.:
            Lhs = np.array([-40.]*len(fr)) 
        if (Hs-sHi)/h1>=1.:
            Lhs = np.array([0.]*len(fr))
        if (Hs-sHi)/h1>1./3. :
            Lhs = -6.17*(1-(Hs-sHi)/h1)*(1-1.37*np.log10(np.sqrt(waveLength*Ws)/Wi))
            if Lhs.any()>0:
                Lhs = np.array([0.]*len(fr))
    
        if (Hr-rHi)/h2<=1./3.:
            Lhr =  np.array([-40.]*len(fr))
        if (Hr-rHi)/h2>=1.:
            Lhr = np.array([0.]*len(fr))
        if (Hr-rHi)/h2>1./3. :
            Lhr = -6.17*(1-(Hr-rHi)/h2)*(1-1.37*np.log10(np.sqrt(waveLength*Wr)/Wi))
            if Lhr.any()>0:
                Lhr = np.array([0.]*len(fr))
    
        minusAcanFlat = 10.5*np.log10(12.46*C1s* rho1**2. * distSR2 * 10.**(0.1*Lhs)/((C3s+Ws)**2.)\
            + 22.24*C1r* rho2**2. * distSR2 * 10.**(0.1*Lhr)/((C3r+Wr)**2.)\
            + 0.05*rho1**2.*rho2**2. * distSR2 * 10.**(0.1*Lhs)*10.**(0.1*Lhr) / ((3.31*h1/np.sqrt(waveLength)+C)*(3.31*h2/np.sqrt(waveLength)+C)))
    #========================================================================================================           
    
        if Hs!=0 or Hr!=0:
            AroofNocan = 0.0
            if Hs!=0 and Hr!=0:
                AroofCan = 5.
            else:
                AroofCan = 2.5
        else:
            AroofNocan = 0.27*minusAbarFlat+2.9
            AroofCan = 0.
        Ainter = np.sqrt(distSR2)/100.0
        if Ainter>5.0:
            Ainter = 5.0
        AbarFlat = -10.*np.log10(10.**(0.1*(minusAbarFlat+AroofNocan)))
        AcanFlat = -10.*np.log10(10.**(0.1*(minusAcanFlat+AroofCan)))
        Adiifr = -10.*np.log10(10.**(-0.1*(AbarFlat))+10.**(-0.1*(AcanFlat)))
        return [Adiifr, Ainter, AbarFlat, AcanFlat]  
    
if __name__=='__main__':
    if len(sys.argv)>=1:
        obj = MultipleReflection()
        obj._execute0()
        plt.xlim([25, 1500])
        fr = [31.5, 63, 125, 250, 500, 1000]     
        plt.grid()
        plt.xticks(fr, fr)
        plt.legend(['Calculated', 'Simulated(FDTD)'], loc='best')
        plt.show()