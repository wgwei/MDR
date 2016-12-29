# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:44:14 2013
calculate the levels in double diffraction case, modified from Pierce paper
Mainly based on the recent study of the high-order diffraction

the diffraction by Pierce, Kawai and previous version could be found in v2
@author: Wei
"""
import numpy as np
import math
import sys, os
from scipy import special
import matplotlib.pylab as plt  # used for check
sys.path.append(r'S:\software\intern\script_wgw\FDTD\FDTD_post_process\beta_11')  # used for check
from post_FDTD import FDTD_2D  # used for check

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
        
def dist2D(p1, p2):
    ''' calculate the distance between p1 and p2'''
    return np.sqrt((p1[0]-p2[0])**2.0+(p1[1]-p2[1])**2.0)    
    
def get_dif_angles(spos, difPoints, rpos, LorR, beta):        
    ''' spos = [x, y]
        difPionts = [[x1, y1], [x2, y2], ...]
        rpos = [x, y]
        LorR = ['L', 'R', 'L', ...] the diffraction is on the left edge or the right edge
        beta = [beta1, beta2, ..], the angles of the diffraction object
    '''
    if  difPoints[0][1]<spos[1]:
        raise 'spos should be lower than the first diffraction edge!'
    else:
        thetasList, thetarList = [], []
        allPts = [spos]+list(difPoints)+[rpos]
        for d, dpt in enumerate(difPoints):
            if beta[d]!= 2*np.pi:
                if LorR[d] =='L':
                    agi = math.atan(abs(dpt[0]-allPts[d][0])/abs(dpt[1]-allPts[d][1]))                    
                    ag2 = math.atan(abs(allPts[d+2][1]-dpt[1])/abs(allPts[d+2][0]-dpt[0]))
                elif LorR[d] =='R':
                    agi = math.atan(abs(dpt[1]-allPts[d][1])/abs(dpt[0]-allPts[d][0]))
                    ag2 = math.atan(abs(allPts[d+2][0]-dpt[0])/abs(allPts[d+2][1]-dpt[1]))      
                thetasList.append(1.5*np.pi-agi)
                thetarList.append(ag2)
            else:
                agi = math.atan(abs(dpt[0]-allPts[d][0])/abs(dpt[1]-allPts[d][1]))
                ag2 = math.atan(abs(allPts[d+2][0]-dpt[0])/abs(allPts[d+2][1]-dpt[1]))      
                thetasList.append(2.*np.pi-agi)
                thetarList.append(ag2)
        return [thetasList, thetarList]

def get_dif_path(difPoints):
    wList = []
    for d in range(len(difPoints)-1):
        wList.append(dist2D(difPoints[d], difPoints[d+1]))
    return wList
                    
class Multiple_DiffractionN():
    def __init__(self):
        pass
    def complex_single_wedge_dif_fun(self, rs, rr, waveLen, theta, theta0, beta):
        [Xplus, Xminus] = self.Xplus_Xminus( rs, rr, waveLen, theta, theta0, beta)
        return np.exp(np.pi/4*1j)/np.sqrt(2)*(self.f_simple(abs(Xplus))+self.f_simple(abs(Xminus)))
        
    def complex_double_edge_dif_fun(self, rs, rr, w, waveLen, thetas, thetar, betas, betar):
        [Xs, Xr] = self.Xs_Xr_double_edge(rs, rr, w, waveLen, thetas, thetar, betas, betar)
        return 1j*(0.37/(0.37+Xs))*(0.37/(0.37+Xr))        
        
    def complex_multiple_dif_fun(self, rs, rr, wList, waveLen, thetasList, thetarList, betaList, adjacentEdgePair=None):
        ''' waveLen = wave length . its a NUMBER instead of an array or a list
            wList is the [w12, w23, w34, ..., wn,n+1] EXCLUDING rr
            thetasList = angle. From the righe diffraction edge to the connecting line
                        between the diffraction point to the source. The source could 
                        be a previous diffraction edge
            thetasList = angle. From the right diffraction edge to the connecting line
                        between the diffraction point to the receiver. The receiver could
                        be a next diffraction edge
            betaList  = angles. Diffraction angles
        '''
        wList = list(wList)+[rr]            
        BiGammai = np.ones(len(wList))
        for i,w in enumerate(wList):                
            BiGammai[i] = np.sqrt(2.*w*(rs+sum(wList[0:i]))/(waveLen*(rs+sum(wList[0:i+1]))))
        nu = np.pi/np.asarray(betaList)
        thetaPlus = np.asarray(thetarList) + np.asarray(thetasList)
        thetaMinus = np.asarray(thetarList) - np.asarray(thetasList)
        MvPlus = (np.cos(nu*np.pi) - np.cos(nu*(thetaPlus)))/(nu*np.sin(nu*np.pi))
        MvMinus = (np.cos(nu*np.pi) - np.cos(nu*(thetaMinus)))/(nu*np.sin(nu*np.pi))
        complexDifFun = 1.
        for i, bg in enumerate(BiGammai):
            complexDifFun *= np.exp(1j*np.pi/4.)/np.sqrt(2)*(0.37/(0.37+bg*abs(MvPlus[i])) + 0.37/(0.37+bg*abs(MvMinus[i])))
        return complexDifFun  
    
    def complex_multiple_dif_fun_kawai(self, rs, rr, wList, waveLen, thetasList, thetarList, betaList, adjacentEdgePair=None):
        ''' waveLen = wave length . its a NUMBER instead of an array or a list
            wList is the [w12, w23, w34, ..., wn,n+1] EXCLUDING rr
            thetasList = angle. From the righe diffraction edge to the connecting line
                        between the diffraction point to the source. The source could 
                        be a previous diffraction edge
            thetasList = angle. From the right diffraction edge to the connecting line
                        between the diffraction point to the receiver. The receiver could
                        be a next diffraction edge
            betaList  = angles. Diffraction angles
        '''
        wList = list(wList)+[rr]            
        BiGammai = np.ones(len(wList))
        for i,w in enumerate(wList):                
            BiGammai[i] = np.sqrt(2.*w*(rs+sum(wList[0:i]))/(waveLen*(rs+sum(wList[0:i+1]))))
        nu = np.pi/np.asarray(betaList)
        thetaPlus = np.asarray(thetarList) + np.asarray(thetasList)
        thetaMinus = np.asarray(thetarList) - np.asarray(thetasList)
        MvPlus = (np.cos(nu*np.pi) - np.cos(nu*(thetaPlus)))/(nu*np.sin(nu*np.pi))
        MvMinus = (np.cos(nu*np.pi) - np.cos(nu*(thetaMinus)))/(nu*np.sin(nu*np.pi))
        complexDifFun = 1.
        for i, bg in enumerate(BiGammai):
            complexDifFun *= np.exp(1j*np.pi/4.)/np.sqrt(2)*(0.37/(0.37+bg*abs(MvPlus[i])) + 0.37/(0.37+bg*abs(MvMinus[i])))
        return complexDifFun  
        
    def Xs_Xr_double_edge(self,rs, rr, w, waveLen, thetas, thetar, betas, betar):
        gammaS = np.sqrt(2.*rs*(rr+w)/(waveLen*(rs+rr+w)))
        nuS = np.pi/betas
        MvS = (np.cos(nuS*np.pi) - np.cos(nuS*(betas-thetas)))/(nuS*np.sin(nuS*np.pi))
        Xs = gammaS*abs(MvS)
        gammaR = np.sqrt(2.*rr*(rs+w)/(waveLen*(rs+rr+w)))
        nuR = np.pi/betar
        MvR = (np.cos(nuR*np.pi) - np.cos(nuR*(betar-thetar)))/(nuR*np.sin(nuR*np.pi))
        Xr = gammaR*abs(MvR)
        B = np.sqrt(w*(rs+rr+w)/((w+rs)*(w+rr)))
        boolX = Xs<Xr
        try:
            if boolX:
                Xs = B*Xs
            else:
                Xr = B*Xr
        except: 
            for n, b in enumerate(boolX):
                if b:
                    Xs[n] = B*Xs[n]
                else:
                    Xr[n] = B*Xr[n]
        return [Xs, Xr]
            
    def Xplus_Xminus(self, rs, rr, waveLen,  theta, theta0, beta):
        ''' calculate X+ and X-'''
        tao = np.sqrt(2.*rs*rr/(waveLen*(rs+rr)))
        nu = np.pi/beta
        Mvplus = (np.cos(nu*np.pi) - np.cos(nu*(theta+theta0)))/(nu*np.sin(nu*np.pi))
        Xplus = tao*Mvplus
        Mvminus = np.cos(nu*np.pi) - np.cos(nu*(theta-theta0))/(nu*np.sin(nu*np.pi))
        Xminus = tao*Mvminus
        return [Xplus, Xminus]
        
    def f_simple(self, X):
        return 0.37/(0.37+X)
    
            
    
if __name__=='__main__':
    if 0:
        # test:
        # single diffraction
        obj = Multiple_DiffractionN()
        rs = 10.
        rr = 10.
        waveLen = 1.
        theta0 = np.pi*5./3.
        beta = np.pi*11./6.
        anglesGrad = np.arange(0, 5./3.*np.pi, 0.1)
        sq = np.ones(len(anglesGrad))
        sqCalc = np.ones(len(anglesGrad))
        for n, theta in enumerate(anglesGrad):
            [xplus, xminus] = obj.Xplus_Xminus(rs, rr, waveLen,  theta, theta0, beta)
            sq[n] = 0.5 * ((f(abs(xplus))+f(abs(xminus)))**2. + (g(abs(xplus))+g(abs(xminus)))**2.)
            complexDf = obj.complex_single_wedge_dif_fun(rs, rr, waveLen,  theta, theta0, beta)
            sqCalc[n] = abs(complexDf)**2.
        plt.plot(anglesGrad, 10.*np.log10(sq), '--', anglesGrad, 10.*np.log10(sqCalc), '-')
        plt.show()   
        
    if 0:
        # test double diffraction
        obj = Multiple_DiffractionN()
        rs = 10.
        rr = 10.
        w = 10.
        waveLen = 1
        thetas = np.pi/4
        theta = np.arange(0, 0.5*np.pi, 0.1)
        betas = 1.5*np.pi
        betar = 1.5*np.pi
        dFL = np.ones(len(theta))
        pDiffr = np.ones(len(theta))
        for n, thetar in enumerate(theta):
            [Xs, Xr] = obj.Xs_Xr_double_edge(rs, rr, w, waveLen, thetas, thetar, betas, betar)
            dF = obj.complex_double_edge_dif_fun(rs, rr, w, waveLen, thetas, thetar, betas, betar)
            dFL[n] = 20.*np.log10(abs(dF))
            pDiffr[n] = (10*np.log10((f(abs(Xs))**2 + g(abs(Xs))**2)*(f(abs(Xr))**2 + g(abs(Xr))**2)))        
        plt.plot(theta, dFL, '-', theta, pDiffr, '--')
        
         # test double diffraction
        obj = Multiple_DiffractionN()
        rs = 10.
        rr = 10.
        wList = [10.]
        waveLen = 1
        thetasList = [np.pi/4+np.pi, 1.5*np.pi]
        theta = np.arange(0, 0.5*np.pi, 0.1)
        betaList = [1.5*np.pi, 1.5*np.pi]
        dFL = np.ones(len(theta))
        pDiffr = np.ones(len(theta))
        for n, thetar in enumerate(theta):
            thetarList = [0, thetar]
            [Xs, Xr] = obj.Xs_Xr_double_edge(rs, rr, wList[0], waveLen, 0.25*np.pi, thetar, betaList[0], betaList[1])
            dF = obj.complex_multiple_dif_fun(rs, rr, wList, waveLen, thetasList, thetarList, betaList, adjacentEdgePair=None)
            dFL[n] = 20.*np.log10(0.5*abs(dF))
            pDiffr[n] = (10*np.log10((f(abs(Xs))**2 + g(abs(Xs))**2)*(f(abs(Xr))**2 + g(abs(Xr))**2)))        
        plt.figure()
        plt.plot(theta, dFL, '-', theta, pDiffr, '--')
        plt.show()
    
    if 0:
        # test multiple diffraction
        dire2 = r'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims'
        simFile = 'BBB-6-12-10-smp64000-wh'
        obj2 = FDTD_2D(dire2, simFile, 'xverloop2')
        [rcvSplFF1, fr] = obj2._simMinusFF('oct')
        objMD = Multiple_DiffractionN()
        waveLen = 340./np.asarray(fr)
        wList = [6., np.sqrt(961)]
        betaList = [1.5*np.pi, 1.5*np.pi, 1.5*np.pi]
        rs = np.sqrt(2.)*5.
        spos = [5.00, 1.00, -0.010]
        rrArray = np.loadtxt(os.path.join(dire2, simFile, 'xverloop2.positions.txt'))
        for r, rpos in enumerate(rrArray):
            d1 = np.sqrt((spos[0]-rpos[0])**2.+(spos[1]-rpos[1])**2.)
            rr = np.sqrt((rpos[0]-47.)**2. + (10.-rpos[0])**2.)
            level = np.ones(len(waveLen))
            d2 = rs+sum(wList)+rr
            for w,wavLambda in enumerate(waveLen):
                thetasList = [1.5*np.pi-math.atan(5./23.), 1.5*np.pi, 1.5*np.pi-math.atan(14./31.)]
                thetarList = [0., math.atan(31./14.), math.atan((rpos[0]-47)/(10.-rpos[1]))]
                complexDifFun = objMD.complex_multiple_dif_fun(rs, rr, wList, wavLambda, thetasList, thetarList, betaList, adjacentEdgePair=None)
                level[w] = 20.*np.log10(0.5*abs(complexDifFun))+10.*np.log10(d1/d2)
            plt.figure()
            plt.semilogx(fr, level, '-', fr, rcvSplFF1[:, r], '--')
            plt.xlim([25, 1500])
        plt.show()            
        
    if 0:
        # test function get_dif_angles(spos, difPoints, rpos, LorR, beta)
        spos = [5.00, 1.00]
        difPoints = [[10., 6.], [26., 12.], [33., 12.], [47., 10.]]
        beta = [1.5*np.pi, 1.5*np.pi, 1.5*np.pi, 1.5*np.pi]            
        LorR = ['L', 'L', 'R', 'R']
        [thetasList, thetarList] = get_dif_angles(spos, difPoints, [49, 2], LorR, beta)
        for tr in thetarList:
            print 'thetar -> %0.1f' %(tr/6.28*360.)
        for ts in thetasList:
            print 'thetas -> %0.1f' %(ts/6.28*360.)
        wList = get_dif_path(difPoints)
        print wList
            
    if 1:
        obj = Multiple_DiffractionN()
        rs = 10.
        rr = 10.
        wList = [10., 12.]
        waveLen = 1
        thetasList = [1.25*np.pi, 0.25*np.pi, 1.2*np.pi]
        thetarList =  [0, 0.3*np.pi, 0.25*np.pi ]
        betaList = [1.5*np.pi, 1.5*np.pi, 1.5*np.pi]
        dF = obj.complex_multiple_dif_fun(rs, rr, wList, waveLen, thetasList, thetarList, betaList, adjacentEdgePair=None)
        print dF
    
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            