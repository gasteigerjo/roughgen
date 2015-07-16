# Created by Hugo Cruz Jimenez, August 2011, KAUST
#def SpecSyn2(N,samp,corr,acf,Rseed):
import scipy.special
import datetime
import time
import random
import math
import numpy as np
from numpy.fft import fft2, ifft2
from numpy.linalg import lstsq
from val import *

def SpecSyn2(*args):
    #function [Y,spar,spec,ierr] = SpecSyn2(N,samp,corr,acf,Rseed)
    #  [Y,spar,spec,ierr] = SpecSyn2(N,samp,corr,'acf',Rseed) 
    #  generates a 2D-random field Y of size (Nz+1 x Nx+1) with
    #  possibly variable spatial sampling in z,x-direction. 
    #  This function simulates anisotropic random fields, i.e. 
    #  the correlation len in both directions can be different 
    #  rectangular grid dimensions are also possible.
    #  corr contains the corr. len ax, az and the Hurstnumber H 
    #  or the fractal dimension D the autocorrelation function to 
    #  be used has to be specified by acf
    #  NOTE: D = 3-H, and larger D (smaller H) yield "rougher" fields
    #
    #  The algorithm is based on the spectral synthesis method by 
    #  Pardo-Iguzquiza, E. and Chica-Olma, M. (1993)
    #  The Fourier integral method: and efficient spectral method for 
    #  simulation of random fields, Mathematical Geology, 25, p177-217.
    #  but extends their method to handle rectangular grids and variable
    #  sampling in the two directions.
    #
    #  INPUT:
    #  N     - grid dimensions [Nz Nx]
    #  samp   - desired sampling, [dz dx] if scalar, then dz = dx
    #  corr  - corr = [az ax]   for 'gs' and 'ex' 
    #     corr = [az ax H] for 'ak' note that 0 <= H <= 1
    #     corr = [D kc]    for 'fr' D is the fractal dimension,
    #      kc: corner wavenumber, spectrum decays linearly for k>kc
    #  acf   - autocorrelation function: 
    #    'gs' or 'GS' - Gaussian
    #    'ex' or 'EX' - Exponential
    #    'ak' or 'AK' - anisotropic vonKarman
    #    'fr' or 'FR' - fractal distribution
    #  Rseed - seeds for random number generators if omitted or empty
    #     Rseed = sum(100*clock) is used (returned in structure spar)
    #     [Rpseed Rsseed] for the phase and small random spectral part
    # 
    #  OUTPUT:
    #  Y    - Y = [z x] random field whose size is determined by the sample
    #    spacing in each direction, the number of points and whether
    #    the pow2-option is given or not. 
    #  spar  - structure with len vectors, sampling in both direction
    #    and other relevant parameters it also contains the random
    #    seed number, useful to reproduce realizations
    #  spec - structure containing the computed power spectrum as well
    #         wavenumber vectors
    #  ierr - 0 when successfully executed 
    #         1 when error in Z,X sampling
    #
    #  Written by Martin Mai (martin@seismo.ifg.ethz.ch) 
    #  originally from 07/16/98, based on SRB-toolbox (Ph. Rio)
    #  last changes 03/01/2000 Nov. 2002;
    # ------------------------------------------------
    
    ierr  = 0              #% error variable
    check = 'n'  	#% set to 'y' if you want to create
      		#% a simple out put plot to check the
      		#% the spectra and the resulting field
    
    #%% check input variables
    if len(args) < 4: print '  Error *** Not enough input arguments ***'
    elif len(args) == 4:
        N=args[0] 
        samp=args[1] 
        corr=args[2] 
        acf=args[3] 
        Rseed = {} 
    elif len(args) == 5:
        N=args[0] 
        samp=args[1] 
        corr=args[2] 
        acf=args[3] 
        Rseed =args[4]
    
    if len(samp) == 1: samp = [samp, samp]
    if len(N) == 1: N = [N, N]
    
    
    #%% error checking on inpur array size and given sampling
    if np.mod(N[1],samp[1]) != 0:
        ierr = 1
        print '** sampling in X does not yield an integer number **'
        print '   of grid points  '
        print '==> BOOM OUT in SpecSyn2<=='
        return

    if np.mod(N[0],samp[0]) != 0:
        ierr = 1
        print '** sampling in Z does not yield an integer number **'
        print '   of grid points ==> abort!'
        print '==> BOOM OUT in SpecSyn2<=='
        return
    
    
    #%% get data values on the correlation len/fractal dimension
    if  acf == 'fr' or acf == 'FR':
        if len(corr) == 2:
            D = corr[0]; kc = corr[1]
        else:
            D = corr[0]; kc = 0.1 
            print('** Corner wavenumber kc not given: set to 0.1 **')
    elif acf == 'ex' or acf == 'EX' or acf == 'gs' or acf == 'GS':
        ax = corr[1]; az = corr[0] 
    elif acf == 'ak' or acf == 'AK':
        ax = corr[1]; az = corr[0]; H = corr[2]
    
    
    #%% set size for spectral synthesis tool that generates
    #%% fields of size (2*rmz+1, 2*rmx+1), i.e. the method 
    #%% requires an ODD number of points in each direction
    nptsX = round(N[1]/samp[1])  #% number of grid-points in X
    nptsZ = round(N[0]/samp[0])  #% number of grid-points in Z
   
    if np.mod(nptsX,2) == 0: rmx = nptsX/2
    elif np.mod(nptsX,2) != 0: rmx = (nptsX-1)/2
    
    if np.mod(nptsZ,2) == 0: rmz = nptsZ/2
    elif np.mod(nptsZ,2) != 0: rmz = (nptsZ-1)/2
    
    
    #%% compose power spectrum for two of the four quadrants
    #%% wavenumber vector in [-pi,pi]
    kx = (rmx/((2*rmx+1)*samp[1]))*np.linspace(-2*np.pi,2*np.pi,2*rmx+1)
    kz = (rmz/((2*rmz+1)*samp[0]))*np.linspace(-2*np.pi,2*np.pi,2*rmz+1)
    rk = [rmz+1,2*rmx+1]
    kr = np.zeros(rk)
    k1 = np.zeros(rk)
  
    for j in range(0,int(2*rmx+1)):        
        for i in range(0,int(rmz+1)):        
            if acf == 'fr' or acf == 'FR':
                kr[i,j] = np.sqrt((kz[i]**2) + (kx[j]**2))
            else:
                kr[i,j] = np.sqrt((az**2)*(kz[i]**2) + (ax**2)*(kx[j]**2))
    
     #%% calculate power spectral density, depending on selected ACF
    if acf == 'gs' or acf == 'GS':
        PS = 0.5*ax*az*np.exp(-0.25*kr**2)
    elif acf == 'ex' or acf == 'EX':
        PS = (ax * az)/(1 + (kr**2))**1.5
    elif acf == 'ak' or acf == 'AK':
        for j in range(0,int(2*rmx+1)):        
            for i in range(0,int(rmz+1)):        
                k1[i,j]=kr[i,j]
        k3=k1.conj().transpose()
        a1,b1=np.shape(k1)
        k2=k3.reshape(a1*b1,1)
        ka = k2.compress((k2>0).flat)
        # In [14]: print k2.compress.__doc__
        # k2.compress(condition, axis=None, out=None)
        #   .flat  also appears if k2 is np.array 
        c1=np.size(ka)
        ka2=ka.reshape(c1,1)
        coef = 4*np.pi*H*ax*az/scipy.special.kv(H,min(ka))
        PS = coef/(1 + (kr**2))**(H+1)
    #    coef = 4*pi*H*ax*az./besselk(H,min(k))
    #     PS = coef./(1 + (kr.**2)).**(H+1)
    elif acf == 'fr' or acf == 'FR':
        decay = 0.5*(8-2*D)
        #% to ensure proper scaling of the power spectrum 
        #% in k-space we cannot allow kr == 0
        if kr.min() == 0:

            a,b=np.shape(kr)
            for i in range(0,a):
                for j in range(0,b):
                    if kr[i,j]==0:
                        p,q= i,j
## NEVER USED ###       k[p,q] = np.mean(kr[p-2:p,q-2:q])
        #% set values below k< kc to constant pi*kc
        sgival(kr,np.pi*kc,np.pi*kc)
        PS = 1/((kr**2)**decay)         	
    # 
    # #%% the IFFT needs the spectrum normalized to max.amp unity
    PS = PS/PS.max()
    AM = np.sqrt(PS)
    # 
    # 
    # #%% compose the random phase
    # #%% initialize random number generator
    if len(Rseed) == 0:
        Rseed = np.zeros(2)
        Rseed[0] = int(time.time())
        Rseed[1] = Rseed[0]*int(time.time())
      
    random.seed(Rseed[0]) 

    # #%% random phase in [0,pi]
    m,n=np.shape(kr)
    PH=np.zeros((m,n))
    for i in range(0,m):
        for j in range(0,n):
            PH[i,j] = 2*np.pi*random.random()
      
    # #%% assemble the random field in FT-domain  
    # # add small random high-wavenumber components
    # #x = (1 + 0.5*randn(size(kr)))
    x = 1
    RAD = AM*x
      
    # #%% set DC-component to different value, if desired
    # #%% NOTE that this only changes the overall 'level' of
    # #%% the field, not its appearance, but has significant
    # #%% impact on the Fourier transform which reflects the
    # #%% DC-value at the smallest wavenumber ("Nugget Effect")
    Neff = 0                               #% "Nugget" alue
    RAD[rmz,rmx] = Neff  	#% "Nugget" effect, zero-mean field
    AM[rmz,rmx]  = RAD[rmz,rmx]
    m1,n1=np.shape(PH)
    Y=np.zeros((m1,n1))
    for i in range(0,m1):
        for j in range(0,n1):
            Y = RAD*np.cos(PH)+1j*RAD*np.sin(PH)
# Matlab and Python need *np.pi/180.0, but here is not used. Why?
  
    # #%% the following takes care of the conjugate symmetry condition
    # #%% in order to ensure the final random field is purely real
    aa=[2*rmz+1,2*rmx+1]
    U = np.zeros(aa)          #% will be conj. sym. field
    Y = np.concatenate((Y , np.conj(np.fliplr(np.flipud(Y[0:rmz,:]  )))))

    for i in range(0,int(rmx)):
        Y[rmz,-i+2*rmx+2-2] = np.conj(Y[rmz,i])
    
    for i in range(0,int(rmz)+1):
        for j in range(0,int(rmx)+1): U[i,j] = Y[i+rmz,j+rmx]

    for i in range(int(rmz)+1,2*int(rmz)+1):
        for j in range(int(rmx)+1,2*int(rmx)+1): U[i,j] = Y[i-rmz-1,j-rmx-1]
   
    for i in range(0,int(rmz)+1):
        for j in range(int(rmx)+1,2*int(rmx)+1): U[i,j] = Y[i+rmz,j-rmx-1]
    
    for i in range(int(rmz)+1,2*int(rmz)+1):
        for j in range(0,int(rmx)+1):
            U[i,j] = Y[i-1-rmz,j+rmx]

#            return U, Y
#    return U[i,j], Y[i-1-rmz,j+rmx]
#            UU=U[i,j]

#        return U[i,j]
    
    # #%% take 2D-inverse FFT to obtain spatial field imaginary parts
    # #%% of order 1e-13 due to machine precision are removed 
    # #%% also, remove mean and scale to unit variance
    Y = np.real(ifft2(U))
    Y = Y/np.std(Y,ddof=1)  			# standard deviation of unity
    if np.mean(Y) < 0: Y = (-1)*Y  # positive mean (though small)

       
    # #%% due to the requirement of [2*rmz+1,2*rmx+1] grid points, the
    # #%% simulated field may be too large by one point in either direction.
    # #%% This is compensated by resampling
    # #Y = resampgrid(Y,[nptsZ nptsX])
    # 
    # 
    # #%% final output structure with len vectors,sample spacing
    # #%% and model parameters
    spar={}
    spar['dim']   = N
    spar['samp']  = samp 
    spar['size']  = np.shape(Y)
    spar['corr']  = corr
    spar['acf']   = acf
    spar['Rseed'] = Rseed
    a1,a2=np.shape(Y)
#    spar['lx'] = np.arange(0,samp[1]*(a2-1)+1,samp[1])
    spar['lx'] = np.arange(0,samp[1]*(a2),samp[1])
#    spar['lx'] = np.arange(0,samp[1]*(a2-1),samp[1])
    spar['lz'] = np.arange(0,samp[0]*(a1),samp[0])
#    spar['lz'] = np.arange(0,samp[0]*(a1-1)+1,samp[0])

    # #% assemble structure with spectra in pos. quadrant
   
    px=[]
    for i in range(0,np.size(kx)):
        if kx[i]>=0:
            px.append(i)

    pz=[]
    for i in range(0,np.size(kz)):
        if kz[i]<=0:
            pz.append(i)

    spec={}
    spec['PD']  = PS[:,px]
    spec['kpx'] = kx[px]
    spec['kpz'] = kz[pz]
    m,n=np.shape(spec['PD'])
    spec['PDx'] = spec['PD'][m-1,]
    spec['PDz'] = spec['PD'][:,0]
      
    print px,pz
#    return px,pz, spec, m
#    return spec
    return Y,spar
      
#        G,spar = SpecSyn2([W, L],samp,corr,acf,SSseed)

    # #%% plot spectrum for error checking
    # if check == 'y':
    # 
    # figure
    # fullpage
    # subplot(221)
    # imagesc(spec.kpx,spec.kpz,log10(spec.PD)) colorbar
    # xlabel('kx') ylabel('kz');
    # title('2D-Power Spectral Density (in log-units)','FontS',12)
    # axis equal axis tight;
    # 
    # 
    # subplot(222)
    # loglog(-spec.kpz,spec.PDz,'r','LineW',2) hold on
    # loglog(spec.kpx,spec.PDx,'b','LineW',2)
    # if acf == 'fr':
    #    loglog(pi*kc*ones(1,10),linspace(1,min(spec.PD(:)),10),'k--')
    #    tstr = ['kc = ',num2str(kc),', dz, dx = ',...
    #   num2str(samp(1)),', ',num2str(samp(2))]
    # else:
    #   tstr = ['az,ax = ',num2str(az),', ',num2str(ax),', dz, dx = ',...
    #   num2str(samp(1)),', ',num2str(samp(2))]
    # end
    # legend('z','x',3)
    # xlabel('wavenumber')
    # ylabel('psd')
    # title(['1D-spectra (',tstr,')'],'FontS',12)
    # axis([0 max([max(spec.kpx) max(-spec.kpz)]) min(spec.PD(:)) 2])
    # axis square
    # hold off
    # 
    # subplot(212)
    # imagesc(spar.lx,spar.lz,Y) colorbar
    # xlabel('X') ylabel('X');
    # title('Resulting Random Field','FontS',12,'FontW','bo')
    # axis equal axis tight
    # 
    # end
