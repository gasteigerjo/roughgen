import numpy as np
import matplotlib.pyplot as plt
from math import sqrt,exp,atan,cos,sin,log10
import random
random.seed('0254887388')
#random.seed()

#######PARAMETERS#######################
#fault Ls*Ld km, minimum roughness considered lambdaMin km
Ls=40.
Ld=20.
L=max(Ls,Ld)

lambdaMin=1.
lambdaMaxS=40.
lambdaMaxD=20.

#for a good discretisation of the largest wavelength (all directions)
#(in any direction)
#we want that 1 << k(lambdaMax)= (here 8)
L=max(Ls,Ld,8*max(lambdaMaxS,lambdaMaxD))

#amplitude to wavelength ratio
alpha=pow(10.,-1.9)
#Hurst index
H=0.8
########################################

N=512*8

Nmax=int(L/lambdaMin)
Pmax=Nmax
Nmin=int(L/lambdaMaxS)
Pmin=int(L/lambdaMaxD)
print(Nmin,Pmin,Nmax)

if max(Nmax,Pmax)>=N/2:
   print("max (Nmax(%d), Pmax (%d)) >= N/2 (%d), increase N" %(Nmax,Pmax,N/2))
   quit()

a = np.zeros(N*N,dtype=complex)
a.shape=(N,N)
beta=2*(H+1.)

for i in range(0,Nmax+1):
   for j in range(0,Pmax+1):

      if max(i,j)==0:
         continue
      k=sqrt(i**2+j**2)

      #dealing with anisotropy in lambdaMax:
      val = pow((1.*i)/Nmin,2)+pow((1.*j)/Pmin,2)
      if val<1.:
         continue
      #same witk lambdaMin:
      val = pow((1.*i)/Nmax,2)+pow((1.*j)/Pmax,2)
      if val>1.:
         continue
      fac=pow(k,-beta*0.5)

      randPhase = random.random()*np.pi*2.
      a[i,j]=fac*np.exp(randPhase*1.j)

      if (i>0):
         randPhase = random.random()*np.pi*2.
         a[N-i,j]=fac*np.exp(randPhase*1.j)
      if (j>0):
         randPhase = random.random()*np.pi*2.
         a[i,N-j]=fac*np.exp(randPhase*1.j)
      if min(i,j)>0:
         randPhase = random.random()*np.pi*2.
         a[N-i,N-j]=fac*np.exp(randPhase*1.j)

a=N*N*a
h= np.real(np.fft.ifft2(a))

#plt.imshow(h)
#plt.colorbar()
#plt.show()

dx=L/N
x=np.arange(-Ls/2,Ls/2+dx,dx)
y=np.arange(0.,Ld+dx,dx)

nx=np.shape(x)[0]
ny=np.shape(y)[0]

X, Y = np.meshgrid(x, y)

h=h[0:ny,0:nx]

plt.pcolormesh(X,Y,h)
plt.colorbar()
plt.show()

#compute rms roughness
hrms=np.std(h)
#print(sqrt(np.sum(np.square(h))/(N*N)))
print hrms

#scale to targeted Hrms
targetHrms=alpha*max(lambdaMaxS, lambdaMaxD)
print("targeted Hrms: %f" %targetHrms)
h=h*targetHrms/hrms
print(np.std(h))

#for the following study
a=a*targetHrms/hrms


plt.pcolormesh(X,Y,h)
plt.colorbar()
plt.show()

fout=open('roughFault_all.dat','w')
Xf=X.flatten()
Yf=Y.flatten()
hf=h.flatten()
for i in range(0,np.shape(hf)[0]):
    fout.write("%f %e %f\n" %(1e3*Xf[i],1e3*hf[i],-1e3*Yf[i]))
fout.close()
