import numpy as np
# pylint: disable=E1101
import matplotlib.pyplot as plt
import random

def genRoughFault(filename, mesh, surfacePts, rseed, length, depth, lambdaMin, alpha, Hurst):
    # Input parameters:
    # filename: Filename for saving .ply file
    # mesh: Boolian value for defining whether to create a mesh
    # surfacePts: Empty array for saving the surface points (for creating an intersection later)
    # rseed: Seed for random number generator
    # length: Length of the fault (in km)
    # depth: Depth of the fault (in km)
    # lambdaMin: Minimal wavelength
    # alpha: Amplitude to wavelength ratio
    # Hurst: Hurst exponent

    random.seed(rseed)

    L = max(length,depth)
    lambdaMax = L

    Ncutoff = int(L/lambdaMin)

    # Factor 2 due to neg & pos. frequencies, refinement, + 1 due to frequency (0, 0)
    refine = 4
    N = refine * Ncutoff * 2 + 1

    freqs = np.zeros(N * N, dtype=complex)
    freqs.shape = (N, N)
    beta = 2 * (Hurst + 1.)

    for kx in range(0, Ncutoff+1):
        for kz in range(0, Ncutoff+1):

            # No frequency of 0 (constant value)
            if max(kx,kz)==0:
                continue

            # Calculate desired amplitude
            k = np.sqrt(kx**2 + kz**2)
            fac = pow(k, -beta * 0.5)

            # Calculate random phases and scale with desired amplitude
            randPhase = random.random()*np.pi*2.
            freqs[kx,kz]=fac*np.exp(randPhase*1.j)

            # Copy conjugates to other quadrants, to get real values in spacial domain
            if kx != 0:
                freqs[N-kx, kz] = np.conjugate(freqs[kx,kz])
            if kz != 0:
                freqs[kx, N-kz] = np.conjugate(freqs[kx,kz])
            if min(kx, kz) != 0:
                freqs[N-kx, N-kz] = freqs[kx,kz]

    Y = np.real(np.fft.ifft2(freqs))

    dx=L/(N - 1) # -1 because of point (0,0)
    x=np.arange(-length/2,length/2+dx,dx)
    z=np.arange(0.,depth+dx,dx)

    X, Z = np.meshgrid(x, z)

    # Crop square to rectangle
    Y=Y[0:len(z),0:len(x)]

    # compute hrms roughness
    hrms = np.std(Y)

    #scale to targeted Hrms
    target_hrms = alpha * lambdaMax
    # print("hrms = {}, targeted hrms = {}".format(hrms, target_hrms))
    Y = Y * target_hrms / hrms
    # print("Corrected hrms: {}".format(np.std(Y)))

    #for the following study
    # freqs=freqs*target_hrms/hrms

    # Show color map
    # plt.pcolormesh(X,Z,Y)
    # plt.colorbar()
    # plt.show()

    # Save surface points for intersection with box
    for i in range(0,len(x)):
        surfacePts.append([x[i],Y[0,i]])

    Xf=X.flatten()
    Yf=Y.flatten()
    Zf=Z.flatten()

    # Write ply-file
    fout=open(filename,'w')
    # Header
    fout.write("ply\n")
    fout.write("format ascii 1.0\n")
    fout.write("element vertex %i\n" % len(Yf))
    fout.write("property float32 x\n")
    fout.write("property float32 y\n")
    fout.write("property float32 z\n")
    # Header for faces
    if mesh:
        fout.write("element face %i\n" % (2 * (len(z) - 1) * (len(x) - 1)))
        fout.write("property list uint8 int32 vertex_index\n")
    fout.write("end_header\n")

    # Vertices
    for i in range(0, len(Yf)):
        if Zf[i] == 0: # prevent -0.000000
            fout.write("%f %f %f\n" %(1e3*Xf[i], 1e3*Yf[i], 0.))
        else:
            fout.write("%f %f %f\n" %(1e3*Xf[i], 1e3*Yf[i], -1e3*Zf[i]))

    # Faces
    if mesh:
        for i in range(0, len(z) - 1):
            for j in range(0, len(x) - 1):
                n = j + i * len(x)
                fout.write("3 %i %i %i\n" %(n, n+len(x), n+1))
                fout.write("3 %i %i %i\n" %(n+1, n+len(x), n+1+len(x)))
    fout.close()

if __name__ == '__main__':
    surfPts = []
    genRoughFault("roughFault.ply", True, surfPts, '0254887388', 40., 20., 1., pow(10.,-1.9), 0.8)
