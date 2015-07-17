from generateRoughFault import genRoughFault
from generateBox import genBox
import os
import subprocess
from operator import sub
from math import floor, ceil

def mergePLY(in1, in2, out):
    fin1=open(in1,'r')
    fin2=open(in2,'r')
    fout=open(out,'w')
    fout.write("ply\n")
    fout.write("format ascii 1.0\n")
    fin1.readline()
    fin1.readline()
    fin2.readline()
    fin2.readline()

    # Header for vertices
    lenV1 = int(fin1.readline().rpartition(' ')[2])
    lenV2 = int(fin2.readline().rpartition(' ')[2])
    fout.write("element vertex %i\n" % (lenV1 + lenV2))
    fout.write("property float32 x\n")
    fout.write("property float32 y\n")
    fout.write("property float32 z\n")
    fin1.readline()
    fin1.readline()
    fin1.readline()
    fin2.readline()
    fin2.readline()
    fin2.readline()

    # Header for faces
    lenF1 = int(fin1.readline().rpartition(' ')[2])
    lenF2 = int(fin2.readline().rpartition(' ')[2])
    fout.write("element face %i\n" % (lenF1 + lenF2))
    fout.write("property list uint8 int32 vertex_index\n")
    fout.write("end_header\n")
    fin1.readline()
    fin1.readline()
    fin2.readline()
    fin2.readline()

    # List of vertices
    for _ in range(0, lenV1):
        fout.write(fin1.readline())
    for _ in range(0, lenV2):
        fout.write(fin2.readline())

    # List of faces
    for _ in range(0, lenF1):
        fout.write(fin1.readline())
    for _ in range(0, lenF2):
        oldStr = fin2.readline()
        newStr = oldStr.split(' ')[0]
        for i in range(0, int(oldStr.split(' ')[0])):
            newStr += ' ' + str(int(oldStr.split(' ')[i+1]) + lenV1)
        fout.write(newStr + "\n")

    fin1.close()
    fin2.close()
    fout.close()

def normCross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    norm = c[0] + c[1] + c[2]
    if norm == 0:
        raise Exception("Cross product is zero! Strange vectors: a={}, b={}".format(a, b))
    else:
        invNorm = 1 / norm
        for i in range(0, 3):
            c[i] = c[i] * invNorm
    return c

def vec_to_str(vec):
    res = str(vec[0])
    for i in range(1, len(vec)):
        res += ' ' + str(vec[i])
    return res

def PLYtoSTL(file_in, file_out):
    f_input = open(file_in,'r')

    f_input.readline()
    f_input.readline()

    # Read vertices header
    lenV = int(f_input.readline().rpartition(' ')[2])
    f_input.readline()
    f_input.readline()
    f_input.readline()

    # Read faces header
    lenF = int(f_input.readline().rpartition(' ')[2])
    f_input.readline()
    f_input.readline()

    # Read in vertices
    vertices = []
    for _ in range(0, lenV):
        vert = f_input.readline().split(' ')
        vertices.append([float(coord) for coord in vert])

    f_output = open(file_out,'w')
    f_output.write("solid rough\n")

    for i in range(0, lenF):
        line = f_input.readline()
        # Get vertice indices
        face = []
        for j in range(0, int(line.split(' ')[0])):
            face.append(int(line.split(' ')[j+1]))
        # Get vertice coordinates
        fVert = []
        for j in range(0, len(face)):
            fVert.append(vertices[face[j]])
        # Calculate face normal
        vec1 = map(sub, fVert[1], fVert[0])
        vec2 = map(sub, fVert[2], fVert[0])
        normal = normCross(vec1, vec2)
        # Write face
        f_output.write("facet normal " + vec_to_str(normal) + "\n")
        f_output.write("\touter loop\n")
        for j in range(0, len(face)):
            f_output.write("\t\tvertex " + vec_to_str(fVert[j]))
        f_output.write("\tendloop\n")
        f_output.write("endfacet\n")

    f_output.write("endsolid rough\n")

    f_input.close()
    f_output.close()

def genRecv(filename, ctrPts, num):
    # Method is asymmetrical! Rounded up intervals on left side (low x).
    fout = open(filename, 'w')
    n_rndUp = (len(ctrPts) - 1) % (num + 1)
    ind = 0
    for i in range(0, num):
        if i < n_rndUp:
            ind += int(ceil((len(ctrPts) - 1) * 1. / (num + 1)))
        else:
            ind += int(floor((len(ctrPts) - 1) * 1. / (num + 1)))
        fout.write(vec_to_str(ctrPts[ind]) + "\n")
    fout.close()

def genMesh(directory, rseed, length, depth, lambdaMin, alpha, Hurst):

    if directory == "":
        directory = "."
    if directory[-1] == "/":
        directory = directory[:-1]

    if not os.path.exists(directory):
        os.makedirs(directory)

    # Generate point clouds
    print("\n--------------- Generating point clouds ---------------")
    surfPts = []
    ctrPts = []
    print("Creating rough fault.")
    genRoughFault("{0}/roughFault.ply".format(directory), True, surfPts, ctrPts, rseed, length, depth, lambdaMin, alpha, Hurst)
    print("Creating box.")
    genBox("{0}/box.ply".format(directory), surfPts)
    print("Creating fault receiver list.")
    genRecv("{0}/Faultreceiverlist.dat".format(directory), ctrPts, 10)

    # Mesh point cloud, create CAD file
    print("\n----------------- Generating CAD file -----------------")
    print("Merging PLY-files.")
    mergePLY("{0}/roughFault.ply".format(directory), "{0}/box.ply".format(directory), "{0}/model.ply".format(directory))
    print("Converting PLY to STL.")
    PLYtoSTL("{0}/model.ply".format(directory), "{0}/model.stl".format(directory))

    # Mesh volume
    print("\n------------------- Generating mesh -------------------")
    subprocess.call("~/PUML/build/bin/pumgen -s simmodsuite -l ~/simmodeler/TUM --vtk {0}/mesh_dmp.vtk --stl ~/roughgen/gen_mesh/meshPar.par --prbfc {0}/model.stl $mpi_ranks {0}/rough_mesh.nc".format(directory), shell=True)

if __name__ == '__main__':
    genMesh("test_mesh", '0254887388', 40., 20., 1., pow(10.,-1.9), 0.8)
