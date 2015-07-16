def genBox(filename, fault_surface):

    # ----- Header -----
    fout=open(filename,'w')
    fout.write("ply\n")
    fout.write("format ascii 1.0\n")
    fout.write("element vertex %i\n" % (10 + len(fault_surface)))
    fout.write("property float32 x\n")
    fout.write("property float32 y\n")
    fout.write("property float32 z\n")
    fout.write("element face %i\n" % (10 + 8 + 2 * (len(fault_surface) - 1)))
    # fout.write("element face %i\n" % (10))
    fout.write("property list uint8 int32 vertex_index\n")
    fout.write("end_header\n")

    # ----- Vertices -----
    # Box
    fout.write("70e3 70e3 0\n")
    fout.write("70e3 -70e3 0\n")
    fout.write("-70e3 -70e3 0\n")
    fout.write("-70e3 70e3 0\n")
    fout.write("70e3 70e3 -75e3\n")
    fout.write("70e3 -70e3 -75e3\n")
    fout.write("-70e3 -70e3 -75e3\n")
    fout.write("-70e3 70e3 -75e3\n")
    # Points for fault intersection
    fout.write( "0 68e3 0\n")
    fout.write( "0 -68e3 0\n")

    # Fault-surface intersection
    for i in range(0, len(fault_surface)):
        fout.write("%f %f %f\n" %(1e3*fault_surface[i][0], 1e3*fault_surface[i][1], 0.0))

    # ----- Faces -----
    # Subsurface box
    fout.write("3 4 5 6\n")
    fout.write("3 4 6 7\n")
    fout.write("3 2 3 7\n")
    fout.write("3 2 7 6\n")
    fout.write("3 1 2 6\n")
    fout.write("3 1 6 5\n")
    fout.write("3 0 1 5\n")
    fout.write("3 0 5 4\n")
    fout.write("3 0 3 7\n")
    fout.write("3 0 7 4\n")

    # Surface plane
    fout.write("3 0 1 %i\n" % (10 + len(fault_surface) - 1))
    fout.write("3 0 %i 8\n" % (10 + len(fault_surface) - 1))
    fout.write("3 0 8 3\n")
    fout.write("3 3 8 10\n")
    fout.write("3 3 10 2\n")
    fout.write("3 1 9 %i\n" % (10 + len(fault_surface) - 1))
    fout.write("3 1 2 9\n")
    fout.write("3 2 10 9\n")
    for i in range(0, len(fault_surface) - 1):
        fout.write("3 9 %i %i\n" % (10 + i, 11 + i))
        fout.write("3 8 %i %i\n" % (11 + i, 10 + i))

    fout.close()
