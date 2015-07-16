#include "MeshSim.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "SimMeshingErrorCodes.h"
#include "SimDiscreteErrorCodes.h"
#include <iostream>

#include <stdio.h>
#include "SimDisplay.h"
using namespace std;

void messageHandler(int type, const char *msg);
void printFaceCoords(pGModel model);

int main(int argc, char *argv[])
{
    if (argc < 3) {
        cerr << "Too few arguments, input or output file missing!" << endl;
        return 1;
    }
    char* inputfile = argv[1];
    char* outputfile = argv[2];

    // You will want to place a try/catch around all Simmetrix calls
    // as errors are thrown.
    try {
        // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
        // pass in the location of a file containing your keys.  For a release
        // product, use Sim_registerKey()
        Sim_readLicenseFile("/work/klicpera/simmodeler/TUM");

        Sim_logOn("createDM.log");
        MS_init();
        SimDiscrete_start(0);

        Sim_setMessageHandler(messageHandler);
        pProgress progress = Progress_new();
        Progress_setDefaultCallback(progress);


        // ------------------------------ Import STL-File ------------------------------

        pMesh mesh = M_new(0,0);
        pDiscreteModel model = 0;
        if(M_importFromSTLFile(mesh, inputfile, progress)) { //check for error
            cerr<<"Error importing file"<<endl;
            M_release(mesh);
            return 1;
        }

        // check the input mesh for intersections
        // this call must occur before the discrete model is created
        if(MS_checkMeshIntersections(mesh,0,progress)) {
            cerr<<"There are intersections in the input mesh"<<endl;
            M_release(mesh);
            return 1;
        }

        // create the Discrete model
        model = DM_createFromMesh(mesh, 1, progress);
        if(!model) { //check for error
            cerr<<"Error creating Discrete model from mesh"<<endl;
            M_release(mesh);
            return 1;
        }

        // define the Discrete model
        DM_findEdgesByFaceNormalsDegrees(model, 70, progress);
        DM_eliminateDanglingEdges(model, progress);
        if(DM_completeTopology(model, progress)) { //check for error
            cerr<<"Error completing Discrete model topology"<<endl;
            M_release(mesh);
            GM_release(model);
            return 1;
        }

        // Since we told the Discrete model to use the input mesh, we release our
        // pointer to it.  It will be fully released when the Discrete model is released.
        M_release(mesh);

        // Print out information about the model
        cout<<"Number of model vertices: "<<GM_numVertices(model)<<endl;
        cout<<"Number of model edges: "<<GM_numEdges(model)<<endl;
        cout<<"Number of model faces: "<<GM_numFaces(model)<<endl;
        cout<<"Number of model regions: "<<GM_numRegions(model)<<endl;


        // ------------------------------ Set boundary conditions ------------------------------

        // Create a new manager. The manager is responsible for creating
        // new attributes, saving/retrieving to/from file, and to look up
        // cases, AttModels etc.
        pAManager attMngr = AMAN_new();

        // Create a case. A case serves as a grouping mechanism. It will (should)
        // contain all the attributes that make up a complete analysis of
        // some kind, although the system has no way of verifying this.
        pACase analysisCase = AMAN_newCase(attMngr,"analysis","",(pModel)model);

        // Now we create an attribute information node. This can be seen
        // as an "attribute generator object", as we will later on create an
        // actual attribute for each geometric face from this attribute
        // information node. We name the attribute T1, and give it the
        // information type "boundaryCondition".
        pAttInfoVoid iSurf = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
        AttNode_setImageClass((pANode)iSurf,"freeSurface");
        pAttInfoVoid iDynRup = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
        AttNode_setImageClass((pANode)iDynRup,"dynamicRupture");
        pAttInfoVoid iAbsorb = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
        AttNode_setImageClass((pANode)iAbsorb,"absorbing");

        // We need to add the Attribute Information Nodes to the case
        AttCase_addNode(analysisCase,(pANode)iSurf);
        AttCase_addNode(analysisCase,(pANode)iDynRup);
        AttCase_addNode(analysisCase,(pANode)iAbsorb);

        // We need to associate the Attribute Information Nodes with the
        // model. To do that we have to create a Model Association for the case
        pModelAssoc aSurf = AttCase_newModelAssoc(analysisCase,(pANode)iSurf);
        pModelAssoc aDynRup = AttCase_newModelAssoc(analysisCase,(pANode)iDynRup);
        pModelAssoc aAbsorb = AttCase_newModelAssoc(analysisCase,(pANode)iAbsorb);

        pGEntity face;
        for (int i = 1; i <= 7; i++) {
            // Get the face
            face = GM_entityByTag(model, 2, i);

            // Add the face to the model association. Note that we passed
            // the Attribute Information Node into the Model Association
            // at the time when the Model Association was created. That prepares
            // the creation of the AttributeVoid on the face as soon as the
            // association process is started
            if (i == 6) {
                AMA_addGEntity(aSurf,face);
            } else if (i == 7) {
                AMA_addGEntity(aDynRup,face);
            } else {
                AMA_addGEntity(aAbsorb,face);
            }
        }

        // printFaceCoords(model);


        // The association process now loops over the case tree and will create
        // actual attributes for each Attribute Information Node it finds in
        // its path. This finalizes the attribute setup process.
        AttCase_associate(analysisCase, progress);

        // ------------------------------ Set meshing parameters ------------------------------

		pACase meshCase = MS_newMeshCase(model);

        // Set global mesh size
        pModelItem modelDomain = GM_domain(model);
        // ( <meshing case>, <entity>, <1=absolute, 2=relative>, <size>, <size expression> )
        MS_setMeshSize(meshCase, modelDomain, 1, 5000, NULL);

        // Set mesh size on fault
        face = GM_entityByTag(model, 2, 7);
        MS_setMeshSize(meshCase, face, 1, 400, NULL);

        // Set gradation relative
        MS_setGlobalSizeGradationRate(meshCase, 0.2);

        // Set target skewness
        MS_setVolumeShapeMetric(meshCase, modelDomain, ShapeMetricType_Skewness, 0.75);

        // ------------------------------ Save file ------------------------------

        // Possible, because discreteModel inherits from NonManifoldModel
        GM_write(model, outputfile, 0, progress);
        GM_release(model);
        Progress_delete(progress);
        SimDiscrete_stop(0);
        MS_exit();
        Sim_logOff();
        Sim_unregisterAllKeys();

    } catch (pSimError err) {
        cerr<<"Simmetrix error caught:"<<endl;
        cerr<<"  Error code: "<<SimError_code(err)<<endl;
        cerr<<"  Error string: "<<SimError_toString(err)<<endl;
        SimError_delete(err);
        return 1;
    } catch (...) {
        cerr<<"Unhandled exception caught"<<endl;
        return 1;
    }
    return 0;
}

void printFaceCoords(pGModel model) {
    GFIter modelFaces;
    pGFace modelFace;
    int ID;
    pPList edgeList;  // Edges bounding a face
    pGEdge thisEdge;
    pSimPolygons poly; // tessellation of a face
    const int maxPolyPoints = 100;
    int polypoint[maxPolyPoints];   // ID of the points of a polygon
    double pntlocation[3];
    double pntnormal[3];
    modelFaces = GM_faceIter(model);
    cout << "\n********Face Information********\n\n";
    while(modelFace=GFIter_next(modelFaces)) { // get the next model face
        ID = GEN_tag(modelFace);
        // edgeList = GF_edges(modelFace);
        // printf("Model Face %d has the following %d edges: ", ID, PList_size(edgeList));
        // void* iter = 0;
        // while(thisEdge=(pGEdge)PList_next(edgeList, &iter)) {
        //       printf("%d ", GEN_tag(thisEdge));
        // }
        // printf("\n");
        // PList_delete(edgeList); // cleanup

        poly = GF_displayRep(modelFace);
        int npolys = SimPolygons_numPolys(poly);
        int npolypnts = SimPolygons_numPoints(poly);
        printf("There are %d polygons and %d points on model face %d, e.g.:\n", npolys, npolypnts, ID);
        int j;
        for (j=0; j<1; j++) { // loop over the polygons
            int myPoints = SimPolygons_polySize(poly, j);
            SimPolygons_poly(poly, j, polypoint);
            printf(" Polygon %d has the following points:", j);
            int k;
            for (k=0; k<myPoints; k++)
                printf(" %d", polypoint[k]);
            printf("\n");
            for (k=0; k<myPoints; k++) {
                int hasnorm = SimPolygons_pointData(poly, polypoint[k], pntlocation, pntnormal);
                printf(" Point %d located at: (%f,%f,%f)\n", polypoint[k], pntlocation[0], pntlocation[1], pntlocation[2]);
            }
        }
        // printf("The points on the face have the following locations and normals: \n");
        // for (j=0; j<npolypnts; j++) { // loop over the points
        //     int hasnorm = SimPolygons_pointData(poly, j, pntlocation, pntnormal);
        //     printf(" Point %d located at: (%f,%f,%f)", j, pntlocation[0], pntlocation[1], pntlocation[2]);
        //     if (hasnorm)
        //         printf("  normal vector: [%f,%f,%f]\n", pntnormal[0], pntnormal[1], pntnormal[2]);
        //     else
        //         printf("  NO NORMAL\n");
        // }
        SimPolygons_delete(poly); // cleanup
        printf("\n");
    }
    GFIter_delete(modelFaces); // cleanup
}

void messageHandler(int type, const char *msg)
{
    switch (type) {
        case Sim_InfoMsg:
            cout<<"Info: "<<msg<<endl;
            break;
        case Sim_DebugMsg:
            cout<<"Debug: "<<msg<<endl;
            break;
        case Sim_WarningMsg:
            cout<<"Warning: "<<msg<<endl;
            break;
        case Sim_ErrorMsg:
            cout<<"Error: "<<msg<<endl;
            break;
    }
    return;
}
