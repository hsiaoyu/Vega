#include "volumetricMeshLoader.h"
#include "corotationalLinearFEM.h"
#include "objMesh.h"
#include "corotationalLinearFEMForceModel.h"
#include "generateSurfaceMesh.h"
#include "generateMassMatrix.h"
#include "implicitBackwardEulerSparse.h"
#include "sceneObjectDeformable.h"
#include <fstream>
#include <sstream>
#include <string>
#include "neoHookeanIsotropicMaterial.h"
#include "isotropicHyperelasticFEMForceModel.h"
#include "isotropicHyperelasticFEM.h"

using namespace std;

int ZConstraints(vector<int> &ConstrainIndex, string FileName);

int main(int argc, char * argv[]){
    if (argc <1){
        printf ("Need to enter config file:");
        return 1;
    }
     
    ifstream ifs;
    ifs.open(argv[1]);
    string str;
    getline(ifs, str);
    VolumetricMesh * volumetricMesh = VolumetricMeshLoader::load(str.c_str());
    if (volumetricMesh == NULL)
        printf("Error: failed to load mesh.\n");
    else
        printf("Success. Number of vertices: %d . Number of elements: %d .\n",
        volumetricMesh->getNumVertices(), volumetricMesh->getNumElements());
    TetMesh * tetMesh;
    if (volumetricMesh->getElementType() == VolumetricMesh::TET)
        tetMesh = (TetMesh*) volumetricMesh; // such down-casting is safe in Vega
    else
    {
       printf("Error: not a tet mesh.\n");
       exit(1);
    }
    
    int r = 3 * tetMesh->getNumVertices(); // total number of DOFs
    
    getline(ifs, str);
    int numConstrainedDOFs, *constrainedDOFs;
    vector<int> constDOFs;
    numConstrainedDOFs = ZConstraints(constDOFs, str);
    constrainedDOFs = &constDOFs[0];
    
    getline(ifs, str);
    ifstream Boundary[2];
    Boundary[0].open(str.c_str());
    getline(ifs, str);
    Boundary[1].open(str.c_str());
    //Get output file name
    getline(ifs, str);
    double force, timestep; // the timestep, in seconds
    double dampingStiffnessCoef; // (primarily) high-frequency damping
    ifs >> dampingStiffnessCoef;
    ifs >> force;
    ifs >> timestep;
    int numTimesteps, FrameRate;
    ifs >> numTimesteps;
    ifs >> FrameRate;

    double * f = new double[r]();
    for(int i=0; i<2; i++){
        int vb;
        while(Boundary[i]>>vb){
            f[3*(vb-1)] = (i==0) ? force : -force ;
        }
    }

//    int enableConstraintResistance = 0;
//    double ConstraintResistance = 0.0;
//    NeoHookeanIsotropicMaterial * neoHookeanModel = new NeoHookeanIsotropicMaterial(tetMesh, enableConstraintResistance, ConstraintResistance);
//    IsotropicHyperelasticFEM * deformableModel = new IsotropicHyperelasticFEM(tetMesh, neoHookeanModel);
//    ForceModel * forceModel = new IsotropicHyperelasticFEMForceModel(deformableModel);

    CorotationalLinearFEM * deformableModel = new CorotationalLinearFEM(tetMesh);
    //create the class to connect the deformable model to the integrator
    ForceModel * forceModel = new CorotationalLinearFEMForceModel(deformableModel);

    SparseMatrix * massMatrix;
    // create consistent (non-lumped) mass matrix
    GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix, true);
    // (tangential) Rayleigh damping
    double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
    // initialize the integrator
    ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef);

    // alocate buffer for external forces


    for(int i=0; i<numTimesteps; i++)
    {
       // important: must always clear forces, as they remain in effect unless changed
       //implicitBackwardEulerSparse->SetExternalForcesToZero();
       if (i == 0) // set some force at the first timestep
       {
         implicitBackwardEulerSparse->SetExternalForces(f);
       }
       implicitBackwardEulerSparse->DoTimestep();
       if (i%FrameRate==0){
          // alocate buffer to read the resulting displacements
          double * u = (double*) malloc (sizeof(double) * r);
          implicitBackwardEulerSparse->GetqState(u);
         
          ofstream ofs; 
          stringstream FileName;
	  FileName << str << "_" << i/FrameRate << ".node"; 
          ofs.open(FileName.str()); 
          ofs << volumetricMesh->getNumVertices() << " 3 0 0 " << endl; 
    //      volumetricMesh->applyDeformation(u);
          for(int j=0; j<volumetricMesh->getNumVertices(); j++){
              Vec3d & v = volumetricMesh->getVertex(j);
              ofs << j+1 << " " << v[0]+u[3*j] << " " << v[1]+u[3*j+1] << " " << v[2]+u[3*j+2] << endl; 
          }
          ofs.close();
          free(u);
       }
    } 
 /*
    SceneObjectDeformable * deformableRenderingMesh = NULL;
    deformableRenderingMesh->SetVertexDeformations(u); 
    deformableRenderingMesh->Render();
  */
    return 0;
}


int ZConstraints(vector<int> &ConstrainIndex, string FileName){
    ifstream Boundary;
    Boundary.open(FileName); // Input the indices at the bottom surface
    int i, sizeC;
    sizeC = 0;
    while(Boundary >> i){
         ConstrainIndex.push_back(3*(i-1)+2); // Constrain only the z component
         sizeC ++;
    }
    return sizeC;
}
