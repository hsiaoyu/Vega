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
#include <vector>
#include <set>
#include "neoHookeanIsotropicMaterial.h"
#include "StVKIsotropicMaterial.h"
#include "isotropicHyperelasticFEMForceModel.h"
#include "isotropicHyperelasticFEM.h"

using namespace std;

double Xmax;

int main(int argc, char * argv[]){
    if (argc <1){
        printf ("Need to enter config file:");
        return 1;
    }
     
    ifstream ifs;
    ifs.open(argv[1]);
    string str, str1;
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
    
   
    string outName;
    getline(ifs, outName);
    double dampingMassCoef; //= 0.01;  (primarily) high-frequency damping
    ifs >> dampingMassCoef;
    double dampingStiffnessCoef; //= 0.01;  (primarily) high-frequency damping
    ifs >> dampingStiffnessCoef;
    double timestep;
    ifs >> timestep;
    int numTimesteps, FrameRate;
    ifs >> numTimesteps;
    ifs >> FrameRate;
        
    int enableConstraintResistance = 0;
    double ConstraintResistance = 0.0;
    
    StVKIsotropicMaterial *stvkModel = new StVKIsotropicMaterial(tetMesh, enableConstraintResistance, ConstraintResistance);
    IsotropicHyperelasticFEM * deformableModel = new IsotropicHyperelasticFEM(tetMesh, stvkModel);
    ForceModel * forceModel = new IsotropicHyperelasticFEMForceModel(deformableModel);
 
    SparseMatrix * massMatrix;
    GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix, true);

    //double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
    //ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, numDynamicConstrainedDOFs, constrainedDynamicDOFs, dampingMassCoef, dampingStiffnessCoef);
    ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, 0, NULL, 0, NULL, dampingMassCoef, dampingStiffnessCoef);
    
     
    for(int i=0; i<numTimesteps; i++)
    {
      //implicitBackwardEulerSparse->SetState(dis);
       double * dis = new double[r]();
       implicitBackwardEulerSparse->SetExternalForcesToZero();
       implicitBackwardEulerSparse->GetqState(dis);
       if (i%FrameRate==0){
          ofstream ofs;
          stringstream FileName; 
          FileName << outName <<  "_" << i/FrameRate << ".node"; 
          ofs.open(FileName.str()); 
          ofs << volumetricMesh->getNumVertices() << " 3 0 0 " << endl; 
          for(int j=0; j<volumetricMesh->getNumVertices(); j++){
              Vec3d & v = volumetricMesh->getVertex(j);
              ofs << j+1 << " " << v[0]+dis[3*j] << " " << v[1]+dis[3*j+1] << " " << v[2]+dis[3*j+2] << endl; 
          }
          ofs.close();
       }

       cout << "timestep : " << i <<" Energy: " << deformableModel->ComputeEnergy(dis) << endl;
       implicitBackwardEulerSparse->DoTimestep();
       delete(dis);
    }
    delete(stvkModel);
    delete(deformableModel);
    delete(forceModel);
    delete(implicitBackwardEulerSparse);
    return 0;
}
