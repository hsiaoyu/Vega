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
#include <time.h>
#include "DisplaceBoundary.h"
#include "neoHookeanIsotropicMaterial.h"
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
    
    Xmax = -1*volumetricMesh->getVertex(0)[0]; //TODO:modify the value based on the mesh
    cout << "Xmax = " << Xmax << endl;
    int r = 3 * tetMesh->getNumVertices(); // total number of DOFs
   
    //Get output file name
    string outName;
    getline(ifs, outName);
    
    //get Boundary for bottom layer
    getline(ifs, str);
    int numConstrainedDOFs;
    set<int> constDOFs;
    numConstrainedDOFs = Constraints(constDOFs, str,3);
    
    //get Boundary left and right
    getline(ifs, str);
    getline(ifs, str1);
    numConstrainedDOFs = Constraints(constDOFs, str,3);
    numConstrainedDOFs = Constraints(constDOFs, str,2);
    numConstrainedDOFs = Constraints(constDOFs, str1,3);
    numConstrainedDOFs = Constraints(constDOFs, str1,2);

    int numDynamicConstrainedDOFs;
    set<int> constDynamicDOFs;
    numDynamicConstrainedDOFs = Constraints(constDynamicDOFs, str,1);
    numDynamicConstrainedDOFs = Constraints(constDynamicDOFs, str1,1);
    
    int constrainedDOFs[numConstrainedDOFs]; 
    set<int>::iterator setIt = constDOFs.begin();
    for(int i=0; i<numConstrainedDOFs; i++){
        constrainedDOFs[i] = *setIt;
        //cout << constrainedDOFs[i]/3 << " " << constrainedDOFs[i]%3 << endl;
        setIt++;
    }
   
    int constrainedDynamicDOFs[numDynamicConstrainedDOFs];
    setIt = constDynamicDOFs.begin();
    for(int i=0; i<numDynamicConstrainedDOFs; i++){
        constrainedDynamicDOFs[i] = *setIt;
    //    cout << constrainedDynamicDOFs[i]/3 << " " << constrainedDynamicDOFs[i]%3 << endl;
        setIt++;
    }

    
    vector<LocalDisplacement> disboundary;
    findBoundary(str, str1, disboundary);

    double dampingStiffnessCoef; //= 0.01;  (primarily) high-frequency damping
    ifs >> dampingStiffnessCoef;
    double DisTot;
    ifs >> DisTot;
    double rampUpTime;
    ifs >> rampUpTime;
    
    double timestep; // the timestep, in seconds
    ifs >> timestep;
    int numTimesteps, FrameRate;
    ifs >> numTimesteps;
    ifs >> FrameRate;
    cout << "Total displacement: " << DisTot <<"Time Step: " << timestep <<  " #steps: " << numTimesteps << " FrameRate: " << FrameRate << endl;

    double * dis = new double[r](); //() initializes all elements with 0
    double delx = (double) (DisTot * timestep / rampUpTime);
    int nRamp = rampUpTime / timestep;
    cout << "delx = " << delx << "  nRamp = " << nRamp << endl; 
    //HardBoundary(delx, volumetricMesh, disboundary);
    //setHardBoundary(disboundary, dis);
    
    SetDisplacement(dis, *volumetricMesh, delx);
    
    //HardBoundary(DisTot, volumetricMesh, disboundary);
    //setHardBoundary(disboundary, dis);
    

    //CorotationalLinearFEM * deformableModel = new CorotationalLinearFEM(tetMesh);
    ////create the class to connect the deformable model to the integrator
    //ForceModel * forceModel = new CorotationalLinearFEMForceModel(deformableModel);

    int enableConstraintResistance = 0;
    double ConstraintResistance = 0.0;
    NeoHookeanIsotropicMaterial * neoHookeanModel = new NeoHookeanIsotropicMaterial(tetMesh, enableConstraintResistance, ConstraintResistance);
    IsotropicHyperelasticFEM * deformableModel = new IsotropicHyperelasticFEM(tetMesh, neoHookeanModel);
    ForceModel * forceModel = new IsotropicHyperelasticFEMForceModel(deformableModel);
    
    SparseMatrix * massMatrix;
    // create consistent (non-lumped) mass matrix
    GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix, true);


    // (tangential) Rayleigh damping
    double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
    // initialize the integrator
    ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, numDynamicConstrainedDOFs, constrainedDynamicDOFs, dampingMassCoef, dampingStiffnessCoef);
    //ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, 0, NULL, dampingMassCoef, dampingStiffnessCoef);
    cout << "Initiating Implicit Solver" << endl;
    // // Set velocity
    // double * vel = new double[r](); //() initializes all elements with 0
    // vector<LocalDisplacement> velboundary;
    // velboundary = disboundary;
    // HardBoundary(DisTot, volumetricMesh, velboundary);
    // setHardBoundary(velboundary, vel);
    // implicitBackwardEulerSparse->SetState(dis, vel);
    implicitBackwardEulerSparse->SetState(dis);
    implicitBackwardEulerSparse->SetExternalForcesToZero();
     
    for(int i=0; i<numTimesteps; i++)
    {
       // important: must always clear forces, as they remain in effect unless changed
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
         // cout << i << "set state" << endl;
         // cout << 0 << " " << vel[0] << " "  << vel[1] << " " << vel[2] <<  endl;
       }
      // cout << i << "  new state" << endl;
      // for(int j=80; j<100; j++)
      //    cout << j+1 << " " << dis[3*j] << endl;

       if(i<nRamp){
          double newDelx;
          newDelx = (i+2) *delx ;
          SetDisplacement(dis,volumetricMesh,delx);
          //HardBoundary(newDelx, volumetricMesh, disboundary);
       }

       
       //setHardBoundary(disboundary, dis);
       //implicitBackwardEulerSparse->SetState(dis);
       //if(i<nRamp){
       //   setHardBoundary(velboundary, vel);
       //   implicitBackwardEulerSparse->SetState(dis,vel);
       //}
       //if(i%FrameRate==0){ 
       //  cout << i << "set state" << endl;
       //  for(int j=0; j<r/3; j++)
       //      cout << j+1 << " " << vel[3*j] << " "  << vel[3*j+1] << " " << vel[3*j+2] <<  endl;
       //}
       implicitBackwardEulerSparse->DoTimestep();
    }

    return 0;
}

//Set X displacement based on the distant to x = 0,
//negative x with a positive displacement and vice versa
//void SetDisplacement(double * u, const  VolumetricMesh &Mesh, double ratio){
void SetDisplacement(double * u, VolumetricMesh &Mesh, double delx){
       srand(time(NULL));
       for(int i=0; i<Mesh.getNumVertices(); i++){
         const Vec3d &v = Mesh.getVertex(i);
	 u[3*i] = -1 * v[0] * delx/Xmax;
         u[3*i+2] += ((double) rand()/ (RAND_MAX)) * 2e-5 - 1e-5;
       }
     
}

void findBoundary(string Filename1, string Filename2, vector<LocalDisplacement> &Bound){
    ifstream Boundary[2];
    Boundary[0].open(Filename1);
    Boundary[1].open(Filename2);

    for(int i=0; i<2; i++){
        int vb;
        while(Boundary[i]>>vb){
              LocalDisplacement temp;
              temp.index = vb-1;
              temp.dis = 0;
              Bound.push_back(temp);
        }
    }

}

void HardBoundary(double delx, VolumetricMesh * Mesh, vector<LocalDisplacement> &Bound){
    for(unsigned int i=0; i<Bound.size(); i++){
           const Vec3d &v = Mesh->getVertex(Bound[i].index);
           Bound[i].dis = (v[0] > 0) ? -delx : delx; 
    }
}

void setHardBoundary(const vector<LocalDisplacement> &hardBound, double * u){
    //cout << "SetHardBoundary" << endl; 
    for(unsigned int i=0; i<hardBound.size(); i++){
         int index = hardBound[i].index;
         double dis = hardBound[i].dis;
         u[3*index] = dis;
    //     cout << index << " " << u[3*index] << endl;
     }
}

int Constraints(set<int> &ConstrainIndex, string FileName, int fCoordinate){
    ifstream Boundary;
    Boundary.open(FileName); // Input the indices at the bottom surface
    int i;
    while(Boundary >> i){
         switch (fCoordinate){
                case 1 : ConstrainIndex.insert(3*(i-1)); 
                         break; // Constrain only the x component
                case 2 : ConstrainIndex.insert(3*(i-1)+1); 
                         break; // Constrain only the y component
                case 3 : ConstrainIndex.insert(3*(i-1)+2); 
                         break; // Constrain only the z component
         }
    }
    return ConstrainIndex.size();
}
