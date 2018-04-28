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
#include <cmath>

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
    
    Xmax = -1*volumetricMesh->getVertex(0)[0]; //TODO:modify the value based on the mesh
    cout << "Xmax = " << Xmax << endl;
   
    //Get output file name
    string outName;
    getline(ifs, outName);
    
    //get Boundary for bottom layer
    getline(ifs, str);
    int numConstrainedDOFs;
    set<int> constDOFs;
    numConstrainedDOFs = Constraints(constDOFs, str,3);
    
    //get Boundary front and back
    getline(ifs, str);
    getline(ifs, str1);
    numConstrainedDOFs = Constraints(constDOFs, str,2);
    numConstrainedDOFs = Constraints(constDOFs, str1,2);

    //get Boundary left and right
    getline(ifs, str);
    getline(ifs, str1);
    //numConstrainedDOFs = Constraints(constDOFs, str,3);
    //numConstrainedDOFs = Constraints(constDOFs, str,2);
    //numConstrainedDOFs = Constraints(constDOFs, str1,3);
    //numConstrainedDOFs = Constraints(constDOFs, str1,2);
    //
    int constrainedDOFs[numConstrainedDOFs]; 
    set<int>::iterator setIt;
    setIt = constDOFs.begin();
    for(int i=0; i<numConstrainedDOFs; i++){
        constrainedDOFs[i] = *setIt;
        //cout << constrainedDOFs[i]/3 << " " << constrainedDOFs[i]%3 << endl;
        setIt++;
    }
   
    int numDynamicConstrainedDOFs;
    set<int> constDynamicDOFs;
    numDynamicConstrainedDOFs = Constraints(constDynamicDOFs, str,1);
    numDynamicConstrainedDOFs = Constraints(constDynamicDOFs, str1,1);
    
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
    double dampingMassCoef; //= 0.01;  (primarily) high-frequency damping
    ifs >> dampingMassCoef;
    double DisTot;
    ifs >> DisTot;
    double delx;
    ifs >> delx;
    int nEqual; // number of equilibruim steps before applying next displacement
    ifs >> nEqual;
    //double tolerance; // number of equilibruim steps before applying next displacement
    //ifs >> tolerance;
    
    double timestep; // the timestep, in seconds
    ifs >> timestep;
    int numTimesteps, FrameRate;
    ifs >> numTimesteps;
    ifs >> FrameRate;
    cout << "Total displacement: " << DisTot <<"Time Step: " << timestep <<  " #steps: " << numTimesteps << " FrameRate: " << FrameRate << endl;

    double * dis = new double[r](); //() initializes all elements with 0
    PerturbZ(dis,r);
    int nRamp = DisTot/delx;

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
    //double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
    // initialize the integrator
    ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, numDynamicConstrainedDOFs, constrainedDynamicDOFs, dampingMassCoef, dampingStiffnessCoef);
    //ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, 0, NULL, dampingMassCoef, dampingStiffnessCoef);
    implicitBackwardEulerSparse->SetExternalForcesToZero();
    implicitBackwardEulerSparse->SetState(dis);
    int nDisInc = 0;
    
    for(int i=0; i<numTimesteps; i++)
    {
       // important: must always clear forces, as they remain in effect unless changed
       if (i%FrameRate==0){
          ofstream ofs;
          stringstream FileName; 
          FileName << outName <<  "_" << i/FrameRate << ".node"; 
          ofs.open(FileName.str()); 
          ofs << volumetricMesh->getNumVertices() << " 3 0 0 " << endl; 
          implicitBackwardEulerSparse->GetqState(dis);
          for(int j=0; j<volumetricMesh->getNumVertices(); j++){
              Vec3d & v = volumetricMesh->getVertex(j);
              ofs << j+1 << " " << v[0]+dis[3*j] << " " << v[1]+dis[3*j+1] << " " << v[2]+dis[3*j+2] << endl; 
          }
          ofs.close();
       }
       
        
       if(nDisInc < nRamp && i%nEqual ==0){
       //if(nDisInc < nRamp ){
       //   implicitBackwardEulerSparse->GetqState(dis);
       //   double * force = new double[r]; //() initializes all elements with 0
       //   double totalForce=0;
       //   deformableModel->ComputeForces(dis,force);
       //   for(int j=0 ; j<r; j++)
       //       totalForce += abs(force[j]);
       //
       //   cout << "Total Force : " << totalForce << endl;

       //   if(totalForce < tolerance)
       //   {
       //      SetDisplacement(dis,*volumetricMesh,delx);
       //      implicitBackwardEulerSparse->SetState(dis);
       //      nDisInc ++;
       //   }
       //   delete(force);
       
          implicitBackwardEulerSparse->GetqState(dis);
          SetDisplacement(dis,*volumetricMesh,delx);
          implicitBackwardEulerSparse->SetState(dis);
          nDisInc ++;
       }
       
       implicitBackwardEulerSparse->DoTimestep();
    }
    
    delete(dis);
    delete(neoHookeanModel); 
    delete(deformableModel);
    delete(forceModel);
    delete(implicitBackwardEulerSparse);
 
    return 0;
}

//Set X displacement based on the distant to x = 0,
//negative x with a positive displacement and vice versa
//void SetDisplacement(double * u, const  VolumetricMesh &Mesh, double ratio){
void SetDisplacement(double * u, VolumetricMesh &Mesh, double delx){
       Xmax = -1*(Mesh.getVertex(0)[0]+u[0]); //TODO:modify the value based on the mesh
       for(int i=0; i<Mesh.getNumVertices(); i++){
         const Vec3d &v = Mesh.getVertex(i);
         double newDx = -1 * (v[0]+u[3*i]) * delx/Xmax;
         u[3*i] += newDx;
       }
}

void PerturbZ(double *u, double r)
{
       srand(time(NULL));
       for(int i=0; i<r/3; i++)
         u[3*i+2] += ((double) rand()/ (RAND_MAX)) * 2e-5 - 1e-5;
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
