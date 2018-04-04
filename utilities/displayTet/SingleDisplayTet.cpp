#include <iostream>
#include <fstream>
#include <igl/viewer/Viewer.h>

using namespace std;
using namespace Eigen;

void displayTet (char* VertsName, char* EleName);

void displayTet (char* VertsName, char* EleName){
//int main(int argc, char* argv[]){
/*    if (argc <1){
        printf ("Need to enter config file:");
        return 1;
    }
     
    ifstream ifs;
    ifs.open(argv[1]);
    string str;
    getline(ifs, str);
    ifstream NodeFile, EleFile;
*/
    //NodeFile.open(str.c_str());
    NodeFile.open(VertsName);
    int nVerts;
    NodeFile >> nVerts;
    cout << nVerts << endl;
    int buffer;
    NodeFile >> buffer; 
    NodeFile >> buffer; 
    NodeFile >> buffer; 
    MatrixXd Verts(nVerts,3);
    for(int i=0; i<nVerts; i++){
        NodeFile >> buffer; 
        NodeFile >> Verts(i,0); 
        NodeFile >> Verts(i,1); 
        NodeFile >> Verts(i,2); 
    }
    cout << Verts << endl;
    

    //getline(ifs, str);
    //EleFile.open(str.c_str());
    EleFile.open(EleName);
    int nEle;
    EleFile >> nEle;
    cout << nEle << endl;
    MatrixXi Elements(nEle,4);
    EleFile >> buffer; 
    EleFile >> buffer; 
    for(int i=0; i< nEle; i++){
        EleFile >> buffer; 
        EleFile >> Elements(i,0); 
        EleFile >> Elements(i,1); 
        EleFile >> Elements(i,2); 
        EleFile >> Elements(i,3); 
        EleFile >> buffer;
    }
    Elements -= MatrixXi::Constant(nEle, 4, 1);
  

/*
    getline(ifs, str);
    EleFile.open(str.c_str());
    int nEle;
    EleFile >> nEle;
    cout << nEle << endl;
    //MatrixXi Elements(nEle,4);
    MatrixXi Elements(nEle,3);
    //EleFile >> buffer; 
    EleFile >> buffer; 
    for(int i=0; i< nEle; i++){
        EleFile >> buffer; 
        EleFile >> Elements(i,0); 
        EleFile >> Elements(i,1); 
        EleFile >> Elements(i,2); 
      //  EleFile >> Elements(i,3); 
        EleFile >> buffer; 
    }
    Elements -= MatrixXi::Constant(nEle, 3, 1);
  */
 
    cout << Elements << endl;
     
    
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(Verts, Elements);
    viewer.data.set_face_based(true);
   // viewer.core.is_animating = true;
   // viewer.core.animation_max_fps = 30.;
   // viewer.callback_pre_draw = &pre_draw;
    viewer.launch();
    
    //return 0;
}
