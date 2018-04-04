#include <iostream>
#include <fstream>
#include "displayTets.h"

using namespace std;
using namespace Eigen;
int currentFrame = 0;
MatrixXd Verts;
MatrixXi Elements;
const char * VertsName;

void displayTet (const char * vertsName, char const* EleName){

 /*
    ifstream NodeFile, EleFile;
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
        //NodeFile >> buffer; 
    }
   
 
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
  */

    VertsName = vertsName;
    stringstream VertsFile;
    VertsFile << VertsName << "_" << 0 << ".node"; 
    Verts2Matrix(VertsFile.str());
   // cout << Verts << endl;
    
    ifstream EleFile;
    EleFile.open(EleName);
    int nEle, buffer;
    EleFile >> nEle;
    //cout << nEle << endl;
    //MatrixXi Elements(nEle,4);
    Elements.resize(nEle,3);
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
  
 
    //cout << Elements << endl;
     
    
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(Verts, Elements);
    viewer.data.set_face_based(true);
   // viewer.core.is_animating = true;
   // viewer.core.animation_max_fps = 30.;
    viewer.callback_key_down = &key_down;
    viewer.launch();
    
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier){
     
/*     switch(key){
           case '1' : currentFrame += 1;
           case '2' : currentFrame += 5;
           case '3' : currentFrame -= 1;
           case '4' : currentFrame -= 5;
     }
*/    
     if (key == '1') 
        currentFrame ++;
     else if (key == '2')
        currentFrame += 5;
     else if (key == '3')
        currentFrame -= 1;
     else if (key == '4')
        currentFrame -= 5;
     //cout << "key " << key << endl;
     stringstream vertsFile;
     vertsFile << VertsName << "_" << currentFrame << ".node";
     //vertsFile << "_" << currentFrame << ".node";
     //cout << vertsFile.str() << endl; 
     Verts2Matrix(vertsFile.str());
     viewer.data.clear();
     viewer.data.set_mesh(Verts,Elements);
     viewer.core.align_camera_center(Verts, Elements);
     cout << "current Frame : " << currentFrame << endl;
     
}


void Verts2Matrix(string NodeName){
    ifstream NodeFile;
    NodeFile.open(NodeName);
    int nVerts;
    NodeFile >> nVerts;
    //cout << nVerts << endl;
    Verts.resize(nVerts,3);
    int buffer;
    NodeFile >> buffer; 
    NodeFile >> buffer; 
    NodeFile >> buffer; 
    for(int i=0; i<nVerts; i++){
        NodeFile >> buffer; 
        NodeFile >> Verts(i,0); 
        NodeFile >> Verts(i,1); 
        NodeFile >> Verts(i,2); 
    }
    NodeFile.close();
}
