#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <igl/viewer/Viewer.h>

void displayTet (const char * VertsName, char const* EleName);
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier);
void Verts2Matrix(std::string NodeName);
