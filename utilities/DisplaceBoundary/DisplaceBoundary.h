#include <iostream>
#include <string>
#include <vector>

typedef struct LocalDisplacement{int index; double dis;} LocalDisplacement;
void SetDisplacement(double * u, VolumetricMesh &Mesh, double delx);
void PerturbZ(double *u, double r);
void findBoundary(std::string Filename1, std::string Filename2, std::vector<LocalDisplacement> &Bound);
void HardBoundary(double delx, VolumetricMesh * Mesh, std::vector<LocalDisplacement> &Bound);
void setHardBoundary(const std::vector<LocalDisplacement> &hardBound, double * u);
int Constraints(std::set<int> &ConstrainIndex, std::string FileName, int fCoordinate);
