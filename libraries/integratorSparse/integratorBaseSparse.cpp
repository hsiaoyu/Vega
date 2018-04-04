/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.1                               *
 *                                                                       *
 * "integrator" library , Copyright (C) 2007 CMU, 2009 MIT, 2016 USC     *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integratorBaseSparse.h"
#include <iostream> //TODO: remove this line
#include <set>

IntegratorBaseSparse::IntegratorBaseSparse(int r, double timestep, SparseMatrix * massMatrix_, ForceModel * forceModel_, int numConstrainedDOFs_, int * constrainedDOFs_, int numDynamicConstrainedDOFs_, int * DynamicConstrainedDOFs_, double dampingMassCoef, double dampingStiffnessCoef): IntegratorBase(r, timestep, dampingMassCoef, dampingStiffnessCoef), massMatrix(massMatrix_), forceModel(forceModel_), numConstrainedDOFs(numConstrainedDOFs_), numDynamicConstrainedDOFs(numDynamicConstrainedDOFs_)
{
  systemSolveTime = 0.0;
  forceAssemblyTime = 0.0;

  constrainedDOFs = (int*) malloc (sizeof(int) * numConstrainedDOFs);
  memcpy(constrainedDOFs, constrainedDOFs_, sizeof(int) * numConstrainedDOFs);
  
  DynamicConstrainedDOFs = (int*) malloc (sizeof(int) * numDynamicConstrainedDOFs);
  memcpy(DynamicConstrainedDOFs, DynamicConstrainedDOFs_, sizeof(int) * numDynamicConstrainedDOFs);
  newDynamicCTIndices = (int*) malloc (sizeof(int) * numDynamicConstrainedDOFs);
  memcpy(newDynamicCTIndices, DynamicConstrainedDOFs_, sizeof(int) * numDynamicConstrainedDOFs);

  numTotConstrainedDOFs = numConstrainedDOFs + numDynamicConstrainedDOFs; 
  TotConstrainedDOFs = (int*) malloc (sizeof(int) * numTotConstrainedDOFs);
  
  ownDampingMatrix = 1;
  SparseMatrixOutline outline(r);
  dampingMatrix = new SparseMatrix(&outline);
  
  tangentStiffnessMatrixOffset = NULL;
  getNewIndices();
  AllConstraints();

}

IntegratorBaseSparse::~IntegratorBaseSparse()
{
  free(constrainedDOFs);
  if (ownDampingMatrix)
    delete(dampingMatrix);
  delete(tangentStiffnessMatrixOffset);
}

void IntegratorBaseSparse::AllConstraints()
{
    std::set<int> temp;
    for(int i=0; i<numConstrainedDOFs ; i++)
        temp.insert(constrainedDOFs[i]);
    for(int i=0; i<numDynamicConstrainedDOFs ; i++)
        temp.insert(DynamicConstrainedDOFs[i]);
    std::set<int>::iterator setIt = temp.begin();
    for(int i=0; i<numTotConstrainedDOFs; i++){
        TotConstrainedDOFs[i] = *setIt;
        setIt ++;
    }
}
//Compute the new indices for dynamic constraints after removing the constant constraints
void IntegratorBaseSparse::getNewIndices()
{
  if(numDynamicConstrainedDOFs<1)
    newDynamicCTIndices = NULL;
  else
  {
    for(int i=0; i<numConstrainedDOFs; i++){
        for(int j=0; j<numDynamicConstrainedDOFs; j++){
            if(constrainedDOFs[i] < DynamicConstrainedDOFs[j])
               newDynamicCTIndices[j]--;
            else if(constrainedDOFs[i] == DynamicConstrainedDOFs[j]){
               printf("Zero constraint and dynamic constraint on the same DOF!!");
               break;
            }
        }
    }
  }
  //for(int i=0; i<numDynamicConstrainedDOFs; i++){
  //    std::cout << "oldDOF: " << DynamicConstrainedDOFs[i] << std::endl;     
  //    std::cout << "newDOF: " << newDynamicCTIndices[i] << std::endl;     
  //}
  
}


void IntegratorBaseSparse::SetDampingMatrix(SparseMatrix * dampingMatrix_)
{
  if (ownDampingMatrix)
    delete(dampingMatrix);

  dampingMatrix = dampingMatrix_;
  ownDampingMatrix = 0;
}

double IntegratorBaseSparse::GetKineticEnergy()
{
  return 0.5 * massMatrix->QuadraticForm(qvel);
}

double IntegratorBaseSparse::GetTotalMass()
{
  return massMatrix->SumEntries();
}

void IntegratorBaseSparse::SetTangentStiffnessMatrixOffset(SparseMatrix * tangentStiffnessMatrixOffset_, int reuseTopology)
{
  if (reuseTopology && (tangentStiffnessMatrixOffset != NULL))
    *tangentStiffnessMatrixOffset = *tangentStiffnessMatrixOffset_;
  else
  {
    delete(tangentStiffnessMatrixOffset);
    tangentStiffnessMatrixOffset = new SparseMatrix(*tangentStiffnessMatrixOffset_);
  }
}

void IntegratorBaseSparse::AddTangentStiffnessMatrixOffset(SparseMatrix * tangentStiffnessMatrixOffset_)
{
  *tangentStiffnessMatrixOffset += *tangentStiffnessMatrixOffset_;
}

void IntegratorBaseSparse::ClearTangentStiffnessMatrixOffset()
{
  delete(tangentStiffnessMatrixOffset);
  tangentStiffnessMatrixOffset = NULL;
}

