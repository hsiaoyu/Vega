/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.1                               *
 *                                                                       *
 * "laplacianMatrix" library , Copyright (C) 2016 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Hongyi Xu, Jernej Barbic                                *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Hongyi Xu, Jernej Barbic                                    *
 *                                                                       *
 * Funding: National Science Foundation                                  *
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

/*
   Computes the Laplacian matrix of a 3D tetrahedral mesh. The
   computed matrices are element-based (as opposed to vertex-based) : 
   they transform scalar quantities defined on the elements. 
*/

#ifndef _LAPLACIANMATRIX_H_
#define _LAPLACIANMATRIX_H_

#include "sparseMatrix.h"
#include "tetMesh.h"

class LaplacianMatrix
{
public:

  // computes the ``umbrella'' discrete Laplacian matrix L on the tets
  // dimension is T x T, where T is the number of tets
  // In each row i, the diagonal entry is the #tets neighboring to tet i. 
  // And the entries corresponding to the neighboring tet indices are set to -1.
  // if biLaplace=1, the routine computes the bilaplace matrix, i.e., L^T L
  static SparseMatrix * ComputeTetLaplacianMatrix(const TetMesh * tetMesh, int biLaplace=0);

  // volume-weighted Laplacian on tets
  static SparseMatrix * ComputeTetVolumeWeightedLaplacianMatrix(const TetMesh * tetMesh);

  // "FEM-style" Laplacian on tets, obtained by conjugating the vertex-based Laplacian with 
  // a volume-weighted matrix converting tet quantities to vertex quantities
  static SparseMatrix * ComputeTetFEMLaplacianMatrix(const TetMesh * tetMesh);
};

#endif

