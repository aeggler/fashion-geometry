//
// Created by Anna Maria Eggler on 17.01.23.
//

#ifndef EXAMPLE_POSTPROCESSING_H
#define EXAMPLE_POSTPROCESSING_H



#include <Eigen/Dense>
#include "seam.h"

using namespace std;
using namespace Eigen;

void smoothBetweenVertices(MatrixXd& currPattern, MatrixXi& Fg_pattern,vector<int>& startAndEnd);
void startRetriangulation(vector<VectorXd>& polylineSelected);


void backTo3Dmapping(MatrixXd& adaptedPattern, MatrixXi& adaptedPattern_faces, MatrixXd& perfectPattern, MatrixXi& perfectPattern_faces ,
                     MatrixXd& perfectPatternIn3d, MatrixXi& perfectPatternIn3d_faces, MatrixXd& adaptedPatternIn3d);

#endif //EXAMPLE_POSTPROCESSING_H
