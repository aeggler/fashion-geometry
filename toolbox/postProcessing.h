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

#endif //EXAMPLE_POSTPROCESSING_H
