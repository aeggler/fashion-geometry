//
// Created by Anna Maria Eggler on 08.12.22.
//

#ifndef EXAMPLE_PATTERNADAPTION_H
#define EXAMPLE_PATTERNADAPTION_H
#include <Eigen/Dense>
#include "patternAdaption.h"
#include "seam.h"
using namespace std;
using namespace Eigen;





void computeTear(MatrixXd& fromPatternFile, MatrixXd&  currPattern, MatrixXi& Fg_pattern,MatrixXi& Fg_pattern_orig, vector<seam*>& seamsList, std::vector<std::vector<int> >& boundaryL, bool& finished);


#endif //EXAMPLE_PATTERNADAPTION_H