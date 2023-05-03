//
// Created by Anna Maria Eggler on 30.03.23.
//

#ifndef EXAMPLE_EXPORTFUNCTIONS_H
#define EXAMPLE_EXPORTFUNCTIONS_H

#endif //EXAMPLE_EXPORTFUNCTIONS_H
#include <Eigen/Dense>
#include <set>

using namespace Eigen;
using namespace std;
void writeMTL(MatrixXd& Ka, MatrixXd& Ks, MatrixXd& Kd, MatrixXd& Vg, MatrixXi& Fg, string garment, string avName, double interp, string dir);

void clipDifference(vector<vector<int>>& boundaryL_adaptedFromPattern,vector<vector<int>>& boundaryL_toPattern,
                    MatrixXd & currPattern, MatrixXd& Vg_to,  vector<vector<VectorXd>>& returnVec);

void computeCols(int num, MatrixXd& cols);
void duplicateInitPattern(MatrixXd& Vg ,MatrixXi& Fg);

void addedSquare(MatrixXi Fg, MatrixXd Vg);
void movePatches();

void   insertToStartEnd(vector<int> &startAndEnd, std::set<int>& cornerset, MatrixXd& currPattern, MatrixXi& Fg_pattern_curr,
                        vector<vector<int>> &bd );