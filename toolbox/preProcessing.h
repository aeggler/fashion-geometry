//
// Created by Anna Maria Eggler on 30.01.23.
//

#ifndef EXAMPLE_PREPROCESSING_H
#define EXAMPLE_PREPROCESSING_H

#endif //EXAMPLE_PREPROCESSING_H
#include <Eigen/Dense>
#include <iostream>
#include <map>


using namespace std;
using namespace Eigen;


void setupCollisionConstraintsCall(Eigen::MatrixXi& collisionVert, vector<int> & pureCollVert, MatrixXd& testMorph_V1left, MatrixXi& testMorph_F1left, MatrixXd& p, int& numVert, double coll_EPS,
                               std::map<int, int> & leftHalfToFullFaceMap, vector<VectorXd> & CleftRight , vector<VectorXd>& NleftRight, VectorXi& closestFaceId,
                               MatrixXd& Vm, MatrixXi& Fm, MatrixXi& Fg);