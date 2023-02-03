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
//void startRetriangulation(vector<VectorXd>& polylineSelected);


void backTo3Dmapping(MatrixXd& adaptedPattern, MatrixXi& adaptedPattern_faces, MatrixXd& perfectPattern, MatrixXi& perfectPattern_faces ,
                     MatrixXd& perfectPatternIn3d, MatrixXi& perfectPatternIn3d_faces, MatrixXd& adaptedPatternIn3d);
void startRetriangulation(vector<VectorXd>& polylineSelected, MatrixXd& V2, MatrixXi& F2 );
void computeAllBetweens(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                        vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                        MatrixXd& currPattern, MatrixXd& Vg_pattern_orig, vector<VectorXd>& polyLineInput, vector<int>& connectedVert );
void computeAllBetweensNew(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                           vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                           MatrixXd& currPattern, MatrixXd& Vg_to, vector<VectorXd>& polyLineInput, vector<int>& connectedVert) ;
void mergeTriagulatedAndPattern(const vector<int> &connectedVert, MatrixXd& Vg_retri, MatrixXi& Fg_retri, MatrixXd& currPattern, MatrixXi& Fg_pattern);

void createHalfAvatarMap(MatrixXd& testMorph_V1, MatrixXi& testMorph_F1,
                         MatrixXd& testMorph_V1left, MatrixXi& testMorph_F1left,
                         MatrixXd& testMorph_V1right, MatrixXi& testMorph_F1right,
                         map<int, int>& leftHalfToFullFaceMap,  map<int, int>& rightHalfToFullFaceMap);

void initialGuessAdaption(MatrixXd& currPattern, MatrixXd& toPattern, MatrixXi& Fg_pattern);
void ensureAngle(MatrixXd& currPattern, MatrixXd& toPattern, MatrixXi& Fg_pattern);
void ensurePairwiseDist(MatrixXd& currPattern, MatrixXd& toPattern, MatrixXi& Fg_pattern);
#endif //EXAMPLE_POSTPROCESSING_H
