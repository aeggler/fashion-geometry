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

void vertex_componentsBasedOnFacet(MatrixXi& Fg_pattern, VectorXi& componentIdPerFace,VectorXi& componentIdPerVert, int n);
void initCollMeshCall( MatrixXd& Vm_left, MatrixXi& Fm_left,
                       MatrixXd& Vm_right, MatrixXi& Fm_right);

void setupCollisionConstraintsCall(Eigen::MatrixXi& collisionVert, vector<int> & pureCollVert,
                                   MatrixXd& testMorph_V1left,
                                   MatrixXi& testMorph_F1left,
                                   MatrixXd& testMorph_V1right,
                                   MatrixXi& testMorph_F1right,
                                   MatrixXd& p, int& numVert, double coll_EPS,
                               std::map<int, int> & leftHalfToFullFaceMap, std::map<int, int> & rightHalfToFullFaceMap, vector<VectorXd> & CleftRight , vector<VectorXd>& NleftRight, VectorXi& closestFaceId,
                               MatrixXd& Vm, MatrixXi& Fm, MatrixXi& Fg);

void createHalfSewingPattern(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, MatrixXd& Vg_pattern_half, MatrixXi& Fg_pattern_half,
                             map<int, int>& halfPatternFaceToFullPatternFace, map<int, int>& fullPatternFaceToHalfPatternFace,
                             map<int, int>& halfPatternVertToFullPatternVert ,
                             map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& insertedIdxToPatternVert,
                             VectorXi& isLeftVertPattern , MatrixXd& rightVert, string skirt);
void preProcessGarment(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, bool insertPlane, int symVert1, int symVert2 ,VectorXd& T, string garment, string garmentEXT);
void smoothLaplacian(MatrixXd& Vg, MatrixXi& Fg);
void splitAndSmooth(MatrixXd& Vg,MatrixXi& Fg,MatrixXd& Vg_pattern,MatrixXi& Fg_pattern,
                    MatrixXd& VgPatternRet,MatrixXi& FgPatternRet,
                    MatrixXd& VgRet, MatrixXi& FgRet, string garment, string garmentEXT );