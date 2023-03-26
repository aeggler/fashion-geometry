//
// Created by Anna Maria Eggler on 17.01.23.
//

#ifndef EXAMPLE_POSTPROCESSING_H
#define EXAMPLE_POSTPROCESSING_H



#include <Eigen/Dense>
#include "seam.h"
#include <set>

using namespace std;
using namespace Eigen;

void smoothBetweenVertices(MatrixXd& currPattern, MatrixXi& Fg_pattern,vector<int>& startAndEnd);
void duplicatePattern( MatrixXd& currPattern, MatrixXi&  Fg_pattern_curr, MatrixXd& addedFabricPatternVg, MatrixXi& addedFabricPatternFg, VectorXd& T_symetry);


void  backTo3Dmapping(MatrixXd& adaptedPattern, MatrixXi& adaptedPattern_faces, MatrixXd& perfectPattern, MatrixXi& perfectPattern_faces ,
                      MatrixXd& perfectPatternIn3d, MatrixXi& perfectPatternIn3d_faces, MatrixXd& adaptedPatternIn3d,  MatrixXi& adaptedPatternIn3d_faces, bool symmetry, string garment );
void startRetriangulation(vector<VectorXd>& polylineSelected, MatrixXd& V2, MatrixXi& F2 );
void computeAllBetweens(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                        vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                        MatrixXd& currPattern, MatrixXd& Vg_pattern_orig, vector<VectorXd>& polyLineInput, vector<vector<int>>& connectedVertVec, vector<int>& patchId, vector<bool>& isAscVec );
void computeAllBetweensNew(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                           vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                           MatrixXd& currPattern, MatrixXd& Vg_to, vector<VectorXd>& polyLineInput, vector<vector<int>>& connectedVertVec, vector<int>& patchId, vector<bool>& isAscVec) ;

void mergeTriagulatedAndPattern(const vector<vector<int>>& connectedVertVec, const vector<int>& patchId,const vector<bool>& isAscVec,
                                MatrixXd& Vg_retri, MatrixXi& Fg_retri, MatrixXd& currPattern, MatrixXi& Fg_pattern, vector<int> & newFaces, string avName, string garment);
void createHalfAvatarMap(MatrixXd& testMorph_V1, MatrixXi& testMorph_F1,
                         MatrixXd& testMorph_V1left, MatrixXi& testMorph_F1left,
                         MatrixXd& testMorph_V1right, MatrixXi& testMorph_F1right,
                         map<int, int>& leftHalfToFullFaceMap,  map<int, int>& rightHalfToFullFaceMap);

void initialGuessAdaption(MatrixXd& currPattern, MatrixXd& mapToVg, MatrixXd& perfectPattern,  MatrixXi& Fg_pattern_curr, MatrixXi& mapToFg, bool symetry,
                          set<int> & cornerSet, map<int, int >& mapCornerToCorner, int origHalfSize, map<int, int>& halfPatternVertToFullPatternVertT, string garment);
void initialGuessAdaptionWithoutT(MatrixXd& currPattern, MatrixXd& mapToVg, MatrixXd& perfectPattern,  MatrixXi& Fg_pattern_curr, MatrixXi& mapToFg, MatrixXi pPFg);
void ensureAngle(MatrixXd& currPattern, MatrixXd& toPattern, MatrixXi& Fg_pattern, MatrixXi& fromPatternFg);
void ensurePairwiseDist(MatrixXd& currPattern, MatrixXd& toPattern, MatrixXi& Fg_pattern);

void createMapCornersToNewCorner(MatrixXd& currPattern,MatrixXd& mapToVg, vector<vector<pair<int, int>>>& cornerPerBoundary,// first is vert id, second ins loop id, but thats bullshit
                                 map<int, int>& mapCornerToCorner, vector<vector<int>>& boundaryL, map<int, int>& halfPatternVertToFullPatternVertT,
                                 map<int, int>& fullPatternVertToHalfPatternVertT, bool symetry, map<int, int>& halfPatternFaceToFullPatternFaceT, map<int, int>& fullPatternFaceToHalfPatternFaceT );

void updateCornerUtils(set<int>& cornerSet, vector<vector<pair<int, int>>>& cornerPerBoundary,
                       map<int, vector<pair<int, int>>>& seamIdPerCorner, map<int, int>& mapCornerToCorner, VectorXd& cornerVertices , map<int, int>& fullToHalfVert);
void updateCornerUtilsInverse(set<int>& cornerSet, vector<vector<pair<int, int>>>& cornerPerBoundary,
                       map<int, vector<pair<int, int>>>& seamIdPerCorner, map<int, int>& mapCornerToCorner, VectorXd& cornerVertices , map<int, int>& fullToHalfVert);

void updateSeamCorner( vector<seam*>& seamsList,  vector<minusOneSeam*> & minusOneSeams, map<int, int>& mapCornerToCorner,
                       vector<vector<int>>& boundaryL);

void stitchSeam(vector<int>& startAndEnd, MatrixXd& currPattern, MatrixXi& Fg_pattern_curr);
void computeAffection(VectorXd& d,double geoDistMax, MatrixXi& Fg_pattern_curr, VectorXd& affectedFaces);

void computeFinalJacobian(MatrixXd& Vg, MatrixXi& Fg,MatrixXd& Vg_gar,
                          MatrixXi& Fg_gar, VectorXd& jacUAdapted, VectorXd& jacVAdapted, VectorXd& jacDiffAdapted);
void stitchAdapted3D(MatrixXd& Vg, MatrixXi& Fg,MatrixXi& Fg_pattern_orig,vector<seam*>& seamsList,map<int, int >& mapCornerToCorner,
                     map<int, int> & halfPatternVertToFullPatternVert, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern );

void computeBoundaryEdges(MatrixXi& changedFitGarF, MatrixXi& boundaryOfToPattern);

#endif //EXAMPLE_POSTPROCESSING_H
