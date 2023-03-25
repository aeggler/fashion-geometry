//
// Created by Anna Maria Eggler on 08.12.22.
//

#ifndef EXAMPLE_PATTERNADAPTION_H
#define EXAMPLE_PATTERNADAPTION_H
#include <Eigen/Dense>
#include "patternAdaption.h"
#include "seam.h"
#include <set>

using namespace std;
using namespace Eigen;

class cutVertEntry{
public:
    int vert;
    int seamType;
    int seamIdInList;
    int patch;
    int leftCorner;
    int rightCorner;
    int cornerInitial;
    bool startCorner;
    bool endCorner;
    bool continuedCorner;
    Vector3d continuedDirection;
    Vector3d leftdirection;
    Vector3d rightdirection;
    Vector3d midVec;
    bool bridgeFlag;
    bool finFlag;
    bool levelOne;
    double stress;
    bool handled;
    int counterPartIdx;
    double stressWithCounter;
    int cutId;
    vector<std::pair<int, int>> dulicatePairs;
    vector<pair<int, int>> boundaryFrac;

    cutVertEntry( int vert, int seamType, int seamIdInList ){
        this-> vert = vert;
        this -> seamType = seamType;
        this -> seamIdInList = seamIdInList;
        handled = false;
        counterPartIdx  = -1;
        midVec = Vector3d::Zero();
        cutId = 1000;
    }
    cutVertEntry ( int vert, int seamType, int seamIdInList, int patch){
        this -> vert = vert;
        this -> seamType = seamType;
        this -> seamIdInList = seamIdInList;
        this -> patch = patch;
        levelOne = true;
        handled = false;
        counterPartIdx = - 1;
        midVec = Vector3d::Zero();
        cutId = 1000;

    }

    bool operator < (cutVertEntry* const & b){
        return stress < b->stress;
    }

};

int computeTear(bool inverseMap,
                MatrixXd& fromPatternFile,
                 MatrixXd& currPattern, // the current vertex positions
                 MatrixXi& Fg_pattern, // faces, can change as we insert new vertices
                 MatrixXd& lengthsOrig, // not to be changed, original faces, matches with Vg and fromPatternFile
                 vector<seam*>& seamsList, // all seams with partner
                 vector<minusOneSeam*> & minusOneSeams, // all boudary seams (seam type -1)
                 std::vector<std::vector<int> >& boundaryL, // the boundary loop, updated after each tear iteration so should be fresh when usig it
                 bool& finished,
                 vector<vector<pair<int, int>>>& cornersPerBoundary, // for each patch the corners on this boudnary loop, saved as vertex id and id in boundary loop (not updated!)
                 map<int, vector<pair<int, int>>>& seamIdPerCorner,// for each corner a map to a pair seamtype and seamid, if type >0 but id<0 its the second side of the seam, access by using (index +1)*(-1)
                 VectorXd& cornerVert,
                 vector<cutVertEntry*>& cutPositions,
                 map<int, pair<int, int>>& releasedVert,
                 std::set<int>& toPattern_boundaryVerticesSet,
                 set<int> & cornerSet,
                 set<int>& handledVerticesSet,
                 bool& prevFinished,
                 const bool & LShapeAllowed,
                 bool& prioInner,
                 bool& prioOuter, double tailor_lazyness,
                 const MatrixXi& mapFromFg,
                 double& setTheresholdlMid,
                 double& setTheresholdBound,
                 map<int, int> & fullPatternVertToHalfPatternVert, map<int, int>& halfPatternVertToFullPatternVert,
                 map<int, int> & halfPatternFaceToFullPatternFace,
                 bool& symetry,
                 set<int>& tipVert,
                 bool midFractureForbidden
                 );



int findWhichEdgeOfFace(int face, int v1, int v2, MatrixXi& Fg);
void updatePositionToIntersection(MatrixXd& p,int next,  const MatrixXd& Vg_bound);
void projectBackOnBoundary(const MatrixXd & mapToVg, MatrixXd& p, const vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,
                      const std::vector<std::vector<int> >& boundaryL_toPattern, const std::vector<std::vector<int> >& boundaryL,
                      map<int, pair<int, int>> & releasedVert, bool inverseMap, map<int,int> & fromtoToVertMapIfSplit,map<int,int> & seamFullHalf);

int  tearFurther(vector<cutVertEntry*>& cutPositions,
                 MatrixXd&  currPattern,
                 MatrixXi& Fg_pattern,
                 vector<seam*>& seamsList,
                 vector<minusOneSeam*>& minusOneSeams,
                 map<int, pair<int, int>> & releasedVert,
                 set<int>& toPattern_boundaryVerticesSet,
                 std::vector<std::vector<int> >& boundaryL,
                 set<int> & cornerSet, set<int>& handledVerticesSet,
                 bool& prevFinished,
                 const bool & preferManySmallCuts,
                 const bool & LShapeAllowed,
                 MatrixXd& patternEdgeLengths_orig,
                 MatrixXd& Vg_pattern_orig,
                 MatrixXi& Fg_pattern_orig,
                 bool& prioInner,
                 bool& prioOuter,
                 double& setTheresholdlMid,
                 double& setTheresholdBound,
                 map<int, int>& fullPatternVertToHalfPatternVert,
                 map<int, int>& halfPatternVertToFullPatternVert,
                 map<int, int> & halfPatternFaceToFullPatternFace);
int tearFurtherVisIdxHelper(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                            map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL,
                            set<int> & cornerSet, set<int>& handledVerticesSet,  bool& prevFinished, const bool & preferManySmallCuts, const bool & LShapeAllowed,
                            MatrixXd& patternEdgeLengths_orig, MatrixXd& Vg_pattern_orig, bool& prioInner,
                            bool& prioOuter );
void smoothCuts(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL, set<int> & cornerSet );
void setUpMap( const std::vector<std::vector<int> >& boundaryL,map<int,int> & fullPatternVertToHalfPatternVert);

void updatePatchId(vector<cutVertEntry*>& cutPositions, const std::vector<std::vector<int> >& boundaryLnew, vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams, map<int, int >& fullPatternVertToHalfPatternVert);
void fitVecToPointSet( MatrixXd& pointVec, VectorXd& vec );
void zipTears(vector<cutVertEntry*>& cutPositions, MatrixXd& Vg, MatrixXi& Fg, MatrixXi& mapFromFg, MatrixXd& mapFromVg, map<int, int>& halfPatternFaceToFullPatternFace, bool inverseMap);

#endif //EXAMPLE_PATTERNADAPTION_H