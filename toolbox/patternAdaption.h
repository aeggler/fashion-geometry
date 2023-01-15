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
    bool bridgeFlag;
    bool finFlag;
    bool levelOne;
    double stress;
    bool handled;
    int counterPartIdx;
    double stressWithCounter;

    cutVertEntry( int vert, int seamType, int seamIdInList ){
        this-> vert = vert;
        this -> seamType = seamType;
        this -> seamIdInList = seamIdInList;
        handled = false;
        counterPartIdx  = -1;
    }
    cutVertEntry ( int vert, int seamType, int seamIdInList, int patch){
        this -> vert = vert;
        this -> seamType = seamType;
        this -> seamIdInList = seamIdInList;
        this -> patch = patch;
        levelOne = true;
        handled = false;
        counterPartIdx = - 1;

    }

    bool operator < (cutVertEntry* const & b){
        return stress < b->stress;
    }

};

void computeTear(MatrixXd& fromPatternFile,
                 MatrixXd& currPattern, // the current vertex positions
                 MatrixXi& Fg_pattern, // faces, can change as we insert new vertices
                 MatrixXi& Fg_pattern_orig, // not to be changed, original faces, matches with Vg and fromPatternFile
                 vector<seam*>& seamsList, // all seams with partner
                 vector<minusOneSeam*> & minusOneSeams, // all boudary seams (seam type -1)
                 std::vector<std::vector<int> >& boundaryL, // the boundary loop, updated after each tear iteration so should be fresh when usig it
                 bool& finished,
                 const std::vector<std::vector<std::pair<int, int>>>& cornersPerBoundary, // for each patch the corners on this boudnary loop, saved as vertex id and id in boundary loop (not updated!)
                 map<int, vector<pair<int, int>>>& seamIdPerCorner,// for each corner a map to a pair seamtype and seamid, if type >0 but id<0 its the second side of the seam, access by using (index +1)*(-1)
                 VectorXd& cornerVert,
                 vector<cutVertEntry*>& cutPositions,
                 map<int, pair<int, int>>& releasedVert,
                 std::set<int>& toPattern_boundaryVerticesSet,
                 set<int> & cornerSet,
                 set<int>& handledVerticesSet,
                 MatrixXd& Vg,
                 bool& prevFinished
                 );


int findWhichEdgeOfFace(int face, int v1, int v2, MatrixXi& Fg);
void updatePositionToIntersection(MatrixXd& p,int next,  const MatrixXd& Vg_bound);
void projectBackOnBoundary(const MatrixXd & Vg_to, MatrixXd& p, const vector<seam*>& seamsList,const vector<minusOneSeam*> & minusOneSeams,
                           const MatrixXi& Fg_pattern,
                           const MatrixXi& Fg_pattern_orig, const std::vector<std::vector<int> >& boundaryL_toPattern,
                           map<int, pair<int, int>>& releasedVert, bool visFlag
                           );

void tearFurther(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                 map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL,
                 set<int> & cornerSet,
                 set<int>& handledVerticesSet,
                 bool& prevFinished,
                 const bool &preferManySmallCuts
);
void smoothCuts(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL, set<int> & cornerSet );

void updatePatchId(vector<cutVertEntry*>& cutPositions, const std::vector<std::vector<int> >& boundaryLnew, vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams);
#endif //EXAMPLE_PATTERNADAPTION_H