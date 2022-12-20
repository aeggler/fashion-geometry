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

class cutVertEntry{
public:
    int vert;
    int seamType;
    int seamIdInList;
    int patch;
    bool startCorner;
    bool endCorner;
    bool continuedCorner;
    Vector3d continuedDirection;
    bool bridgeFlag;
    bool finFlag;
    bool levelOne;

    cutVertEntry( int vert, int seamType, int seamIdInList ){
        this-> vert = vert;
        this -> seamType = seamType;
        this -> seamIdInList = seamIdInList;
    }
    cutVertEntry ( int vert, int seamType, int seamIdInList, int patch){
        this -> vert = vert;
        this -> seamType = seamType;
        this -> seamIdInList = seamIdInList;
        this -> patch = patch;
        levelOne = true;
    }
};

void computeTear(MatrixXd& fromPatternFile, MatrixXd&  currPattern, MatrixXi& Fg_pattern, MatrixXi& Fg_pattern_orig, vector<seam*>& seamsList,
                  vector<minusOneSeam*> & minusOneSeams, std::vector<std::vector<int> >& boundaryL, bool& finished,
                 const std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary, map<int, vector<pair<int, int>>>& seamIdPerCorner,
                 VectorXd& cornerVert,
                 vector<cutVertEntry*>& cutPositions,  map<int, pair<int, int>>& releasedVert);


int findWhichEdgeOfFace(int face, int v1, int v2, MatrixXi& Fg);

void projectBackOnBoundary(const MatrixXd & Vg_to, MatrixXd& p, const vector<seam*>& seamsList,const vector<minusOneSeam*> & minusOneSeams,
                           const MatrixXi& Fg_pattern,
                           const MatrixXi& Fg_pattern_orig, const std::vector<std::vector<int> >& boundaryL_toPattern, const std::vector<std::vector<int> >& boundaryL,
                           map<int, pair<int, int>>& releasedVert
                           );

#endif //EXAMPLE_PATTERNADAPTION_H