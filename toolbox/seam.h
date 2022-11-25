//
// Created by Anna Maria Eggler on 25.11.22.
//

#ifndef EXAMPLE_SEAM_H
#define EXAMPLE_SEAM_H

#include <utility>
#include <map>
#include "Eigen/Dense"

class seam {
private:
    int patch1Id;
    int patch2Id;

    int patch1startCornerId;
    int patch2startCornerId;
    int patch1endCornerId;
    int patch2endCornerId;

    int patch1startBoundaryLoopIdx;
    int patch2startBoundaryLoopIdx;
    int patch1endBoundaryLoopIdx;
    int patch2endBoundaryLoopIdx;

    int length; // to avoid annoying length computation modulo

public:
    seam(int patch1Id, int patch2Id,
         int patch1startCornerId, // boundary loop index [patch1id][patch1startCornerid] to end
         int patch2startCornerId,
         int patch1endCornerId,
         int patch2endCornerId,  int patch1startBoundaryLoopIdx,
    int patch2startBoundaryLoopIdx,
    int patch1endBoundaryLoopIdx,
    int patch2endBoundaryLoopIdx, int length) ;

    std::pair<int, int> getStartAndPatch1();
    std::pair<int, int> getStartAndPatch2();
    int seamLength();

    };
void computeAllSeams(const std::vector<std::vector<int> >& boundaryL, std::map<int,int>& vertexMapPattToGar, std::map<std::pair<int, int>,int>& vertexMapGarAndIdToPatch,
                     std::vector<std::vector<int> >& vfAdj, Eigen::VectorXi& componentIdPerFace, Eigen::VectorXi& componentIdPerVert,
                     Eigen::VectorXd& edgeVertices, std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary, std::vector<seam*>& seamsList
);


#endif //EXAMPLE_SEAM_H
