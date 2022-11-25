//
// Created by Anna Maria Eggler on 25.11.22.
//

#ifndef EXAMPLE_SEAM_H
#define EXAMPLE_SEAM_H

#include <utility>

class seam {
private:
    int patch1Id;
    int patch2Id;

    int patch1startCornerId; // boundary loop index [patch1id][patch1startCornerid] to end
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


#endif //EXAMPLE_SEAM_H
