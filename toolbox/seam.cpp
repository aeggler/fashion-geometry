//
// Created by Anna Maria Eggler on 25.11.22.
//

#include "seam.h"
#include <utility>

seam::seam(int patch1Id, int patch2Id,
int patch1startCornerId, // boundary loop index [patch1id][patch1startCornerid] to end
int patch2startCornerId,
int patch1endCornerId,
int patch2endCornerId,
int patch1startBoundaryLoopIdx,
int patch2startBoundaryLoopIdx,
int patch1endBoundaryLoopIdx,
int patch2endBoundaryLoopIdx, int theirDist
) {
    this -> patch1Id = patch1Id;
    this -> patch2Id = patch2Id;
    this -> patch1startCornerId = patch1startCornerId;
    this -> patch2startCornerId = patch2startCornerId;
    this -> patch1endCornerId = patch1endCornerId;
    this -> patch2endCornerId = patch2endCornerId;
    this -> patch1startBoundaryLoopIdx = patch1startBoundaryLoopIdx;
    this -> patch2startBoundaryLoopIdx = patch2startBoundaryLoopIdx;
    this -> patch1endBoundaryLoopIdx = patch1endBoundaryLoopIdx;
    this -> patch2endBoundaryLoopIdx = patch2endBoundaryLoopIdx;
    this -> length = theirDist;
}

std::pair<int, int> seam::getStartAndPatch1(){
    return std::make_pair(patch1startBoundaryLoopIdx, patch1Id);

}// ATTENTION: start and end correspond, but in order to iterate over the boundary we need to start the iteration from the end
std::pair<int, int> seam::getStartAndPatch2() {
    return std::make_pair(patch2endBoundaryLoopIdx, patch2Id);
}

int seam::seamLength(){
    return length;
}
