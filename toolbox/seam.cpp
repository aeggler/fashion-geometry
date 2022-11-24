//
// Created by Anna Maria Eggler on 25.11.22.
//

#include "seam.h"

seam::seam(int patch1Id,
int patch2Id,
int patch1startCornerId, // boundary loop index [patch1id][patch1startCornerid] to end
int patch2startCornerId,
int patch1endCornerId,
int patch2endCornerId) {
    this->patch1Id = patch1Id;
    this-> patch2Id = patch2Id;
    this-> patch1startCornerId = patch1startCornerId;
    this-> patch2startCornerId = patch2startCornerId;
    this-> patch1endCornerId = patch1endCornerId;
    this-> patch2endCornerId = patch2endCornerId;
}
