//
// Created by Anna Maria Eggler on 25.11.22.
//

#ifndef EXAMPLE_SEAM_H
#define EXAMPLE_SEAM_H


class seam {
private:
    int patch1Id;
    int patch2Id;

    int patch1startCornerId; // boundary loop index [patch1id][patch1startCornerid] to end
    int patch2startCornerId;
    int patch1endCornerId;
    int patch2endCornerId;

};


#endif //EXAMPLE_SEAM_H
