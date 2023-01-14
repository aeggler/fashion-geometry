//
// Created by Anna Maria Eggler on 25.11.22.
//

#ifndef EXAMPLE_SEAM_H
#define EXAMPLE_SEAM_H

#include <utility>
#include <map>
#include "Eigen/Dense"
#include <vector>

class seam {
private:
    int patch1Id;
    int patch2Id;
    int updatedPatchId1; int updatedPatchId2;

    int patch1startCornerId;
    int patch2startCornerId;
    int patch1endCornerId;
    int patch2endCornerId;

    int patch1startBoundaryLoopIdx;//TODO attention they are wrt the original boundary loop index , not any modified inserted versions!!
    int patch2startBoundaryLoopIdx;
    int patch1endBoundaryLoopIdx;
    int patch2endBoundaryLoopIdx;

    int length; // to avoid annoying length computation modulo


public:
    std::map<int, int> duplicates;
    std::map<int, int> duplicates2;// other side of the seam

    seam(int patch1Id, int patch2Id,
         int patch1startCornerId, // boundary loop index [patch1id][patch1startCornerid] to end
         int patch2startCornerId,
         int patch1endCornerId,
         int patch2endCornerId,
         int patch1startBoundaryLoopIdx,
         int patch2startBoundaryLoopIdx,
         int patch1endBoundaryLoopIdx,
         int patch2endBoundaryLoopIdx,
         int length,// TODO     we need both lengths to ensure it's ok even if we insert vertices in one position
         bool inverted
         );

    std::pair<int, int> getStartAndPatch1();

    std::pair<int, int> getStartAndPatch2();
    int getUpdatedPatch1(){
        return updatedPatchId1;
    };
    int getUpdatedPatch2(){
        return updatedPatchId2;
    };
    void updatePatch1(int newPatchId){
        updatedPatchId1 = newPatchId;
    }
    void updatePatch2(int newPatchId){
        updatedPatchId2 = newPatchId;
    }
    int seamLength();
    bool inverted;
    int getPatch1(){
        return patch1Id;
    }
    int getPatch2(){
        return patch2Id;
    }

    int getStart1(){
        return patch1startCornerId;
    }

    int getStart2(){
        return patch2startCornerId;
    }

    std::pair<int, int> getStartAndPatch2ForCorres();
    std::pair<int, int> getEndCornerIds(){
        return std::make_pair(patch1endCornerId, patch2endCornerId);
    }
    void modifyStart(int which, int newVal){
        if(which>0){
            patch1startCornerId= newVal;
        }else{
            patch2startCornerId= newVal;
        }
    }
    int getNextVert1(int currVert, std::vector<int>& boundaryL);
    int getPrevVert1(int currVert, std::vector<int>& boundaryL);
    int getNextVert2(int currVert, std::vector<int>& boundaryL);
    int getPrevVert2(int currVert, std::vector<int>& boundaryL);


};

class minusOneSeam{
private:
    int patchId;
    int updatedPatchId;
    int startVert;
    int endVert;
    int startIdInBoundaryLoop;
    int endIdInBoundaryLoop;
    int len;

public:
    std::map<int, int> duplicates;
    minusOneSeam(int patchId,
    int startVert,
    int endVert,
    int startIdInBoundaryLoop,
    int endIdInBoundaryLoop,
    int len);

    int getNextVert(int currVert, std::vector<int>& boundaryL);
    int getPrevVert(int currVert, std::vector<int>& boundaryL);

    void updatePatch(int newPatchId){
        updatedPatchId = newPatchId;
    }

    int getUpdatedPatch(){
        return updatedPatchId;
    };
    int getPatch(){
        return patchId;
    };
    int getStartVert(){
        return startVert;
    }
    int getStartIdx(){
        return startIdInBoundaryLoop;
    }
    int getEndVert(){
        return endVert;
    }

    int getLength(){
        return len ;
    }
    void modifyStart(int newVal){
        startVert = newVal;
    }
    void modifyEnd(int newVal){
        endVert = newVal;
    }
};


/* Brief: a function that computes all seams of a given patch layout and stores them in a list.
 *
 * Inputs: A list of all boundaries (igl function), a maps between pattern and garment vertices, vertexFace adjacency, for each face and vertex an id which connected component
 * (to access the correct patch vertex), edge vertices that mark beginnings and ends of each seam.
 *
 * Output: edges per boundary and seams list, a list of all seams stored as seam instances
 *
 * */
void computeAllSeams(const std::vector<std::vector<int> >& boundaryL, std::map<int,int>& vertexMapPattToGar,
                     std::map<std::pair<int, int>,int>& vertexMapGarAndIdToPatch,
                     std::vector<std::vector<int> >& vfAdj, Eigen::VectorXi& componentIdPerFace,
                     Eigen::VectorXi& componentIdPerVert,
                     Eigen::VectorXd& cornerVertices, std::vector<std::vector<std::pair<int, int>>>& vertAndLoopIdxPerCornerPerBoundary,
                     std::vector<seam*>& seamsList , std::vector<minusOneSeam*>& minusSeams,
                     std::map<int, std::vector<std::pair<int, int>>>& seamIdPerCorner

);



#endif //EXAMPLE_SEAM_H
