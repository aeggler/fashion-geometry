//
// Created by Anna Maria Eggler on 25.11.22.
//

#include "seam.h"
#include <utility>
#include "iostream"
#include "Eigen/Dense"
#include <map>
#include "adjacency.h"
#include <igl/facet_components.h>
#include <igl/vertex_components.h>

using namespace std;
using namespace Eigen;

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

void computeAllSeams(const std::vector<std::vector<int> >& boundaryL, std::map<int,int>& vertexMapPattToGar, std::map<std::pair<int, int>,int>& vertexMapGarAndIdToPatch,
                      std::vector<std::vector<int> >& vfAdj, Eigen::VectorXi& componentIdPerFace, Eigen::VectorXi& componentIdPerVert,
                      Eigen::VectorXd& edgeVertices, std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary, std::vector<seam*>& seamsList
                     ){
    // we would like a seam to seam mapping where a seam is defined by its two endpoints
    for(int i=0; i<boundaryL.size();  i++){//
        // we walk along the seam and check
        vector<pair<int, int>> edgesForThisBoundary; //pattern id
        vector<int> boundary = boundaryL[i];
        for(int j=0; j<boundary.size(); j++){
            int v0 = boundary[j]; int v0g= vertexMapPattToGar[v0];
            int v1 = boundary[(j+1) % boundary.size()]; int v1g = vertexMapPattToGar[v1];
            int v2 = boundary[(j+2) % boundary.size()]; int v2g = vertexMapPattToGar[v2];

            std::pair<int, int>  facese1, facese0;
            //vi is a pattern id , what is the corresponding 3d ID? it is not the same!
            adjacentFacesToEdge(v0g,v1g, vfAdj, facese0);// in 3D
            adjacentFacesToEdge(v1g, v2g,vfAdj,facese1);// in 3D
            // faces can be -1
            int face1id1 = -1;  int face1id2 = -1;  int face0id1 = -1;  int face0id2 = -1;
            if(facese0.first != -1) face0id1= componentIdPerFace(facese0.first); // component in 2D
            if(facese0.second != -1) face0id2= componentIdPerFace(facese0.second);
            if(facese1.first != -1) face1id1= componentIdPerFace(facese1.first);
            if(facese1.second != -1) face1id2= componentIdPerFace(facese1.second);

            bool samePatch = (face1id1==face0id1 && face1id2 == face0id2);
            bool samePatchCrossover = (face1id1==face0id2 && face1id2 == face0id2);

            if(!(samePatch || samePatchCrossover)){
                // two consecutive edges are not the same patch. We have found a corner
                edgeVertices(v1)= 1;
                edgesForThisBoundary.push_back(make_pair(v1, (j+1) % boundary.size())); //pattern id
                // so far this worked...

                // but one must be the same since they are from the same patch
                if(edgesForThisBoundary.size()>1){
                    // we finish the previous here and add the seam

                    int myPatchId = i;
                    int otherPatchId;
                    if(face0id1== i){
                        otherPatchId = face0id2;

                    }else otherPatchId = face0id1;

                    if(myPatchId>otherPatchId && otherPatchId != -1){
                        //we know the other boundary is processed and we avoid duplicate seams
                        int startId = edgesForThisBoundary[edgesForThisBoundary.size()-2].first; // the previous
                        int startIdInBoundaryIdx = edgesForThisBoundary[edgesForThisBoundary.size()-2].second;
                        int startIdOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[startId], otherPatchId)];
                        int counter=0;
                        while(edgesPerBoundary[otherPatchId][counter].first != startIdOther){
                            counter++;
                        }
                        int startIdOtherInBoundaryIdx = edgesPerBoundary[otherPatchId][counter].second;
                        counter=0;
                        int endIdOther = vertexMapGarAndIdToPatch[std::make_pair(v1g, otherPatchId)];
                        while (edgesPerBoundary[otherPatchId][counter].first != endIdOther){
                            counter++;
                        }
                        int endIdOtherInBoundaryIdx = edgesPerBoundary[otherPatchId][counter].second;
                        int endIdx = (j+1) % boundary.size();
                        int theirDist =  startIdOtherInBoundaryIdx-endIdOtherInBoundaryIdx;
                        if(endIdOtherInBoundaryIdx>startIdOtherInBoundaryIdx){
                            theirDist = startIdOtherInBoundaryIdx  +  (boundaryL[otherPatchId].size()-endIdOtherInBoundaryIdx);
                        }

                        if(abs( (j+1) - startIdInBoundaryIdx) != abs(theirDist)){
                            cout<<boundary.size()<<" size "<<endIdx<<" "<< startIdInBoundaryIdx<<" something wrong about the seams "<<endIdOtherInBoundaryIdx <<" "<< startIdOtherInBoundaryIdx<<endl;
                        }
                        seam* newSeam = new seam (myPatchId, otherPatchId,startId, startIdOther, v1, endIdOther,
                                                  startIdInBoundaryIdx, startIdOtherInBoundaryIdx, endIdx,
                                                  endIdOtherInBoundaryIdx, theirDist
                        );
                        seamsList.push_back(newSeam);
                    }
                }
            }
        }
        edgesPerBoundary.push_back(edgesForThisBoundary);
        // check the first and the last if they are a seam
        if(edgesForThisBoundary.size()>0){
            auto firstcorner = edgesForThisBoundary[0];
            int v0 = boundary[firstcorner.second]; int v0g= vertexMapPattToGar[v0];
            int v1 = boundary[(firstcorner.second - 1) % boundary.size()]; int v1g = vertexMapPattToGar[v1];

            std::pair<int, int>   faces;
            //vi is a pattern id , what is the corresponding 3d ID? it is not the same!
            adjacentFacesToEdge(v0g,v1g, vfAdj, faces);// in 3D
            int facesid0=-1; int facesid1=-1;
            if(faces.first != -1) facesid0 = componentIdPerFace(faces.first);
            if(faces.second != -1) facesid1 = componentIdPerFace(faces.second);

            int otherpatch;
            if(facesid0 == i){
                otherpatch = facesid1;
            }else{
                otherpatch = facesid0;
            }

            if(otherpatch != -1 && otherpatch<i ){

                int startId = edgesForThisBoundary[edgesForThisBoundary.size()-1].first;
                int startIdInBoundaryIdx = edgesForThisBoundary[edgesForThisBoundary.size()-1].second;
                int startIdOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[startId], otherpatch)];
                int counter=0;
                while(edgesPerBoundary[otherpatch][counter].first != startIdOther){
                    counter++;
                }
                int startIdOtherInBoundaryIdx = edgesPerBoundary[otherpatch][counter].second;
                counter=0;
                int endIdOther = vertexMapGarAndIdToPatch[std::make_pair(v0g, otherpatch)];
                while (edgesPerBoundary[otherpatch][counter].first != endIdOther){
                    counter++;
                }
                int endIdOtherInBoundaryIdx = edgesPerBoundary[otherpatch][counter].second;
                int theirDist = firstcorner.second + (boundary.size()-startIdInBoundaryIdx);
                theirDist = theirDist % boundary.size();

                seam* newSeam = new seam (i, otherpatch ,startId, startIdOther, v0, endIdOther,
                                          startIdInBoundaryIdx, startIdOtherInBoundaryIdx, firstcorner.second,
                                          endIdOtherInBoundaryIdx, theirDist
                );
                seamsList.push_back(newSeam);
            }
        }
    }


}
