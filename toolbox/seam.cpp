//
// Created by Anna Maria Eggler on 25.11.22.
//

#include "seam.h"
#include <utility>
#include "iostream"
#include "Eigen/Dense"
#include <map>
#include <set>
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
int patch2endBoundaryLoopIdx, int theirDist, bool invertedFlag
) {
    this -> patch1Id = patch1Id;
    this -> patch2Id = patch2Id;
    this -> updatedPatchId1 = patch1Id;
    this -> updatedPatchId2 = patch2Id;
    this -> patch1startCornerId = patch1startCornerId;
    this -> patch2startCornerId = patch2startCornerId;
    this -> patch1endCornerId = patch1endCornerId;
    this -> patch2endCornerId = patch2endCornerId;
    this -> patch1startCornerIdOld = patch1startCornerId;
    this -> patch2startCornerIdOld = patch2startCornerId;
    this -> patch1endCornerIdOld = patch1endCornerId;
    this -> patch2endCornerIdOld = patch2endCornerId;
    this -> patch1startBoundaryLoopIdx = patch1startBoundaryLoopIdx;
    this -> patch2startBoundaryLoopIdx = patch2startBoundaryLoopIdx;
    this -> patch1endBoundaryLoopIdx = patch1endBoundaryLoopIdx;
    this -> patch2endBoundaryLoopIdx = patch2endBoundaryLoopIdx;
    this -> length = theirDist;
    this -> inverted = invertedFlag;
}

std::pair<int, int> seam::getStartAndPatch1(){
    return std::make_pair(patch1startBoundaryLoopIdx, patch1Id);

}// ATTENTION: start and end correspond, but in order to iterate over the boundary we need to start the iteration from the end
std::pair<int, int> seam::getStartAndPatch2() {
    return std::make_pair(patch2endBoundaryLoopIdx, patch2Id);
}

// ATTENTION: start and end correspond, but in order to iterate over the boundary we need to start the iteration from the end
std::pair<int, int> seam::getStartAndPatch2ForCorres() {
    return std::make_pair(patch2startBoundaryLoopIdx, patch2Id);
}


int seam::seamLength(){
    return length;
}
int seam::getNextVert1(int currVert, std::vector<int>& boundaryL){
    int count=0;
    while (boundaryL[count]!= currVert && count<boundaryL.size()){
        count++;
        count %= boundaryL.size();
    }
    if(boundaryL[count]!=currVert){
        std::cout<<" vertex not found. please stop here"<<std::endl; return -1;
    }
    count++;
    return boundaryL[count % boundaryL.size()];

}
int seam::getPrevVert1(int currVert, std::vector<int>& boundaryL){
    int count=0;
    while (boundaryL[count]!= currVert && count<boundaryL.size()){
        count++;
        count %= boundaryL.size();
    }
    if(boundaryL[count]!=currVert){
        std::cout<<" vertex not found. please stop here"<<std::endl; return -1;
    }
    count--;
    if(count<0) count+= boundaryL.size();
    return boundaryL[count];

}

int seam::getNextVert2(int currVert, std::vector<int>& boundaryL){
    int count=0;
    int nextidx=0;
    while(boundaryL[nextidx]!= currVert && count<boundaryL.size()){
        count++; count %= boundaryL.size();
        nextidx = ( - count);
        if (nextidx < 0) nextidx += boundaryL.size();
        if (this->inverted) nextidx =  count % boundaryL.size();
    }if(boundaryL[nextidx]!= currVert)cout<<" stop here we dont find it!!!"<<endl;
//    cout<<boundaryL[nextidx]<<" found"<<endl;
    count++;
    count %= boundaryL.size();
    nextidx = (-count);
    if (nextidx < 0) nextidx += boundaryL.size();
    if (this->inverted) nextidx =  count % boundaryL.size();
//    cout<<boundaryL[nextidx]<<" returning "<<endl;
    return boundaryL[nextidx];

}
int seam::getPrevVert2(int currVert, std::vector<int>& boundaryL){
    int count=0;
    int nextidx=0;
    while(boundaryL[nextidx]!= currVert && count<boundaryL.size()){
        count++; count %= boundaryL.size();
        nextidx = ( - count);
        if (nextidx < 0) nextidx += boundaryL.size();
        if (this->inverted) nextidx =  count % boundaryL.size();
    }if(boundaryL[nextidx]!= currVert)cout<<" stop here we dont find it!!!"<<endl;
//    cout<<boundaryL[nextidx]<<" found"<<endl;
    count--;
//    count %= boundaryL.size();
    nextidx = (-count);
    if (nextidx < 0) nextidx += boundaryL.size();
    if (this->inverted) nextidx =  count % boundaryL.size();
//    cout<<boundaryL[nextidx]<<" returning "<<endl;
    return boundaryL[nextidx];

}
void push_to_Map(int whichKindOfSeam, int& idx, int& startId,int& startIdOther, map<int, vector<pair<int, int>>>& seamIdPerCorner){
    pair<int, int> whatSeamAndIdFirst = make_pair(whichKindOfSeam, idx);
    if(seamIdPerCorner.find(startId) == seamIdPerCorner.end()){
        vector<pair<int, int>> seamInfo ;
        seamInfo.push_back(whatSeamAndIdFirst);
        seamIdPerCorner[startId] = seamInfo;
    }else {
//        cout<<"found existing"<<endl;
        seamIdPerCorner[startId].push_back(whatSeamAndIdFirst);
    }

    if(whichKindOfSeam<0)return;

    pair<int, int> whatSeamAndIdSecond = make_pair(whichKindOfSeam,  -1* idx -1);
    if(seamIdPerCorner.find(startIdOther) == seamIdPerCorner.end()){
        vector<pair<int, int>> seamInfo ;
        seamInfo.push_back(whatSeamAndIdSecond);
        seamIdPerCorner[startIdOther] = seamInfo;
    }else{
//        cout<<"found existing second"<<endl;
        seamIdPerCorner[startIdOther].push_back(whatSeamAndIdSecond);
    }
}
void computeAllSeams(const std::vector<std::vector<int> >& boundaryL, std::map<int,int>& vertexMapPattToGar, std::map<std::pair<int, int>,int>& vertexMapGarAndIdToPatch,
                     std::vector<std::vector<int> >& vfAdj, Eigen::VectorXi& componentIdPerFace, Eigen::VectorXi& componentIdPerVert,
                     Eigen::VectorXd& cornerVertices, std::vector<std::vector<std::pair<int, int>>>& vertAndLoopIdxPerCornerPerBoundary, std::vector<seam*>& seamsList, std::vector<minusOneSeam*>& minusSeams,
                     map<int, vector<pair<int, int>>>& seamIdPerCorner, string garment
                     ){

    // we would like a seam to seam mapping where a seam is defined by its two endpoints
    Eigen::VectorXi isBoundaryVertexVec= Eigen::VectorXi::Zero(componentIdPerVert.rows());
    set<int> additionalCorners;
    if(garment == "top"){
        additionalCorners.insert(958);
        additionalCorners.insert( 726);
        additionalCorners.insert( 354);
        additionalCorners.insert( 267);
        additionalCorners.insert( 35);
        additionalCorners.insert( 1045);
    }else if (garment == "skirt_no2"){
        additionalCorners.insert(369);
        additionalCorners.insert(31);

        additionalCorners.insert(772);
        additionalCorners.insert(1110);

        additionalCorners.insert(1043);
        additionalCorners.insert(302);

        additionalCorners.insert(296);
        additionalCorners.insert(1037);
    }

//
    for(int i=0; i< boundaryL.size(); i++){
        for(int j=0; j < boundaryL[i].size(); j++){
            isBoundaryVertexVec(boundaryL[i][j]) = 1;
        }

    }
    for(int i=0; i< boundaryL.size(); i++){
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
            bool samePatchCrossover = (face1id1==face0id2 && face1id2 == face0id1);

            if(!(samePatch || samePatchCrossover) || (additionalCorners.find(v1) != additionalCorners.end())) {
                // two consecutive edges are not the same patch. We have found a corner
                cornerVertices(v1) = 1;
                edgesForThisBoundary.push_back(make_pair(v1, (j + 1) % boundary.size())); //pattern id
                // so far this worked...
            }

        }
        vertAndLoopIdxPerCornerPerBoundary.push_back(edgesForThisBoundary);
    }

    for(int i=0; i < boundaryL.size();  i++){//
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
            bool samePatchCrossover = (face1id1==face0id2 && face1id2 == face0id1);

            if(!(samePatch || samePatchCrossover)|| (additionalCorners.find(v1) != additionalCorners.end())){

                edgesForThisBoundary.push_back(make_pair(v1, (j + 1) % boundary.size())); //pattern id
                // but one must be the same since they are from the same patch

                if(edgesForThisBoundary.size()>1){
                    // we finish the previous here and add the seam

                    int myPatchId = i;
                    int otherPatchId;
                    if(face0id1== i){
                        otherPatchId = face0id2;

                    }else otherPatchId = face0id1;

                    if(myPatchId==otherPatchId){

                        int startId = edgesForThisBoundary[edgesForThisBoundary.size()-2].first;
                        int startIdInBoundaryIdx = edgesForThisBoundary[edgesForThisBoundary.size()-2].second;

                        int startIdOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[startId], otherPatchId)];
                        int maxID = componentIdPerVert.maxCoeff();
                        if( startIdOther == startId &&
                        vertexMapGarAndIdToPatch.find(make_pair(vertexMapPattToGar[startId], maxID+1+ otherPatchId)) != vertexMapGarAndIdToPatch.end()) startIdOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[startId], maxID+1+ otherPatchId)];

                        if(vertexMapPattToGar[startIdOther] != vertexMapPattToGar[startId]){
                            cout<<startId<<" "<<startIdInBoundaryIdx<<" "<<startIdOther<<" ------------------------the beginnings of a same patch pair do not match ------------------"<<endl;
                            cout<<(vertexMapGarAndIdToPatch.find(make_pair(vertexMapPattToGar[startId], maxID+1+ otherPatchId)) == vertexMapGarAndIdToPatch.end())<<endl;
                        }

                        // iterate over edges for this boundary since we only look for indices smaller than start id in boundary index and upt to there its same as edgesPerBoundary[i]
                        int counter=0;
                        int startIdOtherInBoundaryIdx;
                        while(edgesForThisBoundary[counter].first != startIdOther && counter <= startIdInBoundaryIdx){//edgesForThisBoundary.size()
                            counter++;
                        }
                        if(edgesForThisBoundary[counter].first != startIdOther) {
//                            cout<<" should not be there yet, dropping it "<<endl<<endl;
                            continue;
                        }
                        startIdOtherInBoundaryIdx = edgesForThisBoundary[counter].second;

                        counter=0;
                        int endIdOtherInBoundaryIdx;
                        int endIdOther = vertexMapGarAndIdToPatch[std::make_pair(v1g, otherPatchId)];
                        if(endIdOther == v1 &&
                                vertexMapGarAndIdToPatch.find(make_pair(vertexMapPattToGar[v1], maxID+1+ otherPatchId)) != vertexMapGarAndIdToPatch.end()
                        ) endIdOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[v1], maxID+1+ otherPatchId)];
                        if(vertexMapPattToGar[endIdOther] != vertexMapPattToGar[v1]){
                            cout<<"--------------------------- the ends of a same patch pair do not match ------------------"<<endl;
                        }
                        while (vertAndLoopIdxPerCornerPerBoundary[i][counter].first != endIdOther && counter < vertAndLoopIdxPerCornerPerBoundary[i].size()){
                             counter++;
                        }
                        if(vertAndLoopIdxPerCornerPerBoundary[i][counter].first != endIdOther){
                             cout<<" something is fishy, drop it "<<endl; continue;
                        }

                        endIdOtherInBoundaryIdx = vertAndLoopIdxPerCornerPerBoundary[i][counter].second;

                        int endIdx = (j+1) ;

                        int theirDist = startIdOtherInBoundaryIdx- endIdOtherInBoundaryIdx;if(theirDist<0) theirDist+= boundaryL[otherPatchId].size();
                        int mydist = endIdx - startIdInBoundaryIdx;
                        /* trial */
//                        cout<<myPatchId<<" "<<otherPatchId<<" "<<startId<<" "<< startIdOther<<" "<< v1<<" "<< endIdOther<<" "<<
//                                startIdInBoundaryIdx<<" "<< startIdOtherInBoundaryIdx<<" end1:"<< endIdx<<" end2:"<<
//                                endIdOtherInBoundaryIdx<<" dist:"<< theirDist<<endl;
//                        cout<<boundaryL[myPatchId][(endIdx-startIdInBoundaryIdx)/2+ startIdInBoundaryIdx]<<"trial"<<endl;
////nno!!
//                        int newmid = boundaryL[myPatchId][(endIdx-startIdInBoundaryIdx)/2+ startIdInBoundaryIdx];
//                        cornerVertices[newmid]= 1;
//                        v1= newmid;
//                        endIdOther = newmid;
//                        endIdx = (endIdx-startIdInBoundaryIdx)/2+ startIdInBoundaryIdx;
//                        endIdOtherInBoundaryIdx= (endIdx-startIdInBoundaryIdx)/2+ startIdInBoundaryIdx;
//                        theirDist/=2;
//                        edgesForThisBoundary.push_back(make_pair(newmid,(endIdx-startIdInBoundaryIdx)/2+ startIdInBoundaryIdx ));

                        /*end trial */
                        seam* newSeam = new seam (myPatchId, otherPatchId,startId, startIdOther, v1, endIdOther,
                                                  startIdInBoundaryIdx, startIdOtherInBoundaryIdx, endIdx,
                                                  endIdOtherInBoundaryIdx, theirDist, false
                        );

//                        cout<<" we use end tiral "<<otherPatchId <<" "<<startId<<endl<<endl<<endl;

//                        cout<<"seam no "<<  seamsList.size()<< " with start "<<startId<<" "<<startIdOther<<endl;
                        int seamListsize = seamsList.size();
                        if(startIdOther>=startId){
                            push_to_Map(1, seamListsize, startId, startIdOther, seamIdPerCorner);
                            seamsList.push_back(newSeam);
                        }

//                        cout<<seamIdPerCorner.size()<<" the size of the map int first"<<endl;


                    }
                    else if(myPatchId > otherPatchId && otherPatchId != -1){
                        //we know the other boundary is processed and we avoid duplicate seams
                        int startId = edgesForThisBoundary[edgesForThisBoundary.size()-2].first; // the previous
                        int startIdInBoundaryIdx = edgesForThisBoundary[edgesForThisBoundary.size()-2].second;
                        // corner can be duplicated, use the one before the corner
                        int nextbef = (startIdInBoundaryIdx+1) % boundaryL[myPatchId].size();
                        int idBeforeStart = boundary[nextbef];
                        int idBeforeStartOther;// = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[idBeforeStart], otherPatchId)];
                        if(vertexMapGarAndIdToPatch.find(make_pair(vertexMapPattToGar[idBeforeStart], otherPatchId))== vertexMapGarAndIdToPatch.end()) {
                            nextbef = startIdInBoundaryIdx-1;
                            if(nextbef < 0) nextbef +=  boundaryL[myPatchId].size();
                            idBeforeStart = boundary[nextbef];
                            idBeforeStartOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[idBeforeStart], otherPatchId)];
                            cout<<" new before start other fail!!!!!!!!!!!!!!!!!!!!!!!! "<< idBeforeStartOther<< endl;
                        }else{
                            idBeforeStartOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[idBeforeStart], otherPatchId)];
                        }
                        if(isBoundaryVertexVec(idBeforeStartOther)==0){
                            cout<<" we detected a pocket"<<endl; continue;
                            //TODO
                        }
                        // TODO SPECIAL CASE FOR SEAM LENGTH 2

                        int counter=0;
                        while(boundaryL[otherPatchId][counter]!= idBeforeStartOther && counter < boundaryL[otherPatchId].size()){
                            counter++;
                        }
                        if(boundaryL[otherPatchId][counter]!= idBeforeStartOther){
                            cout<<myPatchId<<"issue detected, not found"<<otherPatchId<<" "<<idBeforeStartOther<<endl;
                            cout<<componentIdPerVert[idBeforeStartOther]<<" iis it on this patch "<<otherPatchId<<endl;
                            for (auto it: boundaryL[otherPatchId]){
//                                cout<<it<<" item"<<endl;
                            }
                        // if partner not boundary vertex then issue solved we oonly look at thiis one
                        }

                        bool invertedflag = false;
                        int startIdOtherInBoundaryIdx =   (counter + 1)% boundaryL[otherPatchId].size();
                        if(vertexMapPattToGar[startId] != vertexMapPattToGar[boundaryL[otherPatchId][(counter + 1)% boundaryL[otherPatchId].size()]] ){
                            startIdOtherInBoundaryIdx = (counter - 1)% boundaryL[otherPatchId].size();
//                            cout<<"--------------inverted"<<endl;
                            invertedflag = true;
                        }

                        int startIdOther = (vertexMapPattToGar[startId] == vertexMapPattToGar[boundaryL[otherPatchId][(counter + 1)% boundaryL[otherPatchId].size()]] ) ?
                                boundaryL[otherPatchId][(counter+1)%boundaryL[otherPatchId].size() ] :  boundaryL[otherPatchId][counter-1]  ;
                        if(vertexMapPattToGar[startIdOther]!=vertexMapPattToGar[startId]){
                            cout<<"-----------they start dont match we have a problem ------------- "<<otherPatchId<<endl;
                            cout<<    boundaryL[otherPatchId][(counter+1)%boundaryL[otherPatchId].size() ] <<" ---------- "<<  boundaryL[otherPatchId][counter-1]<<endl;
                        }

                        counter=0;
                        int beforeEnd = v0g;
                        int otherBeforeEnd = vertexMapGarAndIdToPatch[std::make_pair(v0g, otherPatchId)];
                        if(isBoundaryVertexVec(otherBeforeEnd)==0){
                            cout<<" we detected a pocket "<<myPatchId<<endl; continue;
                            //TODO
                        }
                        while(boundaryL[otherPatchId][counter] != otherBeforeEnd && counter < boundaryL[otherPatchId].size()){
                            counter++;
                        }
                        int endIdOtherInBoundaryIdx;
                        if(counter==0) counter = boundaryL[otherPatchId].size();


                        int endIdOther =  (vertexMapPattToGar[v1]== vertexMapPattToGar[boundaryL[otherPatchId][(counter-1)]])?
                                boundaryL[otherPatchId][counter-1]:  boundaryL[otherPatchId][counter+1] ;//vertexMapGarAndIdToPatch[std::make_pair(v1g, otherPatchId)];

                        endIdOtherInBoundaryIdx = (vertexMapPattToGar[v1]== vertexMapPattToGar[boundaryL[otherPatchId][(counter-1)]])? counter-1: counter+1;

                        int endIdx = (j+1) % boundary.size();
//                        if(endIdx - startIdInBoundaryIdx<=2)cout<<" detected small with diff patches: id "<<seamsList.size()<<endl;
//
                        int theirDist =  startIdOtherInBoundaryIdx-endIdOtherInBoundaryIdx;
                        if(invertedflag) theirDist = endIdOtherInBoundaryIdx- startIdOtherInBoundaryIdx;
                        if(endIdOtherInBoundaryIdx > startIdOtherInBoundaryIdx && !invertedflag){
                            theirDist = startIdOtherInBoundaryIdx  +  (boundaryL[otherPatchId].size()-endIdOtherInBoundaryIdx);
                        }
                        if( endIdOtherInBoundaryIdx< startIdOtherInBoundaryIdx && invertedflag){
                            theirDist = endIdOtherInBoundaryIdx + (boundaryL[otherPatchId].size()-startIdOtherInBoundaryIdx);
                        }
                        if(vertexMapPattToGar[endIdOther]!=vertexMapPattToGar[v1]){
                            cout<<"-----------the ends dont match we have a problem ------------- "<<otherPatchId<<endl;
                        }

                        if(abs( (j+1) - startIdInBoundaryIdx) != abs(theirDist)){
                            cout<<boundary.size()<<" size "<<endIdx<<" "<< startIdInBoundaryIdx<<" something wrong about the seams "<<endIdOtherInBoundaryIdx <<" "<< startIdOtherInBoundaryIdx<<endl;
                            cout<<otherPatchId<<" other and my Id "<<myPatchId<<endl;
                            cout<<startId<<" "<< startIdOther<<" "<< v1<<" "<< endIdOther<<" the absolute found indices "<<endl;
                            cout<< endIdOtherInBoundaryIdx- startIdOtherInBoundaryIdx<<" inverted dist and mine "<< abs( (j+1) - startIdInBoundaryIdx)<<endl;
                        }

                        seam* newSeam = new seam (myPatchId, otherPatchId,startId, startIdOther, v1, endIdOther,
                                                      startIdInBoundaryIdx, startIdOtherInBoundaryIdx, endIdx,
                                                      endIdOtherInBoundaryIdx, theirDist, invertedflag);
//                        if(endIdx - startIdInBoundaryIdx<=2)cout<<endl<<"-----------------"<<" detected small with same patches: id "<<seamsList.size()<<" size: "<<theirDist<<" and patch id: "<<myPatchId<<endl;
//                        cout<<" we use a function "<<otherPatchId <<" "<<startId<<endl<<endl<<endl;

//                        cout<<"seam no "<<  seamsList.size()<< " with start "<<startId<<" "<<startIdOther<<endl;
                        int seamListsize = seamsList.size();
                        int one=1;
                        push_to_Map(one, seamListsize, startId, startIdOther, seamIdPerCorner);
//                        cout<<" ret "<< seamIdPerCorner[startId][0].first<<" "<<  seamIdPerCorner[startId][0].second<<endl;
//                        cout<<seamIdPerCorner.size()<<" the size of the map in second "<<endl;
                        seamsList.push_back(newSeam);

                    }else if(otherPatchId != -1 ){//todo case if it is just a pocket with smaller id than it' s parent
                         }else if (otherPatchId==-1){
                        // we found a global boundary
                        int startId = edgesForThisBoundary[edgesForThisBoundary.size()-2].first;
                        int startIdInBoundaryIdx = edgesForThisBoundary[edgesForThisBoundary.size()-2].second;
                        int endIdx = (j+1);
                        int mydist = endIdx - startIdInBoundaryIdx;
                        int endId = v1;
                        int len = endIdx - startIdInBoundaryIdx;
                        if(len<0) len+= boundaryL[i].size();
                        minusOneSeam* m0 = new minusOneSeam (myPatchId, startId, endId, startIdInBoundaryIdx ,endIdx, len);
//                        cout<<"negative seam no "<<  minusSeams.size()<< " with start "<<startId<<endl;
                        int seamListsize = minusSeams.size();
                        int one=-1; int zero=0;
                        push_to_Map(one, seamListsize, startId, zero, seamIdPerCorner);
//                        cout<<seamIdPerCorner.size()<<" the size of the map in third "<<endl;
//                        cout<<" ret "<< seamIdPerCorner[startId][0].first<<" "<<  seamIdPerCorner[startId][0].second<<endl;
//                        seamIdPerCorner(startId) = -minusSeams.size();
                        minusSeams.push_back(m0);

                    }
                }
            }
        }
//        edgesPerBoundary.push_back(edgesForThisBoundary);
        // check the first and the last if they are a seam
        if(edgesForThisBoundary.size()>0){//    todo check code above this is not updated

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
// todo check code above this is not updated
            if(otherpatch != -1 && otherpatch<i ){

                int startId = edgesForThisBoundary[edgesForThisBoundary.size()-1].first;
                int startIdInBoundaryIdx = edgesForThisBoundary[edgesForThisBoundary.size()-1].second;

                int startIdOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[startId], otherpatch)];
                bool otherExists = false;
                if(vertexMapGarAndIdToPatch.find(make_pair(vertexMapPattToGar[startId], otherpatch+ boundaryL.size()))!= vertexMapGarAndIdToPatch.end()){
                     otherExists = true;
                }
                int counter=0;
                while(vertAndLoopIdxPerCornerPerBoundary[otherpatch][counter].first != startIdOther){
                    counter++;
                }
                int startIdOtherInBoundaryIdx = vertAndLoopIdxPerCornerPerBoundary[otherpatch][counter].second;
                if(otherExists){

                    int prevMe = boundaryL[i][(startIdInBoundaryIdx+1)%( boundaryL[i]).size()];
                    int temp = startIdInBoundaryIdx-1; if(temp<0) temp+= boundaryL[i].size();
                    int nextMe = boundaryL[i][temp];
                    int prevOther = boundaryL[otherpatch][(startIdOtherInBoundaryIdx+1) % boundaryL[otherpatch].size() ];
                    temp = startIdOtherInBoundaryIdx -1; if(temp<0) temp+=  boundaryL[otherpatch].size();
                    int nextOther = boundaryL[otherpatch][temp];
                    if(vertexMapPattToGar[prevMe]==vertexMapPattToGar[prevOther] ||
                       vertexMapPattToGar[nextMe]==vertexMapPattToGar[nextOther]||
                       vertexMapPattToGar[prevMe]==vertexMapPattToGar[nextOther]||
                       vertexMapPattToGar[nextMe]==vertexMapPattToGar[prevOther]
                            ){

                    }else{
                        startIdOther = vertexMapGarAndIdToPatch[make_pair(vertexMapPattToGar[startId], otherpatch+ boundaryL.size())];
                        counter=0;
                        while(vertAndLoopIdxPerCornerPerBoundary[otherpatch][counter].first != startIdOther){
                            counter++;
                        }
                        startIdOtherInBoundaryIdx = vertAndLoopIdxPerCornerPerBoundary[otherpatch][counter].second;
                    }
                }
                counter=0;
                int endIdOther = vertexMapGarAndIdToPatch[std::make_pair(v0g, otherpatch)];
                while (vertAndLoopIdxPerCornerPerBoundary[otherpatch][counter].first != endIdOther){
                    counter++;
                }
                int endIdOtherInBoundaryIdx = vertAndLoopIdxPerCornerPerBoundary[otherpatch][counter].second;
                int theirDist = firstcorner.second + (boundary.size()-startIdInBoundaryIdx);
                theirDist = theirDist % boundary.size();

                seam* newSeam = new seam (i, otherpatch ,startId, startIdOther, v0, endIdOther,
                                          startIdInBoundaryIdx, startIdOtherInBoundaryIdx, firstcorner.second,
                                          endIdOtherInBoundaryIdx, theirDist, false
                );

                int seamListsize = seamsList.size();
                int one=1; int zero=0;
                push_to_Map(one, seamListsize, startId, startIdOther, seamIdPerCorner);
                seamsList.push_back(newSeam);

            }
            else if (otherpatch == -1){
                int startId = edgesForThisBoundary[edgesForThisBoundary.size()-1].first;
                int startIdInBoundaryIdx = edgesForThisBoundary[edgesForThisBoundary.size()-1].second;

                int endIdx = edgesForThisBoundary[0].second;
                int mydist =endIdx+ boundaryL[i].size() - startIdInBoundaryIdx;
                mydist = mydist % boundaryL[i].size();
                int endId = v0;
                minusOneSeam* m0 = new minusOneSeam (i, startId, endId, startIdInBoundaryIdx ,endIdx, mydist);

                if(boundaryL[i][startIdInBoundaryIdx] != startId || boundaryL[i][endIdx] !=endId){
                    cout<<" somehting is wrong in the seams computation of the minus one boundaries "<<endl;
                }

                int seamListsize = minusSeams.size();
                int one=-1; int zero=0;
                push_to_Map(one, seamListsize, startId, zero, seamIdPerCorner);
//                cout<<seamIdPerCorner.size()<<" the size of the map in fifth "<<endl;
                minusSeams.push_back(m0);



            }
        }
    }
//    cout<<seamsList.size()<<endl;
//    cout<<minusSeams.size()<<endl;

    for(int j=0; j<seamsList.size(); j++){
        seam* firstSeam = seamsList[j];
//        cout<<endl <<j<<" seam "<<endl;
        auto stP1 = firstSeam-> getStartAndPatch1();
        auto stP2 = firstSeam-> getStartAndPatch2ForCorres();
//        cout<<boundaryL[stP1.second][stP1.first]<<" starts "<<boundaryL[stP2.second][stP2.first]<<endl;
        int len = firstSeam -> seamLength();
//        cout<<len<<" len"<<endl;
        int boundLen1 = boundaryL[stP1.second].size();
//        cout<<boundLen1<<" bound "<<stP1.second<<" len"<<endl;
        int boundLen2 = boundaryL[stP2.second].size();
//        cout<<boundLen2<<" bound "<<stP2.second<<" len"<<endl;

    }

}

minusOneSeam::minusOneSeam(int patchId, int startVert, int endVert, int startIdInBoundaryLoop, int endIdInBoundaryLoop,
                           int len) {
    this-> patchId = patchId;
    this-> updatedPatchId = patchId;
    this -> startVert = startVert;
    this -> endVert = endVert;
    this -> startVertOld = startVert;
    this -> endVertOld = endVert;
    this -> startIdInBoundaryLoop = startIdInBoundaryLoop;
    this -> endIdInBoundaryLoop = endIdInBoundaryLoop;
    this -> len = len;

}
int minusOneSeam::getNextVert(int currVert, std::vector<int>& boundaryL){
    int count=0;
    while (boundaryL[count]!= currVert && count<boundaryL.size()){
        count++;
        count %= boundaryL.size();
    }
    if(boundaryL[count]!=currVert){
        std::cout<<" vertex not found. please stop here"<<std::endl; return -1;
    }
    count++;
    return boundaryL[count % boundaryL.size()];

}
int minusOneSeam::getPrevVert(int currVert, std::vector<int>& boundaryL){
    int count=0;
    while (boundaryL[count]!= currVert && count<boundaryL.size()){
        count++;
        count %= boundaryL.size();
    }
    if(boundaryL[count]!=currVert){
        std::cout<<" vertex not found. please stop here"<<std::endl; return -1;
    }
    count--;
    if (count<0) count+=  boundaryL.size();
    return boundaryL[count];

}