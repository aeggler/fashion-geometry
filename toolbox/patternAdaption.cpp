//
// Created by Anna Maria Eggler on 08.12.22.
//
#include <Eigen/Dense>
#include "patternAdaption.h"
#include <iostream>
#include <queue>
#include "seam.h"
#include <igl/edge_lengths.h>
#include <igl/adjacency_list.h>
#include "igl/boundary_loop.h"
#include "igl/is_vertex_manifold.h"
#include "/Library/gurobi1000/macos_universal2/include/gurobi_c++.h"
#include "adjacency.h"
#include <igl/HalfEdgeIterator.h>
#include "MathFunctions.h"
#include <igl/barycentric_coordinates.h>
#include <igl/barycentric_interpolation.h>

using namespace std;
using namespace Eigen;

MatrixXd lengthsOrig;
map<int, cutVertEntry *> cveStartPositionsSet;
set<int> cutThroughCornerVertices;
int findWhichEdgeOfFace(int face, int v1, int v2, MatrixXi& Fg){
    int faceidxv1, faceidxv2;
    for(int j=0; j<3; j++){
        if(Fg(face, j) == v1) faceidxv1 = j;
        if(Fg(face, j) == v2) faceidxv2 = j;

    }
    int whichEdge = -1;
    //[1,2],[2,0],[0,1] is the order of igl edge lengths
    if(  (faceidxv1==1 && faceidxv2==2)   ||  (faceidxv2==1 && faceidxv1==2)){
        whichEdge = 0;
    }  else if ( (faceidxv1==0 && faceidxv2==2)   ||  (faceidxv2==0 && faceidxv1==2) ){
        whichEdge=1;
    }else {
        whichEdge = 2;
    }
    return whichEdge;
}

void findIndxInBoundaryloop(vector<int> & boundary, int& target,int& idx){
    idx = 0;
    while (boundary[idx] != target && idx < boundary.size()){
        idx++;
    }
    if(boundary[idx] != target){
        cout<<idx<<" ERROR IN LP WE HAVE NOT FOUND THE INDEX of "<<target<<" "<<boundary.size()<<endl;
        for(int i=0; i< boundary.size(); i++){
            cout<<boundary[i]<<" ";
        }
    }
}
bool canBeSplit(int vertId,  vector<vector<int> >& vfAdj){
    return vfAdj[vertId].size()>1;
}

void updateWithNewId(MatrixXi& Fg, int vert,int rightFaceId, int newVertIdx ){
    if(Fg(rightFaceId , 1) == vert ){
        Fg(rightFaceId , 1) = newVertIdx;
    }
    else if(Fg(rightFaceId  , 0)== vert ){
        Fg(rightFaceId , 0) = newVertIdx;
    }
    else {
        Fg(rightFaceId , 2) = newVertIdx;
    }
}
void addToDulicateList( cutVertEntry*& cve, vector<seam*>& seamsList,  vector<minusOneSeam*> & minusOneSeams, int newVertIdx, bool global){
    if(global){
        if(cve->seamType<0){
            minusOneSeams[cve->seamIdInList]->duplicatesGlob[cve->vert]= newVertIdx;
        }else{
            if(cve->seamIdInList >= 0){
                seamsList[cve-> seamIdInList]-> duplicatesGlob[cve->vert]=newVertIdx;
            }else{
                seamsList[(cve-> seamIdInList+1)*(-1)]->duplicatesGlob2[cve->vert]= newVertIdx;
            }
        }
    }else{
        if(cve->seamType<0){
            minusOneSeams[cve->seamIdInList]->duplicates[cve->vert]= newVertIdx;
        }else{
            if(cve->seamIdInList >= 0){
                seamsList[cve-> seamIdInList]-> duplicates[cve->vert]=newVertIdx;
            }else{
                seamsList[(cve-> seamIdInList+1)*(-1)]->duplicates2[cve->vert]= newVertIdx;
            }
        }
    }

}
bool isRight(Vector3d a,Vector3d b,Vector3d c ){
    // right gets  newVertIdx
    if(((b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0))) ==0 )cout<<"zero"<<endl;
//    cout<<((b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0)))<<" is right if smaller 0 "<<endl;
    // TODO HANDLE THIS CASE PROPERLY!!
    return ((b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0))) < 0;

}
void computeMidVecBasedOnPCA(VectorXd& midVec, vector<vector<int>>& vvAdj, vector<vector<int>>& vfAdj, MatrixXd& Vg ,
                             MatrixXd& lengthsCurr, MatrixXd& lengthsOrig, MatrixXi& Fg_pattern,int& vert, VectorXd& ws, MatrixXd& dirs, double& lenMid){
    int n = vvAdj[vert].size();
    ws.resize(n);
    midVec.resize(3);
    dirs.resize(n, 3);
    std::pair<int, int>  faces;

    for(int i=0; i<vvAdj[vert].size(); i++) {
        int otherVert = vvAdj[vert][i];
        dirs.row(i) = (Vg.row(otherVert) - Vg.row(vert)).normalized();
        adjacentFacesToEdge(vert, otherVert, vfAdj, faces );
        //face area weighting needed?

        if(faces.first!= -1){
            int whichEdgeLeft = findWhichEdgeOfFace(faces.first, vert, otherVert, Fg_pattern);
            ws(i) += lengthsCurr(faces.first, whichEdgeLeft)/lengthsOrig(faces.first, whichEdgeLeft);
        }
        else if(faces.second!= -1){
            int whichEdgeRight = findWhichEdgeOfFace(faces.second, vert, otherVert, Fg_pattern);
            ws(i) += lengthsCurr(faces.second, whichEdgeRight)/lengthsOrig(faces.second, whichEdgeRight);
        }
        dirs.row(i) *= ws(i);
        dirs.row(i) += Vg.row(vert);
        cout<<dirs.row(i)<<" direction , stress "<<ws(i)<<" for " <<otherVert<<endl;
    }
    fitVecToPointSet( dirs, midVec );
    midVec(2) = 0;// Vg(vert, 2); // just to get the 3rd constant dimension right
    lenMid =  midVec.norm();
    midVec = midVec.normalized();
    cout<<" the newly computed midvec of "<<vert<<" is "<<midVec.transpose()<<endl;
}

void computeMidVecBasedOnStress(VectorXd& midVec, vector<vector<int>>& vfAdj, MatrixXd& Vg ,MatrixXd& Vg_orig ,
                           MatrixXi& Fg, int& vert, double& lenMid){
    // new approach: look at the stress vectors of all adjacent faces, take the average
    Vector3d uVecs, vVecs;
    for(int i=0; i<vfAdj[vert].size(); i++){
        int face = vfAdj[vert][i];
        Vector3d v0 = Vg_orig.row(Fg(face, 0)).transpose();
        Vector3d v1 = Vg_orig.row(Fg(face, 1)).transpose();
        Vector3d v2 = Vg_orig.row(Fg(face, 2)).transpose();

        Vector3d bary = (v0 + v1 + v2)/3;
        Vector3d u = bary; u(0)+= 1;
        Vector3d v = bary; v(1)+= 1;
        MatrixXd ubary, vbary;
        MatrixXd input(1, 3);

        input.row(0)= u;
        igl::barycentric_coordinates( input, v0.transpose(), v1.transpose(), v2.transpose(), ubary);
        input.row(0)= v;
        igl::barycentric_coordinates( input, v0.transpose(), v1.transpose(), v2.transpose(), vbary);

        Vector3d v0new = Vg.row(Fg(face, 0)).transpose();
        Vector3d v1new = Vg.row(Fg(face, 1)).transpose();
        Vector3d v2new = Vg.row(Fg(face, 2)).transpose();
        Vector3d barynew = (v0new + v1new + v2new)/3;

        Vector3d uVec = (ubary(0) * v0new + ubary(1) * v1new + ubary(2) * v2new) -  barynew;
        Vector3d vVec = (vbary(0) * v0new + vbary(1) * v1new + vbary(2) * v2new) -  barynew;
        uVecs += uVec;
        vVecs += vVec;

    }
    uVecs/= vfAdj[vert].size();
    vVecs/= vfAdj[vert].size();
    VectorXd no;
    if(uVecs.norm() > vVecs.norm()){
        lenMid = uVecs.norm();
        no = uVecs.normalized();
    }else{
        lenMid = vVecs.norm();
        no = vVecs.normalized();
    }
    // rotate counterclockwise 90 deg. if wrong direction will be fixed after!
    midVec.resize(3); // = no;  //
    midVec(2) = 0;
    midVec (0)= -no(1);// -y
    midVec (1)= no(0);// x

}
bool checkIfTearIsUseful(int vert, Vector3d& cutDirection,  vector<vector<int>>& vvAdj,  vector<vector<int>>& vfAdj, MatrixXd& Vg ,
                         MatrixXd& lengthsCurr,  MatrixXd& lengthsOrig, MatrixXi& Fg_pattern, VectorXd& ws, bool preComputed ){
    // check if it makes sense, i.e. releases stress cutting in the direction
    // if the dot product for at least one adjacent vertex of the one we are going to cut is large enough we allow to cut further
    // todo 3.1. maybe a better option is not to ignore ones with wrong direction but to actually clamp the stess

    double thereshold = 1.051;// todo fix this!! it cuts pretty much always
    double thW = 1.1; double thDot = 0.3;
    bool flag = false;
    if(!preComputed) ws.resize(vvAdj[vert].size());

    std::pair<int, int>  faces;
    for(int i=0; i<vvAdj[vert].size(); i++){
        int otherVert = vvAdj[vert][i];
        double w=0;
        if(!preComputed){
            adjacentFacesToEdge(vert, otherVert, vfAdj, faces );
            if(faces.first!= -1){
                int whichEdgeLeft = findWhichEdgeOfFace(faces.first, vert, otherVert, Fg_pattern);
                w += lengthsCurr(faces.first, whichEdgeLeft)/lengthsOrig(faces.first, whichEdgeLeft);
                cout<<otherVert<<" Length now:"<<lengthsCurr(faces.first, whichEdgeLeft)<<" and orig"<<lengthsOrig(faces.first, whichEdgeLeft)<<endl;
            }else if(faces.second!= -1){
                int whichEdgeRight = findWhichEdgeOfFace(faces.second, vert, otherVert, Fg_pattern);
                w += lengthsCurr(faces.second, whichEdgeRight)/lengthsOrig(faces.second, whichEdgeRight);
                cout<<otherVert<<" Length now:"<< lengthsCurr(faces.second, whichEdgeRight) <<" and orig"<< lengthsOrig(faces.second, whichEdgeRight) <<endl;

            }
            ws(i) = w;
        }else{
            w=ws(i);
        }

        Vector3d vecDir = Vg.row(otherVert);
        vecDir -= Vg.row(vert);
        vecDir = vecDir.normalized();

        if(w > 2){
            w = 2;
        }
        auto dotdir = (1-abs(vecDir.dot(cutDirection.normalized())));
        if(dotdir * w > thereshold || (w > thW && dotdir > thDot ) ){
            flag= true;
        }

    }
    if(!flag){
        cout<<cutDirection.transpose()<<" direction"<<endl;
        // do everythiing again but now print
        for(int i=0; i<vvAdj[vert].size(); i++){
            int otherVert = vvAdj[vert][i];
            double w = ws(i);

            Vector3d vecDir = Vg.row(otherVert);
            vecDir -= Vg.row(vert);
            vecDir = vecDir.normalized();

          cout<<otherVert<<" RESULT thereshold "<< (1-abs(vecDir.dot(cutDirection.normalized()))) <<" "<< w <<" = "<<(1-abs(vecDir.dot(cutDirection.normalized()))) * w <<endl;
        }

    }
    return flag;
}
map<int, vector<int>> releasedVertNew;
map<int, vector<int>> cornerToSeams;
void addToMapIfNotExisting( int key, int i){
    if(cornerToSeams.find(key) != cornerToSeams.end()){
        cornerToSeams[key].push_back(i);
    }else{
        vector<int> vec;
        vec.push_back(i);
        cornerToSeams[key] = vec;
    }

}
void getBoundaryNextVert(const bool& startCorner ,const int seamType, const int seamIdInList, const int minusOne, const int plusOne, int& next ){

    if(startCorner){// does not matter if it is a starter or not

        if(seamType>0){
            if(seamIdInList>=0){
                next = minusOne;
            }else{
                next = plusOne;
            }
        }else{
            cout<<"the next one is "<<minusOne<<endl;

            next = minusOne;
        }

    }else{
        if(seamType > 0){
            if(seamIdInList >= 0){
                next = plusOne;
            }else{
                next = minusOne;
            }
        }else{
            next = plusOne;
        }
    }
}

void computeOrientationOfFaces(VectorXd& signs,vector<int> vfAdj, MatrixXi& Fg, MatrixXd& Vg){
    signs.resize(vfAdj.size());
    for(int i=0; i<vfAdj.size(); i++){
        int idx1, idx2, idx3;
        idx1 = Fg(vfAdj[i], 0);
        idx2 = Fg(vfAdj[i], 1);
        idx3 = Fg(vfAdj[i], 2);
        double part1 = (Vg(idx2, 1) - Vg(idx1, 1)) * (Vg(idx2, 0) - Vg(idx1, 0));
        double part2 = (Vg(idx3, 2) - Vg(idx1,1)) * (Vg(idx3, 0) - Vg(idx1, 0));

        signs(i) = ((part1-part2) );
    }
}
//void insertFracPart(int vert, vector<vector<int> >& vfAdj, MatrixXd& Vg, MatrixXd& Vg_pattern_orig,  MatrixXi& Fg, bool& firstLayer,  <pair<pair<int, VectorXd>, int>>& frac ){
//    VectorXd midVec;
//    double lenMidVec;
//    computeMidVecBasedOnStress(midVec, vfAdj, Vg, Vg_pattern_orig, Fg, vert, lenMidVec);
//
//}
//void preComputePath(cutVertEntry*& cve,
//                    MatrixXd& Vg, // this is the current pattern we modify
//                    MatrixXi& Fg, // this will be modified and have entries that are not in the original pattern
//                    vector<vector<int> >& vfAdj,
//                    std::vector<std::vector<int> >& boundaryL,
//                    vector<seam*>& seamsList,
//                    vector<minusOneSeam*> & minusOneSeams,
//                    MatrixXd& Vg_pattern_orig){
//    // contains a pair of face and bary coords and which vert is closest, i.e at which vec we would have to continue
//    queue frac <pair<pair<int, VectorXd>, int>>;
//    insertFracPart( cve->vert, vfAdj, Vg, Vg_pattern_orig, Fg, true, frac );
//    while(get<1>frac ! in boundary){
//        //we insert another part
//        int currVert = lastElementOfFrac;
//        insertFracePart(currVert, vfAdj, Vg, Vg_pattern_orig, Fg, true, frac);
//
//    }
//
//}

void splitVertexFromCVE( cutVertEntry*& cve,
                         MatrixXd& Vg, // this is the current pattern we modify
                         MatrixXi& Fg, // this will be modified and have entries that are not in the original pattern
                         vector<vector<int> >& vfAdj,
                         std::vector<std::vector<int> >& boundaryL,
                         vector<seam*>& seamsList,
                         vector<minusOneSeam*> & minusOneSeams,
                         map<int, pair<int, int>> & releasedVert,
                         set<int>& toPattern_boundaryVerticesSet,
                         MatrixXd& lengthsCurr,
                         set<int> & cornerSet,
                         set<int>& handledVerticesSet ,
                         const bool & LShapeAllowed,
                         MatrixXd& Vg_pattern_orig){

    cout<<" patch "<<cve->patch<<" of "<<boundaryL.size()<<" L Shape allowed? "<<LShapeAllowed <<endl;

    if(cve-> finFlag) {
        cout<<cve->vert<<" done already "<<endl;
        return;
    }
    cve->handled = true;

// todo als the end of a cut through can be released
    if(handledVerticesSet.find(cve->vert) != handledVerticesSet.end()){
        cout<<"Handled by other seams already."<<endl;

        // if it is a corner (and cut from both sides(,) or a cut position itself, mabe we can still cut
        bool cornerOrStart = ((cveStartPositionsSet.find(cve->vert) != cveStartPositionsSet.end()) || cornerSet.find(cve->vert) != cornerSet.end()
                || cutThroughCornerVertices.find(cve->vert) != cutThroughCornerVertices.end());
        cout<<(cutThroughCornerVertices.find(cve->vert) != cutThroughCornerVertices.end())<<" in this set? "<<endl;
        if( cornerOrStart && LShapeAllowed){
            cout<<"But it is a corner and we allow L Shapes. Go on."<<endl;

        }else{
            cout<<"And really it is in the middle. Stop ."<<endl;
            cve->finFlag = true;
            return ;
        }

    }

    // has to be the original boundary Loop, else the patch id might not correspond with the new
    // solution: mapBoundaryToPatch with new patches
    auto boundary = boundaryL[cve->patch];

    double eps = 0.01;
    vector<vector<int>> vvAdj;
    igl::adjacency_list(Fg,vvAdj);
    vector<int> adjacentFaces = vfAdj[cve->vert];

    int newVertIdx = Vg.rows();
    MatrixXd newVg (newVertIdx + 1, 3);
    newVg.block(0,0, newVertIdx, 3)= Vg;
    handledVerticesSet.insert(newVertIdx);
    handledVerticesSet.insert(cve-> vert);

    int idx = 0;
    int plusOneId, minusOneId;

    // this has to be on the new boundary
    while (boundary[idx] != cve->vert && idx < boundary.size()){
        idx++;
    }
    if(boundary[idx] != cve->vert){
        cout<<"ERROR VERTEX NOT FOUND "<<boundary.size()<<" patch "<<cve->patch<<" of "<<boundaryL.size() <<endl;
        cout<<cve->vert<<" the vert "<<endl;
        for (int j = 0; j<boundary.size(); j++){
            cout<<boundary[j]<<" ";
        }
    }


    plusOneId = (idx + 1) % boundary.size();
    minusOneId = (idx -1);
    if(minusOneId<0) minusOneId+= boundary.size();

    int leftFaceId = adjacentFaceToEdge(boundaryL[cve->patch][plusOneId], cve-> vert, -1, vfAdj);
    int rightFaceId = adjacentFaceToEdge(boundaryL[cve->patch][minusOneId], cve-> vert, -1, vfAdj);

    if(leftFaceId ==-1 || rightFaceId ==-1){
        cout<<endl<<boundaryL[cve->patch][plusOneId]<<" something went wrong, we have no neighbor faces. "<<idx<<" "<<boundaryL[cve->patch][minusOneId]<<endl;

        if(cve->seamType>0){
            seam *helper;
            if(cve->seamIdInList>=0){
                helper = seamsList[cve->seamIdInList];
                if(leftFaceId == -1) plusOneId = helper->duplicates[plusOneId];
                if(rightFaceId == -1) minusOneId = helper->duplicates[minusOneId];
            }else{
                helper = seamsList[(-1)*cve->seamIdInList+1];
                if(leftFaceId == -1) plusOneId = helper->duplicates2[plusOneId];
                if(rightFaceId == -1) minusOneId = helper->duplicates2[minusOneId];
            }

        }else{
            minusOneSeam*helper = minusOneSeams[cve->seamIdInList];
            if(leftFaceId == -1) plusOneId = helper->duplicates[plusOneId];
            if(rightFaceId == -1) minusOneId = helper->duplicates[minusOneId];
        }
//       leftFaceId =  adjacentFaceToEdge(boundaryL[cve->patch][plusOneId], cve-> vert, -1, vfAdj );
        leftFaceId = adjacentFaceToEdge(boundaryL[cve->patch][plusOneId], cve-> vert, -1, vfAdj);
        rightFaceId = adjacentFaceToEdge(boundaryL[cve->patch][minusOneId], cve-> vert, -1, vfAdj);
        if(leftFaceId ==-1 || rightFaceId ==-1){cout<<"still no success with "<<plusOneId<<" "<<minusOneId<<endl;   }
        cve->finFlag=true;
        return;
    }

    double thereshold = 1.05;

    Vector3d toLeft = Vg.row(boundaryL[cve->patch][plusOneId]) - Vg.row(cve->vert);
    Vector3d toRight = Vg.row(boundaryL[cve->patch][minusOneId])- Vg.row(cve->vert);

    // if it's a corner all get the new index
    if(cve->startCorner|| cve-> endCorner ){
        cout<<" start or end corner case"<<endl;
        int nextVertOnBoundary;
        getBoundaryNextVert(cve-> startCorner , cve-> seamType, cve-> seamIdInList,
                            boundaryL[cve->patch][minusOneId], boundaryL[cve->patch][plusOneId],nextVertOnBoundary );

        Vector3d cutDirection = Vg.row(nextVertOnBoundary) - Vg.row(cve->vert); //cve-> continuedDirection;
        if (cornerSet.find(cve->vert) != cornerSet.end()&& cve->vert != cve-> cornerInitial){
            cutDirection = (toLeft(0) == cutDirection(0) && toLeft(1) == cutDirection(1)) ? toRight : toLeft;

        }
        VectorXd ws;
        vector<int> fn = vfAdj[cve->vert]; // the face neighbors
        bool tearIsUseful = false;
       for(int faceIdx =0; faceIdx <fn.size(); faceIdx++){
           int adjFace = fn[faceIdx];
           VectorXd faceBary = (Vg_pattern_orig.row(Fg(adjFace, 0))+
                                Vg_pattern_orig.row(Fg(adjFace, 1))+
                                Vg_pattern_orig.row(Fg(adjFace, 2))).transpose();
           faceBary /= 3;

           //90* angle to cut direrction is where we measure the stress
           faceBary(0) -= cutDirection(1);
           faceBary(1) += cutDirection(0); // perp to midvec measure stress
           MatrixXd distB;
           MatrixXd input(1, 3);
           input.row(0)= faceBary;
           // is the face elongated in cut direction compared to its original shape ? (=> rest shape)
           igl::barycentric_coordinates(input, Vg_pattern_orig.row(Fg(adjFace, 0)), Vg_pattern_orig.row(Fg(adjFace, 1)),
                                        Vg_pattern_orig.row(Fg(adjFace, 2)), distB);
           double lenThen = cutDirection.norm();

           VectorXd distNow = (Vg.row(Fg(adjFace, 0)) * distB(0)+
                               Vg.row(Fg(adjFace, 1)) * distB(1) +
                               Vg.row(Fg(adjFace, 2))* distB(2)).transpose();
           distNow -= (Vg.row(Fg(adjFace, 0)) * (1./3)+
                       Vg.row(Fg(adjFace, 1)) *  (1./3) +
                       Vg.row(Fg(adjFace, 2))*  (1./3) ).transpose();
           double lenNow = distNow.norm();
           cout<<"Face: "<<adjFace<<" "<<lenNow<<" dist now and then "<<lenThen<<" ratio is "<<lenNow/lenThen<<endl;
//        cout<<"pos now "<<distNow.transpose() <<" and normed"<<distNow.normalized().transpose() <<endl;//it has to be for any of the adjacent ones, not just this single one

           double checkIfTearIsUsefulThereshold = 1.051;
           if(lenNow/lenThen >checkIfTearIsUsefulThereshold ){
               tearIsUseful= true;
           }
           /*END TEST IF USEFUL*/

       }

        if(! checkIfTearIsUseful(cve-> vert, cutDirection, vvAdj, vfAdj, Vg, lengthsCurr, lengthsOrig, Fg, ws, false)){

            cout<<"stopping now bc tearing is not useful anymore. but we have the new check so let it decide  "<<endl;
//            cve-> finFlag = true;
//            return;
        }
        if(!tearIsUseful){
            cout<<"stopping now because of adj face stress condition  "<<endl;

            cve->finFlag = true;
            return;
        } cout<<" Confirmed tear is useful! "<<endl;

        // set of extra cases
        pair<int, int> valPair = make_pair(cve->seamType, cve ->seamIdInList);
        int seamComp;// get the seam id, m -1 seams are stored negatively
        if(cve->seamType>0){
            if(cve->seamIdInList>=0){
                seamComp = cve ->seamIdInList;
            }else{
                seamComp = (cve ->seamIdInList+1)*(-1);
            }
        }else{
            seamComp = (-1)* (cve ->seamIdInList+1);
        }
        // this is the seam in which we make the decision to release,(!) not from this seam but from the subsequent
        releasedVert[cve->vert]= valPair;
        cout<<"released from seam "<<valPair.first<<" "<<valPair.second<<endl;
        // now we need to find the seam from which we release.
        // each corner is adjacent to two seams. If one is the one we make the decision from, then the other is the one from which it is released
        if(cornerToSeams[cve->cornerInitial][0] == seamComp ){
            // if we allow releasing from two seams this has to be a vec
            if(releasedVertNew.find(cve->vert)!= releasedVertNew.end()){
                releasedVertNew[cve->vert].push_back( cornerToSeams[cve->cornerInitial][1]);
            }else{
                vector<int> temp;
                temp.push_back(cornerToSeams[cve->cornerInitial][1]);
                releasedVertNew[cve->vert] = temp;
            }
        }else if (cornerToSeams[cve->cornerInitial][1] == seamComp ){
            // if we allow releasing from two seams this has to be a vec
            if(releasedVertNew.find(cve->vert)!= releasedVertNew.end()){
                releasedVertNew[cve->vert].push_back( cornerToSeams[cve->cornerInitial][0]);
            }else{
                vector<int> temp;
                temp.push_back(cornerToSeams[cve->cornerInitial][0]);
                releasedVertNew[cve->vert] = temp;
            }


        }else{
            cout<<cornerToSeams[cve->cornerInitial][0]<<" we have a problem, it is not found "<<cornerToSeams[cve->cornerInitial][1]<<" but what we have is "
            <<cve->seamType<<" "<< cve ->seamIdInList<<endl;
            cout<<endl<<endl<<"------------"<<endl<<endl<<"----------"<<endl<<endl;
        }

        newVg.row(newVertIdx) = Vg.row(cve->vert);

        cve->continuedCorner = true;
        cve->finFlag = (cornerSet.find(cve->vert) != cornerSet.end()&& cve->vert != cve-> cornerInitial); //if it is a corner we are done
        int nextVertComp;
        getBoundaryNextVert(cve-> startCorner ,cve-> seamType, cve-> seamIdInList, boundary[minusOneId], boundary[ plusOneId], nextVertComp );
        cve->vert = nextVertComp;
        Vg.resize(Vg.rows()+1, 3);
        Vg= newVg;
        return;

    }
    // if it's a bridge there is no next and we set fin flag
    if(cve-> bridgeFlag){
        // todo I guess this is obsolete
        // this is the final cut, we are done after: should look like  ><
        cout<<"cutting the bridge, but better dont for now , it messes up the patches "<<endl;
       // updateWithNewId(Fg, cve->vert, rightFaceId, newVertIdx);

        cve-> finFlag = true;
        return;
    }

    Vector3d midVec;

    double lenMidVec;
    VectorXd ws, newMidVec, stressMidVec;
    MatrixXd dirs;

    computeMidVecBasedOnStress(stressMidVec, vfAdj, Vg, Vg_pattern_orig, Fg, cve->vert, lenMidVec );

    computeMidVecBasedOnPCA(newMidVec, vvAdj, vfAdj, Vg, lengthsCurr, lengthsOrig, Fg, cve->vert, ws, dirs, lenMidVec );
    midVec = (-1) * stressMidVec; // newMidVec;
    // todo sketchy, for whatever reason it breaks without this. does it do the transposing?
    cout<<" stress and pca midvec "<<midVec.transpose()<<endl;//<<stressMidVec.transpose()

    Vector3d cutDirection = midVec;
    if(!checkIfTearIsUseful(cve-> vert, cutDirection, vvAdj, vfAdj, Vg, lengthsCurr, lengthsOrig, Fg, ws, true)){
        cout<<"stopping now with new condition, but test and go on  "<<endl;
        //        cve-> finFlag = true;
//        return;
    }

    //  we have the midvec direction, but don't know for sure if adding or subtract. Take the direction that has longer distance to the existing boundary to not go backwards.
    // This is a heuristic.
    Vector3d A, Btemp;
    if(cve->levelOne){
        A = Vg.row(boundary[plusOneId]);
        Btemp = Vg.row(cve -> vert);
    }else{
        A = Vg.row(cve -> leftCorner);
        Btemp = Vg.row(cve -> rightCorner);
    }

    //(AB,AM), where M(X,Y) is the query point:
    //https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
    Vector3d C = Vg.row(cve -> vert).transpose()+ midVec;
    auto dist1 = abs((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0)));
    dist1 /= sqrt((Btemp(0)-A(0) )*(Btemp(0)-A(0) ) + (Btemp(1)-A(1))*(Btemp(1)-A(1)));

    C = Vg.row(cve -> vert).transpose() + newMidVec;
    auto dist2 = abs((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0)));
    dist2 /= sqrt((Btemp(0)-A(0) )*(Btemp(0)-A(0) ) + (Btemp(1)-A(1))*(Btemp(1)-A(1)));

    // we have a vector from the pca, but do we need to change the sing?
    if(!cve-> levelOne){
        if(dist2 > dist1 ){
            midVec = newMidVec;
            cout<<"WE CHANGED THE SIGN OF THE MIDVEC!"<<endl;
        }
    }else{
        cout<<"initial!"<<endl;
        // insert the corner

        // if levelone take the side of the other interior vertices
        Vector3d C = Vg.row(cve -> vert).transpose()+ midVec;
        //(AB,AM), where M(X,Y) is the query point:
        //https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
        auto sign1 =(((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0))) > 0);
        // take the side of interior vertices
        for(int j=0; j<vvAdj[cve -> vert].size(); j++) {
            int testInt = vvAdj[cve->vert][j];
            if ((testInt == boundary[minusOneId] || testInt == boundary[plusOneId]))continue;
            // it is an interior
            C = Vg.row(testInt);
            auto sign2 = (((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0))) > 0);
//            cout<<"side of interior vert "<<sign2<<endl;
            if(sign1!= sign2){
                midVec = newMidVec;
                cout<<"WE CHANGED THE SIGN OF THE MIDVEC!"<<endl;
            }
            break;
        }
    }
    // the new midvec is secured now. But to ensure a smooth line we want to include the previous midvec too
    double prior = 0.5;
    if(cve->midVec != Vector3d::Zero()){
        cout<<"using a prior of "<<prior<<" * "<<cve->midVec.transpose()<<endl;
        midVec = (1-prior)* midVec + prior * cve->midVec;
    }
    cve-> midVec = midVec;

    /*-----------------begin trial of snapping to pca direction ------------------*/

    Vector2d mid2 = midVec.transpose().leftCols(2);
    int insertIdx = -1;
    int intersectingFace = -1 ;
    Vector2d intersectPosition;
    VectorXd newPos = Vg.row(cve->vert);
    bool tearIsUseful = false;
    for(int i=0; i < vfAdj[cve->vert].size(); i++){
        int adjFace = vfAdj[cve->vert][i];
        int k=0;
        if(Fg(adjFace, k) == cve->vert) k++;

        int j = k+1;
        if(Fg(adjFace, j) == cve->vert) j++;


        /* TEST IF IT IS ACTUALLY USEFUL*/
        VectorXd distThen = (Vg_pattern_orig.row(Fg(adjFace, 0))+
                             Vg_pattern_orig.row(Fg(adjFace, 1))+
                             Vg_pattern_orig.row(Fg(adjFace, 2))).transpose();
        distThen /= 3;
        distThen(0) -= mid2(1);
        distThen(1) += mid2(0); // perp to midvec measure stress
        MatrixXd distB;
        MatrixXd input(1, 3);
        input.row(0)= distThen;
        igl::barycentric_coordinates(input, Vg_pattern_orig.row(Fg(adjFace, 0)), Vg_pattern_orig.row(Fg(adjFace, 1)),
                                     Vg_pattern_orig.row(Fg(adjFace, 2)), distB);
        double lenThen = mid2.norm();
        VectorXd distNow = (Vg.row(Fg(adjFace, 0)) * distB(0)+
                            Vg.row(Fg(adjFace, 1)) * distB(1) +
                            Vg.row(Fg(adjFace, 2))* distB(2)).transpose();
        distNow -= (Vg.row(Fg(adjFace, 0)) * (1./3)+
                    Vg.row(Fg(adjFace, 1)) *  (1./3) +
                    Vg.row(Fg(adjFace, 2))*  (1./3) ).transpose();
        double lenNow = distNow.norm();
        cout<<"Face: "<<adjFace<<" "<<lenNow<<" dist now and then "<<lenThen<<" ratio is "<<lenNow/lenThen<<endl;
//        cout<<"pos now "<<distNow.transpose() <<" and normed"<<distNow.normalized().transpose() <<endl;//it has to be for any of the adjacent ones, not just this single one

        double checkIfTearIsUsefulThereshold = 1.051;
        if(lenNow/lenThen >checkIfTearIsUsefulThereshold ){
            tearIsUseful= true;
        }
        /*END TEST IF USEFUL*/

        if(raySegmentIntersection(Vg.row(cve->vert).leftCols(2), Vg.row(Fg(adjFace, k)).leftCols(2),
                                                                            Vg.row(Fg(adjFace, j)).leftCols(2), mid2,
                                                                            lengthsOrig(adjFace, 0)*100, intersectPosition)){

            intersectingFace = adjFace;
            newPos(0)= intersectPosition(0);
            newPos(1)= intersectPosition(1);
            newPos(2) = Vg(Fg(adjFace, k), 2);
            double dist1 = (Vg.row(Fg(adjFace, k))- newPos).squaredNorm();
            double dist2 = (Vg.row(Fg(adjFace, j))- newPos).squaredNorm();

            insertIdx = (dist1<dist2) ? Fg(adjFace, k): Fg(adjFace, j);
            cout<< insertIdx <<" the insert idx is found "<<endl;
//            break; // it intersects only one
        }

    }
    if(!tearIsUseful){
            cve->finFlag = true;
            return;
    }

    // todo no face found, maybe one is exactly on the intersection? might be! in that case we have to search again...
    if(intersectingFace == -1 || insertIdx ==-1){
        cout<<"ERRRROR WE FOUND NO FACE!!! "<<endl;
        cout<< mid2.transpose()<<" "<<cve->vert<<endl;

        for(int i=0; i < vfAdj[cve->vert].size(); i++){
            int adjFace = vfAdj[cve->vert][i];
            int k=0;
            if(Fg(adjFace, k) == cve->vert) k++;

            int j = k+1;
            if(Fg(adjFace, j) == cve->vert) j++;
            cout<<adjFace<<" adj face;"<<endl;
            cout<<"to "<<Fg(adjFace, j)<<": "<< (Vg.row(Fg(adjFace, j)).leftCols(2)-Vg.row(cve->vert).leftCols(2)) <<endl;
            cout<<"to "<<Fg(adjFace, k)<<": "<< (Vg.row(Fg(adjFace, k)).leftCols(2)-Vg.row(cve->vert).leftCols(2)) <<endl;


            if(raySegmentIntersection(Vg.row(cve->vert).leftCols(2), Vg.row(Fg(adjFace, k)).leftCols(2),
                                      Vg.row(Fg(adjFace, j)).leftCols(2), mid2,
                                      lengthsOrig(adjFace, i)*10, intersectPosition)){
                cout<<Fg.row(adjFace)<<", face id "<<adjFace<<" k="<<k<<" j="<<j<<endl;
                intersectingFace = adjFace;
                newPos(0)= intersectPosition(0);
                newPos(1)= intersectPosition(1);
                newPos(2) = Vg(Fg(adjFace, k), 2);
                double dist1 = (Vg.row(Fg(adjFace, k))- newPos).squaredNorm();
                double dist2 = (Vg.row(Fg(adjFace, j))- newPos).squaredNorm();

                insertIdx = (dist1<dist2) ? Fg(adjFace, k): Fg(adjFace, j);
                cout<< insertIdx <<" the insert idx is found "<<endl;
                break; // it intersects only one
            }

        }
    }
    VectorXd signs, signsAfter;
    // preparation, calc orientation to check if flipped after
    computeOrientationOfFaces(signs, vfAdj[insertIdx], Fg, Vg);

    MatrixXd insertIdxInBary;
    MatrixXd input(1, 3);
    input.row(0)= newPos;
    igl::barycentric_coordinates(input, Vg.row(Fg(intersectingFace, 0)), Vg.row(Fg(intersectingFace, 1)),
                                 Vg.row(Fg(intersectingFace, 2)), insertIdxInBary);

    newVg.row(insertIdx) = newPos;
    cout<<newPos.transpose()<<" new pos"<<endl;

    VectorXd updatedRestShapeVertPos = insertIdxInBary(0) * Vg_pattern_orig.row(Fg(intersectingFace, 0)) ;
    updatedRestShapeVertPos += insertIdxInBary(1) * Vg_pattern_orig.row(Fg(intersectingFace, 1)) ;
    updatedRestShapeVertPos += insertIdxInBary(2) * Vg_pattern_orig.row(Fg(intersectingFace, 2) );
    //  upadte all original edge lengths -> in main
    cout<<  Vg_pattern_orig.row(insertIdx)<<" before restshape "<<endl;
    Vg_pattern_orig.row(insertIdx) = updatedRestShapeVertPos;
    cout<<  Vg_pattern_orig.row(insertIdx)<<" after restshape"<<endl;

    computeOrientationOfFaces(signsAfter, vfAdj[insertIdx], Fg, Vg);
    if(signsAfter != signs){
        for(int i=0; i<vfAdj[insertIdx].size(); i++){
            if(signs(i) != signsAfter(i) ) {
                cout <<" ERROR WE FLIPPED SOMETING in face "<< vfAdj[insertIdx][i]<<endl;
                cout << endl;
            }// mess with the original edge length for computation . this is not clean II thing
            // todo
        }
    }

    /*-----------------end trial of snapping to pca direction ------------------*/

    MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(Fg, TT, TTi);
    std::pair<int, int>  faces;
    adjacentFacesToEdge(cve->vert, insertIdx, vfAdj, faces );
    cout<<insertIdx<<" the insert idx "<<faces.first<<" and second "<<faces.second<<endl;
    int helperWhich = -1;
    if(Fg(faces.first, 0) == cve-> vert){
        helperWhich = 0;
    }else if (Fg(faces.first, 1) == cve-> vert){
        helperWhich = 1;
    }else{
        helperWhich = 2;
    }

    igl::HalfEdgeIterator <MatrixXi, MatrixXi, MatrixXi>hei (Fg, TT, TTi, faces.first, helperWhich, false);
    igl::HalfEdgeIterator <MatrixXi, MatrixXi, MatrixXi>hei2 (Fg, TT, TTi, faces.second, helperWhich, false);
    int lastF = hei.Fi();
    int lastV= hei.Vi();
    bool nextEnd = hei.NextFE();
    if( hei.Fi()!= faces.second){
        // then we skip the first and directly go from there
        // we need to change the position
        updateWithNewId(Fg, cve->vert, lastF, newVertIdx);
    }
    if(0!= nextEnd) {
        updateWithNewId(Fg, cve->vert, hei.Fi(), newVertIdx);
        while (hei.NextFE() != 0) {
            updateWithNewId(Fg, cve->vert, hei.Fi(), newVertIdx);

        }
    }

    newVg.row(newVertIdx) = Vg.row(cve->vert);

    if(releasedVertNew.find(cve->vert) != releasedVertNew.end() ){
        // if that one was released and we cut here (only if l shapes are allowed), we have to release the duplicate too!
        releasedVertNew[newVertIdx] = releasedVertNew[cve-> vert];
        pair<int, int> valPair = make_pair(cve->seamType, cve ->seamIdInList);
        releasedVert[newVertIdx] = valPair;

    }
    cveStartPositionsSet[newVertIdx] =cve;
    cout<<insertIdx<<" the inserted index"<<endl;// TODO CASE IT IS NOT RIGHT NEITHER LEFT BUT ACTUALLY ON!!

    Eigen::MatrixXi B;
    bool isManifold = igl::is_vertex_manifold( Fg, B);
    if(B(insertIdx, 0)!= 1){
        cout<<" ATTENTION WE ARE CREATING A NON MANIFOLD MESH! CUT THROUGH IMMEDIATELY"<<endl;
        int newnewVertIdx = newVertIdx+1;
        if(Fg(faces.first, 0) == insertIdx){
            helperWhich = 0;
        }else if (Fg(faces.first, 1) == insertIdx){
            helperWhich = 1;
        }else{
            helperWhich = 2;
        }
        igl::HalfEdgeIterator <MatrixXi, MatrixXi, MatrixXi>hei (Fg, TT, TTi, faces.first, helperWhich, false);

        int lastF = hei.Fi();
        int lastV= hei.Vi();
        bool nextEnd = hei.NextFE();
        if( hei.Fi()!= faces.second){
            // then we skip the first and directly go from there
            // we need to change the position
            cout<<"updating face last F  "<<lastF<<endl;
            updateWithNewId(Fg, insertIdx, lastF, newnewVertIdx);
        }
        if(0!= nextEnd) {
            cout<<"updating face  0 neq next end "<<hei.Fi()<<endl;
            updateWithNewId(Fg, insertIdx, hei.Fi(), newnewVertIdx);
            while (hei.NextFE() != 0) {

                cout<<"updating face "<<hei.Fi()<<endl;
                updateWithNewId(Fg, insertIdx, hei.Fi(), newnewVertIdx);

            }
        }

        MatrixXd newnewVg(newVg.rows()+1, 3);
        newnewVg.block(0,0, newVg.rows(), 3)= newVg;
        newnewVg.row(newnewVertIdx)= newVg.row(insertIdx);//+ (eps * toRight).transpose();

        handledVerticesSet.insert(newnewVertIdx);
        cutThroughCornerVertices.insert(newnewVertIdx);
        cout<<"inserted "<<newnewVertIdx<<" and "<<insertIdx<<endl;
        cutThroughCornerVertices.insert(insertIdx);
//        cout<<newnewVg.row(insertIdx)<<" old insert idx"<<endl;
//        cout<<newnewVg.row(newnewVertIdx)<<" new insert idx duplicate "<<endl;

        newVg.resize(newnewVg.rows(), 3);
        newVg = newnewVg;
        cout<<" fin operation"<<endl;

        cve->finFlag= true;

    }

    Vg.resize(newVg.rows(), 3);
    Vg = newVg;

    if(toPattern_boundaryVerticesSet.find(insertIdx)!= toPattern_boundaryVerticesSet.end()){
            cout<<"we are nearly done here, it's cut through! "<<endl;
            cve-> vert = insertIdx;
            cve -> bridgeFlag = true;
    }

    // only the first level i.e. boundary duplicate has to be projected, hence only this one is to be added
    if(cve->levelOne) {
        addToDulicateList(cve,seamsList, minusOneSeams, newVertIdx, false );
//        cout<<minusOneSeams[cve->seamIdInList]->duplicates[cve->vert]<<" check if new is in duplicate list "<<endl ;
        toPattern_boundaryVerticesSet.insert(newVertIdx);// the duplicate is also on the boundary, hence insert it
        cve->leftdirection = toLeft;
        cve->rightdirection = toRight;
        cve-> leftCorner =  cve-> vert;
        cve-> rightCorner = newVertIdx;
        cout<<" inserted"<<endl;
    }
    addToDulicateList(cve,seamsList, minusOneSeams, newVertIdx, true );


    cve-> vert = insertIdx;
    cve->levelOne = false;

}
void splitCounterPart(vector<cutVertEntry*>& cutPositions, int idxOfCVE,  cutVertEntry*& cve, MatrixXd& Vg, MatrixXi& Fg, vector<vector<int> >& vfAdj,
                      std::vector<std::vector<int> >& boundaryL,  vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams,
                      map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  MatrixXd& lengthsCurr,
                      MatrixXi& Fg_pattern, set<int> & cornerSet,  set<int>& handledVerticesSet, const bool & LShapeAllowed,MatrixXd& Vg_pattern_orig){
    //idx of cve has higher stress, hence we have searched before already, just need to increment
    //if not found we have handled it already
    if(cve->seamType == -1 ) {
        // there is no counterpart
        cout<<"NO counterpart"<<endl; return ;

    }
    int counterID  = (cve->seamIdInList+1)*(-1);
    if(!cve->startCorner && !cve->endCorner){
        // then there is only one other of this seam -> search it
        int idx = -1;
        for(int i=idxOfCVE; i<cutPositions.size(); i++){
            if(cutPositions[i]->seamType == 1 &&cutPositions[i]->seamIdInList == counterID){
                idx= i;
            }
        }
        if(idx == -1) return;
        splitVertexFromCVE( cutPositions[idx], Vg, Fg, vfAdj,
                            boundaryL, seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,
                            lengthsCurr, cornerSet, handledVerticesSet, LShapeAllowed, Vg_pattern_orig);
        return;
    }

    int lookFor = -1;
    seam* currSeam;
    if(counterID<0){
        // then we are > 0
        currSeam = seamsList[cve->seamIdInList];

    }else{
        // we are < 0
        currSeam = seamsList[(cve->seamIdInList+1)*(-1)];
    }
    if(cve->startCorner){
        // we need to find another startcorner
        if(cve->vert == currSeam-> getStart1()){
            lookFor = currSeam-> getStart2();
        }else if (cve->vert == currSeam-> getStart2()) {
            lookFor = currSeam->getStart1();
        }else{
            cout<<"we have a problem, should find cve as corner but seems like its not "<<cve-> vert<<endl; return;
        }
    }
      // find the mapped vertex of this one
    else {
        if (cve->vert == currSeam->getEndCornerIds().first) {
            lookFor = currSeam->getEndCornerIds().second;
        } else if (cve->vert == currSeam->getEndCornerIds().second) {
            lookFor = currSeam->getEndCornerIds().first;
        } else {
            cout << "we have a problem, should find cve as corner but seems like its not " << cve->vert << endl;return;
        }
    }

    int idx = -1;
    for(int i = idxOfCVE; i<cutPositions.size(); i++){
        if(cutPositions[i]->seamType == 1 && cutPositions[i]->vert == lookFor){
            idx= i;
        }
    }
    if(idx == -1) {
        cout<<"not found, problem with cutting the corresponding, did not find "<<lookFor<<endl;
        return;
    }
    splitVertexFromCVE( cutPositions[idx], Vg, Fg, vfAdj,
                        boundaryL, seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet, lengthsCurr,
                        cornerSet, handledVerticesSet, LShapeAllowed, Vg_pattern_orig);
    return;
}
/*
|
.___b

.___a
|
 if we open the cut horizontally at a, we want to open at b as well since the two together form a seam (though the common seam has it's cut in the middle)
 */
 int openParallelPosition(int& cornerInitial, int& seamType, vector<seam*>& seamsList, vector <cutVertEntry*>& cutPositions){
     // if a paralell cut position exists we find and open it . But attention, it is not guaranteed the position exists. But we don't enforce it (might be worth a thought tough)
    if(seamType < 0){
        return -1; // no parallel to open
    }
    seam* currSeam = seamsList[seamType];
    // figure out if the initial was pos or neg side, the one we are searching is the other side
    int searchedVert;
    if(currSeam->getStart1() == cornerInitial){
        searchedVert = currSeam->getStart2();
    }else if(currSeam->getStart2() == cornerInitial) {
        searchedVert = currSeam->getStart1();
    }else if(currSeam->getEndCornerIds().first == cornerInitial) {
        searchedVert = currSeam-> getEndCornerIds().second;
    }else if(currSeam->getEndCornerIds().second == cornerInitial) {
        searchedVert = currSeam-> getEndCornerIds().first;
    }else{
//        cout<<" partner not found, we have a huge problem in opening the parallel positions "<<seamType<<endl;
        return -1;
    }

    int size =  cornerToSeams[searchedVert].size();
    if( size != 2) {
        cout<<searchedVert<<" the size is not 2, it's "<<size<<endl;
    }
    int otherSeamId;
    if(cornerToSeams[searchedVert][0] == seamType){

        otherSeamId = cornerToSeams[searchedVert][1];

    }else{
        otherSeamId = cornerToSeams[searchedVert][0];
    }

    // identify the index of the paralell cut postion by comparing the vertex with the one we search
    for(int i=0; i<cutPositions.size(); i++){
        if(cutPositions[i]->vert == searchedVert || cutPositions[i]->cornerInitial == searchedVert){
            return i;

        }
    }
    return -1;

}
void findCorrespondingCounterCutPosition(vector<cutVertEntry*>& cutPositions, int idxOfCVE, cutVertEntry*& cve, MatrixXd& Vg, MatrixXi& Fg, vector<vector<int> >& vfAdj,
                      std::vector<std::vector<int> >& boundaryL,  vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams,
                      map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  MatrixXd& lengthsCurr,
                      MatrixXi& Fg_pattern, set<int> & cornerSet,  set<int>& handledVerticesSet){
    //idx of cve has higher stress, hence we have searched before already, just need to increment
    //if not found we have handled it already
    if(cve->seamType == -1 ) {
        // there is no counterpart
        cve->stressWithCounter = cve->stress;
        return ;

    }
    if(cve-> counterPartIdx != -1){
//        cout<<"found counter already, it's "<< cve-> vert <<" and "<< cutPositions[cve->counterPartIdx]->vert<<endl;
        return;
    }

    int counterID  = (cve->seamIdInList+1)*(-1);
    if(!cve->startCorner && !cve->endCorner){
        // then there is only one other of this seam -> search it
        int idx = -1;
        for(int i = 0; i < cutPositions.size(); i++){
            if(cutPositions[i]-> seamType == 1 && cutPositions[i]-> seamIdInList == counterID){
                idx= i;
            }
        }
        if(idx == -1) {
            cout<<"Problem: found a middle vertex but no counterpart to it "<<cve->seamIdInList<<" with  vertex "<<cve->vert<<endl;
            return;
        }
        cve-> counterPartIdx = idx;
        cutPositions[idx]-> counterPartIdx = idxOfCVE;

        return;
    }

    int lookFor = -1;
    seam* currSeam;
    currSeam = (counterID<0) ? seamsList[cve->seamIdInList] : seamsList[(cve->seamIdInList+1)*(-1)];

    if(cve->startCorner){
        // we need to find another startcorner
        if(cve->vert == currSeam-> getStart1()){
            lookFor = currSeam-> getStart2();
        }else if (cve->vert == currSeam-> getStart2()) {
            lookFor = currSeam->getStart1();
        }else{
            cout<<"we have a problem, should find cve as corner but seems like its not "<<cve-> vert<<endl; return;
        }
    }
        // find the mapped vertex of this one
    else {
        if (cve->vert == currSeam->getEndCornerIds().first) {
            lookFor = currSeam->getEndCornerIds().second;
        } else if (cve->vert == currSeam->getEndCornerIds().second) {
            lookFor = currSeam->getEndCornerIds().first;
        } else {
            cout << "we have a problem, should find cve as corner but seems like its not " << cve->vert << endl;return;
        }
    }

    int idx = -1;
    for(int i = 0; i<cutPositions.size(); i++){
        if(cutPositions[i]->seamType == 1 && cutPositions[i]->vert == lookFor && cutPositions[i]-> seamIdInList == counterID){
            idx= i;
        }
    }
    if(idx == -1) {
        cout<<cve->vert<<" not found, problem with cutting the corresponding, did not find "<<lookFor<<endl;
        return;
    }
    cve-> counterPartIdx = idx;
    cutPositions[idx]-> counterPartIdx = idxOfCVE;

    double stressSum = cve->stress + cutPositions[idx]->stress;
    cve-> stressWithCounter = stressSum;
    cutPositions[idx]->stressWithCounter = stressSum;

    return;
}

int addoncount=0;
void addVarToModel (int vert, int prevVert, int nextVert, vector<vector<int>> & vfAdj, bool isConstrained, int& varCount, GRBVar* & cutVar,
                      MatrixXi& Fg_pattern,MatrixXd& lengthsOrig, MatrixXd& lengthsCurr, map <int, cutVertEntry*> & mapVarIdToVertId,
                      int seamType, int seamIdInList, double tailor_lazyness, bool corner, map<pair<int,int>, int>& mapVertAndSeamToVar, int counterpart, GRBModel& model ){

    double w_init = 0;
    int count = 0;
    if(nextVert != -1){
        count++;
        int faceIdx = adjacentFaceToEdge(vert, nextVert, -1, vfAdj );
        int whichEdge = findWhichEdgeOfFace(faceIdx, vert, nextVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        w_init += lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);
    }
    if(prevVert != -1){
        count++;
        int faceIdx = adjacentFaceToEdge(vert, prevVert, -1, vfAdj );
        int whichEdge = findWhichEdgeOfFace(faceIdx, vert, prevVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        w_init += lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);
    }
//TODO SPECIAL CASE IF PREV -1
//
    if(count>1) w_init/=count;

    if(isConstrained){
        try {
            cutVar[varCount].set(GRB_DoubleAttr_Obj, 0);
        }catch(GRBException e){
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
            cout<<varCount<<endl;
        }
    }else{

        int tailorLazyFactor = 1; if(corner) tailorLazyFactor*= tailor_lazyness;
        try{
            cutVar[varCount].set(GRB_DoubleAttr_Obj, tailorLazyFactor * 1/w_init);

        }catch(GRBException e){
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
            cout<<varCount<<endl;

        }
    }

    cutVertEntry* cve = new cutVertEntry(vert, seamType, seamIdInList);
    mapVarIdToVertId[varCount]= cve;

    if(seamType ==-1 ){
        varCount++;
        return ;
    }
    int seamIdToCompare = seamIdInList;
    if(seamIdInList<0){
        seamIdToCompare = (seamIdInList+1 )*(-1);
    }
    int thisVar = varCount;
    if(mapVertAndSeamToVar.find(make_pair(counterpart, seamIdToCompare))!= mapVertAndSeamToVar.end()){
        // if the counterpart exists already, we use another variable for the XNOR constraint
        //https://yetanothermathprogrammingconsultant.blogspot.com/2022/06/xnor-as-linear-inequalities.html
        int otherVar = mapVertAndSeamToVar[make_pair(counterpart, seamIdToCompare)];
        // the helper, appear in the objective with higher weight since we need consistency over perfect cuts
        cutVar[varCount+1].set(GRB_DoubleAttr_Obj, 10);

        model.addConstr(1- cutVar[varCount+1] <=  cutVar[ varCount] + cutVar[ otherVar] );
        model.addConstr(1- cutVar[varCount+1] >=  cutVar[ otherVar] - cutVar[ varCount] );
        model.addConstr(1- cutVar[varCount+1] >=  cutVar[ varCount] - cutVar[otherVar ] );
        model.addConstr(1- cutVar[varCount+1] <= 2 - cutVar[ varCount] - cutVar[ otherVar] );
        varCount++;
        addoncount++;

    }
    else{
        mapVertAndSeamToVar[make_pair(vert, seamIdToCompare)] = thisVar;
    }
    varCount++;


}

double computeSeamLength(pair<int, int>& seamId, vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams, std::vector<std::vector<int> >& boundaryL, MatrixXd& Vg ){
    double len=0;
    int type = seamId.first;
    int seamIdx = seamId.second;
    if (type<0){
        minusOneSeam* currSeam = minusOneSeams[seamIdx];
        int vertStart = currSeam->getStartVert();
        int vertEnd = currSeam->getEndVert();
        vector<int> boundary = boundaryL[currSeam->getPatch()];
        int startIdx;
        findIndxInBoundaryloop( boundary, vertStart , startIdx);

        while(boundary[startIdx] != vertEnd){
            int next = (startIdx+1)% boundary.size();

            len += (Vg.row(startIdx)-Vg.row(next)).norm();
            startIdx = next;
        }

    }else{
        if(seamIdx >=0){
            seam* currSeam = seamsList[seamIdx];
            int vertStart = currSeam->getStart1();
            int vertEnd = currSeam->getEndCornerIds().first;
            vector<int> boundary = boundaryL[currSeam->getPatch1()];
            int startIdx;
            findIndxInBoundaryloop( boundary, vertStart , startIdx);

            while(boundary[startIdx] != vertEnd){
                int next = (startIdx+1)% boundary.size();

                len += (Vg.row(startIdx)-Vg.row(next)).norm();
                startIdx = next;
            }

        }else{
            seam* currSeam = seamsList[(seamIdx+1)*(-1)];
            int vertStart = currSeam->getStart2();
            int vertEnd = currSeam->getEndCornerIds().second;
            vector<int> boundary = boundaryL[currSeam->getPatch2()];
            int startIdx;
            findIndxInBoundaryloop( boundary, vertStart , startIdx);
            int count=0;
            while(boundary[startIdx] != vertEnd){
                count++;
                int next = (startIdx-1);
                if(next<0) next+= boundary.size();
                if(currSeam->inverted) next = (startIdx + count) % boundary.size();

                len += (Vg.row(startIdx)-Vg.row(next)).norm();
                startIdx = next;
            }

        }

    }


    return len;
}

void getPrevAndNextVertAndStress(int seamType, int seamId, int vert, int & prevVert, int & nextVert, double & prevStress, double & nextStress,
                                 vector<seam*>& seamsList, const vector<minusOneSeam*>& minusOneSeams, std::vector<std::vector<int> >& boundaryL, MatrixXi& Fg_pattern,
                                 MatrixXd& lengthsOrig, MatrixXd& lengthsCurr, vector<vector<int>> & vfAdj, bool inverted ){
    int patch = -1;
    for(int i=0; i<boundaryL.size(); i++){
        for(int j=0; j<boundaryL[i].size(); j++){
            if(boundaryL[i][j]==vert){
                patch = i;
                break;
            }
        }
        if(patch != -1)break;
    }
    if(patch == -1) cout<< "no patch found ERROR"<<endl;

    if(seamType == -1 ){

        nextVert = minusOneSeams[seamId]->getNextVert(vert, boundaryL[patch]);
        prevVert= minusOneSeams[seamId]->getPrevVert(vert, boundaryL[patch]);
        int faceIdx = adjacentFaceToEdge(vert, nextVert, -1, vfAdj );
        int whichEdge = findWhichEdgeOfFace(faceIdx, vert, nextVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        nextStress = lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);

        faceIdx = adjacentFaceToEdge(vert, prevVert, -1, vfAdj );
        whichEdge = findWhichEdgeOfFace(faceIdx, vert, prevVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        prevStress = lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);
//

    }else if(seamId>=0|| inverted){
        if(inverted){cout<<"inverted issue"<<endl;
        seamId = (seamId+1)*(-1);}
//        int patch  = seamsList[seamId]->getUpdatedPatch1();
        int patchsize =  boundaryL[patch].size();
        nextVert = seamsList[seamId]->getNextVert1(vert, boundaryL[patch]);
        prevVert = seamsList[seamId]->getPrevVert1(vert, boundaryL[patch]);

        int faceIdx = adjacentFaceToEdge(vert, nextVert, -1, vfAdj );
        int whichEdge = findWhichEdgeOfFace(faceIdx, vert, nextVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        nextStress = lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);
        faceIdx = adjacentFaceToEdge(vert, prevVert, -1, vfAdj );
        whichEdge = findWhichEdgeOfFace(faceIdx, vert, prevVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        prevStress = lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);

    }else{
//        int patch = seamsList[(seamId+1)*(-1)] -> getUpdatedPatch2();
        int patchsize =  boundaryL[patch].size();
        int count=0;
        nextVert = seamsList[(seamId +1)*(-1)]->getPrevVert2(vert, boundaryL[patch]);
        prevVert = seamsList[(seamId +1)*(-1)]->getNextVert2(vert, boundaryL[patch]);
        int faceIdx = adjacentFaceToEdge(vert, nextVert, -1, vfAdj );
        int whichEdge = findWhichEdgeOfFace(faceIdx, vert, nextVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        nextStress = lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);
        faceIdx = adjacentFaceToEdge(vert, prevVert, -1, vfAdj );
        whichEdge = findWhichEdgeOfFace(faceIdx, vert, prevVert, Fg_pattern);
        // in case we stretch this is lower 1, the greater it is the smaller it gets, hence 1/w
        prevStress = lengthsCurr(faceIdx, whichEdge) / lengthsOrig(faceIdx, whichEdge);

    }
}

void setLP(std::vector<std::vector<int> >& boundaryL , vector<vector<int>> & vfAdj, MatrixXi& Fg_pattern,
        MatrixXd& lengthsOrig, MatrixXd& lengthsCurr,const std::vector<std::vector<std::pair<int, int>>>& cornersPerBoundary,
        map<int, vector<pair<int, int>>>& seamIdPerCorner, vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,
        double tailor_lazyness, double minConstrained, vector <cutVertEntry*>& cutPositions, VectorXd& cornerVert,
        MatrixXd& currPattern, const bool & LShapeAllowed){

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();
    // Create an empty model
    GRBModel model = GRBModel(env);
    int varCount =0;
    int numVar=0; // whole boundary and all corners duplicate
    for(int i=0; i<cornersPerBoundary.size(); i++){
        numVar += cornersPerBoundary[i].size();
        numVar += boundaryL[i].size();
    }
    // with the constraints we have for each counterpart another variable
    if(numVar%2 != 0 )cout<<"----------------------------------not divisible by 2, the number of constraints is wrong------------------"<<endl<<endl<<endl;
    numVar += (numVar/2);// 574 without mapping -> 287 added ones ? we don't need them all bc there are no common ones for th -1s
    map<int, int> trackCornerIds; // keeps track of the var Id per corner - to understand if it exists already or not
    map <int, cutVertEntry*>  mapVarIdToVertId;
    map <pair<int, int>, int> mapVertAndSeamToVar;

    // a map to track the gurobi id of a corner.  whenever we set a corner we add the corner id, so we know once we set the second and connect them
    GRBVar* cutVar = model.addVars(numVar, GRB_BINARY);

    for(int i = 0; i < cornersPerBoundary.size(); i++) {
        int boundSize = boundaryL[i].size();
        int startVarPatch = varCount;
        for (int j = 0; j < cornersPerBoundary[i].size(); j++) {

            // first is the absolute index, second the index wrt the boundary loop
            auto cornerPair = cornersPerBoundary[i][j];
            if(seamIdPerCorner.find(cornerPair.first) == seamIdPerCorner.end()) continue;

            vector<pair<int, int>> seamId = seamIdPerCorner[cornerPair.first];// all seams that start at this corner, this can be max 2
            cout<<seamId.size()<<" corner "<<cornerPair.first<<" size of seams per corner here"<<endl;
            if(seamId.size()>2) cout<<" something is veryy odd!! we have more than two seams for a corner. impossible."<<endl;

            for(int si = 0; si < seamId.size(); si++) {
                cout<<"corner "<<cornerPair.first<<" info: "<<seamId[si].first<<" "<<seamId[si].second<<endl;

                GRBLinExpr innerSumConstr = 0; // interior sum
                GRBLinExpr lSumConstr = 0; // left sum
                GRBLinExpr rSumConstr = 0; // right sum
                int startVarOfThis = varCount;
                int count = 0;
                double distStartToEnd = computeSeamLength(seamId[si], seamsList, minusOneSeams, boundaryL, currPattern);
                //todo
                double widthThereshold = 50;
                // check for direction
                // first if it is a seam or a -1 seam
                // second gives the direction, if less than 0 take i -1 in negative direction
                if ((seamId[si].first >= 0 && seamId[si].second >= 0) || seamId[si].first < 0) {
                    int idx;
                    int vert;// first corner
                    int end;// last corner
                    int length;

                    pair<int, int> startAndPatchOther;
                    int idxOther;
                    int vertOther;
                    int counterPart = -1;
                    int boundSizeOther;


                    if (seamId[si].first >= 0) {
                        seam *seam = seamsList[seamId[si].second];
                        //gives an index and it is updated -> is it right? to check
                        auto startAndPatch = seam->getStartAndPatch1();
                        startAndPatchOther = seam-> getStartAndPatch2ForCorres();
                        if (startAndPatch.second != i)
                            cout << " now in th e+1 seams the patch does not match where we are in the loop "<< startAndPatch.second << endl;
                        int startVal = seam->getStart1();
                        findIndxInBoundaryloop(boundaryL[i], startVal, idx);
                        startVal = seam->getStart2();
                        findIndxInBoundaryloop(boundaryL[startAndPatchOther.second], startVal, idxOther);
//                        idx = startAndPatch.first;
//                        idxOther = startAndPatchOther.first;
                        vert = boundaryL[startAndPatch.second][idx];
                        vertOther =  boundaryL[startAndPatchOther.second][idxOther];
                        boundSizeOther = boundaryL[startAndPatchOther.second].size();
                        end = seam->getEndCornerIds().first;
                        length = seam->seamLength();
//                        safe falsch!!

                    } else if (seamId[si].first < 0) {
                        minusOneSeam *currSeam = minusOneSeams[seamId[si].second];
                        int patch = currSeam->getPatch();

                        if (patch != i)
                            cout << " now in th e-1 seams the patch does not match where we are in the loop , would be patch "<<patch << endl;
                        int startVal = currSeam->getStartVert();
                        findIndxInBoundaryloop(boundaryL[patch], startVal, idx);
//                        idx = currSeam->getStartIdx();
                        vert = boundaryL[i][idx];
                        end = currSeam->getEndVert();
                        length = currSeam->getLength();
                    }


                    // check if the corners exist already. If so then connect with corner, else add the indices
                    // todo for L cutting allowance
                    if(trackCornerIds.find(vert) != trackCornerIds.end()){
                        model.addConstr(cutVar[ trackCornerIds[vert]] + cutVar[varCount] <= 1);
                    }else{
                        trackCornerIds[vert] = varCount;
                    }

                    int prevVert = -1;
                    while (vert != end) {
                        // might need a map for patch and vert id to seam adn first or Second
                        double relId = ((double) count) / length;
                        bool isConstrained = true;
                        bool corner = (count == 0);


                        if (relId == 0 ||
                        (relId > minConstrained && relId < (1 - minConstrained)  && distStartToEnd > widthThereshold)
                        ) {
                            isConstrained = false;
                        }else{
                            model.addConstr(cutVar[varCount] == 0);
                        }
                        int currVar = varCount;

                        addVarToModel(vert, prevVert, boundaryL[i][(idx + 1) % boundSize], vfAdj, isConstrained, varCount, cutVar,
                                      Fg_pattern, lengthsOrig, lengthsCurr, mapVarIdToVertId, seamId[si].first, seamId[si].second,
                                      tailor_lazyness,corner,mapVertAndSeamToVar, vertOther, model );

                        if (!isConstrained) {
                            lSumConstr += cutVar[currVar];
                            if (relId != 0) {
                                rSumConstr += cutVar[currVar];
                                innerSumConstr += cutVar[currVar];
                            }
                        }
                        idx = (idx + 1) % boundSize;
                        prevVert = vert;
                        vert = boundaryL[i][idx];
                        count++;
                        if(seamId[si].first >= 0){
                            idxOther =  (startAndPatchOther.first - ( count)) % boundSizeOther;
                            if (idxOther < 0) idxOther += boundSizeOther;
                            if (seamsList[seamId[si].second]->inverted) idxOther = (startAndPatchOther.first + (count)) % boundSize;
                            vertOther = boundaryL[startAndPatchOther.second][idxOther];
                        }
                    }// handle the last, add it to the right sum

                    //  make sure duplicate corner is chosen only once, and if we found the last consider its duplicate with start of this patch
                    if(trackCornerIds.find(end) != trackCornerIds.end()){
                        model.addConstr(cutVar[ trackCornerIds[end]] + cutVar[varCount] <= 1);
                    }else{
                        trackCornerIds[end] = varCount;
                    }
                    rSumConstr += cutVar[varCount];
                    addVarToModel(vert, prevVert , -1, vfAdj, false, varCount, cutVar, Fg_pattern,
                                  lengthsOrig, lengthsCurr, mapVarIdToVertId, seamId[si].first, seamId[si].second, tailor_lazyness, true, mapVertAndSeamToVar, vertOther, model );
                } else {
                    // iterate in inverse direction
                    seam *seam = seamsList[(seamId[si].second ) * -1 -1];
//                    cout<<" seam id "<<(seamId[si].second ) * -1 -1<<endl;
                    auto startAndPatch = seam -> getStartAndPatch2ForCorres();
                    auto startAndPatchOther = seam -> getStartAndPatch1();
                    int idx = startAndPatch.first;

                    if (startAndPatch.second != i) cout << " something is wrong, it should be in patch " << i << " but the seam is from patch  " << startAndPatch.second << endl;
                    int vert = boundaryL[startAndPatch.second][idx];
                    int otherVert = boundaryL[startAndPatchOther.second][startAndPatchOther.first];
                    int end = seam->getEndCornerIds().second;

                    int nextidx = (startAndPatch.first - (1 + count)) % boundSize;
                    if (nextidx < 0) nextidx += boundSize;
                    if (seam->inverted) nextidx = (startAndPatch.first + (1 + count)) % boundSize;
                    int nextvert = boundaryL[startAndPatch.second][nextidx];

                    // check if vert is in map, if so then add constraint to connect it
                    if(trackCornerIds.find(vert) != trackCornerIds.end()){
                        model.addConstr(cutVar[ trackCornerIds[vert]] + cutVar[varCount] <= 1);
                    }else{
                        trackCornerIds[vert] = varCount;
                    }
                    int prevVert=-1;
                    while (vert != end) {

                        double relId = ((double) count) / seam->seamLength();
                        bool isConstrained = true;
                        bool corner = false;
                        if(count ==0) corner = true;
                        int currVar = varCount;
                        if (relId == 0 ||
                            (relId > minConstrained && relId < (1 - minConstrained)  && distStartToEnd > widthThereshold)

                        ) {
                            isConstrained = false;
                        }else{
                            model.addConstr(cutVar[varCount] == 0);
                        }
                        addVarToModel(vert, prevVert, nextvert, vfAdj, isConstrained, varCount, cutVar, Fg_pattern, lengthsOrig,
                                      lengthsCurr, mapVarIdToVertId, seamId[si].first, seamId[si].second, tailor_lazyness, corner, mapVertAndSeamToVar, otherVert, model );
                        if (!isConstrained) {
                            lSumConstr += cutVar[currVar];
                            if (relId != 0) {
                                rSumConstr += cutVar[currVar];
                                innerSumConstr += cutVar[currVar];
                            }
                        }
                        prevVert = vert;
                        vert = nextvert;
                        count++;
                        nextidx = (startAndPatch.first - (1 + count)) % boundSize;
                        if (nextidx < 0) nextidx += boundSize;
                        if (seam->inverted) nextidx = (startAndPatch.first + (1 + count)) % boundSize;
                        nextvert = boundaryL[startAndPatch.second][nextidx];
                        otherVert =  boundaryL[startAndPatchOther.second][(startAndPatchOther.first + count) %  boundaryL[startAndPatchOther.second].size() ];

                    }// handle the least element and set constraints
                    if(trackCornerIds.find(end) != trackCornerIds.end()){
                        model.addConstr(cutVar[ trackCornerIds[end]] + cutVar[varCount] <= 1);
                    }else{
                        trackCornerIds[end] = varCount;
                    }
                    rSumConstr += cutVar[varCount];


                    addVarToModel(vert, prevVert, -1, vfAdj, false, varCount,
                                  cutVar, Fg_pattern, lengthsOrig, lengthsCurr, mapVarIdToVertId,
                                  seamId[si].first, seamId[si].second, tailor_lazyness, true, mapVertAndSeamToVar, otherVert, model);

                }

                // set lin expr
                model.addConstr(innerSumConstr <= 1);
                model.addConstr(lSumConstr <= 1);
                model.addConstr(rSumConstr <= 1);
            }
        }

    }
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    model.optimize();


    for(int i =0; i< varCount; i++){
        if( cutVar[i].get(GRB_DoubleAttr_X ) >0.99){
            if(mapVarIdToVertId.find(i)== mapVarIdToVertId.end()){
                continue;// it is one of the XNOR helper variables,noting associated with them
            }
            int vert = mapVarIdToVertId[i]->vert;
            int seamType = mapVarIdToVertId[i]->seamType;
            int seamId = mapVarIdToVertId[i]->seamIdInList;
            int patch ;
            if(mapVarIdToVertId[i]->seamType > 0 ){
                if(mapVarIdToVertId[i]->seamIdInList>=0){
                    patch = seamsList[mapVarIdToVertId[i]->seamIdInList]->getStartAndPatch1().second;

                }else{
                    patch = seamsList[(mapVarIdToVertId[i]->seamIdInList +1)*(-1)]->getStartAndPatch2().second;
                }

            }else{
                patch = minusOneSeams[mapVarIdToVertId[i]->seamIdInList]->getPatch();
            }
            cout<<patch<<" patch, Id=  "<< mapVarIdToVertId[i]->vert<<" from seam "<<mapVarIdToVertId[i]->seamType<<" and id "<<mapVarIdToVertId[i]->seamIdInList <<" ";// ALRIGHT BUT NOW WE NEED A SEAM ID

            cutVertEntry* cve = new cutVertEntry ( vert, seamType, seamId, patch);


            cveStartPositionsSet[vert] = cve;
            cve-> leftCorner =  -1;
            cve-> rightCorner = -1;
            int nextVert, prevVert;
            double nextStress, prevStress ;

            bool inverted = false;

            if(seamId<0) inverted = seamsList[(-1)*(seamId+1)]->inverted;
             getPrevAndNextVertAndStress( seamType, seamId, vert, prevVert, nextVert, prevStress, nextStress,
                                              seamsList, minusOneSeams, boundaryL, Fg_pattern, lengthsOrig, lengthsCurr, vfAdj, inverted );
            if(vert == 788|| vert == 783){
                cout<<vert<<": "<<prevVert<<" prev, next "<<nextVert<<" "<<prevStress<<" stresses "<<nextStress<<endl;
            }
            if(cornerVert[vert]==1){
                cout<<"corner "<<endl;
                cve->cornerInitial = vert;
//                // left or right corner?
                int firstInSeam;
                if(seamType == -1 ){
                    firstInSeam = minusOneSeams[seamId]->getStartVert() ;
                }else if(seamId >= 0){
                    firstInSeam = seamsList[mapVarIdToVertId[i]->seamIdInList]->getStart1();
                }else{
                    firstInSeam = seamsList[(mapVarIdToVertId[i]->seamIdInList+1)*(-1)]->getStart2();
                }
                if(firstInSeam == vert){
                    cve -> startCorner = true;
                    cve->stress = nextStress;
                    //experiment
                    if(seamId<0 && !inverted) cve->stress = prevStress;

//                    cve->continuedDirection = Vg.row(nextVert)- Vg.row(vert);
                }else{
                    cve-> endCorner = true;
                    cve->stress = prevStress;
                    //experiment
                    if(seamId<0 &&  !inverted) cve->stress = nextStress;
//                    cve->continuedDirection = Vg.row(prevVert)- Vg.row(vert);
                }

            }else{
                cve->stress = (nextStress + prevStress)/2;
            }
            cout<<cve->stress<<" the stress there "<<endl<<endl;
            cutPositions.push_back(cve);
        }
    }
    cout<<cutPositions.size()<<" size"<<endl<<endl;
}

void updateStress(vector<cutVertEntry*>& cutPositions, vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                  std::vector<std::vector<int> >& boundaryL, MatrixXi& Fg_pattern, vector<vector<int>> & vfAdj, MatrixXd& lengthsCurr, bool& prioInner,
                  bool& prioOuter){
    cout<<"Updating the stress"<<endl;
    for(int i=0; i< cutPositions.size(); i++){
        cutVertEntry* cve = cutPositions[i];
        int nextVert, prevVert;
        double nextStress, prevStress;
        bool inverted = false;
        if(cve->seamIdInList < 0){
            inverted = seamsList[(cve->seamIdInList + 1)*(-1)]->inverted;
        }
        getPrevAndNextVertAndStress(cve -> seamType, cve -> seamIdInList, cve -> vert, prevVert, nextVert, prevStress, nextStress,
                                     seamsList, minusOneSeams, boundaryL, Fg_pattern, lengthsOrig, lengthsCurr, vfAdj, inverted );

                // it is a start corner and it still is!
        if(cve->startCorner && (cve->vert == cve->cornerInitial)){
            cve->stress = nextStress;
            if(cve->seamIdInList < 0 && ! inverted ) cve->stress = prevStress;
        }else if (cve->endCorner && (cve->vert == cve->cornerInitial)){
            cve->stress = prevStress;
            if(cve->seamIdInList < 0 && !inverted ) cve->stress = nextStress;

        }else{
            cve -> stress = (prevStress + nextStress)/2;
        }
        if(cve->startCorner || cve->endCorner){
            // check what kind of seam the other side is
            int thisSeam;
            if(cve->seamType == -1){
                thisSeam = (-1) * cve->seamIdInList -1;
            }else if (cve->seamIdInList<0) {
                thisSeam = (-1) * cve->seamIdInList -1;
            }else{
               thisSeam =  cve->seamIdInList;
            }

            int otherSeam = (cornerToSeams[cve->cornerInitial][0] == thisSeam) ? cornerToSeams[cve->cornerInitial][1] : cornerToSeams[cve->cornerInitial][0];
            if(prioOuter && otherSeam < 0){
                cve->stress +=1;
            }
            if(prioInner && otherSeam >= 0){
                cve->stress +=1;
            }

        }else{
            if(prioOuter && cve->seamType==-1){
                cve->stress +=1;
            }
            if(prioInner && cve->seamType ==1){
                cve->stress +=1;
            }
        }

    }

}
int tearFurtherVisIdxHelper(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL,
                set<int> & cornerSet, set<int>& handledVerticesSet,  bool& prevFinished, const bool & preferManySmallCuts, const bool & LShapeAllowed,
                MatrixXd& patternEdgeLengths_orig, MatrixXd& Vg_pattern_orig, bool& prioInner,
                bool& prioOuter ){
    cout<<endl<<endl<<"-----------------------"<<endl<<endl;

    int returnPosition = -1;
    //when releasing the boundary it can turn into a non manifold mesh. not sure if this causes further problems
    Eigen::MatrixXi B;
    bool isManifold = igl::is_vertex_manifold( Fg_pattern, B);

    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    lengthsOrig = patternEdgeLengths_orig;

    updateStress( cutPositions, seamsList, minusOneSeams, boundaryL,  Fg_pattern, vfAdj, lengthsCurr,  prioInner, prioOuter);

    if(prevFinished || preferManySmallCuts){
        // if we want many small cuts we sort always and there is no need to finish a seam before handling the next one!
        cout<<"It's time to sort again"<<endl;
        sort(cutPositions.begin(), cutPositions.end(), []( cutVertEntry* &a,  cutVertEntry* &b) { return a->stress > b-> stress; });


    }

    int  count = 0 ;
    for(int i = 0; i < 1; i++){

        int currVert = cutPositions[count]->vert;
        bool parallelFinFlag = true;
        int parallel;
        int thisSeam = cutPositions[count]->seamIdInList;
        if(thisSeam<0) thisSeam = thisSeam*(-1)-1;

        if(cutPositions[count]->finFlag ){
                    i--;
                    count ++;
        }else{

            cout<<endl<< cutPositions[count]->vert<<" vertex up next handling with i= "<<count<<" /"<<cutPositions.size()-1<<endl;
            cout<<currPattern.row(cutPositions[count]->vert)<<"final position in rest shape"<<endl ;

            returnPosition = cutPositions[count] ->vert;
            return returnPosition;
        }

    }
}

int tearFurther(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,
                vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                 map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL,
                 set<int> & cornerSet, set<int>& handledVerticesSet,  bool& prevFinished, const bool & preferManySmallCuts, const bool & LShapeAllowed,
                 MatrixXd& patternEdgeLengths_orig, MatrixXd& Vg_pattern_orig, bool& prioInner,
                bool& prioOuter ){
    cout<<endl<<endl<<"-----------------------"<<endl<<endl;

            int returnPosition = -1;
    //when releasing the boundary it can turn into a non manifold mesh. not sure if this causes further problems
    Eigen::MatrixXi B;
    bool isManifold = igl::is_vertex_manifold( Fg_pattern, B);

    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    lengthsOrig = patternEdgeLengths_orig;

    updateStress( cutPositions, seamsList, minusOneSeams, boundaryL,  Fg_pattern, vfAdj, lengthsCurr,  prioInner, prioOuter);

    if(prevFinished || preferManySmallCuts){
        // if we want many small cuts we sort always and there is no need to finish a seam before handling the next one!
        cout<<"It's time to sort again"<<endl;
        sort(cutPositions.begin(), cutPositions.end(), []( cutVertEntry* &a,  cutVertEntry* &b) { return a->stress > b-> stress; });


    }

    int  count = 0 ;
    for(int i = 0; i < 1; i++){
        if(count == cutPositions.size()){
            cout<<endl<<"* * * * * * * * * * * * * * * * * * *"<<endl;
            cout<<"All done. Stop now. "<<endl;
            cout<<"* * * * * * * * * * * * * * * * * * * *"<<endl;
            return -1;
        }
        int currVert = cutPositions[count]->vert;
        bool parallelFinFlag = true;
        int parallel;
        int thisSeam = cutPositions[count]->seamIdInList;
        if(thisSeam<0) thisSeam = thisSeam*(-1)-1;

        if(cutPositions[count]->finFlag ){
            // if it is a corner and it has been released
            if((cutPositions[count]->startCorner || cutPositions[count]->endCorner) &&
                releasedVertNew.find( cutPositions[count]->cornerInitial) != releasedVertNew.end() ){
                int seamPotentiallyReleasedFrom = (cornerToSeams[cutPositions[count]-> cornerInitial][0] == thisSeam) ? cornerToSeams[cutPositions[count]-> cornerInitial][1] : cornerToSeams[cutPositions[count]-> cornerInitial][0];
                if(releasedVertNew[ cutPositions[count]->cornerInitial].end() != std::find( releasedVertNew[ cutPositions[count]->cornerInitial].begin(), releasedVertNew[ cutPositions[count]->cornerInitial].end(), seamPotentiallyReleasedFrom)){
                    parallel = openParallelPosition(cutPositions[count]-> cornerInitial, seamPotentiallyReleasedFrom, seamsList, cutPositions);//releasedVertNew[cutPositions[count]->cornerInitial]
                    if(parallel >= 0) parallelFinFlag = cutPositions[parallel]->finFlag;
                    if(parallelFinFlag){
                        i--;
                        count ++;
                    }else{
                        returnPosition = cutPositions[parallel] ->vert;
                        splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                                           minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, cornerSet,
                                           handledVerticesSet, LShapeAllowed, Vg_pattern_orig);
                        prevFinished = cutPositions[parallel]->finFlag;
                    }
                }else{
                    i--;
                    count ++;
                }

            }else{
                i--;
                count ++;
            }

        }else{

            cout<<endl<< cutPositions[count]->vert<<" vertex up next handling with i= "<<count<<" /"<<cutPositions.size()-1<<endl;
            returnPosition = cutPositions[count] ->vert;
            splitVertexFromCVE(cutPositions[count], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList, minusOneSeams, releasedVert,
                                   toPattern_boundaryVerticesSet,lengthsCurr, cornerSet, handledVerticesSet, LShapeAllowed, Vg_pattern_orig );
            cout<<"finished ? "<<cutPositions[count]->finFlag<<endl<<endl;
            prevFinished = cutPositions[count] -> finFlag;

            parallel = -1;
            // if it is a corner and it has been released
            if(cutPositions[count]->startCorner || cutPositions[count]->endCorner){
                // once we finished cutting one side, check if it was a side opening. If so we can go on with the other side
                //this is the released seam, so the other one is the one we have in common
                int seamPotentiallyReleasedFrom = (cornerToSeams[cutPositions[count]-> cornerInitial][0] == thisSeam) ? cornerToSeams[cutPositions[count]-> cornerInitial][1] : cornerToSeams[cutPositions[count]-> cornerInitial][0];
                parallel = openParallelPosition(cutPositions[count]-> cornerInitial, seamPotentiallyReleasedFrom, seamsList, cutPositions);
                // open the other side of a released seam .
                // attention this is not the other side of the seam but the same 3D corner of a different patch.
                // note that there is no guarantee there is a cut position. If not, we do not enforce it.
            }

            if(parallel<0) {
                cout<<"no proper parallel found"<<endl;
                continue;
            }else if(!cutPositions[parallel]->finFlag){
                cout<<"split parallel "<<endl;

                splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                                   minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, cornerSet,
                                   handledVerticesSet, LShapeAllowed, Vg_pattern_orig);
                cout<<"finished  p ? "<<cutPositions[parallel]->finFlag<<endl;
                prevFinished = (prevFinished && cutPositions[parallel] -> finFlag);

            }else{
                cout<<"parallel already finished"<<endl;
            }

        }

    }
    cout<<prevFinished<<" --------------------"<<endl;

    return returnPosition;
}

void smoothSingleCut( cutVertEntry*& cve, MatrixXd& Vg, std::vector<std::vector<int> >& boundaryL){
    auto boundary = boundaryL[cve->patch];
    int startCut;

    if(cve->startCorner || cve->endCorner){
        cout<<endl<<" side case"<<endl;
        int idx = 0;
        while(boundary[idx] != cve-> cornerInitial && boundary[idx] != cve->vert ){
            idx++;
        }
        int endCut;
        Vector3d R, Q;
        if(boundary[idx]== cve->cornerInitial){
             cout<<"if"<<endl;
            startCut = cve->cornerInitial;
            endCut = cve->vert;
            R = Vg.row(startCut);
            Q =Vg.row(cve->vert);
        }else{
            cout<<"e"<<endl;
            startCut = cve->vert;
            endCut = cve->cornerInitial;
            Q = Vg.row(startCut);
            R =Vg.row(cve->cornerInitial);

        }
        cout<<" start "<<startCut<<" end "<<endCut<<endl ;

//        Vector3d R = Vg.row(startCut);
//        Vector3d Q =Vg.row(cve->vert);

        while(boundary[idx] != endCut){
            cout<<boundary[idx]<<" one intermediate side ";
//        updatePositionToIntersection(Vg, boundary[idx], line);
            double t = (R-Q).dot(Q-Vg.row(boundary[idx]).transpose())/((R-Q).dot(R-Q));
            cout<<t<<" the dist"<<endl;
            Vg.row(boundary[idx])= Q-t*(R-Q);

            idx++;
        }

        return;
    }// we have a normal cut
    if(cve->leftCorner==-1 || cve->rightCorner ==-1){
        cout<< "there is nothing to smooth" <<endl;
        return;
    }
    cout<<endl<<"starting search with endpoints "<< cve -> rightCorner<<" and "<< cve -> leftCorner<<" and vertex "<<cve->vert<<endl;


    int idx = 0;

    while(boundary[idx] != cve-> leftCorner && boundary[idx]!= cve->rightCorner ){
        idx++;
    }
    int endCut;
    if(boundary[idx]== cve-> rightCorner){
        startCut = cve -> rightCorner;
        endCut = cve -> leftCorner;
    }else{
        startCut = cve -> leftCorner;
        endCut = cve -> rightCorner;
    }
    Vector3d R = Vg.row(startCut);
    Vector3d Q =Vg.row(cve->vert);

    while(boundary[idx] != cve->vert){
        cout<<boundary[idx]<<" one intermediate l ";
//        updatePositionToIntersection(Vg, boundary[idx], line);
        double t = (R-Q).dot(Q-Vg.row(boundary[idx]).transpose())/((R-Q).dot(R-Q));
        cout<<t<<" the dist"<<endl;
        Vg.row(boundary[idx])= Q-t*(R-Q);

        idx++;
    }
    idx++;

    R = Vg.row(endCut);
    while(boundary[idx] != endCut){
        cout<<boundary[idx]<<" one intermediate r ";
        double t = (R-Q).dot(Q-Vg.row(boundary[idx]).transpose())/((R-Q).dot(R-Q));
        cout<<t<<" the dist"<<endl;
        Vg.row(boundary[idx])= Q-t*(R-Q);

        idx++;
    }

}

void smoothCuts(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
            map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL, set<int> & cornerSet ){

    for(int i=0; i<8; i++){
        // increment from right to vert and then vert to left
        smoothSingleCut(cutPositions[i] , currPattern, boundaryL);
    }
}

int computeTear(Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern_curr, MatrixXd& patternlengthsOrig,
                 vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams, std::vector<std::vector<int> >& boundaryL, bool & finished,
                 const std::vector<std::vector<std::pair<int, int>>>& cornersPerBoundary, map<int, vector<pair<int, int>>>& seamIdPerCorner,
                 VectorXd& cornerVert, vector<cutVertEntry*>& cutPositions, map<int, pair<int, int>> & releasedVert,
                 set<int>& toPattern_boundaryVerticesSet, set<int> & cornerSet, set<int>& handledVerticesSet ,
                 bool& prevFinished,
                 const bool & LShapeAllowed,
                 bool& prioInner,
                 bool& prioOuter )
                 {
    cout<<boundaryL.size()<<" boundaryL size"<<endl;
    std::vector<std::vector<int> > boundaryLnew;
    igl::boundary_loop(Fg_pattern_curr, boundaryLnew);
    boundaryL.clear();
    boundaryL = boundaryLnew;
    cout<<boundaryL.size()<<" boundaryL size after"<<endl;


                     vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern_curr, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern_curr, lengthsCurr);

    lengthsOrig =  patternlengthsOrig;
//    igl::edge_lengths(fromPattern, Fg_pattern_orig, lengthsOrig);

    for(int i =0; i<seamsList.size(); i++){
        seam* currSeam = seamsList[i];
        int start1 = currSeam -> getStart1();
        int start2 = currSeam ->getStart2();
        int end1 = currSeam-> getEndCornerIds().first;
        int end2 = currSeam-> getEndCornerIds().second;

        addToMapIfNotExisting(start1, i);
        addToMapIfNotExisting(start2, i);
        addToMapIfNotExisting(end1, i);
        addToMapIfNotExisting(end2, i);

    }
    for(int i =0; i<minusOneSeams.size(); i++){
       minusOneSeam* currSeam = minusOneSeams[i];
       int start = currSeam-> getStartVert();
       int end = currSeam -> getEndVert();
       addToMapIfNotExisting( start, -i-1);
       addToMapIfNotExisting( end, -i-1);

    }

    // first we need to know where to tear, set up LP for this
    // information we need: stress. For stress, we need lengths old and ned
    // for lengths old and new we need which edge of face
    // for which edge of face we need face
    // for face we need vertices
    // for vertices we need boundary loop
    double tailor_lazyness = 1;
    double minConstrained = 0.25;
    setLP(boundaryL, vfAdj, Fg_pattern_curr, lengthsOrig, lengthsCurr, cornersPerBoundary, seamIdPerCorner,
          seamsList, minusOneSeams, tailor_lazyness, minConstrained, cutPositions,
          cornerVert, currPattern, LShapeAllowed);// Vg= currpattern

    cout<<"finished set lp "<<endl;
    //update with the preferences before sorting
    updateStress( cutPositions, seamsList, minusOneSeams, boundaryL,  Fg_pattern_curr, vfAdj, lengthsCurr,  prioInner, prioOuter);
    cout<<"update stress"<<endl;


    //  here we need to sort and check if handled already
    sort(cutPositions.begin(), cutPositions.end(), []( cutVertEntry* &a,  cutVertEntry* &b) { return a->stress > b-> stress; });


    for(int i = 0; i < cutPositions.size(); i++){
        findCorrespondingCounterCutPosition(cutPositions, i, cutPositions[i], currPattern, Fg_pattern_curr, vfAdj, boundaryL,
                                            seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr,Fg_pattern_curr, cornerSet, handledVerticesSet );

    }

    int count=0;
    prevFinished = false;
    int returnPosition = -1;
    // we cut the first one
    for(int i = 0; i < 1; i++){// cutPositions.size(); i++){
        if(cutPositions[count]->handled){
            cout<<"handled case, go to next"<<endl;
            i--;
            count ++;
            continue;
        }
        cout<<cutPositions[count]->stress<<" curr stress "<<cutPositions[count]->stressWithCounter<<endl;
        int currVert = cutPositions[count]->vert;
        returnPosition = cutPositions[count]->vert;

        splitVertexFromCVE(cutPositions[count], currPattern, Fg_pattern_curr, vfAdj, boundaryL, seamsList,
                           minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr,
                           cornerSet, handledVerticesSet, LShapeAllowed, fromPattern);
        prevFinished = cutPositions[count] -> finFlag;

        if( releasedVertNew.find(currVert) != releasedVertNew.end()){
            // dann können müssen wir ja auch die passende andere seite des cuts öffnen
            // achtung, das ist nicht die corresponding seam
            // es kann aber noch nicht von einer anderen Seite geöffnett sein. daher nimm idx [0] von den released vert
            int parallel = openParallelPosition(cutPositions[count]-> cornerInitial, releasedVertNew[currVert][0], seamsList, cutPositions);
            if(parallel<0) {cout<<"no proper parallel found"<<endl;continue;}

            splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern_curr, vfAdj, boundaryL, seamsList,
                               minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr,
                               cornerSet, handledVerticesSet, LShapeAllowed, fromPattern);
            prevFinished = (prevFinished && cutPositions[parallel] -> finFlag);

        }


//        int counterPart = cutPositions[count]->counterPartIdx;
//        if(counterPart >0){
//            cout<<endl<<" and corresponding "<<cutPositions[count]->counterPartIdx<<" "<<cutPositions[cutPositions[count]->counterPartIdx]->vert<<" with stress "<<cutPositions[cutPositions[count]->counterPartIdx]->stress<<endl;
//            splitVertexFromCVE(cutPositions[counterPart], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
//                             minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, Fg_pattern, cornerSet, handledVerticesSet);
//
//        }else {
//            cout<<"-1 seam, no counter to split"<<endl;
//        }

    }
//    for(int i=0; i<cutPositions.size(); i++){
//        cout<<cutPositions[i]->cornerInitial<<endl;
//    }
cout<<"fin compute tear "<<endl ;
    return returnPosition;
}
void updatePositionToIntersection(MatrixXd& p,int next, const MatrixXd& Vg_bound, bool shouldBeLeft){

    // we know that the first two indices of the face define the edge we want to intersect with
    // derive where QR and P meet = t, https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line
    double minDist = std::numeric_limits<double>::max();
    Vector3d minDistTarget;
    int mini=-1;double mint;
    for(int i=0; i< Vg_bound.rows()-1; i++){
        Vector3d R = Vg_bound.row(i);
        Vector3d Q =  Vg_bound.row(i+1);

        double t = (R-Q).dot(p.row(next).transpose()-Q)/((R-Q).dot(R-Q));
        double tbefore = t;

        t = max(0., t);
        t = min(1., t);
        Vector3d targetPos = Q+t*(R-Q);

        double dist = (targetPos - p.row(next).transpose()).norm();
        if( dist < minDist){
            mint = tbefore;
            minDist = dist;
            minDistTarget = targetPos;
            mini = i;
        }
    }
    double stiffness = 0.8; //todo
    p.row(next) += stiffness * (minDistTarget.transpose()-p.row(next));

}

/*we build a new structure to find the closest position on the original boundary
        *     v0*-------------------*v1*---------------------*v2
        *           - v_(len+1)-          - v_(len+1+1) -
        *
        *           where v_(len+i) are defined as (vi + v(i+1))/2 + eps* n w
        *           where n is the normal
        *
*/

/*   Idea:
     iterate over every seam, create a triangle mesh of the original boundary of that seam
     then for each boundary vertex, project it's current position to the boundary
     update the current position towards the projected
     */

void fillMatrixWithBoundaryVert(const vector<int>& boundary, const int& start, const int& end, const MatrixXd& mapToVg, MatrixXd& Vg_seam1to, bool inverted ){

    int countLen = 2;
    int startIdx = 0;
    int boundLen = boundary.size();

    while(boundary[startIdx] != start && startIdx <= boundLen){
        startIdx++;
        startIdx = startIdx % boundLen ;
    }
    if(boundary[startIdx] != start ) cout<<"START NOT FOUND ERROR "<<start<<endl;
    int endIdx = (inverted)? startIdx-1: startIdx+1;

    while(boundary[endIdx] != end && endIdx != startIdx ){// messy?
        endIdx = (inverted)? endIdx-1 : endIdx+1;
        if(endIdx <0) endIdx+= boundLen;
        endIdx = endIdx % boundLen;
        countLen++;
    }

    if(boundary[endIdx] != end ) {
        cout<<"END NOT FOUND ERROR "<<end<<endl;
        for(int i=0; i<boundary.size(); i++){
            cout<<boundary[i]<<endl;
        }
    }

    Vg_seam1to= MatrixXd(countLen , 3);
    for(int i=0; i < countLen; i++){
        int bid = (inverted) ? (startIdx - i)  :  ((startIdx + i) % boundLen);
        if(bid < 0){
            bid += boundLen;
        }
        int v1_oneSide = boundary[bid ];
        Vg_seam1to.row(i) = mapToVg.row(v1_oneSide);

    }

}
void projectBackOnBoundary(const MatrixXd & mapToVg, MatrixXd& p, const vector<seam*>& seamsList
                           ,const vector<minusOneSeam*> & minusOneSeams,
                           const std::vector<std::vector<int> >& boundaryL_toPattern, const std::vector<std::vector<int> >& boundaryL,
                           map<int, pair<int, int>> & releasedVert ,bool visFlag){

    int numSeams = seamsList.size();
    int count =0;
//    if(releasedVert.size()>0){
//        auto iter = releasedVertNew.begin();
//        cout<<"Released vert: "<<endl;
//        while (iter != releasedVertNew.end()) {
//            cout << count<<" [" << iter->first << ","
//                 << iter->second[0] << "]\n";
//            ++iter;
//            count++;
//        }
//    }


    for (int j = 0; j<numSeams; j++){
        seam* currSeam  = seamsList[j];
        auto stP1= currSeam-> getStartAndPatch1();
        auto stP2 =  currSeam -> getStartAndPatch2ForCorres(); // attention this is with respect to the original pattern
//        int len  = currSeam -> seamLength();
//        int boundLen1 = boundaryL_toPattern[stP1.second].size();
//        int boundLen2 = boundaryL_toPattern[stP2.second].size();

        // build the structure for closest search
        MatrixXd Vg_seam1to, Vg_seam2to;
        fillMatrixWithBoundaryVert(boundaryL_toPattern[currSeam->getStartAndPatch1().second], currSeam-> patch1startCornerIdOld , currSeam-> patch1endCornerIdOld, mapToVg, Vg_seam1to, false );
        fillMatrixWithBoundaryVert(boundaryL_toPattern[currSeam->getStartAndPatch2().second], currSeam-> patch2startCornerIdOld , currSeam-> patch2endCornerIdOld, mapToVg, Vg_seam2to, true );

        bool shoulBeLeft =true; // for 2 case

        // for each interior (=not corner) vertex of the new boundary we need to find the closest position on the polyline and map it there
        // todo never ever cut the corner or change the corner index
        pair<int, int> ends = currSeam->getEndCornerIds();

        int bsize = boundaryL[stP1.second].size();
        int next = currSeam -> getStart1();
        int nextIdx = 0 ;
        while(boundaryL[stP1.second][nextIdx] != next && nextIdx<bsize){
            nextIdx ++;
        }
        if(boundaryL[stP1.second][nextIdx] != next){
            cout<<"PROJECTION ERROR we dont find the index "<<endl;
        }

        while( next!= ends.first ){
            // it is not released, project on boundary
            if(releasedVert.find(next) == releasedVert.end()){
                updatePositionToIntersection( p, next,Vg_seam1to, true);
            }
            // else it is released somehow. But from which seam? If it is released from another seam then pull it to this boundary still
            else if( std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end()){
               updatePositionToIntersection( p, next,Vg_seam1to, true);
            }
            nextIdx++;
            next = boundaryL[stP1.second][(nextIdx) % bsize];
        }
        // the last corner. Again if it is constrained from another side pull it to boundary, else ignore since handled by corner
        if(releasedVert.find(next) == releasedVert.end()){
            updatePositionToIntersection( p, next,Vg_seam1to, true);

        } else if (std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end())){
            updatePositionToIntersection( p, next,Vg_seam1to, true);
        }

        /**********  second side  *********/
//        int i2 = 0;
        bsize = boundaryL[stP2.second].size();
        nextIdx = 0;
        next = currSeam -> getStart2();
        while (boundaryL[stP2.second][nextIdx] != next && nextIdx < bsize ){
            nextIdx ++;
        }
        if(boundaryL[stP2.second][nextIdx] != next){
            cout<<"PROJECTION ERROR 2 we dont find the index "<<next<<" should be on patch "<<stP2.second<<endl;
        }

        while( next!= ends.second ){
            // general case an interior vertex , if it is not constrained pull it to boundary
            if(releasedVert.find(next) == releasedVert.end() ){
                updatePositionToIntersection( p, next,Vg_seam2to, shoulBeLeft);
            } else if( std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end()){
                updatePositionToIntersection( p, next,Vg_seam2to, shoulBeLeft);
            }
            nextIdx -=1;// (nextidx - i2) % (bsize);

            if(nextIdx < 0) {nextIdx += bsize;}
            next = boundaryL[stP2.second][nextIdx];
        }
        if(releasedVert.find(next) == releasedVert.end()){
            updatePositionToIntersection( p, next,Vg_seam2to, shoulBeLeft);
        }else if (std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end()){
            updatePositionToIntersection( p, next,Vg_seam2to, shoulBeLeft);
        }

        // also project all duplicates of interior cut vertices
        for(const auto & addedVert : currSeam->duplicates){
            if(releasedVert.find(addedVert.second) == releasedVert.end() ){
                updatePositionToIntersection(p, addedVert.second, Vg_seam1to, true);
            } else if( std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end()){
                updatePositionToIntersection(p, addedVert.second, Vg_seam1to, true);
            }
        }

        for(const auto & addedVert : currSeam->duplicates2){

            if(releasedVert.find(addedVert.second) == releasedVert.end() ){
                updatePositionToIntersection(p, addedVert.second, Vg_seam2to, shoulBeLeft);
            } else if( std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end()){
                updatePositionToIntersection(p, addedVert.second, Vg_seam2to, shoulBeLeft);
            }
        }

    }
    for(int j = 0; j < minusOneSeams.size(); j++){

        minusOneSeam* currSeam  = minusOneSeams[j];
        int patch = currSeam -> getPatch();
        int startVert = currSeam -> getStartVert();
        int endVert = currSeam -> getEndVert();
        int boundLen = boundaryL[patch].size();

        // build the structure for closest search
        MatrixXd Vg_seamto;
        fillMatrixWithBoundaryVert(boundaryL_toPattern[patch], currSeam->startVertOld, currSeam->endVertOld,  mapToVg, Vg_seamto, false );

       int startidx = 0;
       while ( boundaryL[patch][(startidx)] != startVert && startidx < boundLen+1){
           startidx++;
       }
       if( boundaryL[patch][(startidx)] != startVert){
           cout<<   "ERROR IN -1 SEAM , WE CANNOT FIND THE START VERT"<<endl ;
       }

        int next = startVert;int counter=0;
        while(next != endVert && counter < 1100){
            counter++;
            // general case, it is not released hence pull it to the boundary
            if(releasedVert.find(next) == releasedVert.end()){
                updatePositionToIntersection( p, next,Vg_seamto , true);
            }else if(std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(),  (-1)*(j+1)) == releasedVertNew[next].end()){
                // it is released but not from this seam,thus it has to stay on the projection
                updatePositionToIntersection( p, next,Vg_seamto, true);
            }
            startidx++;
            startidx = startidx % boundLen;
            next = boundaryL[patch][ startidx];
        }
        // it is released for another side hence we have to pull it to our side
        if(releasedVert.find(next) != releasedVert.end() && std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(),  (-1)*(j+1)) == releasedVertNew[next].end()){
            updatePositionToIntersection( p, next,Vg_seamto, true);

        }
        // also map all projections
        for(const auto & addedVert : currSeam -> duplicates){

            if(releasedVert.find(addedVert.second) == releasedVert.end() ){
                updatePositionToIntersection(p, addedVert.second, Vg_seamto, true);
            } else if( std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end()){
                updatePositionToIntersection(p, addedVert.second, Vg_seamto, true);
            }
//            updatePositionToIntersection(p, addedVert.second, Vg_seamto, true);

        }

    }

}

void updatePatchId(vector<cutVertEntry*>& cutPositions, const std::vector<std::vector<int> >& boundaryLnew, vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams){
    cout<<" The number of patches has changed. Create a mapping from one to the other and update all cve "<<endl;
    map<int, int> mapVertToNewPatch;
    for(int i = 0; i < boundaryLnew.size(); i++){
        for(int j = 0; j < boundaryLnew[i].size(); j++){
            mapVertToNewPatch[boundaryLnew[i][j]] = i;
        }
    }
    for(int i = 0; i < cutPositions.size(); i++){
        cutPositions[i] -> patch = mapVertToNewPatch[cutPositions[i] -> vert];
    }
//  ATTENTTION THE SEAM ID OF THE PATCH IS NOT UPDATED !!!
//  some  SEAMS ARE SPLIT BETWEE TWO PATCHES (THE NEW PATCH) AND THUS THE SEAMID IS NOT THE SAME FOR START END END
    for(int i=0; i < seamsList.size(); i++){
        seamsList[i]-> updatePatch1(mapVertToNewPatch[seamsList[i]-> getStart1()] );
        seamsList[i]-> updatePatch2(mapVertToNewPatch[seamsList[i]-> getStart2()] );
        seamsList[i]->seamSplit1 = (mapVertToNewPatch[seamsList[i]-> getStart1()] == mapVertToNewPatch[seamsList[i]-> getEndCornerIds().first]);
        seamsList[i]->seamSplit2 = (mapVertToNewPatch[seamsList[i]-> getStart2()] == mapVertToNewPatch[seamsList[i]-> getEndCornerIds().second]);


    }

    for(int i=0; i < minusOneSeams.size(); i++){
        minusOneSeams[i]-> updatePatch(mapVertToNewPatch[minusOneSeams[i]-> getStartVert()]);
        minusOneSeams[i]->seamSplit = (mapVertToNewPatch[minusOneSeams[i]-> getStartVert()] == mapVertToNewPatch[minusOneSeams[i]-> getEndVert()]);

    }

}

void computeCovarianceMatrix( MatrixXd& pointVec, VectorXd& barycenter, Matrix2d& m){
    int n;
    n = pointVec.rows();
    // first compute the barycenter
    barycenter = VectorXd::Zero(pointVec.cols()-1);
    for(int i = 0; i < n; i++){
        barycenter += pointVec.row(i).leftCols(2).transpose();
    }
    barycenter /= n;

    // compute covariance matrix
    m.resize(2, 2);
    m.setZero();
    MatrixXd p(2,1);

    for(int i = 0; i < n; i++){
        p = (pointVec.row(i).leftCols(2) - barycenter.transpose()).transpose();
        m += p*p.transpose();

    }

}
void fitVecToPointSet( MatrixXd& pointVec, VectorXd& vec ){

    Eigen::Matrix2d covMat ;
    VectorXd b;
    computeCovarianceMatrix(pointVec, b, covMat);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(covMat);
    Eigen::VectorXd eval = eig.eigenvalues();
    Eigen::Matrix2d evec = eig.eigenvectors();
    eval = eval.cwiseAbs();
    int minInd;
    eval.minCoeff(&minInd);

    vec(0) = evec(0,minInd);
    vec(1) = evec(1,minInd);
    vec = vec.normalized();
    vec *= eval;
//    vec += b;


}
