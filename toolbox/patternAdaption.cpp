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
#include <igl/writeOBJ.h>
using namespace std;
using namespace Eigen;

MatrixXd lengthsOrig;
map<int, cutVertEntry *> cveStartPositionsSet;
set<int> cutThroughCornerVertices;
double boundThereshold ;
double middleThereshold ;
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
            minusOneSeams[cve->seamIdInList]->duplicatesGlob[cve->vert] = newVertIdx;
        }else{
            if(cve->seamIdInList >= 0){
                seamsList[cve-> seamIdInList]-> duplicatesGlob[cve->vert] = newVertIdx;
            }else{
                seamsList[(cve-> seamIdInList+1)*(-1)]->duplicatesGlob2[cve->vert] = newVertIdx;
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

void computeMidVecBasedOnStress(VectorXd& midVec, vector<vector<int>>& vfAdj, MatrixXd& Vg ,MatrixXd& Vg_orig ,
                           MatrixXi& Fg, const MatrixXi& Fg_orig, int& vert, double& lenMid, map<int, int>& halfPatternVertToFullPatternVert, map<int, int>& fullPatternVertToHalfPatternVert,
                           map<int, int>& halfPatternFaceToFullPatternFace ){
    // new approach: look at the stress vectors of all adjacent faces, take the average
    Vector3d uVecs, vVecs;
    for(int i=0; i<vfAdj[vert].size(); i++){
        int face = vfAdj[vert][i];
        int faceOrig =halfPatternFaceToFullPatternFace[face];
        Vector3d v0 = Vg_orig.row(Fg_orig(faceOrig, 0)).transpose();
        Vector3d v1 = Vg_orig.row(Fg_orig(faceOrig, 1)).transpose();
        Vector3d v2 = Vg_orig.row(Fg_orig(faceOrig, 2)).transpose();

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

map<int, vector<int>> releasedVertNew;
map<int, vector<int>> cornerToSeams;

map<int, int> patchMapToHalf, patchMapToHalfInverse;bool inverseMapping, sym;
int addoncount=0;
MatrixXd uperFace, vperFace;

void addToMapIfNotExisting( int key, int i){
    if(cornerToSeams.find(key) != cornerToSeams.end()){
        cornerToSeams[key].push_back(i);
    }else{
        vector<int> vec;
        vec.push_back(i);
        cornerToSeams[key] = vec;
    }

}
void getBoundaryPrevVert(const bool& startCorner ,const int seamType, const int seamIdInList, const int minusOne, const int plusOne, int& next){
    if(startCorner){// does not matter if it is a starter or not

        if(seamType>0){
            if(seamIdInList>=0){
                next = plusOne;
            }else{
                next = minusOne;
            }
        }else{
            cout<<"the next one is "<<minusOne<<endl;

            next = plusOne;
        }

    }else{
        if(seamType > 0){
            if(seamIdInList >= 0){
                next = minusOne;
            }else{
                next = plusOne;
            }
        }else{
            next = minusOne;
        }
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

void findCorrectSeamAndAddToDuplicates(vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams, int vertIdx,
                                       int duplVertIdx, vector<int> boundary, int patch, map<int, int> & fullPatternVertToHalfPatternVert ){
    cout<<"finding Corrrect dupl"<<endl;
    bool found = false;
    int seam; bool inv = false;
    for(int i=0; i<seamsList.size(); i++){
        auto currSeam = seamsList[i];
        if(patchMapToHalf[currSeam->getPatch1()] == patch){
            int startVert = fullPatternVertToHalfPatternVert[currSeam->getStart1()];
            int end = fullPatternVertToHalfPatternVert[currSeam->getEndCornerIds().first];
            int idxOfStart = 0;
            while (boundary[idxOfStart]!= startVert){
                idxOfStart++;
            }
            while(boundary[idxOfStart] != end){

                if(boundary[idxOfStart] == vertIdx){
                    seamsList[i]-> duplicates[vertIdx] = duplVertIdx;
                    cout<<"located in seam "<<i<<endl;
                    return;
                }
                idxOfStart ++; idxOfStart  = idxOfStart % boundary.size();
            }
            if(boundary[idxOfStart] == vertIdx){
                seamsList[i]-> duplicates[vertIdx] = duplVertIdx;
                return;
            }


        }
        if(patchMapToHalf[currSeam->getPatch2()] == patch){
            int startVert = fullPatternVertToHalfPatternVert[currSeam->getStart2()];
            int end = fullPatternVertToHalfPatternVert[currSeam->getEndCornerIds().second];
            int idxOfStart = 0;
            while (boundary[idxOfStart]!= startVert){
                idxOfStart++;
            }
            while(boundary[idxOfStart] != end){

                if(boundary[idxOfStart] == vertIdx){
                    seamsList[i]-> duplicates2[vertIdx] = duplVertIdx;
                    cout<<"located in seam "<<i<<endl;

                    return;

                }
                idxOfStart --;if( idxOfStart<0)   idxOfStart += boundary.size();
            }
            if(boundary[idxOfStart] == vertIdx){
                seamsList[i]-> duplicates2[vertIdx] = duplVertIdx;
                return;
            }

        }

    }
    for(int i = 0; i<minusOneSeams.size(); i++){
        auto currSeam = minusOneSeams[i];
        if(currSeam->getPatch() != patch )continue;
        int startVert = fullPatternVertToHalfPatternVert[currSeam->getStartVert()];
        int end = fullPatternVertToHalfPatternVert[currSeam->getEndVert()];
        int idxOfStart = 0;
        while (boundary[idxOfStart]!= startVert){
            idxOfStart++;
        }
        while(boundary[idxOfStart] != end){
            if(boundary[idxOfStart] == vertIdx){
                minusOneSeams[i]-> duplicates[vertIdx] = duplVertIdx;
                return;
            }
            idxOfStart ++;
            idxOfStart = idxOfStart % boundary.size();
        }
        if(boundary[idxOfStart] == vertIdx){
            minusOneSeams[i]-> duplicates[vertIdx] = duplVertIdx;
            return;
        }


    }


    cout<<"No seam found on which we could locate "<<vertIdx<<endl;

}

void splitVertexFromCVE( cutVertEntry*& cve,
                         MatrixXd& Vg, // this is the current pattern we modify
                         MatrixXi& Fg, // this will be modified and have entries that are not in the original pattern
                         vector<vector<int> >& vfAdj,
                         std::vector<std::vector<int> >& boundaryL,
                         vector<seam*>& seamsList,
                         vector<minusOneSeam*> & minusOneSeams,
                         map<int, pair<int, int>> & releasedVert,
                         set<int>& toPattern_boundaryVerticesSet,
                         set<int> & cornerSet,
                         set<int>& handledVerticesSet ,
                         const bool & LShapeAllowed,
                         MatrixXd& Vg_pattern_orig,const MatrixXi& Fg_pattern_orig,
                         map<int, int> & fullPatternVertToHalfPatternVert,
                         map<int, int> & halfPatternVertToFullPatternVert,
                         map<int, int> & halfPatternFaceToFullPatternFace
                         ){

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
        bool cornerOrStart = ((cveStartPositionsSet.find(cve->vert) != cveStartPositionsSet.end()) || cornerSet.find(halfPatternVertToFullPatternVert[cve->vert]) != cornerSet.end()
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
            minusOneSeam* helper = minusOneSeams[cve->seamIdInList];
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

    Vector3d toLeft = Vg.row(boundaryL[cve->patch][plusOneId]) - Vg.row(cve->vert);
    Vector3d toRight = Vg.row(boundaryL[cve->patch][minusOneId])- Vg.row(cve->vert);

    // if it's a corner all get the new index
    if(cve->startCorner|| cve-> endCorner ){
        cout<<" start or end corner case"<<endl;
        int nextVertOnBoundary;
        getBoundaryNextVert(cve-> startCorner , cve-> seamType, cve-> seamIdInList,
                            boundaryL[cve->patch][minusOneId], boundaryL[cve->patch][plusOneId],nextVertOnBoundary );

        if(cve->vert!= cve->cornerInitial &&
        (cornerSet.find(cve->vert) != cornerSet.end()||
        cornerSet.find((-1)*cve->vert) != cornerSet.end()||
        cveStartPositionsSet.find(cve->vert) != cveStartPositionsSet.end() )
        ){
            getBoundaryPrevVert(cve-> startCorner , cve-> seamType, cve-> seamIdInList,
                                boundaryL[cve->patch][minusOneId], boundaryL[cve->patch][plusOneId],nextVertOnBoundary );
            cout<<" we should take the previous vert to get the right direction!"<<endl;
        }

        Vector3d cutDirection = Vg.row(nextVertOnBoundary) - Vg.row(cve->vert); //cve-> continuedDirection;
        cout<<"next vert on boundary "<<nextVertOnBoundary<<" "<<cutDirection.transpose()<<endl;

        VectorXd ws;
        vector<int> fn = vfAdj[cve->vert]; // the face neighbors
        bool tearIsUseful = false;
        for(auto faceIdx : vfAdj[cve->vert]){
            auto cd = cutDirection.normalized();
            cd(0)= -cutDirection.normalized()(1);
            cd(1)= cutDirection.normalized()(0);

            auto dot = uperFace.row(faceIdx).normalized().transpose().dot( cd);
            double actU = abs(dot);
            auto dotv = vperFace.row(faceIdx).normalized().transpose().dot( cd);
            double actV = abs(dotv);

           double w = actU * uperFace.row(faceIdx).norm() + vperFace.row(faceIdx).norm() * actV;
            cout<<w <<" = w, "<<actU<<" , "<<actV<< " contribution,  u norm "<<uperFace.row(faceIdx).norm()<<" ,v norm"<<vperFace.row(faceIdx).norm()<<", direction "<<cutDirection.transpose()<<endl;
            if(w > boundThereshold && (uperFace.row(faceIdx).norm()>1 || vperFace.row(faceIdx).norm()>1) ){
                tearIsUseful= true;
            }
           /*END TEST IF USEFUL*/
       }

        if(!tearIsUseful){
            cout<<"STOP now because of adj face stress condition  "<<endl;

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
        // now we need to find the seam from which we release.
        // each corner is adjacent to two seams. If one is the one we make the decision from, then the other is the one from which it is released

        int currComp ;
        if(cornerToSeams[cve->cornerInitial][0] == seamComp ){
            currComp = cornerToSeams[cve->cornerInitial][1];
            // if we allow releasing from two seams this has to be a vec
            if(releasedVertNew.find(cve->vert)!= releasedVertNew.end()){
                releasedVertNew[cve->vert].push_back( cornerToSeams[cve->cornerInitial][1]);

            }else{
                vector<int> temp;
                temp.push_back(cornerToSeams[cve->cornerInitial][1]);
                releasedVertNew[cve->vert] = temp;
            }
//            if(endInsert) releasedVertNew[cve->vert].push_back( furtherSeam);

        }else if (cornerToSeams[cve->cornerInitial][1] == seamComp ){
            currComp = cornerToSeams[cve->cornerInitial][0];

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

//        newVg.row(newVertIdx) = Vg.row(cve->vert);

        cve->continuedCorner = true;
        int searchVert= (halfPatternVertToFullPatternVert.find(cve->vert) == halfPatternVertToFullPatternVert.end()) ? (-1)*(cve->vert) : halfPatternVertToFullPatternVert[cve->vert];
        cve->finFlag = (cornerSet.find(searchVert) != cornerSet.end() && cve->vert != cve-> cornerInitial); //if it is a corner we are done
        int nextVertComp;
        getBoundaryNextVert(cve-> startCorner ,cve-> seamType, cve-> seamIdInList, boundary[minusOneId], boundary[ plusOneId], nextVertComp );
        cve->vert = nextVertComp;
//        Vg.resize(Vg.rows()+1, 3);
//        Vg= newVg;

        return;

    }
//todo enforce 90* more to make it more rigid!!

    Vector3d midVec;
    double lenMidVec;
    VectorXd  stressMidVec;

    computeMidVecBasedOnStress(stressMidVec, vfAdj, Vg, Vg_pattern_orig, Fg, Fg_pattern_orig, cve->vert, lenMidVec, halfPatternVertToFullPatternVert, fullPatternVertToHalfPatternVert, halfPatternFaceToFullPatternFace );

    midVec = (-1) * stressMidVec; // newMidVec;
    // todo sketchy, for whatever reason it breaks without this. does it do the transposing?
    cout<<" stress midvec "<<midVec.transpose()<<endl;//<<stressMidVec.transpose()

    Vector3d cutDirection = midVec;

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

    C = Vg.row(cve -> vert).transpose() + stressMidVec;
    auto dist2 = abs((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0)));
    dist2 /= sqrt((Btemp(0)-A(0) )*(Btemp(0)-A(0) ) + (Btemp(1)-A(1))*(Btemp(1)-A(1)));

    // we have a vector from the pca, but do we need to change the sing?
    if(!cve-> levelOne){
        if(dist2 > dist1 ){
            midVec = stressMidVec;
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
                midVec = stressMidVec;
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
        VectorXd distThen = (Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[adjFace], 0))+
                             Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[adjFace], 1))+
                             Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[adjFace], 2))).transpose();
//        cout<<Fg_pattern_orig.row(halfPatternFaceToFullPatternFace[adjFace])<<" dis then:"<<(distThen/3).transpose()<<" "<<mid2.transpose()<<endl;

        distThen /= 3;
        distThen(0) -= mid2(1);
        distThen(1) += mid2(0); // perp to midvec measure stress
        MatrixXd distB;
        MatrixXd input(1, 3);
        input.row(0)= distThen;
        igl::barycentric_coordinates(input, Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[adjFace], 0)), Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[adjFace], 1)),
                                     Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[adjFace], 2)), distB);
        double lenThen = mid2.norm();

        VectorXd distNow = (Vg.row(Fg(adjFace, 0)) * distB(0)+
                            Vg.row(Fg(adjFace, 1)) * distB(1) +
                            Vg.row(Fg(adjFace, 2))* distB(2)).transpose();
//        cout<<distNow.transpose()<<" newPos, and bary "<< (Vg.row(Fg(adjFace, 0)) * (1./3)+ Vg.row(Fg(adjFace, 1)) *  (1./3) +  Vg.row(Fg(adjFace, 2))*  (1./3) )<<endl;
        distNow -= (Vg.row(Fg(adjFace, 0)) * (1./3)+
                    Vg.row(Fg(adjFace, 1)) *  (1./3) +
                    Vg.row(Fg(adjFace, 2))*  (1./3) ).transpose();
        double lenNow = distNow.norm();
//        cout<<"new face verts "<<Fg.row(adjFace)<<endl;
        cout<<"Face: "<<adjFace<<" or "<<halfPatternFaceToFullPatternFace[adjFace]<<" "<<lenNow<<" dist now and then "<<lenThen<<" ratio is "<<lenNow/lenThen<<" th:"<<middleThereshold<<endl;

        if(lenNow/lenThen > middleThereshold ){
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
            cout<<"Not useful. STOP"<<endl;
            cve->finFlag = true;
            return;
    }
    MatrixXd newVg (newVertIdx + 1, 3);
    newVg.block(0,0, newVertIdx, 3)= Vg;
    newVg.row(newVertIdx) = Vg.row(cve->vert);

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

    VectorXd updatedRestShapeVertPos = insertIdxInBary(0) * Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[intersectingFace], 0)) ;
    updatedRestShapeVertPos += insertIdxInBary(1) * Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[intersectingFace], 1)) ;
    updatedRestShapeVertPos += insertIdxInBary(2) * Vg_pattern_orig.row(Fg_pattern_orig(halfPatternFaceToFullPatternFace[intersectingFace], 2) );
    //  update all original edge lengths -> in main
    Vg_pattern_orig.row(halfPatternVertToFullPatternVert[insertIdx]) = updatedRestShapeVertPos;

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
    handledVerticesSet.insert(newVertIdx);
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
        newVg.resize(newnewVg.rows(), 3);
        newVg = newnewVg;
        cout<<" fin operation"<<endl;
        findCorrectSeamAndAddToDuplicates(seamsList, minusOneSeams, insertIdx, newnewVertIdx, boundaryL[cve->patch], cve->patch, fullPatternVertToHalfPatternVert );

            // if it is on the boundary anyways, find on which seam it is (the original) and add it to the duplicates of theat seam

        cve->finFlag= true;

    }

    Vg.resize(newVg.rows(), 3);
    Vg = newVg;

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
/*
|
.___b

.___a
|
 if we open the cut horizontally at a, we want to open at b as well since the two together form a seam (though the common seam has it's cut in the middle)
 */
int openParallelPositionInvSym(int& cornerInitial, int& seamType, vector<seam*>& seamsList, vector <cutVertEntry*>& cutPositions,
                         map<int, int> & fullPatternVertToHalfPatternVert, map<int, int> & halfPatternVertToFullPatternVert ){
    // if a paralell cut position exists we find and open it . But attention, it is not guaranteed the position exists. But we don't enforce it (might be worth a thought tough)
    // what is a parallel position? If we tear along the boundary we effectivelyy open another seam (the one we release from). This behaviour should be mirrored on the other side of the released from seam
    if(seamType < 0){
        return -1; // no parallel to open
    }
    seam* currSeam = seamsList[seamType];
    // figure out if the initial was pos or neg side, the one we are searching is the other side
    int searchedVert;
    int cI = (-1)*cornerInitial ; if(halfPatternVertToFullPatternVert.find(cornerInitial) != halfPatternVertToFullPatternVert.end())cI = halfPatternVertToFullPatternVert[cornerInitial];
    if(currSeam->getStart1() == halfPatternVertToFullPatternVert[cornerInitial]){
        searchedVert = currSeam->getStart2();
    }else if(currSeam->getStart2() == halfPatternVertToFullPatternVert[ cornerInitial]) {
        searchedVert = currSeam->getStart1();
    }else if(currSeam->getEndCornerIds().first == halfPatternVertToFullPatternVert[cornerInitial] ) {
        searchedVert = currSeam-> getEndCornerIds().second;
    }else if(currSeam->getEndCornerIds().second == halfPatternVertToFullPatternVert[cornerInitial] ) {
        searchedVert = currSeam-> getEndCornerIds().first;
    }else{
//        cout<<" partner not found, we have a huge problem in opening the parallel positions "<<seamType<<endl;
        return -1;
    }
    if(fullPatternVertToHalfPatternVert.find(searchedVert) == fullPatternVertToHalfPatternVert.end() ){
        cout<<" parallell is on other half. Not here. "<<endl;
        return -1;
    }
    searchedVert = fullPatternVertToHalfPatternVert[searchedVert];

//    int size =  cornerToSeams[searchedVert].size();
//    if( size != 2) {
//        cout<<searchedVert<<" the size is not 2, it's "<<size<<endl;
//    }
    int otherSeamId;
    if(cornerToSeams[searchedVert][0] == seamType){

        otherSeamId = cornerToSeams[searchedVert][1];

    }else{
        otherSeamId = cornerToSeams[searchedVert][0];
    }

    // identify the index of the paralell cut postion by comparing the vertex with the one we search
    for(int i=0; i<cutPositions.size(); i++){
        int cI = searchedVert ; if(halfPatternVertToFullPatternVert.find(searchedVert) != halfPatternVertToFullPatternVert.end()) cI = halfPatternVertToFullPatternVert[searchedVert];
        if(cutPositions[i]->vert == halfPatternVertToFullPatternVert[searchedVert] || cutPositions[i]->cornerInitial == halfPatternVertToFullPatternVert[searchedVert]){
            if(cutPositions[i]->seamIdInList == otherSeamId){
                return i;
            }
        }
    }
    return -1;

}
 int openParallelPosition(int& cornerInitial, int& seamType, vector<seam*>& seamsList, vector <cutVertEntry*>& cutPositions,
        map<int, int> & fullPatternVertToHalfPatternVert, map<int, int> & halfPatternVertToFullPatternVert ){
     // if a paralell cut position exists we find and open it . But attention, it is not guaranteed the position exists. But we don't enforce it (might be worth a thought tough)
    // what is a parallel position? If we tear along the boundary we effectivelyy open another seam (the one we release from). This behaviour should be mirrored on the other side of the released from seam
     if(seamType < 0){
        return -1; // no parallel to open
    }
    seam* currSeam = seamsList[seamType];
    // figure out if the initial was pos or neg side, the one we are searching is the other side
    int searchedVert;
    if(currSeam->getStart1() == halfPatternVertToFullPatternVert[cornerInitial]){
        searchedVert = currSeam->getStart2();
    }else if(currSeam->getStart2() == halfPatternVertToFullPatternVert[ cornerInitial]) {
        searchedVert = currSeam->getStart1();
    }else if(currSeam->getEndCornerIds().first == halfPatternVertToFullPatternVert[cornerInitial] ) {
        searchedVert = currSeam-> getEndCornerIds().second;
    }else if(currSeam->getEndCornerIds().second == halfPatternVertToFullPatternVert[cornerInitial] ) {
        searchedVert = currSeam-> getEndCornerIds().first;
    }else{
//        cout<<" partner not found, we have a huge problem in opening the parallel positions "<<seamType<<endl;
        return -1;
    }
    cout<<currSeam->getStart1()<<" "<<currSeam->getEndCornerIds().first<<"; "<< currSeam->getStart2()<<" "<<currSeam->getEndCornerIds().second<<endl;
    cout<<searchedVert<<" searched vert for corner initial "<<cornerInitial ;
    if(fullPatternVertToHalfPatternVert.find(searchedVert) == fullPatternVertToHalfPatternVert.end() ){
        cout<<" parallell is on other half. Not here. "<<endl;
        return -1;
    }
    searchedVert = fullPatternVertToHalfPatternVert[searchedVert];
    cout<<" updated to "<<searchedVert;

//    int size =  cornerToSeams[searchedVert].size();
//    if( size != 2) {
//        cout<<searchedVert<<" the size is not 2, it's "<<size<<endl;
//    }
    int otherSeamId;
    if(cornerToSeams[searchedVert][0] == seamType){

        otherSeamId = cornerToSeams[searchedVert][1];

    }else{
        otherSeamId = cornerToSeams[searchedVert][0];
    }

    // identify the index of the paralell cut postion by comparing the vertex with the one we search
    for(int i=0; i<cutPositions.size(); i++){
        if(cutPositions[i]->vert == halfPatternVertToFullPatternVert[searchedVert] || cutPositions[i]->cornerInitial == halfPatternVertToFullPatternVert[searchedVert]){
            if(cutPositions[i]->seamIdInList == otherSeamId){
                return i;
            }
        }
    }
    return -1;

}
void findCorrespondingCounterCutPosition(vector<cutVertEntry*>& cutPositions, int idxOfCVE, cutVertEntry*& cve, MatrixXd& Vg, MatrixXi& Fg, vector<vector<int> >& vfAdj,
                      std::vector<std::vector<int> >& boundaryL,  vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams,
                      map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet, set<int>& handledVerticesSet, map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& halfPatternVertToFullPatternVert){
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
        if(counterID > 0 && patchMapToHalf.find(seamsList[counterID]->getPatch1()) == patchMapToHalf.end() ){
            cout<<cve->vert<<" the corresponding seam is in the other half. We don't use it"<<endl;
            return;
        }
        if(counterID < 0 && patchMapToHalf.find(seamsList[cve->seamIdInList]->getPatch2()) == patchMapToHalf.end() ){
            cout<<cve->vert<<" the corresponding seam is in the other half. We don't use it"<<endl;
            return;
        }

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
        if(halfPatternVertToFullPatternVert[cve->vert] == currSeam-> getStart1()){
            lookFor = currSeam-> getStart2();
        }else if (halfPatternVertToFullPatternVert[cve->vert] == currSeam-> getStart2()) {
            lookFor = currSeam->getStart1();
        }else{
            cout<<"we have a problem, should find cve as corner but seems like its not "<<cve-> vert<<endl; return;
        }
    }
        // find the mapped vertex of this one
    else {
        if (halfPatternVertToFullPatternVert[cve->vert] == currSeam->getEndCornerIds().first) {
            lookFor = currSeam->getEndCornerIds().second;
        } else if (halfPatternVertToFullPatternVert[cve->vert] == currSeam->getEndCornerIds().second) {
            lookFor = currSeam->getEndCornerIds().first;
        } else {
            cout << "we have a problem, should find cve as corner but seems like its not " << cve->vert << endl;return;
        }
    }
    if(fullPatternVertToHalfPatternVert.find(lookFor)== fullPatternVertToHalfPatternVert.end()){
        return;
    }
    int idx = -1;
    for(int i = 0; i < cutPositions.size(); i++){

        if(cutPositions[i]->seamType == 1 && halfPatternVertToFullPatternVert[cutPositions[i]->vert] == lookFor && cutPositions[i]-> seamIdInList == counterID){
            idx= i;
        }
    }
    if(idx == -1) {
        cout<<cve->vert<<" aka "<< halfPatternVertToFullPatternVert[cve->vert]<<" not found, problem with cutting the corresponding, did not find "<<lookFor<<endl;
        return;
    }
    cve-> counterPartIdx = idx;
    cutPositions[idx]-> counterPartIdx = idxOfCVE;

    double stressSum = cve->stress + cutPositions[idx]->stress;
    cve-> stressWithCounter = stressSum;
    cutPositions[idx]->stressWithCounter = stressSum;

    return;
}

 void computePerFaceUV( const MatrixXi& Fg_pattern_curr, const MatrixXi&  mapFromFg, const MatrixXd& mapFromVg,const MatrixXd& currPattern, map<int, int>& halfPatternFaceToFullPatternFace, bool inverseMap){
     uperFace.resize(Fg_pattern_curr.rows(), 3);
     vperFace.resize(Fg_pattern_curr.rows(), 3);
     for(int i=0; i<Fg_pattern_curr.rows(); i++){
         // we use bary to left and bary to top +1 of the old and compute in bary coords
         // get the same of the new, check the new length
         //start form from pattern, that is the rest shape of the shape we have
         int idx = (inverseMap) ? i : halfPatternFaceToFullPatternFace[i];
         VectorXd v0 = mapFromVg.row(mapFromFg(idx, 0)).transpose();
         VectorXd v1 = mapFromVg.row(mapFromFg(idx, 1)).transpose();
         VectorXd v2 = mapFromVg.row(mapFromFg(idx, 2)).transpose();

         VectorXd bary = (v0 + v1 + v2)/3;
         VectorXd u = bary; u(0)+= 1;
         VectorXd v = bary; v(1)+= 1;
         VectorXd ubary, vbary;

         igl::barycentric_coordinates(u.transpose(),v0.transpose(),v1.transpose(),v2.transpose(),ubary);
         igl::barycentric_coordinates(v.transpose(),v0.transpose(),v1.transpose(),v2.transpose(),vbary);

         VectorXd v0new = currPattern.row(Fg_pattern_curr(i, 0)).transpose();
         VectorXd v1new = currPattern.row(Fg_pattern_curr(i, 1)).transpose();
         VectorXd v2new = currPattern.row(Fg_pattern_curr(i, 2)).transpose();
         VectorXd startPerEdge = ((v0new + v1new + v2new)/3).transpose();

         uperFace.row(i) = ((ubary(0) * v0new + ubary(1) * v1new + ubary(2) * v2new) - startPerEdge).transpose();
         vperFace.row(i) = ((vbary(0) * v0new + vbary(1) * v1new + vbary(2) * v2new) - startPerEdge).transpose();

     }
 }
void addVarToModel (bool inverseMap, int vert, int prevVert, int nextVert, vector<vector<int>> & vfAdj, bool isConstrained, int& varCount, GRBVar* & cutVar,
                      MatrixXi& Fg_pattern, map <int, cutVertEntry*> & mapVarIdToVertId,
                      int seamType, int seamIdInList, double tailor_lazyness, bool corner, map<pair<int,int>, int>& mapVertAndSeamToVar, int counterpart, GRBModel& model, const MatrixXd& currPattern ){

    // idea: when inverse mapping allow only cuts from the boundary, but L shapes should be ok !
    double w_init = 0;
    int count = 0;
    VectorXd nextDir= VectorXd::Zero(3);
    if(nextVert != -1){
        nextDir += currPattern.row(vert) - currPattern.row(nextVert);
    }if(prevVert != -1){
        nextDir += currPattern.row(prevVert) - currPattern.row(vert);
    }
    for(auto faceIdx : vfAdj[vert]){
        count++;

        auto dot = uperFace.row(faceIdx).normalized().transpose().dot( nextDir.normalized());
        double actU = abs(dot);
        auto dotv = vperFace.row(faceIdx).normalized().transpose().dot( nextDir.normalized());
        double actV = abs(dotv);

        w_init += uperFace.row(faceIdx).norm() * (actU);
        w_init += vperFace.row(faceIdx).norm() * (actV);

    }
    w_init/= count;
    if(!inverseMap) {
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
//            cout<<tailorLazyFactor<<" "<<vert <<" taylor lazy factor and weight "<<w_init<< " "<<count<<endl;
            try{
                cutVar[varCount].set(GRB_DoubleAttr_Obj, tailorLazyFactor * w_init);

            }catch(GRBException e){
                cout << "Error code = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
                cout<<varCount<<endl;

            }
        }
    }else {
        if (!corner) {
            try {
                cutVar[varCount].set(GRB_DoubleAttr_Obj, 0);
            } catch (GRBException e) {
                cout << "Error code = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
                cout << varCount << endl;
            }
        } else {

            int tailorLazyFactor = 1;
            if (corner) tailorLazyFactor *= tailor_lazyness;
//            cout <<vert<<" "<< tailorLazyFactor << " taylor lazy factor and weight " << w_init << " " << varCount << endl;
            try {
                cutVar[varCount].set(GRB_DoubleAttr_Obj, tailorLazyFactor * w_init);

            } catch (GRBException e) {
                cout << "Error code = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
                cout << varCount << endl;

            }
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

double computeSeamLength(pair<int, int>& seamId, vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,
                         vector<vector<int>>& boundaryL, MatrixXd& Vg, map<int, int>& fullPatternVertToHalfPatternVert ){

    double len=0;
    int type = seamId.first;
    int seamIdx = seamId.second;
    if (type<0){
        minusOneSeam* currSeam = minusOneSeams[seamIdx];
        int vertStart = currSeam->getStartVert();
        int vertEnd = currSeam->getEndVert();
        if(fullPatternVertToHalfPatternVert.find(vertStart)== fullPatternVertToHalfPatternVert.end() || fullPatternVertToHalfPatternVert.find(vertEnd)== fullPatternVertToHalfPatternVert.end()){
            return 0;
        }else{
            vertStart = fullPatternVertToHalfPatternVert[vertStart];
            vertEnd = fullPatternVertToHalfPatternVert[vertEnd];
        }
        vector<int> boundary = boundaryL[patchMapToHalf[currSeam->getPatch()]];
        int startIdx;
        findIndxInBoundaryloop( boundary, vertStart , startIdx);

        while(boundary[startIdx] != vertEnd){
            int next = (startIdx + 1)% boundary.size();

            len += (Vg.row(boundary[startIdx])-Vg.row(boundary[next])).norm();
            startIdx = next;
        }

    }else{
        if(seamIdx >=0){
            seam* currSeam = seamsList[seamIdx];
            int vertStart = currSeam->getStart1();
            int vertEnd = currSeam->getEndCornerIds().first;
            if(fullPatternVertToHalfPatternVert.find(vertStart)== fullPatternVertToHalfPatternVert.end() || fullPatternVertToHalfPatternVert.find(vertEnd)== fullPatternVertToHalfPatternVert.end()){
                return 0;
            }else{
                vertStart = fullPatternVertToHalfPatternVert[vertStart];
                vertEnd = fullPatternVertToHalfPatternVert[vertEnd];
            }
            vector<int> boundary = boundaryL[patchMapToHalf[currSeam->getPatch1()]];
            int startIdx;
            findIndxInBoundaryloop( boundary, vertStart , startIdx);

            while(boundary[startIdx] != vertEnd){
                int next = (startIdx+1)% boundary.size();

                len += (Vg.row(boundary[startIdx])-Vg.row(boundary[next])).norm();
                startIdx = next;
            }

        }else{
            seam* currSeam = seamsList[(seamIdx+1)*(-1)];
            int vertStart = currSeam->getStart2();
            int vertEnd = currSeam->getEndCornerIds().second;
            if(fullPatternVertToHalfPatternVert.find(vertStart)== fullPatternVertToHalfPatternVert.end() || fullPatternVertToHalfPatternVert.find(vertEnd)== fullPatternVertToHalfPatternVert.end()){
                return 0;
            }else{
                vertStart = fullPatternVertToHalfPatternVert[vertStart];
                vertEnd = fullPatternVertToHalfPatternVert[vertEnd];
            }
            vector<int> boundary = boundaryL[patchMapToHalf[currSeam->getPatch2()]];
            int startIdx;
//            cout<<"The patch "<<patchMapToHalf[currSeam->getPatch2()]<<" "<<currSeam->getStart2()<<" "<<vertStart<<endl;
            findIndxInBoundaryloop( boundary, vertStart , startIdx);
//            cout<<"found index "<<startIdx<<endl;
            int count=0;
            while(boundary[startIdx] != vertEnd){
                count++;
                int next = (startIdx-1);
                if(next<0) next+= boundary.size();
                if(currSeam->inverted) cout<<"iiNv issue"<<endl; //next = (startIdx + count) % boundary.size();
                if(seamIdx== -8) cout<< (Vg.row(boundary[startIdx])-Vg.row(boundary[next])).norm()<<endl;

                len += (Vg.row(boundary[startIdx])-Vg.row(boundary[next])).norm();
                startIdx = next;
            }

        }

    }


    return len;
}

void getStressAtVert(int seamType, int seamId, int vert, int & prevVert, int & nextVert, double & Stress,
                                 vector<seam*>& seamsList, const vector<minusOneSeam*>& minusOneSeams, std::vector<std::vector<int> >& boundaryL, MatrixXi& Fg_pattern,
                                MatrixXd& currPattern, vector<vector<int>> & vfAdj, bool inverted, map<int, int>& fullPatternVertToHalfPatternVert ){
    int patch = -1;
    int idx= -1; nextVert = -1; prevVert = -1;
    // uggly but needed after patch update
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
        // ensure it is not the end
        if(fullPatternVertToHalfPatternVert[ minusOneSeams[seamId]->getEndVert()] != vert)nextVert = minusOneSeams[seamId]->getNextVert(vert, boundaryL[patch]);
        if(fullPatternVertToHalfPatternVert[ minusOneSeams[seamId]->getStartVert()] != vert)prevVert= minusOneSeams[seamId]->getPrevVert(vert, boundaryL[patch]);
        cout<<patch<<" Computing stress for "<<vert<<" with adj "<<prevVert<<" "<<nextVert<<" of Seam"<<seamType<<seamId<<endl;


    }else if(seamId>=0){
        if(inverted){
            cout<<"inverted issue, patch should be 2"<<endl;
            seamId = (seamId+1)*(-1);
        }
        if(fullPatternVertToHalfPatternVert[ seamsList[seamId ]->getEndCornerIds().first] != vert) nextVert = seamsList[seamId]->getNextVert1(vert, boundaryL[patch]);
        if(fullPatternVertToHalfPatternVert[ seamsList[seamId ]->getStart1()] != vert) prevVert = seamsList[seamId]->getPrevVert1(vert, boundaryL[patch]);
        cout<<patch<<" Computing stress for "<<vert<<" with adj "<<prevVert<<" "<<nextVert<<" of Seam"<<seamType<<seamId<<endl;

    }else{
        if(fullPatternVertToHalfPatternVert[seamsList[(seamId +1)*(-1) ]-> getStart2() ] != vert) nextVert = seamsList[(seamId +1)*(-1)]->getPrevVert2(vert, boundaryL[patch]);
        if(fullPatternVertToHalfPatternVert[ seamsList[(seamId +1)*(-1) ]-> getEndCornerIds().second ] != vert) prevVert = seamsList[(seamId +1)*(-1)]->getNextVert2(vert, boundaryL[patch]);
        cout<<patch<<" Computing stress for "<<vert<<" with adj "<<prevVert<<" "<<nextVert<<" of Seam"<<seamType<<seamId<<endl;

    }

    double w_init = 0;
    int count = 0;
    VectorXd nextDir= VectorXd::Zero(3);
    if(nextVert != -1){
        nextDir += currPattern.row(vert) - currPattern.row(nextVert);
    }
    if(prevVert != -1){
        nextDir += currPattern.row(prevVert) - currPattern.row(vert);
    }


    for(auto faceIdx : vfAdj[vert]){
        count++;
        auto dot = uperFace.row(faceIdx).normalized().transpose().dot( nextDir.normalized());
        double actU = abs(dot);
        auto dotv = vperFace.row(faceIdx).normalized().transpose().dot( nextDir.normalized());
        double actV = abs(dotv);

        w_init += uperFace.row(faceIdx).norm() * (actU);
        w_init += vperFace.row(faceIdx).norm() * (actV);
        cout<< uperFace.row(faceIdx).norm()<<" * "<<actU<<" + "<<vperFace.row(faceIdx).norm()<<" * "<<actV<<" = ";


    }
    w_init/= count;
    cout<<w_init<<endl;
    Stress = w_init ;

}

void setLP(bool inverseMap, std::vector<std::vector<int> >& boundaryL , vector<vector<int>> & vfAdj, MatrixXi& Fg_pattern
           ,const std::vector<std::vector<std::pair<int, int>>>& cornersPerBoundary,
        map<int, vector<pair<int, int>>>& seamIdPerCorner, vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,
        double tailor_lazyness, double minConstrained, vector <cutVertEntry*>& cutPositions, VectorXd& cornerVert,
        MatrixXd& currPattern, const bool & LShapeAllowed, const MatrixXd & fromPattern, const MatrixXi mapFromFg, map<int, int>& fullPatternVertToHalfPatternVert,
        map<int, int>& halfPatternVertToFullPatternVert,map<int, int>& halfPatternFaceToFullPatternFace ){
    computePerFaceUV(Fg_pattern, mapFromFg, fromPattern, currPattern,halfPatternFaceToFullPatternFace, inverseMap );
    for(auto it:fullPatternVertToHalfPatternVert ){
        if(it.first<0){
            cout<<"it is in!!"<<it.first<<" "<<it.second<<endl;
        }
    }
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
    numVar += (numVar/2);// 574 without mapping -> 287 added ones ? we don't need them all bc there are no common ones for th -1s
    map<int, int> trackCornerIds; // keeps track of the var Id per corner - to understand if it exists already or not
    map <int, cutVertEntry*>  mapVarIdToVertId;
    map <pair<int, int>, int> mapVertAndSeamToVar;
//    if(inverseMap && (fullPatternVertToHalfPatternVert != halfPatternVertToFullPatternVert)) numVar *= 2; // we have symetry, account for additional vars?
    cout<<"num var "<<numVar<<endl<<endl;
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
            cout<<endl<<seamId.size()<<" corner "<<cornerPair.first<<" size of seams per corner here"<<endl;
            if(seamId.size()>2) cout<<" something is veryy odd!! we have more than two seams for a corner. impossible."<<endl;

            for(int si = 0; si < seamId.size(); si++) {
                cout<<"corner "<<cornerPair.first<<" info: "<<seamId[si].first<<" "<<seamId[si].second<<endl;

                GRBLinExpr innerSumConstr = 0; // interior sum
                GRBLinExpr lSumConstr = 0; // left sum
                GRBLinExpr rSumConstr = 0; // right sum
                GRBLinExpr totalSum =0;
                int startVarOfThis = varCount;
                int count = 0;
                double distStartToEnd = computeSeamLength(seamId[si], seamsList, minusOneSeams, boundaryL, currPattern, fullPatternVertToHalfPatternVert);

                double widthThereshold = 50;
                // check for direction
                // first if it is a seam or a -1 seam
                // second gives the direction, if less than 0 take i -1 in negative direction
                if ((seamId[si].first >= 0 && seamId[si].second >= 0) || seamId[si].first < 0) {
                    int idx;
                    int vert;// first corner
                    int end, endOther;// last corner
                    int length;
                    bool inHalf, inHalfOther; int patch1, patch2;
                    pair<int, int> startAndPatchOther;
                    int idxOther,idxOtherOrig ;
                    int vertOther = -1;
                    int boundSizeOther;

                    if (seamId[si].first >= 0) {
                        seam *seam = seamsList[seamId[si].second];
                        //gives an index and it is updated -> is it right? to check
                        auto startAndPatch = seam->getStartAndPatch1();
                        startAndPatchOther = seam-> getStartAndPatch2ForCorres();

                        int startVal = seam->getStart1();
                        inHalf = (fullPatternVertToHalfPatternVert.find(startVal) != fullPatternVertToHalfPatternVert.end());
                        if(inHalf){
                            cout<<"in half"<<endl;
                            startVal = fullPatternVertToHalfPatternVert[startVal];
                            findIndxInBoundaryloop(boundaryL[i], startVal, idx);
                            patch1 =  patchMapToHalf[startAndPatch.second];
                            vert = boundaryL[patch1][idx];
                            end = fullPatternVertToHalfPatternVert[seam->getEndCornerIds().first];
                        }else{
                            continue;
                        }

                        startVal = seam->getStart2();
                        inHalfOther = (fullPatternVertToHalfPatternVert.find(startVal) != fullPatternVertToHalfPatternVert.end());
                        if(inHalfOther){
                            startVal = fullPatternVertToHalfPatternVert[startVal];
                            cout<<" start other "<<startVal<<endl;
                            patch2 =  patchMapToHalf[startAndPatchOther.second];
                            findIndxInBoundaryloop(boundaryL[patch2], startVal, idxOther);idxOtherOrig = idxOther;
                            vertOther = boundaryL[patch2][idxOther];
                            boundSizeOther = boundaryL[patch2].size();
                            endOther = fullPatternVertToHalfPatternVert[seam->getEndCornerIds().second];
                            cout<<" end other "<<endOther<<endl;


                        }

                        length = seam->seamLength();// heuristic for constrained area

                    } else if (seamId[si].first < 0) {
                        minusOneSeam *currSeam = minusOneSeams[seamId[si].second];
                        bool used = (patchMapToHalf.find(currSeam->getPatch()) != patchMapToHalf.end());
                        if(!used) continue;

                        int patch = patchMapToHalf[currSeam->getPatch()];

                        if (patch != i)
                            cout << " now in th e-1 seams the patch does not match where we are in the loop , would be patch "<<patch << endl;
                        int startVal = fullPatternVertToHalfPatternVert[ currSeam->getStartVert()];
                        findIndxInBoundaryloop(boundaryL[patch], startVal, idx);
//                        idx = currSeam->getStartIdx();
                        vert = boundaryL[i][idx];
                        end =fullPatternVertToHalfPatternVert[ currSeam->getEndVert()];
                        length = currSeam->getLength();
                    }

                    // check if the corners exist already. If so then connect with corner, else add the indices
                    // todo for L cutting allowance
                    cout<<vertOther<<" counterpart "<<vert<<" vert and end "<<end<<endl;
                    if(!inverseMap){
                        if(trackCornerIds.find(vert) != trackCornerIds.end()){
                            model.addConstr(cutVar[ trackCornerIds[vert]] + cutVar[varCount] <= 1);
                        }else{
                            trackCornerIds[vert] = varCount;
                        }
                    }

                    int prevVert = -1;
                    while (vert != end) {// first part always exits
                        // might need a map for patch and vert id to seam adn first or Second
                        if(end == 1463) cout<<" curr at "<<vert<<endl;
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
                        addVarToModel(inverseMap, vert, prevVert, boundaryL[i][(idx + 1) % boundSize], vfAdj, isConstrained, varCount, cutVar,
                                      Fg_pattern, mapVarIdToVertId, seamId[si].first, seamId[si].second,
                                      tailor_lazyness,corner,mapVertAndSeamToVar, vertOther, model, currPattern );
                        if(end == 1463) cout<<" finished add of  "<<vert<<endl;

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
                        if(seamId[si].first >= 0 && inHalfOther){
                            idxOther =  (idxOtherOrig - ( count));
                            if(idxOther>= boundSizeOther) idxOther= idxOther % boundSizeOther;
                            if (idxOther < 0) idxOther += boundSizeOther;
                            if (seamsList[seamId[si].second]->inverted) idxOther = (idxOtherOrig + (count)) % boundSize;
                            vertOther = boundaryL[patch2][idxOther];
                        }
                    }// handle the last, add it to the right sum

                    //  make sure duplicate corner is chosen only once, and if we found the last consider its duplicate with start of this patch
                    if(!inverseMap) {
                        if (trackCornerIds.find(end) != trackCornerIds.end()) {
                            model.addConstr(cutVar[trackCornerIds[end]] + cutVar[varCount] <= 1);
                        } else {
                            trackCornerIds[end] = varCount;
                        }
                    }
                    rSumConstr += cutVar[varCount];
                    totalSum = lSumConstr + cutVar[varCount];
                    if(inverseMap && seamId[si].first >= 0 && inHalfOther){
                        // make sure the last counterpart is actually the least one!
                        vertOther = endOther;
                    }
                    addVarToModel(inverseMap, vert, prevVert , -1, vfAdj, false, varCount, cutVar, Fg_pattern, mapVarIdToVertId, seamId[si].first, seamId[si].second,
                                  tailor_lazyness, true, mapVertAndSeamToVar, vertOther, model, currPattern );
                } else {
//                    cout<<" case inverse direction"<<endl;
                    // iterate in inverse direction
                    seam *seam = seamsList[(seamId[si].second ) * -1 -1];
                    auto startAndPatch = seam -> getStartAndPatch2ForCorres();
                    auto startAndPatchOther = seam -> getStartAndPatch1();
                    int startVal = seam->getStart2();
                    bool inHalf = (fullPatternVertToHalfPatternVert.find(startVal) != fullPatternVertToHalfPatternVert.end());
                    int idx, idxOrig, patch1, vert, end, endOther;
                    if(inHalf){
                            startVal = fullPatternVertToHalfPatternVert[startVal];
                            findIndxInBoundaryloop(boundaryL[i], startVal, idx);
                            idxOrig = idx;
                            patch1 =  patchMapToHalf[startAndPatch.second];
                            vert = boundaryL[patch1][idx];
//                            cout<<startVal<<" start value and start index "<<idx<<endl;
                            end = fullPatternVertToHalfPatternVert[seam->getEndCornerIds().second];
                    }


                    if (patch1 != i) cout << " something is wrong, it should be in patch " << i << " but the seam is from patch  " << startAndPatch.second << endl;

                    int startValOther = seam->getStart1();
                    bool inHalfOther = (fullPatternVertToHalfPatternVert.find(startValOther) != fullPatternVertToHalfPatternVert.end());
                    int idx2, idx2Orig, patch2, otherVert, end2;
                    otherVert = -1;
                    if(inHalfOther){
                        startValOther = fullPatternVertToHalfPatternVert[startValOther];
                        patch2 =  patchMapToHalf[startAndPatchOther.second];

                        findIndxInBoundaryloop(boundaryL[patch2], startValOther, idx2);
                        idx2Orig = idx2;
                        otherVert = boundaryL[patch2][idx2];
                        end2 = fullPatternVertToHalfPatternVert[seam->getEndCornerIds().first];
                    }

                    int nextidx = (idxOrig - (1 + count)) % boundSize;
                    if (nextidx < 0) nextidx += boundSize;
                    if (seam->inverted) nextidx = (idxOrig + (1 + count)) % boundSize;
                    int nextvert = boundaryL[patch1][nextidx];

                    // check if vert is in map, if so then add constraint to connect it
                    if(!inverseMap){
                        if(trackCornerIds.find(vert) != trackCornerIds.end()){
                            model.addConstr(cutVar[ trackCornerIds[vert]] + cutVar[varCount] <= 1);
                        }else{
                            trackCornerIds[vert] = varCount;
                        }
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
                        addVarToModel(inverseMap, vert, prevVert, nextvert, vfAdj, isConstrained, varCount, cutVar, Fg_pattern
                                      , mapVarIdToVertId, seamId[si].first, seamId[si].second, tailor_lazyness, corner,
                                      mapVertAndSeamToVar, otherVert, model, currPattern );
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
                        nextidx = (idxOrig - (1 + count)) % boundSize;
                        if (nextidx < 0) nextidx += boundSize;
                        if (seam->inverted) nextidx = (idxOrig + (1 + count)) % boundSize;
                        nextvert = boundaryL[patch1][nextidx];
                        if(inHalfOther) otherVert =  boundaryL[patch2][(idx2Orig + count) %  boundaryL[patch2].size() ];

                    }// handle the least element and set constraints
                    if(!inverseMap) {
                        if (trackCornerIds.find(end) != trackCornerIds.end()) {
                            model.addConstr(cutVar[trackCornerIds[end]] + cutVar[varCount] <= 1);
                        } else {
                            trackCornerIds[end] = varCount;
                        }
                    }
                    rSumConstr += cutVar[varCount];
                    totalSum = lSumConstr + cutVar[varCount];
                    if(inverseMap && inHalfOther){
                        otherVert = end2;
                    }
                    addVarToModel(inverseMap, vert, prevVert, -1, vfAdj, false, varCount,
                                  cutVar, Fg_pattern, mapVarIdToVertId,
                                  seamId[si].first, seamId[si].second, tailor_lazyness, true,
                                  mapVertAndSeamToVar, otherVert, model, currPattern);

                }

                // set lin expr
                model.addConstr(innerSumConstr <= 1);
                model.addConstr(lSumConstr <= 1);
                model.addConstr(rSumConstr <= 1);
//                model.addConstr(totalSum >= 1);
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
//                    todo this need to be mapped, no? are we sure it exists
                    patch = seamsList[mapVarIdToVertId[i]->seamIdInList]->getStartAndPatch1().second;
                    if(patchMapToHalf.find(patch)== patchMapToHalf.end()){
                        cout<<" it is not in a good patch!!! "<<endl;
                    }else{
                        patch = patchMapToHalf[patch];
                    }

                }else{
                    patch = seamsList[(mapVarIdToVertId[i]->seamIdInList +1)*(-1)]->getStartAndPatch2().second;
                    if(patchMapToHalf.find(patch)== patchMapToHalf.end()){
                        cout<<" it is not in a good patch!!! "<<endl;
                    }else{
                        patch = patchMapToHalf[patch];
                    }
                }

            }else{
                patch = minusOneSeams[mapVarIdToVertId[i]->seamIdInList]->getPatch();
                if(patchMapToHalf.find(patch)== patchMapToHalf.end()){
                    cout<<" it is not in a good patch!!! "<<endl;
                }else{
                    patch = patchMapToHalf[patch];
                }
            }
            cout<<patch<<" patch, Id=  "<< mapVarIdToVertId[i]->vert<<" from seam "<<mapVarIdToVertId[i]->seamType<<" and id "<<mapVarIdToVertId[i]->seamIdInList <<" ";// ALRIGHT BUT NOW WE NEED A SEAM ID

            cutVertEntry* cve = new cutVertEntry ( vert, seamType, seamId, patch);


            cveStartPositionsSet[vert] = cve;
            cve-> leftCorner =  -1;
            cve-> rightCorner = -1;
            int nextVert, prevVert;
            double Stress ;

            bool inverted = false;

            if(seamId<0) inverted = seamsList[(-1)*(seamId+1)]->inverted;
            getStressAtVert( seamType, seamId, vert, prevVert, nextVert, Stress,
                     seamsList, minusOneSeams,  boundaryL, Fg_pattern, currPattern,vfAdj, inverted , fullPatternVertToHalfPatternVert);
            int searchedVert = ( halfPatternVertToFullPatternVert.find(vert) == halfPatternVertToFullPatternVert.end()) ? (-1)*vert : vert;
//            for(int ll=0; ll<cornerVert.size(); ll++){if(cornerVert[ll]==0) cout<<ll<<" "; }
            if( inverseMap|| cornerVert[halfPatternVertToFullPatternVert[vert]]==1 ){
                cout<<"corner   TODO CHECK IF STARTER OR END! "<<endl;
                cve->cornerInitial = vert;
//                // left or right corner?
                int firstInSeam;
                if(seamType == -1 ){
                    firstInSeam = fullPatternVertToHalfPatternVert[ minusOneSeams[seamId]->getStartVert()] ;
                }else if(seamId >= 0){
                    firstInSeam = fullPatternVertToHalfPatternVert[seamsList[mapVarIdToVertId[i]->seamIdInList]->getStart1()];
                }else{
                    firstInSeam = fullPatternVertToHalfPatternVert[seamsList[(mapVarIdToVertId[i]->seamIdInList+1)*(-1)]->getStart2()];
                }
                if(firstInSeam == vert){
                    cve -> startCorner = true;
                    cve->stress = Stress;

                }else{
                    cve-> endCorner = true;
                    cve->stress = Stress;
                }

            }else{
                cve->stress = Stress ;
            }
            cout<<cve->stress<<" the stress there "<<endl<<endl;
            cutPositions.push_back(cve);
        }
    }
    cout<<cutPositions.size()<<" size"<<endl<<endl;
}

void updateStress(vector<cutVertEntry*>& cutPositions, vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                  std::vector<std::vector<int> >& boundaryL, MatrixXi& Fg_pattern, vector<vector<int>> & vfAdj , bool& prioInner,
                  bool& prioOuter, MatrixXd& currPattern, map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& halfPatternVertToFullPatternVert){
    for(int i=0; i< cutPositions.size(); i++){
        cutVertEntry* cve = cutPositions[i];
        int nextVert, prevVert;
        double  Stress;
        bool inverted = false;
        if(cve->seamIdInList < 0){
            inverted = seamsList[(cve->seamIdInList + 1)*(-1)]->inverted;
        }

        getStressAtVert(cve -> seamType, cve -> seamIdInList, cve -> vert, prevVert, nextVert, Stress,
                                     seamsList, minusOneSeams, boundaryL, Fg_pattern, currPattern, vfAdj, inverted, fullPatternVertToHalfPatternVert );

        cve -> stress = Stress;
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

int tearFurther(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,
                vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                 map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL,
                 set<int> & cornerSet, set<int>& handledVerticesSet,  bool& prevFinished, const bool & preferManySmallCuts, const bool & LShapeAllowed,
                 MatrixXd& patternEdgeLengths_orig, MatrixXd& Vg_pattern_orig, MatrixXi& Fg_pattern_orig, bool& prioInner,
                bool& prioOuter, double& setTheresholdlMid, double& setTheresholdBound, map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& halfPatternVertToFullPatternVert,
                map<int, int> & halfPatternFaceToFullPatternFace){
    middleThereshold = setTheresholdlMid;
    boundThereshold = setTheresholdBound;

    cout<<endl<<endl<<"-----------------------"<<endl<<endl;

    int returnPosition = -1;
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    computePerFaceUV(Fg_pattern, Fg_pattern_orig, Vg_pattern_orig, currPattern, halfPatternFaceToFullPatternFace, inverseMapping);

//
    if(prevFinished || preferManySmallCuts){
        // if we want many small cuts we sort always and there is no need to finish a seam before handling the next one!
        cout<<"It's time to sort again"<<endl;
        updateStress( cutPositions, seamsList, minusOneSeams, boundaryL,  Fg_pattern, vfAdj,  prioInner, prioOuter, currPattern, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert);

        sort(cutPositions.begin(), cutPositions.end(), []( cutVertEntry* &a,  cutVertEntry* &b) { return a->stress > b-> stress; });
//        for(int i = 0; i < cutPositions.size(); i++){
//            findCorrespondingCounterCutPosition(cutPositions, i, cutPositions[i], currPattern, Fg_pattern, vfAdj, boundaryL,
//                                                seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet, handledVerticesSet,
//                                                fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert );
//        }
    }
//
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
//
        if(cutPositions[count]->finFlag ){
//            // if it is a corner and it has been released
            if((cutPositions[count]->startCorner || cutPositions[count]->endCorner) &&
                releasedVertNew.find( cutPositions[count]->cornerInitial) != releasedVertNew.end() ){
                int seamPotentiallyReleasedFrom = (cornerToSeams[cutPositions[count]-> cornerInitial][0] == thisSeam) ? cornerToSeams[cutPositions[count]-> cornerInitial][1] : cornerToSeams[cutPositions[count]-> cornerInitial][0];
                if(releasedVertNew[ cutPositions[count]->cornerInitial].end() != std::find( releasedVertNew[ cutPositions[count]->cornerInitial].begin(), releasedVertNew[ cutPositions[count]->cornerInitial].end(), seamPotentiallyReleasedFrom)){
                    if(sym && inverseMapping) {
                        parallel = openParallelPositionInvSym(cutPositions[count]->cornerInitial,
                                                              seamPotentiallyReleasedFrom, seamsList, cutPositions,
                                                              fullPatternVertToHalfPatternVert,
                                                              halfPatternVertToFullPatternVert);
                    }else parallel = openParallelPosition(cutPositions[count]-> cornerInitial, seamPotentiallyReleasedFrom,
                                                          seamsList, cutPositions,fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert );//releasedVertNew[cutPositions[count]->cornerInitial]
                    if(parallel >= 0) parallelFinFlag = cutPositions[parallel]->finFlag;
                    if(parallelFinFlag){
                        i--;
                        count ++;
                    }else{
                        returnPosition = cutPositions[parallel] ->vert;
                        cout<<"only handling paralell:"<<endl;
                        splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                                           minusOneSeams, releasedVert, toPattern_boundaryVerticesSet, cornerSet,
                                           handledVerticesSet, LShapeAllowed, Vg_pattern_orig, Fg_pattern_orig,
                                           fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace);
                        prevFinished = cutPositions[parallel]->finFlag;
                    }
                }else{
                    i--;
                    count ++;
                }
//
            }else{
                i--;
                count ++;
            }
//
        }else{
//
            cout<<endl<<endl << cutPositions[count]->vert<<" vertex up next from seam  "<<cutPositions[count]->seamType<<" "<<cutPositions[count]->seamIdInList<<" the stress here is "<<cutPositions[count]->stress<<endl;
            returnPosition = cutPositions[count] ->vert;
            splitVertexFromCVE(cutPositions[count], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList, minusOneSeams, releasedVert,
                                   toPattern_boundaryVerticesSet, cornerSet, handledVerticesSet, LShapeAllowed, Vg_pattern_orig, Fg_pattern_orig,
                               fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace);
            cout<<"finished ? "<<cutPositions[count]->finFlag<<endl<<endl;
            prevFinished = cutPositions[count] -> finFlag;

            parallel = -1;
            // if it is a corner and it has been released
            if(cutPositions[count]->startCorner || cutPositions[count]->endCorner){
//                // once we finished cutting one side, check if it was a side opening. If so we can go on with the other side
//                //this is the released seam, so the other one is the one we have in common
                int seamPotentiallyReleasedFrom = (cornerToSeams[cutPositions[count]-> cornerInitial][0] == thisSeam) ? cornerToSeams[cutPositions[count]-> cornerInitial][1] : cornerToSeams[cutPositions[count]-> cornerInitial][0];
                if(inverseMapping && sym) {
                    parallel = openParallelPositionInvSym(cutPositions[count]->cornerInitial,
                                                          seamPotentiallyReleasedFrom, seamsList, cutPositions,
                                                          fullPatternVertToHalfPatternVert,
                                                          halfPatternVertToFullPatternVert);
                }else {
                    // open the other side of a released seam .
//                // attention this is not the other side of the seam but the same 3D corner of a different patch.
                    parallel = openParallelPosition(cutPositions[count]->cornerInitial, seamPotentiallyReleasedFrom,
                                                    seamsList, cutPositions, fullPatternVertToHalfPatternVert,
                                                    halfPatternVertToFullPatternVert);
                }
            }
            if(parallel<0) {
                cout<<"no proper parallel found, only avalilable for boundary cuts "<<endl;
                continue;
            }else if(!cutPositions[parallel]->finFlag){
                cout<<"split parallel "<<endl;

                splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                                   minusOneSeams, releasedVert, toPattern_boundaryVerticesSet, cornerSet,
                                   handledVerticesSet, LShapeAllowed, Vg_pattern_orig, Fg_pattern_orig, fullPatternVertToHalfPatternVert,
                                   halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace);
                cout<<"finished  p ? "<<cutPositions[parallel]->finFlag<<endl;
                prevFinished = (prevFinished && cutPositions[parallel] -> finFlag);

            }else{
                cout<<"parallel already finished"<<endl;
            }
        }

    }
    cout<<prevFinished<<" --------------------"<<endl;
    cout<<currPattern.rows()<<" in end  "<<endl;

    return returnPosition;
}

int computeTear(bool inverseMap, Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern_curr, MatrixXd& patternlengthsOrig,
                 vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams, std::vector<std::vector<int> >& boundaryL, bool & finished,
                  std::vector<std::vector<std::pair<int, int>>>& cornersPerBoundary, map<int, vector<pair<int, int>>>& seamIdPerCorner,
                 VectorXd& cornerVert, vector<cutVertEntry*>& cutPositions, map<int, pair<int, int>> & releasedVert,
                 set<int>& toPattern_boundaryVerticesSet, set<int> & cornerSet, set<int>& handledVerticesSet ,
                 bool& prevFinished,const bool & LShapeAllowed,bool& prioInner,
                 bool& prioOuter, double tailor_lazyness, const MatrixXi& mapFromFg, double& setTheresholdlMid, double& setTheresholdBound,
                 map<int, int> & fullPatternVertToHalfPatternVert, map<int, int>& halfPatternVertToFullPatternVert,
                 map<int, int> & halfPatternFaceToFullPatternFace,
                 bool& symetry )
                 {
     inverseMapping = inverseMap;
     sym = symetry;
    middleThereshold = setTheresholdlMid;
    boundThereshold = setTheresholdBound;

    std::vector<std::vector<int> > boundaryLnew;
    igl::boundary_loop(Fg_pattern_curr, boundaryLnew);
    boundaryL.clear();
    boundaryL = boundaryLnew;
    vector<vector<pair<int, int>>> cornersPerBoundaryHalf;
    for(int i=0; i<cornersPerBoundary.size(); i++){
       for(int j=0; j<cornersPerBoundary[i].size(); j++){
          int ver = cornersPerBoundary[i][j].first;
//          cout<<i<<" with vert "<<ver<<endl;
            if(ver<0 && fullPatternVertToHalfPatternVert.find(ver) == fullPatternVertToHalfPatternVert.end()){
                fullPatternVertToHalfPatternVert[ver]= -1*ver;
            }
          if(fullPatternVertToHalfPatternVert.find(ver) == fullPatternVertToHalfPatternVert.end()){
              continue;
          }
           int newVer = fullPatternVertToHalfPatternVert[ver];
          for(int ii=0; ii < boundaryL.size(); ii++){
              for(int jj= 0; jj<boundaryL[ii].size(); jj++){
                  if(boundaryL[ii][jj] == newVer){
                      patchMapToHalf[i]= ii;

                  }
              }
          }
       }
    }
    for(int i=0; i<cornersPerBoundary.size(); i++){
        if(patchMapToHalf.find(i)!= patchMapToHalf.end()){
            cornersPerBoundaryHalf.push_back(cornersPerBoundary[i]);
        }
    }
    if(symetry){
        cornersPerBoundary.clear();
        cornersPerBoundary = cornersPerBoundaryHalf;

    }
    cout<<"************ all the half corners *******************"<<endl;

    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern_curr, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern_curr, lengthsCurr);

    lengthsOrig =  patternlengthsOrig;

    for(int i =0; i<seamsList.size(); i++){
        seam* currSeam = seamsList[i];
        int start1 = currSeam -> getStart1();
        int start2 = currSeam ->getStart2();
        int end1 = currSeam-> getEndCornerIds().first;
        int end2 = currSeam-> getEndCornerIds().second;

        if(patchMapToHalf.find(currSeam->getPatch1()) != patchMapToHalf.end()){
            addToMapIfNotExisting(fullPatternVertToHalfPatternVert[start1], i);

        }
        if(patchMapToHalf.find(currSeam->getPatch2()) != patchMapToHalf.end()){
            addToMapIfNotExisting(fullPatternVertToHalfPatternVert[start2], i);
        }
        if(patchMapToHalf.find(currSeam->getPatch1()) != patchMapToHalf.end()){
            addToMapIfNotExisting(fullPatternVertToHalfPatternVert[end1], i);

        }
        if(patchMapToHalf.find(currSeam->getPatch2()) != patchMapToHalf.end()){
            addToMapIfNotExisting(fullPatternVertToHalfPatternVert[end2], i);

        }
    }
    for(int i =0; i<minusOneSeams.size(); i++){
       minusOneSeam* currSeam = minusOneSeams[i];
       int start = currSeam-> getStartVert();
       int end = currSeam -> getEndVert();

        if(patchMapToHalf.find(currSeam->getPatch()) != patchMapToHalf.end()){
            addToMapIfNotExisting( fullPatternVertToHalfPatternVert[start], -i-1);
            addToMapIfNotExisting( fullPatternVertToHalfPatternVert[end], -i-1);
        }
    }

    // first we need to know where to tear, set up LP for this
    // information we need: stress. For stress, we need lengths old and ned
    // for lengths old and new we need which edge of face
    // for which edge of face we need face
    // for face we need vertices
    // for vertices we need boundary loop

    double minConstrained = 0.25;
    setLP(inverseMap, boundaryL, vfAdj, Fg_pattern_curr, cornersPerBoundary, seamIdPerCorner,
          seamsList, minusOneSeams, tailor_lazyness, minConstrained, cutPositions,
          cornerVert, currPattern, LShapeAllowed, fromPattern, mapFromFg, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace);

    cout<<"finished set lp "<<endl;
    //update with the preferences before sorting
    updateStress( cutPositions, seamsList, minusOneSeams, boundaryL,  Fg_pattern_curr, vfAdj,  prioInner, prioOuter, currPattern, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert);

    // sort and check if handled already
    sort(cutPositions.begin(), cutPositions.end(), []( cutVertEntry* &a,  cutVertEntry* &b) { return a->stress > b-> stress; });


//    for(int i = 0; i < cutPositions.size(); i++){
//        findCorrespondingCounterCutPosition(cutPositions, i, cutPositions[i], currPattern, Fg_pattern_curr, vfAdj, boundaryL,
//                                            seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet, handledVerticesSet, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert);
//
//    }
    cout<<" found all counter positions "<<endl;
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
        cout<<cutPositions[count]->stress<<" curr stress "<<cutPositions[count]->vert<<endl;
        int currVert = cutPositions[count]->vert;
        returnPosition = cutPositions[count]->vert;

        splitVertexFromCVE(cutPositions[count], currPattern, Fg_pattern_curr, vfAdj, boundaryL, seamsList,
                           minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,
                           cornerSet, handledVerticesSet, LShapeAllowed, fromPattern, mapFromFg, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace);

        prevFinished = cutPositions[count] -> finFlag;

        if( releasedVertNew.find(currVert) != releasedVertNew.end()){
            // dann können müssen wir ja auch die passende andere seite des cuts öffnen
            // achtung, das ist nicht die corresponding seam
            // es kann aber noch nicht von einer anderen Seite geöffnett sein. daher nimm idx [0] von den released vert
            int parallel;
            if(symetry && inverseMap){
                parallel = openParallelPositionInvSym(cutPositions[count]-> cornerInitial, releasedVertNew[currVert][0], seamsList, cutPositions, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert);

            }else parallel = openParallelPosition(cutPositions[count]-> cornerInitial, releasedVertNew[currVert][0], seamsList, cutPositions, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert);
            if(parallel<0) {cout<<"no proper parallel found"<<endl;continue;}else{
                cout<<endl<<"Working on parallel: "<<parallel<<endl<<endl;
            }

            splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern_curr, vfAdj, boundaryL, seamsList,
                               minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,
                               cornerSet, handledVerticesSet, LShapeAllowed, fromPattern, mapFromFg,
                               fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace);
            prevFinished = (prevFinished && cutPositions[parallel] -> finFlag);

        }

    }

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
// todo this could be done only once ;-)
void fillMatrixWithBoundaryVert(const vector<int>& boundary, const int& start, const int& end, const MatrixXd& mapToVg, MatrixXd& Vg_seam1to, bool inverted, const map<int, int>& fullPatternVertToHalfPatternVert ){

    int countLen = 2;
    int startIdx = 0;
    int boundLen = boundary.size();

    while(boundary[startIdx] != start && startIdx <= boundLen){
        startIdx++;
        startIdx = startIdx % boundLen ;
    }
    if(boundary[startIdx] != start ) cout<<"START NOT FOUND ERROR "<<start<<endl;
    int endIdx = (inverted)? startIdx-1: startIdx+1;
    if(endIdx <0) endIdx += boundLen;

    while(boundary[endIdx] != end && endIdx != startIdx  ){// messy? && (fullPatternVertToHalfPatternVert.find(boundary[endIdx]) != fullPatternVertToHalfPatternVert.end())
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
// set up with to Pattern
void setUpMap( const std::vector<std::vector<int> >& boundaryL,map<int,int> & fullPatternVertToHalfPatternVert){
    int count = 0;
    for(int i=0; i<boundaryL.size(); i++){
        int found = false;
        for(int j=0; j< boundaryL[i].size(); j++){
            if(found) break;
            if(fullPatternVertToHalfPatternVert.find(boundaryL[i][j])!= fullPatternVertToHalfPatternVert.end()){
                patchMapToHalfInverse[i]= count;
                count ++;
                found = true;
            }
        }

    }
    cout<<(fullPatternVertToHalfPatternVert.find(1512)== fullPatternVertToHalfPatternVert.end())<<" found? "<<endl;
}
void projectBackOnBoundary(const MatrixXd & mapToVg, MatrixXd& p, const vector<seam*>& seamsList
                           ,const vector<minusOneSeam*> & minusOneSeams,
                           const std::vector<std::vector<int> >& boundaryL_toPattern, const std::vector<std::vector<int> >& boundaryL,
                           map<int, pair<int, int>> & releasedVert ,bool inverseMap, map<int,int> & fromtoToVertMapIfSplit,map<int,int> & seamFullHalf ){

    int numSeams = seamsList.size();
    int count =0;
    for (int j = 0; j<numSeams; j++){

        seam* currSeam  = seamsList[j];
        auto stP1= currSeam-> getStartAndPatch1();
        auto stP2 =  currSeam -> getStartAndPatch2ForCorres(); // attention this is with respect to the original pattern
//
        MatrixXd Vg_seam1to, Vg_seam2to;
        // build the structure for closest search
        if(seamFullHalf.find(currSeam-> patch1startCornerIdOld) != seamFullHalf.end()  ||
                seamFullHalf.find(currSeam-> patch1endCornerIdOld) != seamFullHalf.end()){
            int toPatch = currSeam->getStartAndPatch1().second;
            int start = currSeam-> patch1startCornerIdOld;
            int end = currSeam-> patch1endCornerIdOld;
            fillMatrixWithBoundaryVert(boundaryL_toPattern[toPatch], start ,end, mapToVg, Vg_seam1to, false , fromtoToVertMapIfSplit);

        }
        if(seamFullHalf.find(currSeam-> patch2startCornerIdOld) != seamFullHalf.end()  ||
                seamFullHalf.find(currSeam-> patch2endCornerIdOld) != seamFullHalf.end() ){

            int toPatch = currSeam->getStartAndPatch2().second;
            int start =  currSeam-> patch2startCornerIdOld;
            int end =  currSeam-> patch2endCornerIdOld;
            fillMatrixWithBoundaryVert(boundaryL_toPattern[toPatch], start, end, mapToVg, Vg_seam2to, true, fromtoToVertMapIfSplit );

        }

        bool shoulBeLeft =true; // for 2 case
        // for each interior (=not corner) vertex of the new boundary we need to find the closest position on the polyline and map it there
        pair<int, int> ends = currSeam->getEndCornerIds();


        if(seamFullHalf.find(currSeam->getStart1()) != seamFullHalf.end() ||
                seamFullHalf.find(ends.first) != seamFullHalf.end() ){// only if at least one of them exists on the half pattern iit makes sense to iterate

            int patchUsed = patchMapToHalfInverse[stP1.second];
            if(!inverseMap) patchUsed = stP1.second;
            int bsize = boundaryL[patchUsed].size();
            int next = (inverseMap)? seamFullHalf[currSeam -> getStart1()] : currSeam -> getStart1();
            int nextIdx = 0 ;
            while(boundaryL[patchUsed][nextIdx] != next && nextIdx < bsize){
                nextIdx ++;
            }
            if(boundaryL[patchUsed][nextIdx] != next){
                cout<<patchUsed<<" patch, ,searched for "<<next<<endl;
                cout<<"PROJECTION ERROR we dont find the index "<<endl;
                for(auto it:boundaryL[patchUsed] ){cout<<it<<" "; }cout<<endl;
            }
//
            int endsFirst = (inverseMap)? seamFullHalf[ends.first] : ends.first;
            int searchNExt;
            while( next!= endsFirst ){
                // it is not released, project on boundary
                 searchNExt = (inverseMap) ? next : seamFullHalf[next];
                if(releasedVert.find(searchNExt) == releasedVert.end()){
                    updatePositionToIntersection( p, searchNExt,Vg_seam1to, true);
                }
                    // else it is released somehow. But from which seam? If it is released from another seam then pull it to this boundary still
                else if( std::find(releasedVertNew[searchNExt].begin(),
                                   releasedVertNew[searchNExt].end(), j) == releasedVertNew[searchNExt].end()){
                    updatePositionToIntersection( p, searchNExt,Vg_seam1to, true);
                }
                nextIdx++;
                next = boundaryL[patchUsed][(nextIdx) % bsize];
            }
            searchNExt = (inverseMap) ? next : seamFullHalf[next];

            // the last corner. Again if it is constrained from another side pull it to boundary, else ignore since handled by corner
            if(releasedVert.find(searchNExt) == releasedVert.end()){
                updatePositionToIntersection( p, searchNExt,Vg_seam1to, true);

            } else if (std::find(releasedVertNew[searchNExt].begin(), releasedVertNew[searchNExt].end(), j) == releasedVertNew[searchNExt].end()){
                updatePositionToIntersection( p, searchNExt,Vg_seam1to, true);
            }
        }

        /**********  second side  *********/
        if(seamFullHalf.find(currSeam->getStart2()) != seamFullHalf.end() ||
                seamFullHalf.find(ends.second) != seamFullHalf.end() ){
            int patchUsed = (inverseMap)? patchMapToHalfInverse[stP2.second] : stP2.second;
            int bsize = boundaryL[patchUsed].size();
            int nextIdx = 0;
            int next = (inverseMap)? seamFullHalf[currSeam -> getStart2()]: currSeam -> getStart2();
            while (boundaryL[patchUsed][nextIdx] != next && nextIdx < bsize ){
                nextIdx ++;
            }
            if(boundaryL[patchUsed][nextIdx] != next){
                cout<<"PROJECTION ERROR 2 we dont find the index "<<next<<" should be on patch "<<stP2.second<<endl;
            }
            int endSecond = (inverseMap) ? seamFullHalf[ends.second]: ends.second;
            int nextSeach;
            while( next!= endSecond ){
                // general case an interior vertex , if it is not constrained pull it to boundary
                 nextSeach = (inverseMap)? next: seamFullHalf[next];
                if(releasedVert.find(nextSeach) == releasedVert.end() ){
                    updatePositionToIntersection( p, nextSeach,Vg_seam2to, shoulBeLeft);
                } else if( std::find(releasedVertNew[nextSeach].begin(),
                                     releasedVertNew[nextSeach].end(), j) == releasedVertNew[nextSeach].end()){
                    updatePositionToIntersection( p,nextSeach ,Vg_seam2to, shoulBeLeft);
                }
                nextIdx -=1;// (nextidx - i2) % (bsize);

                if(nextIdx < 0) {nextIdx += bsize;}
                next = boundaryL[patchUsed][nextIdx];
            }
            nextSeach = (inverseMap)? next: seamFullHalf[next];
            if(releasedVert.find(nextSeach) == releasedVert.end()){
                updatePositionToIntersection( p, nextSeach,Vg_seam2to, shoulBeLeft);
            }else
                if (std::find(releasedVertNew[nextSeach].begin(),
                                releasedVertNew[nextSeach].end(), j) == releasedVertNew[nextSeach].end()){
                updatePositionToIntersection( p, nextSeach,Vg_seam2to, shoulBeLeft);
            }

        }

        // also project all duplicates of interior cut vertices
        for(const auto & addedVert : currSeam->duplicates){
            int avs = addedVert.second;
            if(releasedVert.find(avs) == releasedVert.end() ){
                updatePositionToIntersection(p, avs, Vg_seam1to, true);
            } else if( std::find(releasedVertNew[avs].begin(), releasedVertNew[avs].end(), j) == releasedVertNew[avs].end()){
                updatePositionToIntersection(p, avs, Vg_seam1to, true);
            }
        }

        for(const auto & addedVert : currSeam->duplicates2){

            int avs = addedVert.second;
            if(releasedVert.find(addedVert.second) == releasedVert.end() ){
                updatePositionToIntersection(p, addedVert.second, Vg_seam2to, shoulBeLeft);
            } else if( std::find(releasedVertNew[avs].begin(), releasedVertNew[avs].end(), j) == releasedVertNew[avs].end()){
                updatePositionToIntersection(p, avs, Vg_seam2to, shoulBeLeft);
            }
        }
    }
    for(int j = 0; j < minusOneSeams.size(); j++){//

        minusOneSeam* currSeam  = minusOneSeams[j];
        int patch = currSeam -> getPatch();
        int startVert = currSeam -> getStartVert();
        int endVert =  currSeam -> getEndVert();

        // build the structure for closest search
        MatrixXd Vg_seamto;
        if(seamFullHalf.find(currSeam->startVertOld) != seamFullHalf.end() ||
                seamFullHalf.find(currSeam->endVertOld) != seamFullHalf.end()) {

            fillMatrixWithBoundaryVert(boundaryL_toPattern[patch], currSeam->startVertOld, currSeam->endVertOld,
                                       mapToVg, Vg_seamto, false, fromtoToVertMapIfSplit );

            int patchUsed = (inverseMap)? patchMapToHalfInverse[patch] : patch;//
            int boundLen = boundaryL[patchUsed].size();
            int startidx = 0;
            int next = (inverseMap)? seamFullHalf[startVert]: startVert;
            endVert =  (inverseMap)? seamFullHalf[ currSeam -> getEndVert()] : currSeam -> getEndVert();
            while (boundaryL[patchUsed][startidx] != next){
                startidx++;
            }
            int nextSearch ;
            while(next != endVert ){
                // general case, it is not released hence pull it to the boundary
                nextSearch= (inverseMap) ? next : seamFullHalf[next];
                if(releasedVert.find(nextSearch) == releasedVert.end()){
                    updatePositionToIntersection( p,nextSearch,Vg_seamto , true);
                }else if(std::find(releasedVertNew[nextSearch].begin(), releasedVertNew[nextSearch].end(),
                                   (-1)*(j+1)) == releasedVertNew[nextSearch].end()){
                    // it is released but not from this seam,thus it has to stay on the projection
                    updatePositionToIntersection( p, nextSearch,Vg_seamto, true);
                }
                startidx++;
                startidx = startidx % boundLen;
                next = boundaryL[patchUsed][ startidx];
            }// it is released for another side hence we have to pull it to our side
            nextSearch= (inverseMap) ? next : seamFullHalf[next];
            if(releasedVert.find(nextSearch) != releasedVert.end() && std::find(releasedVertNew[nextSearch].begin(),releasedVertNew[nextSearch].end(),  (-1)*(j+1)) == releasedVertNew[nextSearch].end()){
                updatePositionToIntersection( p, nextSearch,Vg_seamto, true);

            }


            // also map all projections
            for(const auto & addedVert : currSeam -> duplicates){
                int avs = addedVert.second;
                if(releasedVert.find(avs) == releasedVert.end() ){
                    updatePositionToIntersection(p, avs, Vg_seamto, true);
                } else if( std::find(releasedVertNew[next].begin(), releasedVertNew[next].end(), j) == releasedVertNew[next].end()){
                    updatePositionToIntersection(p, avs, Vg_seamto, true);
                }

            }
        }




    }
}

void updatePatchId(vector<cutVertEntry*>& cutPositions, const std::vector<std::vector<int> >& boundaryLnew, vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams,
                   map<int, int >& fullPatternVertToHalfPatternVert){
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
        if(fullPatternVertToHalfPatternVert.find(seamsList[i]-> getStart1()) != fullPatternVertToHalfPatternVert.end()){
            seamsList[i]-> updatePatch1(mapVertToNewPatch[fullPatternVertToHalfPatternVert[seamsList[i]-> getStart1()]] );
            seamsList[i]->seamSplit1 = (mapVertToNewPatch[fullPatternVertToHalfPatternVert[seamsList[i]-> getStart1()]] ==
                    mapVertToNewPatch[fullPatternVertToHalfPatternVert[seamsList[i]-> getEndCornerIds().first]]);

        }
        if(fullPatternVertToHalfPatternVert.find(seamsList[i]-> getStart2()) != fullPatternVertToHalfPatternVert.end()){
            seamsList[i]-> updatePatch2(mapVertToNewPatch[fullPatternVertToHalfPatternVert[seamsList[i]-> getStart2()]] );
            seamsList[i]->seamSplit2 = (mapVertToNewPatch[fullPatternVertToHalfPatternVert[seamsList[i]-> getStart2()]] ==
                    mapVertToNewPatch[fullPatternVertToHalfPatternVert[seamsList[i]-> getEndCornerIds().second]]);

        }
    }

    for(int i=0; i < minusOneSeams.size(); i++){
        if(fullPatternVertToHalfPatternVert.find(minusOneSeams[i]-> getStartVert()) != fullPatternVertToHalfPatternVert.end()){
            minusOneSeams[i]-> updatePatch(mapVertToNewPatch[fullPatternVertToHalfPatternVert[minusOneSeams[i]-> getStartVert()]]);
            minusOneSeams[i]->seamSplit = (mapVertToNewPatch[fullPatternVertToHalfPatternVert[minusOneSeams[i]-> getStartVert()]] ==
                                           mapVertToNewPatch[fullPatternVertToHalfPatternVert[minusOneSeams[i]-> getEndVert()]]);

        }

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
