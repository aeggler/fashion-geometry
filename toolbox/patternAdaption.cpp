//
// Created by Anna Maria Eggler on 08.12.22.
//
#include <Eigen/Dense>
#include "patternAdaption.h"
#include <iostream>
#include "seam.h"
#include <igl/edge_lengths.h>
#include <igl/adjacency_list.h>
#include "igl/boundary_loop.h"
//
#include "/Library/gurobi1000/macos_universal2/include/gurobi_c++.h"
#include "adjacency.h"
#include <igl/HalfEdgeIterator.h>

using namespace std;
using namespace Eigen;

class Node{
public:
    int vertId;
    Node *next;
    Node *prev;
    int faceId;
    int edgeOfFaceId;
    int duplicate;
    double stretchLeft;
};

MatrixXd lengthsOrig;
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

// push item to the back of the list
void push(Node** tail_ref, int vert_Id, int dupl, int face_id ,int whichEdge,  double stretch_Left, int& listLength){

    Node* new_node = new Node();
    new_node -> prev = (*tail_ref);
    new_node -> next = NULL;

    new_node -> vertId = vert_Id;
    new_node -> duplicate = dupl;
    new_node -> faceId = face_id;
    new_node -> edgeOfFaceId = whichEdge;
    new_node -> stretchLeft = stretch_Left;

    if((*tail_ref)!= NULL){
        (*tail_ref)-> next = new_node;
    }
    (*tail_ref) = new_node;

    if(new_node-> prev != NULL){
//        new_node -> stretchRight  = new_node-> prev->stretchLeft;
//        new_node -> faceIdRight = new_node -> prev -> faceId;
//        new_node-> edgeOfFaceIdRight = new_node -> prev -> edgeOfFaceId;
    }
    listLength++;
}
bool canBeSplit(int vertId,  vector<vector<int> >& vfAdj){
    return vfAdj[vertId].size()>1;

}
int printList(Node * node, int length,  vector<vector<int> >& vfAdj){
    // this is the thereshold
    double maxStretch = 0.0001;
    int idxMax = -1;
    int currIdx = 0;
//    while (node != NULL){
    while(currIdx < length){
//        cout<<node->vertId<<" "<< node->stretchLeft <<" and face to the left  "<<node->faceId<<" and id "<<node-> next->vertId <<endl;
        if(node->stretchLeft + node->prev->stretchLeft > maxStretch && canBeSplit(node->vertId, vfAdj)){
            maxStretch = node->stretchLeft + node->prev->stretchLeft;
            idxMax = currIdx;
        }
        currIdx++;
        node = node -> next;
    }
    if (idxMax == -1){
        cout<<" nothing to tear! "<<endl;
    }
    return idxMax;
}
void printDataOfIdx(Node * node, int whichTear){
    int currIdx = 0;
    while (currIdx<whichTear){
        node = node-> next;
        currIdx++;
    }
    cout<<node->vertId<<" "<< node->stretchLeft + node->prev->stretchLeft<<endl;
    cout<<node -> faceId<<" "<< node -> edgeOfFaceId<<" "<< node -> duplicate<<" the metadata"<<endl;
    cout<<  node -> prev->faceId <<" "<< node-> prev->edgeOfFaceId <<" the right side"<<endl;
}
Node getDataOfIdx(Node * node, int whichTear){
    int currIdx = 0;
    while (currIdx<whichTear){
        node = node-> next;
        currIdx++;
    }
    return *node;

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
void addToDulicateList( cutVertEntry*& cve, vector<seam*>& seamsList,  vector<minusOneSeam*> & minusOneSeams, int newVertIdx){

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
bool isRight(Vector3d a,Vector3d b,Vector3d c ){
    // right gets  newVertIdx
    if(((b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0))) ==0 )cout<<"zero"<<endl;
    cout<<((b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0)))<<" is right if smaller 0 "<<endl;
    // TODO HANDLE THIS CASE PROPERLY!!
    return ((b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0))) < 0;

}

void splitVertexFromCVE( cutVertEntry*& cve, MatrixXd& Vg, MatrixXi& Fg, vector<vector<int> >& vfAdj,
                         std::vector<std::vector<int> >& boundaryL,  vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams,
                         map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  MatrixXd& lengthsCurr, MatrixXi& Fg_pattern, set<int> & cornerSet){
    if(cve-> finFlag) {
        cout<<cve->vert<<" done already "<<endl;
        return;
    }

    cout<<"handling vert: "<<cve->vert<<endl;
    auto boundary = boundaryL[cve->patch];

    if(cve->vert == 449){
        for(int j=0; j<boundary.size(); j++){
            cout<<boundary[j]<<" ";
        }
    }
    double eps = 0.01;
    vector<vector<int>> vvAdj;
    igl::adjacency_list(Fg,vvAdj);
    vector<int> adjacentFaces = vfAdj[cve->vert];

    int newVertIdx = Vg.rows();
    MatrixXd newVg (newVertIdx+1, 3);
    newVg.block(0,0, newVertIdx, 3)= Vg;
    cout<<"searching index "<<endl;

    int idx = 0;
    int leftId, rightId;
    while (boundary[idx] != cve->vert){
        idx++;
    }
    cout<<"found index "<<idx<<endl;
    leftId = (idx+1) % boundary.size();
    rightId = (idx -1 );
    cout<<leftId<<" left and right idx "<<rightId<<endl;


    cout<<boundary[leftId]<<" left and right "<<boundary[rightId]<<endl;
    if(rightId<0) rightId += boundary.size();
    int leftFaceId = adjacentFaceToEdge(boundaryL[cve->patch][leftId], cve-> vert, -1, vfAdj );
    int rightFaceId = adjacentFaceToEdge(boundaryL[cve->patch][rightId], cve-> vert, -1, vfAdj );
    if(leftFaceId ==-1 || rightFaceId ==-1){
        cve->finFlag=true;
        cout<<" something went wrong, we have no neighbor faces. Maybe cut too fast for simulation"<<endl;
        return;
    }

    cout<<leftFaceId<<" left and right face index "<<rightFaceId<<endl;
    int whichEdgeLeft = findWhichEdgeOfFace(leftFaceId, cve->vert, boundaryL[cve->patch][leftId], Fg_pattern);
    int whichEdgeRight= findWhichEdgeOfFace(rightFaceId, cve->vert, boundaryL[cve->patch][rightId], Fg_pattern);

    double w_init = lengthsCurr(leftFaceId, whichEdgeLeft)/lengthsOrig(leftFaceId, whichEdgeLeft);
    w_init +=  lengthsCurr(rightFaceId, whichEdgeRight)/lengthsOrig(rightFaceId, whichEdgeRight);
    w_init /= 2;
    double thereshold = 1.005;
    if(w_init < thereshold){
        cve-> finFlag = true;
        cout<<"stopping now "<<endl;
        return;
    }


    Vector3d toLeft = Vg.row(boundaryL[cve->patch][leftId])- Vg.row(cve->vert);
    Vector3d toRight = Vg.row(boundaryL[cve->patch][rightId])- Vg.row(cve->vert);

    // if it's a corner all get the new index
    if(cve->startCorner|| cve-> endCorner ){
        cout<<" extra case"<<endl;
        // set of extra cases
        pair<int, int> valPair = make_pair(cve->seamType, cve ->seamIdInList);
        releasedVert[cve->vert]= valPair;
        newVg.row(newVertIdx)= Vg.row(cve->vert);

        cve->continuedCorner = true;
        if(cve->startCorner){// does not matter if it is a starter or not
            cout<<"starter"<<endl;
//            newVg.row(cve->vert) += (eps * toLeft).transpose();
            if(cve-> seamIdInList>=0){
                // we count pos hence the previous is negative
                cve->vert = boundary[rightId];
            }else{
                cve->vert = boundary[leftId];
            }


        }else{
            cout<<"ending"<<endl;
            if(cve-> seamIdInList>=0){
                // we count pos hence the next +1
                cve->vert = boundary[leftId];
            }else{
                cve->vert = boundary[rightId] ;
            }
//            if first stop here

        }
        cve->finFlag = (cornerSet.find(cve->vert) != cornerSet.end()); //if it is a corner we are done
        Vg.resize(Vg.rows()+1, 3);
        Vg= newVg;
        cout<<" Vg.rows "<<Vg.rows()<<endl;
        return;

    }
    // if it's a bridge there is no next and we set fin flag
    if(cve-> bridgeFlag){
        // this is the final cut, we are done after: should look like  ><
        cout<<"cutting the bridge, but better dont for now , it messes up the patches "<<endl;
       // updateWithNewId(Fg, cve->vert, rightFaceId, newVertIdx);

        cve-> finFlag = true;
        return;
    }

//    https://github.com/libigl/libigl/blob/main/include/igl/internal_angles.cpp
//    auto s = toLeft.normalized().cross(toRight.normalized()).norm();
//    auto c = toLeft.normalized().dot(toRight.normalized());
//    double angle = atan2(s, c);
//    angle = 2*M_PI - angle;
//    double deg = angle*180/M_PI;
//    double cosbeta = cos(angle/2);
//    double sinbeta = sin(angle /2);
//    // we rotate midVec by half to get to the middle
//    Vector3d leftRot = toLeft.normalized();
//    leftRot(0)= cosbeta * toLeft.normalized()(0) - sinbeta * toLeft.normalized()(1);
//    leftRot(1)=  sinbeta * toLeft.normalized()(0) + cosbeta * toLeft.normalized()(1);

    Vector3d midVec;// = leftRot;
    Vector3d midVect = Vg.row(boundaryL[cve->patch][leftId])- Vg.row(boundaryL[cve->patch][rightId]);
    if(!cve->levelOne){
        midVect = cve->leftdirection - cve->rightdirection;
    }
    midVect= midVect.normalized();
    midVec(0)= -midVect(1);
    midVec(1) = midVect(0);

    cout<<" midvec "<<midVec.transpose()<<endl;


    // all other normal cases
    double dist = std::numeric_limits<double>::max();
    int idxofClosest = -1;
    for(int i=0; i<vvAdj[cve -> vert ].size(); i++){
        int adjVert = vvAdj[cve -> vert][i];
        if(adjVert== rightId || adjVert == leftId) continue; // we want a middle one

        Vector3d edgeVec = Vg.row(adjVert)- Vg.row(cve -> vert);
        edgeVec= edgeVec.normalized();
        // both have unit distance, so as a measure we can take the distance form another
        if((midVec - edgeVec).norm()<dist){
            idxofClosest= i;
            dist = (midVec - edgeVec).norm();
        }
    }
    // is the new to be inserted vertex
    int insertIdx =  vvAdj[cve -> vert][idxofClosest];
    MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(Fg, TT, TTi);
    std::pair<int, int>  faces;
    adjacentFacesToEdge(cve->vert, insertIdx, vfAdj, faces );
    cout<<insertIdx<<" the insert idx "<<faces.first<<" and second "<<faces.second<<endl;
    int helperWhich = -1;
    if(Fg(faces.first, 0)== cve-> vert){
        helperWhich = 0;
    }else if (Fg(faces.first, 1)== cve-> vert){
        helperWhich = 1;
    }else{
        helperWhich = 2;
    }
    cout<<Fg.row(faces.first)<<" found and identified index "<<helperWhich<<endl;
//    igl::HalfEdgeIterator<MatrixXi, MatrixXi, MatrixXi>;
    igl::HalfEdgeIterator <MatrixXi, MatrixXi, MatrixXi>hei (Fg, TT, TTi, faces.first, helperWhich, false);
    igl::HalfEdgeIterator <MatrixXi, MatrixXi, MatrixXi>hei2 (Fg, TT, TTi, faces.second, helperWhich, false);
    int lastF = hei.Fi();
    int lastV= hei.Vi();
    bool nextEnd = hei.NextFE();
    if( hei.Fi()!= faces.second){
        // then we skip the first and directly go from there
        // we need to change the position
        cout<<"change position of old "<<lastF<<" "<<lastV<<endl;
        updateWithNewId(Fg, cve->vert, lastF, newVertIdx);
    }
    if(0!= nextEnd) {
        cout << "change position of " << hei.Fi()<<" "<<hei.Vi() << endl;
        updateWithNewId(Fg, cve->vert, hei.Fi(), newVertIdx);
        while (hei.NextFE() != 0) {
            cout << " while change position of " << hei.Fi() <<" "<<hei.Vi()<< endl;
            updateWithNewId(Fg, cve->vert, hei.Fi(), newVertIdx);

        }
        cout << " vertex of hei " << hei.Vi() << endl;// whil ehei.next != 0
//    cout<<hei.NextFE()<<hei.Fi()<<" the next, and after "<<hei.NextFE()<<hei.Fi()<<" "<<hei.NextFE()<<hei.Fi()<<endl;
    }
    cout<<Fg.row(156)<<", row 156"<<endl;
    cout<<Fg.row(782)<<", row 782"<<endl;
    cout<<Fg.row(417)<<", row 417"<<endl;
    newVg.row(newVertIdx) = Vg.row(cve->vert);


    cout<<insertIdx<<" the inserted index"<<endl;// TODO CASE IT IS NOT RIGHT NEITHER LEFT BUT ACTUALLY ON!!
//    bool compareRight = isRight(Vg.row(cve->vert), Vg.row(cve-> vert)+midVec.transpose(), Vg.row(boundaryL[cve->patch][rightId]) ); ;//isRight(Vg.row(boundaryL[cve->patch][rightId]), Vg.row(cve-> vert)+midVec.transpose(), Vg.row(boundaryL[cve->patch][rightId]) );
////
////    updateWithNewId(Fg, cve->vert, rightFaceId, newVertIdx);
////
////    // adapt the position
//    if(cve->levelOne){
//        if(isRight(Vg.row(cve->vert), Vg.row(cve-> vert)+midVec.transpose(),  Vg.row(cve->vert) + (eps * toRight).transpose()))
//        {
//            newVg.row(newVertIdx) = Vg.row(cve->vert) + (eps* 10 * toRight).transpose();
//            newVg.row(cve->vert) = Vg.row(cve->vert) + (eps* 10 * toLeft).transpose();
//
//        }else{
//            newVg.row(cve->vert)= Vg.row(cve->vert) + (eps * 10* toRight).transpose();
//            newVg.row(newVertIdx) = Vg.row(cve->vert) + (eps * 10* toLeft).transpose();
//        }
//
//    }else{
//        if(isRight(Vg.row(cve->vert), Vg.row(cve-> vert)+midVec.transpose(),  Vg.row(cve->vert) + (eps * toRight).transpose())){
//            newVg.row(newVertIdx) = Vg.row(cve->vert) + (eps * 10* cve -> rightdirection).transpose();
//            newVg.row(cve->vert) = Vg.row(cve->vert) + (eps * 10* cve -> leftdirection).transpose();
//        }else{
//            newVg.row( cve->vert)= Vg.row(cve->vert) + (eps * 10* cve -> rightdirection).transpose();
//            newVg.row(newVertIdx) = Vg.row(cve->vert) + (eps* 10 * cve -> leftdirection).transpose();
//        }
//
//
//
//    }
//cout<< newVg.row(newVertIdx).transpose()<<" and "<< newVg.row(cve->vert).transpose()<<endl;



    Vg.resize(Vg.rows()+1, 3);
    Vg= newVg;

    // for each adjacent face we update it to the original or new vertex
//    for(int i=0; i< adjacentFaces.size(); i++){
////        if(adjacentFaces[i] == leftFaceId || adjacentFaces[i] == rightFaceId ) continue; // they are handled separately
//        int testVert =-1;
//        // we take one edge and check it's side
//        if(Fg(adjacentFaces[i], 0)!= cve-> vert && Fg(adjacentFaces[i], 0)!= insertIdx ){
//            testVert = Fg(adjacentFaces[i], 0);
//        }else if( Fg(adjacentFaces[i], 1)!= cve-> vert && Fg(adjacentFaces[i], 1)!= insertIdx ){
//            testVert = Fg(adjacentFaces[i], 1);
//        }
//        else {
//            testVert = Fg(adjacentFaces[i], 2);
//        }
//        if(testVert == -1 )cout<<" no suitable test vert found, something is wrong "<<endl;
//        cout<<testVert<<" vert from face "<<adjacentFaces[i]<<endl;
////        bool checkRight = isRight(Vg.row(cve->vert), Vg.row(cve-> vert)+midVec.transpose(), Vg.row(testVert) );
//        bool checkRight = isRight(Vg.row(cve->vert), Vg.row(cve-> vert)+midVec.transpose(), Vg.row(testVert) );
//
//
////        if (checkRight == compareRight){
//        if(checkRight){
//            cout<<"is updated"<<endl;
//            updateWithNewId(Fg, cve->vert, adjacentFaces[i], newVertIdx);
//        }
//        if(adjacentFaces[i] == 62|| adjacentFaces[i] == 776 || adjacentFaces[i] ==781){
//            cout<<" face "<<adjacentFaces[i] <<" : "<<Fg.row(adjacentFaces[i])<<endl;
//        }
//    }
    cout<<" right "<< Fg.row(rightFaceId)<<endl;
    cout<<"left "<<Fg.row(leftFaceId)<<endl;


    if(toPattern_boundaryVerticesSet.find(insertIdx)!= toPattern_boundaryVerticesSet.end()){
            cout<<"we are nearly done here, it's cut through! "<<endl;
            cve-> vert = insertIdx;
            cve -> bridgeFlag = true;
    }

    // only the first level i.e. boundary duplicate has to be projected, hence only this one is to be added
    if(cve->levelOne) {
        addToDulicateList(cve,seamsList, minusOneSeams, newVertIdx );
        toPattern_boundaryVerticesSet.insert(newVertIdx);// the duplicate is also on the boundary, hence insert it
        cve->leftdirection = toLeft;
        cve->rightdirection = toRight;
        cout<<" inserted"<<endl;
    }
    cve-> vert = insertIdx;
    cve->levelOne = false;
    cout<<"fin"<<endl<<endl;



}
void splitVertex(Node** head, int & listLength, int  whichTear, MatrixXd& Vg, MatrixXi& Fg, vector<vector<int> >& vfAdj, MatrixXd& lengthsOrig, MatrixXd& lengthsCurr){
    cout<<"splitting "<<endl;

    Node* node = (*head);
    int currIdx = 0;
    while (currIdx<whichTear){
        node = node-> next;
        currIdx++;
    }
    whichTear = node->vertId;

    // to which vert is the vertex adjacent? from there we get all the adjacent edges
    vector<vector<int>> vvAdj;
    igl::adjacency_list(Fg,vvAdj);
    vector<int> adjacentFaces = vfAdj[whichTear];
    if(adjacentFaces.size()<2){
        // there is nothing to split, only left and right vertex are adjacent
        cout<<" no splitting here"<<endl;
        return;
    }

    int newVertIdx = Vg.rows();
    MatrixXd newVg (newVertIdx+1, 3);
    newVg.block(0,0, newVertIdx, 3)= Vg;

    node -> duplicate = newVertIdx;
    // we need to know which is closest to normal in order to insert this and new vertex and in right order
    int leftId = node -> next -> vertId;
    int rightId = node -> prev -> vertId;

    Vector3d toLeft = Vg.row(leftId)- Vg.row(whichTear);
    Vector3d toRight = Vg.row(rightId) - Vg.row(whichTear);
    Vector3d midVec = (toLeft.normalized()+ toRight.normalized())/2;
    midVec(0)= -toLeft.normalized()(1);
    midVec(1)= toLeft.normalized()(0);
//    cout<<midVec.transpose()<<" mid vector"<<endl;//TODO THIS IS a heuristic
//    double angle = acos((toLeft.normalized()).dot(toRight.normalized()));
//    double deg = angle*180/M_PI;

    // first: find closest to mid, this is the new inserted voundary point
    double dist = std::numeric_limits<double>::max();
    int idxofClosest = -1;
    for(int i=0; i<vvAdj[whichTear].size(); i++){
        int adjVert = vvAdj[whichTear][i];
        if(adjVert== rightId || adjVert == leftId) continue; // we want a middle one
        Vector3d edgeVec = Vg.row(adjVert)- Vg.row(whichTear);
        edgeVec= edgeVec.normalized();
        // both have unit distance, so as a measure we can take the distance form another
        if((midVec - edgeVec).norm()<dist){
            idxofClosest= i;
            dist = (midVec - edgeVec).norm();
        }

    }
    // is the new to be inserted vertex
    int insertIdx =  vvAdj[whichTear][idxofClosest];
    cout<<insertIdx<<" the inserted index"<<endl;
    // TODO if this is a boundary vertex we disconnect the mesh and are finished

    double x1 = Vg(whichTear, 0); double y1 = Vg(whichTear, 1);
    double x2 = Vg(insertIdx, 0); double y2 = Vg(insertIdx, 1);

    // it should be sufficient to check one edge per face since it cannot happen that a face has edges on both sides
    // if it coudl then it is (a) not stressed or (b) we found no splitnormal

    // right gets  newVertIdx
    double dRight = (Vg(rightId, 0) - x1) * (y2-y1) - (Vg(rightId, 1) - y1) * (x2 - x1);
    // defines the side of right, left is the opposite
    bool rightDSmaller = (dRight<0);

    if(Fg(node->prev->faceId , 1) == whichTear){
        Fg(node->prev->faceId, 1) = newVertIdx;
    }
    else if(Fg(node-> prev->faceId , 0)== whichTear){
        Fg(node->prev->faceId, 0) = newVertIdx;
    }
    else {
        Fg(node->prev->faceId, 2) = newVertIdx;
    }

    // adapt the position
    double eps = 0.01;
    cout<<"orig "<<Vg.row(whichTear)<<endl;
    cout<<toRight.transpose()<<" right"<<endl;
    cout<<toLeft.transpose()<<" left"<<endl;
    newVg.row(newVertIdx)= Vg.row(whichTear)+ (eps * toRight).transpose();
    newVg.row(whichTear) = Vg.row(whichTear) + (eps * toLeft).transpose();
    Vg.resize(Vg.rows()+1, 3);
    Vg= newVg;

    // adapt the faces to have the right id
    for(int i=0; i<adjacentFaces.size(); i++){
        // we take one edge and check it's side
        if(adjacentFaces[i]== leftId || adjacentFaces[i]== rightId) continue; // they are handled seperately
        int testVert =-1;
        if(Fg(adjacentFaces[i], 0)!= whichTear &&
        Fg(adjacentFaces[i], 0)!= newVertIdx &&
        Fg(adjacentFaces[i], 0)!= leftId &&
        Fg(adjacentFaces[i], 0)!= rightId &&
        Fg(adjacentFaces[i], 0)!= insertIdx
        ){
            testVert = Fg(adjacentFaces[i], 0);
        }else if(
                Fg(adjacentFaces[i], 1)!= whichTear &&
                Fg(adjacentFaces[i], 1)!= newVertIdx &&
        Fg(adjacentFaces[i], 1)!= leftId &&
        Fg(adjacentFaces[i], 1)!= rightId &&
        Fg(adjacentFaces[i], 1)!= insertIdx ){
            testVert = Fg(adjacentFaces[i], 1);
        }
        else if(Fg(adjacentFaces[i], 2)!= whichTear &&
                Fg(adjacentFaces[i], 2)!= newVertIdx &&
                 Fg(adjacentFaces[i], 2)!= leftId &&
                 Fg(adjacentFaces[i], 2)!= rightId &&
                 Fg(adjacentFaces[i], 2)!= insertIdx ){
            testVert = Fg(adjacentFaces[i], 2);
        }
        else{
            cout<<"face "<<adjacentFaces[i]<<" becomes skipped and remains "<< Fg.row(adjacentFaces[i])<<endl;
            continue;
        }

        double d = (Vg(testVert, 0) - x1)*(y2-y1)-(Vg(testVert, 1) - y1)*(x2 - x1);

        // if it is one of the same side as the one we call right
        if (rightDSmaller == (d<0)){
            if(adjacentFaces[i]==5319)cout<<" in loop"<<endl;
            // same side as what we call right
            if(Fg(adjacentFaces[i], 1)== whichTear && testVert!= -1  ){
                Fg(adjacentFaces[i], 1) = newVertIdx;
            }
            else if(Fg(adjacentFaces[i], 0)== whichTear && testVert!= -1){
                Fg(adjacentFaces[i], 0) = newVertIdx;
            }
            else if(testVert!= -1 && Fg(adjacentFaces[i], 2)== whichTear){
                Fg(adjacentFaces[i], 2) = newVertIdx;
            }

        }
        cout<<"face "<<adjacentFaces[i]<<" becomes "<< Fg.row(adjacentFaces[i])<<endl;
    }

    // update vertexFace
    createVertexFaceAdjacencyList(Fg, vfAdj);
    int faceOfMid = adjacentFaceToEdge( insertIdx, whichTear, -1, vfAdj );
    int faceOfRight = adjacentFaceToEdge( insertIdx, newVertIdx, -1, vfAdj );

    // now insert into the linked list
// node is the left side, insert idx follows right and then the new one
    Node* middle_node = new Node();
    middle_node->vertId = insertIdx;
    middle_node->faceId = faceOfMid;
    middle_node -> duplicate = -1;

    Node* right_noed = new Node();
    right_noed-> vertId = newVertIdx;
    right_noed->faceId = faceOfRight; //right_noed -> faceIdRight = node->faceIdRight;
    right_noed-> duplicate = whichTear;

    middle_node-> next = node;
    middle_node->prev = right_noed;
    right_noed-> next = middle_node;
    right_noed -> prev = node->prev;
    node->prev ->next = right_noed;
    node ->prev = middle_node;

    listLength+=2;
    // find middle node edge of face id
    // find right node edge of face id

    middle_node->edgeOfFaceId = findWhichEdgeOfFace(faceOfMid, insertIdx, whichTear, Fg);
//    cout<< middle_node->edgeOfFaceId<<" edge left "<<endl;
    right_noed->edgeOfFaceId = findWhichEdgeOfFace(faceOfRight, newVertIdx, whichTear, Fg);
//    cout<< right_noed->edgeOfFaceId<<" edge right "<<endl;

    double origLength = lengthsOrig(faceOfMid, middle_node->edgeOfFaceId);
    double newLength = lengthsCurr ( faceOfMid, middle_node->edgeOfFaceId);

    middle_node-> stretchLeft = newLength/origLength;// assuming the new one is stretched it is certainly longer

    origLength = lengthsOrig(faceOfRight, right_noed->edgeOfFaceId);
    newLength = lengthsCurr (faceOfRight, right_noed->edgeOfFaceId);

    right_noed-> stretchLeft = newLength/origLength;// assuming the new one is stretched it is certainly longer

}

int addoncount=0;
void addVarToModel (int vert, int nextVert, vector<vector<int>> & vfAdj, bool isConstrained, int& varCount, GRBVar* & cutVar,
                      MatrixXi& Fg_pattern,MatrixXd& lengthsOrig, MatrixXd& lengthsCurr, map <int, cutVertEntry*> & mapVarIdToVertId,
                      int seamType, int seamIdInList, double tailor_lazyness, bool corner, map<pair<int,int>, int>& mapVertAndSeamToVar, int counterpart, GRBModel& model ){

    int faceIdx = adjacentFaceToEdge(vert, nextVert, -1, vfAdj );
    int whichEdge = findWhichEdgeOfFace(faceIdx, vert, nextVert, Fg_pattern);
    double w_init = lengthsCurr(faceIdx, whichEdge)/lengthsOrig(faceIdx, whichEdge);

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

void setLP(std::vector<std::vector<int> >& boundaryL , vector<vector<int>> & vfAdj, MatrixXi& Fg_pattern,
MatrixXd& lengthsOrig, MatrixXd& lengthsCurr,const std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary, map<int, vector<pair<int, int>>>& seamIdPerCorner,
           vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,
            int tearPatch, int tearVert, int tearVertInBoundaryIndex, double tailor_lazyness, double minConstrained, vector <cutVertEntry*>& cutPositions,
            VectorXd& cornerVert){

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();
    // Create an empty model
    GRBModel model = GRBModel(env);
    int varCount =0;
    int numVar=0; // whole boundary and all corners duplicate
    for(int i=0; i<edgesPerBoundary.size(); i++){
        numVar += edgesPerBoundary[i].size();
        numVar += boundaryL[i].size();
    }
    // with the constraints we have for each counterpart another variable
    if(numVar%2 != 0 )cout<<"----------------------------------not divisible by 2, the number of constraints is wrong------------------"<<endl<<endl<<endl;
    numVar += (numVar/2);// 574 without mapping -> 287 added ones ? we don't need them all bc there are no common ones for th -1s
    map<int, int> trackCornerIds;
    map <int, cutVertEntry*>  mapVarIdToVertId;
    map <pair<int, int>, int> mapVertAndSeamToVar;

    // a map to track the gurobi id of a corner.  whenever we set a corner we add the corner id, so we know once we set the second and connect them
    GRBVar* cutVar = model.addVars(numVar, GRB_BINARY);
    //edgesPerBoundary.size();
    for(int i = 0; i < edgesPerBoundary.size(); i++) {
        int boundSize = boundaryL[i].size();
        int startVarPatch = varCount;
        for (int j = 0; j < edgesPerBoundary[i].size(); j++) {
            // first is the absolute index, second the index wrt the boundary loop


            auto cornerPair = edgesPerBoundary[i][j];
            if(seamIdPerCorner.find(cornerPair.first)== seamIdPerCorner.end()) continue;

            vector<pair<int, int>> seamId = seamIdPerCorner[cornerPair.first];
            cout<<seamId.size()<<" the size ";
            if(seamId.size()>2) cout<<" something is veryy odd!! we have more than two seams for a corner. impossible."<<endl;
            for(int si = 0; si<seamId.size(); si++) {
                cout<<" info "<<seamId[si].first<<" "<<seamId[si].second<<endl;

                GRBLinExpr innerSumConstr = 0; // interior sum
                GRBLinExpr lSumConstr = 0; // left sum
                GRBLinExpr rSumConstr = 0; // right sum
                int startVarOfThis = varCount;
                int count = 0;
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
                        auto startAndPatch = seam->getStartAndPatch1();
                        startAndPatchOther = seam-> getStartAndPatch2ForCorres();
                        if (startAndPatch.second != i)
                            cout << " now in th e+1 seams the patch does not match where we are in the loop "<< startAndPatch.second << endl;
                        idx = startAndPatch.first;
                        idxOther = startAndPatchOther.first;
                        vert = boundaryL[startAndPatch.second][idx];
                        vertOther =  boundaryL[startAndPatchOther.second][idxOther];
                        boundSizeOther = boundaryL[startAndPatchOther.second].size();
                        end = seam->getEndCornerIds().first;
                        length = seam->seamLength();

                    } else if (seamId[si].first < 0) {
                        minusOneSeam *currSeam = minusOneSeams[seamId[si].second];
                        int patch = currSeam->getPatch();

                        if (patch != i)
                            cout << " now in th e-1 seams the patch does not match where we are in the loop , would be patch "<<patch << endl;
                        idx = currSeam->getStartIdx();
                        vert = boundaryL[i][idx];
                        end = currSeam->getEndVert();
                        length = currSeam->getLength();
                    }
                    // check if the corners exist already. If so then connect with corner, else add the indices
                    if(trackCornerIds.find(vert) != trackCornerIds.end()){
                        model.addConstr(cutVar[ trackCornerIds[vert]] + cutVar[varCount] <= 1);
                    }else{
                        trackCornerIds[vert] = varCount;
                    }
                    while (vert != end) {
                        // might need a map for patch and vert id to seam adn first or Second
                        double relId = ((double) count) / length;
                        bool isConstrained = true;
                        bool corner = false;
                        if(count ==0) corner = true;
                        if (relId == 0 || relId > minConstrained && relId < (1 - minConstrained)) {
                            isConstrained = false;
                        }else{
                            model.addConstr(cutVar[varCount] == 0);
                        }
                        int currVar = varCount;

                        addVarToModel(vert, boundaryL[i][(idx + 1) % boundSize], vfAdj, isConstrained, varCount, cutVar,
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
                    addVarToModel(vert, boundaryL[i][(idx + 1) % boundSize], vfAdj, false, varCount, cutVar, Fg_pattern,
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

                    while (vert != end) {

                        double relId = ((double) count) / seam->seamLength();
                        bool isConstrained = true;
                        bool corner = false;
                        if(count ==0) corner = true;
                        int currVar = varCount;
                        if (relId == 0 || relId > minConstrained && relId < (1 - minConstrained)) {
                            isConstrained = false;
                        }else{
                            model.addConstr(cutVar[varCount] == 0);
                        }
                        addVarToModel(vert, nextvert, vfAdj, isConstrained, varCount, cutVar, Fg_pattern, lengthsOrig,
                                      lengthsCurr, mapVarIdToVertId, seamId[si].first, seamId[si].second, tailor_lazyness, corner, mapVertAndSeamToVar, otherVert, model );
                        if (!isConstrained) {
                            lSumConstr += cutVar[currVar];
                            if (relId != 0) {
                                rSumConstr += cutVar[currVar];
                                innerSumConstr += cutVar[currVar];
                            }
                        }

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
                    addVarToModel(vert, boundaryL[startAndPatch.second][(idx + 1) % boundSize], vfAdj, false, varCount,
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
            cout<<patch<<" chosen Id  "<< mapVarIdToVertId[i]->vert<<" from which kind of seam "<<mapVarIdToVertId[i]->seamType<<" and id "<<mapVarIdToVertId[i]->seamIdInList <<endl;// ALRIGHT BUT NOW WE NEED A SEAM ID

            cutVertEntry* cve = new cutVertEntry ( vert, seamType, seamId, patch);
            if(cornerVert[vert]==1){
                cout<<"corner"<<endl;
//                // left or right corner?
                int firstInSeam;
                if(seamType == -1 ){
                    firstInSeam = minusOneSeams[seamId]->getStartVert() ;
                }else if(seamId>=0){
                    firstInSeam = seamsList[mapVarIdToVertId[i]->seamIdInList]->getStartAndPatch1().first;
                }else{
                    firstInSeam = seamsList[(mapVarIdToVertId[i]->seamIdInList+1)*(-1)]->getStartAndPatch2ForCorres().first;

                }
                if(firstInSeam==vert){
                    cve -> startCorner = true;
                }else{
                    cve-> endCorner = true;
                }
            }
            cutPositions.push_back(cve);

        }

    }


}

void tearFurther(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                 map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL, set<int> & cornerSet
){
    cout<<"-----------------------"<<endl<<endl;
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);

    for(int i = 0; i < 5; i++){// cutPositions.size(); i++){
        cout<<endl<< cutPositions[i]->vert<<" vertex up next handling"<<endl;
        splitVertexFromCVE(cutPositions[i], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList, minusOneSeams, releasedVert,
                           toPattern_boundaryVerticesSet,lengthsCurr, Fg_pattern, cornerSet );
        cout<<"back in loop "<<cutPositions.size()<<endl;
        cout<< cutPositions[4]->vert<<" vertex up next handling end"<<endl;

    }
    cout<<"--------------------"<<endl;
}
void computeTear(Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern, MatrixXi& Fg_pattern_orig,
                 vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams, std::vector<std::vector<int> >& boundaryL, bool & finished,
                 const std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary, map<int, vector<pair<int, int>>>& seamIdPerCorner,
                 VectorXd& cornerVert, vector<cutVertEntry*>& cutPositions,  map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet, set<int> & cornerSet )
                 {

    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    igl::edge_lengths(fromPattern, Fg_pattern_orig, lengthsOrig);

    // first we need to know where to tear, set up LP for this
    // information we need: stress. For stress, we need lengths old and ned
    // for lengths old and new we need which edge of face
    // for which edge of face we need face
    // for face we need vertices
    // for vertices we need boundary loop
    int tearPatch, tearVert, tearVertInBoundaryIndex;
    double tailor_lazyness = 1;
    double minConstrained = 0.25;
    setLP( boundaryL, vfAdj, Fg_pattern, lengthsOrig, lengthsCurr, edgesPerBoundary, seamIdPerCorner,
           seamsList, minusOneSeams, tearPatch, tearVert, tearVertInBoundaryIndex, tailor_lazyness, minConstrained, cutPositions,
           cornerVert);


    // we cut the first one
    for(int i = 0; i < 5; i++){// cutPositions.size(); i++){
        splitVertexFromCVE(cutPositions[i], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                           minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, Fg_pattern, cornerSet);
    }

cout<<"fin compute tear "<<endl ;
}

void updatePositionToIntersection(MatrixXd& p,int next,  const MatrixXd& Vg_bound){

    // we know that the first two indices of the face define the edge we want to intersect with
    // derive where QR and P meet = t, https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line
    double minDist = std::numeric_limits<double>::max();
    Vector3d minDistTarget;
    for(int i=1; i< Vg_bound.rows(); i++){
        Vector3d R = Vg_bound.row(i-1);
        Vector3d Q =  Vg_bound.row(i);

        double t = (R-Q).dot(Q-p.row(next).transpose())/((R-Q).dot(R-Q));
        t = min(1., t);
        t = max(0., t);
        Vector3d targetPos = Q-t*(R-Q);

        double dist = (targetPos - p.row(next).transpose()).norm();
        if( dist < minDist){
            minDist = dist;
            minDistTarget = targetPos;
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
void projectBackOnBoundary(const MatrixXd & Vg_to, MatrixXd& p, const vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,  const MatrixXi& Fg_pattern,
                           const MatrixXi& Fg_pattern_orig, const std::vector<std::vector<int> >& boundaryL_toPattern, const std::vector<std::vector<int> >& boundaryL,
                           map<int, pair<int, int>> & releasedVert){
    // add duplicates to the mapping

    // Idea:
    // iterate over every seam, create a triangle mesh of the original boundary of that seam
    // then for each boundary vertex, project it's current position to the boundary
    // update the current position towards the projected
    int numSeams = seamsList.size();
    for (int j = 0; j<numSeams; j++){

        seam* currSeam  = seamsList[j];
        auto stP1= currSeam-> getStartAndPatch1();
        auto stP2 =  currSeam -> getStartAndPatch2ForCorres(); // attention this is with respect to the original pattern
        int len  = currSeam -> seamLength();
        int boundLen1 = boundaryL_toPattern[stP1.second].size();
        int boundLen2 = boundaryL_toPattern[stP2.second].size();

        // build the structure for closest search
        MatrixXd Vg_seam1to(len+1 , 3);
        MatrixXd Vg_seam2to(len+1 , 3);


        for(int i = 0; i<= len ; i++){
            int v1_oneSide = boundaryL_toPattern[stP1.second][(stP1.first+i) % boundLen1];
            int v1_otherSide_idx = (stP2.first-i) % boundLen2;
            if(v1_otherSide_idx < 0) {
                v1_otherSide_idx +=boundLen2;
            }
            if(seamsList[j]->inverted) v1_otherSide_idx = (stP2.first + i) % boundLen2;
            int v1_otherSide = boundaryL_toPattern[stP2.second][v1_otherSide_idx];

            Vg_seam1to.row(i) = Vg_to.row(v1_oneSide);
            Vg_seam2to.row(i) = Vg_to.row(v1_otherSide);

        }
        // for each interior (=not corner) vertex of the new boundary we need to check the closest position
        // todo never ever cut the corner or change the corner index
        auto ends = currSeam->getEndCornerIds();
        int i1=0;
        int bsize = boundaryL_toPattern[stP1.second].size();
        int next = boundaryL_toPattern[stP1.second][(stP1.first+i1)% bsize];
//        cout<<next<<" searching for end "<<ends.first<<endl;
        pair<int, int> compPair = make_pair(1,j );

        while( next!= ends.first ){
            if(i1 == 0){
                if(releasedVert.find(next) != releasedVert.end() && releasedVert[next]== compPair){
                    // it is released for another side hence we have to pull it to our side
                    updatePositionToIntersection( p, next,Vg_seam1to);
                }
            }

           else if(releasedVert.find(next) == releasedVert.end() ){
               // general case an interior vertex , if it is not constrained pull it to boundary
                updatePositionToIntersection( p, next,Vg_seam1to);
            }else if(releasedVert[next]== compPair){
               // if it is constrained but from another side pull it to boundary
                updatePositionToIntersection( p, next,Vg_seam1to);
            }

            i1++;
            next = boundaryL_toPattern[stP1.second][(stP1.first+i1)% bsize];
        }
        // the last corner. Again if it is constrained from another side pull it to boundary, else ignore since handled by corner
        if(releasedVert.find(next) != releasedVert.end() && releasedVert[next] == compPair ){
            updatePositionToIntersection( p, next,Vg_seam1to);
        }

        int sizeOneSide = i1+1; // account for  0and last

        int i2 = 0;
        bsize = boundaryL_toPattern[stP2.second].size();
        int nextidx = (stP2.first- i2) % bsize;
        if(nextidx < 0) nextidx += bsize;
        next = boundaryL_toPattern[stP2.second][nextidx];
//        cout<<next<<"and now "<<ends.second<<endl;
        compPair = make_pair(1,-j-1 );

        while( next!= ends.second ){

            if(i2 == 0){
                if(releasedVert.find(next) != releasedVert.end() && releasedVert[next] == compPair){
                    // it is released for another side hence we have to pull it to our side
                    updatePositionToIntersection( p, next,Vg_seam2to);
                }
            }
            else if(releasedVert.find(next) == releasedVert.end()){
                updatePositionToIntersection(p, next,Vg_seam2to);
            }else if(releasedVert[next] == compPair){
                updatePositionToIntersection(p, next,Vg_seam2to);

            }

            i2++;
            nextidx = (stP2.first - i2) % (bsize);

            if(nextidx < 0) {nextidx += bsize;}
            if(seamsList[j]->inverted) nextidx = (stP2.first + i2) % bsize;
            next = boundaryL_toPattern[stP2.second][nextidx];
        }

        if(releasedVert.find(next) != releasedVert.end() && releasedVert[next] == compPair){
            // last is released for another side hence we have to pull it to our side
            updatePositionToIntersection( p, next,Vg_seam2to);
        }
        // also project all duplicates
        for(const auto & addedVert : currSeam->duplicates){
            updatePositionToIntersection(p, addedVert.second, Vg_seam1to);
        }
        for(const auto & addedVert : currSeam->duplicates2){
            updatePositionToIntersection(p, addedVert.second, Vg_seam2to);
        }

    }
    for(int j=0; j<minusOneSeams.size(); j++){
//
        minusOneSeam* currSeam  = minusOneSeams[j];
        int patch = currSeam -> getPatch();
        int startidx = currSeam -> getStartIdx();
        int endVert = currSeam -> getEndVert();
        int len = currSeam -> getLength();
        int boundLen = boundaryL_toPattern[patch].size();

        // build the structure for closest search
        MatrixXd Vg_seamto(len+1 , 3);

        for(int i = 0; i<= len ; i++){
            int v1 = boundaryL_toPattern[patch][(startidx+i)% boundLen];
            Vg_seamto.row(i) = Vg_to.row(v1);
        }

        int i1=1;
        int next = boundaryL_toPattern[patch][(startidx+i1) % boundLen];
        pair<int, int> compPair = make_pair(-1,j );
        while(next != endVert){
            if(i1==0){
                if(releasedVert.find(next) != releasedVert.end() && releasedVert[next] == compPair){
                    // it is released for another side hence we have to pull it to our side
                    updatePositionToIntersection( p, next,Vg_seamto);
                }
            }
            else if(releasedVert.find(next)== releasedVert.end()){
                updatePositionToIntersection( p, next,Vg_seamto);
            }else if(releasedVert[next] == compPair){
                updatePositionToIntersection( p, next,Vg_seamto);

            }
            i1++;
            next = boundaryL_toPattern[patch][( startidx+i1 ) % boundLen];
        }

        if(releasedVert.find(next) != releasedVert.end() && releasedVert[next] == compPair){
            // it is released for another side hence we have to pull it to our side
            updatePositionToIntersection( p, next,Vg_seamto);
        }
        // also map all projections
        for(const auto & addedVert : currSeam->duplicates){
            updatePositionToIntersection(p, addedVert.second, Vg_seamto);

        }

    }


}
