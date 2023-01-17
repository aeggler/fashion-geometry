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
#include "igl/is_vertex_manifold.h"
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
void computeMidVecBasedOnPCA(VectorXd& midVec, vector<vector<int>>& vvAdj, vector<vector<int>>& vfAdj, MatrixXd& Vg ,
                             MatrixXd& lengthsCurr, MatrixXd& lengthsOrig, MatrixXi& Fg_pattern,int& vert, VectorXd& ws, MatrixXd& dirs){
    int n = vvAdj[vert].size();
    ws.resize(n);
    midVec.resize(3);
    dirs.resize(n, 3);
    std::pair<int, int>  faces;

    for(int i=0; i<vvAdj[vert].size(); i++) {
        int otherVert = vvAdj[vert][i];
        dirs.row(i) = (Vg.row(otherVert) - Vg.row(vert)).normalized();
        adjacentFacesToEdge(vert, otherVert, vfAdj, faces );

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
    }
    fitVecToPointSet( dirs, midVec );
    midVec(2) = 0;// Vg(vert, 2); // just to get the 3rd constant dimension right
    midVec = midVec.normalized();
    cout<<" the newly computed midvec of "<<vert<<" is "<<midVec.transpose()<<endl;
}

bool checkIfTearIsUseful(int vert, Vector3d& cutDirection,  vector<vector<int>>& vvAdj,  vector<vector<int>>& vfAdj, MatrixXd& Vg ,
                         MatrixXd& lengthsCurr,  MatrixXd& lengthsOrig, MatrixXi& Fg_pattern, VectorXd& ws, bool preComputed ){
    // check if it makes sense, i.e. releases stress cutting in the direction
    // if the dot product for at least one adjacent vertex of the one we are going to cut is large enough we allow to cut further
    // todo 3.1. maybe a better option is not to ignore ones with wrong direction but to actually clamp the stess

    double thereshold = 1.051;// todo fix this!! it cuts pretty mich always
    double thW = 1.1; double thDot = 0.5;
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
            }else if(faces.second!= -1){
                int whichEdgeRight = findWhichEdgeOfFace(faces.second, vert, otherVert, Fg_pattern);
                w += lengthsCurr(faces.second, whichEdgeRight)/lengthsOrig(faces.second, whichEdgeRight);
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
map<int, int> releasedVertNew;
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
                         set<int>& handledVerticesSet){

    cout<<" patch "<<cve->patch<<" of "<<boundaryL.size() <<endl;

    if(cve-> finFlag) {
        cout<<cve->vert<<" done already "<<endl;
        return;
    }
    cve->handled = true;
// todo this prohibits any kind of cross cutting! maybe not always desirable ! removing?
    if(handledVerticesSet.find(cve->vert) != handledVerticesSet.end()){
        cout<<"Handled by other seams already. Stop here"<<endl;
        cve->finFlag = true;
        return ;
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
        Vector3d cutDirection = cve-> continuedDirection;
        VectorXd ws;
        if(! checkIfTearIsUseful(cve-> vert, cutDirection, vvAdj, vfAdj, Vg, lengthsCurr, lengthsOrig, Fg, ws, false)){
            cve-> finFlag = true;
            cout<<"stopping now bc tearing is not useful anymore "<<endl;

            return;
        }

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
        if(cornerToSeams[cve->cornerInitial][0] == seamComp ){
            releasedVertNew[cve->vert]= cornerToSeams[cve->cornerInitial][1];

        }else if (cornerToSeams[cve->cornerInitial][1] == seamComp ){
            releasedVertNew[cve->vert]= cornerToSeams[cve->cornerInitial][0];
        }else{
            cout<<cornerToSeams[cve->cornerInitial][0]<<" we have a problem, it is not found "<<cornerToSeams[cve->cornerInitial][1]<<" but what we have is "
            <<cve->seamType<<" "<< cve ->seamIdInList<<endl;
            cout<<endl<<endl<<"------------"<<endl<<endl<<"----------"<<endl<<endl;
        }


        //todo 4.1.

        newVg.row(newVertIdx) = Vg.row(cve->vert);

        cve->continuedCorner = true;
        cve->finFlag = (cornerSet.find(cve->vert) != cornerSet.end()&& cve->vert != cve-> cornerInitial); //if it is a corner we are done

        if(cve->startCorner){// does not matter if it is a starter or not
            cout<<"starter"<<endl;

            if(cve->seamType>0){
                if(cve->seamIdInList>=0){
                    cve->vert = boundary[minusOneId];
                }else{
                    cve->vert = boundary[plusOneId];
                }
            }else{
                cout<<"the next one is "<<boundary[minusOneId];

                cve->vert = boundary[minusOneId];
            }

        }else{
//            cout << cve->vert <<" ending " << boundary[plusOneId] << "  FIGURE OUT WHICH ONE ATM ITS WRONG " << boundary[minusOneId] << " " << endl;
//            cout << "for seam type " << cve->seamType <<" "<< cve-> seamIdInList << endl;
//            seam* helper= seamsList[8];
//            cout<<helper->getStart1()<<" "<< helper->getEndCornerIds().first <<endl;

           if(cve->seamType>0){
                if(cve->seamIdInList>=0){
                    cve->vert = boundary[plusOneId];
                }else{
                    cve->vert = boundary[minusOneId];
                }
            }else{
                cve->vert = boundary[plusOneId];
            }
        }

        Vg.resize(Vg.rows()+1, 3);
        Vg= newVg;
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

    Vector3d midVec;// = leftRot;
    Vector3d midVect = Vg.row(boundaryL[cve->patch][plusOneId]) - Vg.row(boundaryL[cve->patch][minusOneId]);
    if(!cve->levelOne){
        midVect = cve->leftdirection - cve->rightdirection;
    }
    midVect= midVect.normalized();
    midVec(0)= -midVect(1);
    midVec(1) = midVect(0);

    VectorXd newMidVec;
    VectorXd ws; MatrixXd dirs;
    cout<<"starting midvec pca calculation"<<endl;
    computeMidVecBasedOnPCA(newMidVec, vvAdj, vfAdj, Vg, lengthsCurr, lengthsOrig, Fg, cve->vert, ws, dirs );
    midVec = (-1) * newMidVec;
    // todo sketchy, for whatever reason it breaks without this. does it do the transposing?
    cout<<" midvec "<<midVec.transpose()<<endl;

    Vector3d cutDirection = midVec;
    cutDirection(0)= -cutDirection(1);
    cutDirection(1) = midVec(0);
    if(!checkIfTearIsUseful(cve-> vert, cutDirection, vvAdj, vfAdj, Vg, lengthsCurr, lengthsOrig, Fg, ws, true)){
        cve-> finFlag = true;
        cout<<"stopping now with new condition "<<endl;
        return;
    }

    //  we have the midvec direction, but don't know for sure if adding or subtract. Take tthe direction that has longer distance to the existing boundary to not go backwards.
    // This is a heuristic.
    Vector3d A, Btemp;
    if(cve->levelOne){
        A = Vg.row(boundary[plusOneId]);
        Btemp = Vg.row(cve -> vert);
    }else{
        A = Vg.row(cve -> leftCorner);
        Btemp = Vg.row(cve -> rightCorner);
    }


    Vector3d C = Vg.row(cve -> vert).transpose()+ midVec;
    //(AB,AM), where M(X,Y) is the query point:
    //https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
    auto sign1 = abs((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0)));
    sign1 /= sqrt((Btemp(0)-A(0) )*(Btemp(0)-A(0) ) + (Btemp(1)-A(1))*(Btemp(1)-A(1)));
    C = Vg.row(cve -> vert).transpose() + newMidVec;
    auto sign2 = abs((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0)));
    sign2 /= sqrt((Btemp(0)-A(0) )*(Btemp(0)-A(0) ) + (Btemp(1)-A(1))*(Btemp(1)-A(1)));



    for(int i=0; i<vvAdj[cve -> vert].size(); i++) {
        int adjVert = vvAdj[cve->vert][i];
        if(! (adjVert ==boundary[ minusOneId] || adjVert == boundary[plusOneId]) )continue; // we want a border
        //todo if levelone take the side of the other interior vertices
        // the direction we want has a longer distance to the border
        Vector3d edgeVec = Vg.row(adjVert)- Vg.row(cve -> vert);
        edgeVec= edgeVec.normalized();
        // both have unit distance, so as a measure we can take the distance from another
        double posVec = (midVec - edgeVec).norm();
        double negVec = (newMidVec - edgeVec).norm();
        cout<<posVec<<" pos and neg "<<negVec<<endl;
        cout<<sign1<<" sign 1 and 2 "<<sign2<<endl;
        // take the longer sign. If they are the same (= initial), take the one on the same side as interior vertices
        double eps = 0.0001;
        if(sign2 > sign1){
            midVec = newMidVec;
            cout<<"WE CHANGED THE SIGN OF THE MIDVEC!"<<endl;
        }else if ((sign1- sign2) < eps){
            cout<<"initial!"<<endl;

            Vector3d C = Vg.row(cve -> vert).transpose()+ midVec;
            //(AB,AM), where M(X,Y) is the query point:
            //https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
            auto sign1 =(((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0))) > 0);
            cout<<"side of midvec "<<sign1<<endl;
            // take the side of interior vertices
            for(int j=0; j<vvAdj[cve -> vert].size(); j++) {
                int testInt = vvAdj[cve->vert][j];
                if ((testInt == boundary[minusOneId] || testInt == boundary[plusOneId]))continue;
                // it is an interior
                C = Vg.row(testInt);
                auto sign2 = (((Btemp(0)-A(0) )*(A(1)-C(1))-(Btemp(1)-A(1))*(A(0)-C(0))) > 0);
                cout<<"side of interior vert "<<sign2<<endl;
                if(sign1!= sign2){
                    midVec = newMidVec;
                    cout<<"WE CHANGED THE SIGN OF THE MIDVEC!"<<endl;
                }
                break;
            }

        }
//        if(negVec > posVec){
//            midVec = newMidVec;
//            cout<<"WE CHANGED THE SIGN OF THE MIDVEC!"<<endl;
//        }
        break;

    }

    // all other normal cases
    double dist = std::numeric_limits<double>::max();
    int idxofClosest = -1;
    for(int i=0; i<vvAdj[cve -> vert ].size(); i++){
        int adjVert = vvAdj[cve -> vert][i];
        if(adjVert== boundary[minusOneId ]|| adjVert == boundary[plusOneId]) continue; // we want a middle one

        Vector3d edgeVec = Vg.row(adjVert)- Vg.row(cve -> vert);
        edgeVec= edgeVec.normalized();
        // both have unit distance, so as a measure we can take the distance from another
        if((midVec - edgeVec).norm() < dist){
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

    cout<<insertIdx<<" the inserted index"<<endl;// TODO CASE IT IS NOT RIGHT NEITHER LEFT BUT ACTUALLY ON!!
    //todo handle the case where this becomes a non manifold vertex and cut through immediately. but it destroys all patch ids
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

//        cout<<newnewVertIdx<<" newnew and new "<<newVertIdx<<", check with matrix sizes "<< newVg.rows()<<endl;
        MatrixXd newnewVg(newVg.rows()+1, 3);
        newnewVg.block(0,0, newVg.rows(), 3)= newVg;
        newnewVg.row(newnewVertIdx)= newVg.row(insertIdx);//+ (eps * toRight).transpose();

        handledVerticesSet.insert(newnewVertIdx);
//        cout<<newnewVg.row(insertIdx)<<" old insert idx"<<endl;
//        cout<<newnewVg.row(newnewVertIdx)<<" new insert idx duplicate "<<endl;

        newVg.resize(newnewVg.rows(), 3);
        newVg = newnewVg;
        cout<<" fin operation"<<endl;
//        cout<<newVg.row(newnewVertIdx)<<" newnew"<<endl;
//        cout<<newVg.row(insertIdx)<<" insert"<<endl;
//        cout<<newVg.row(newVertIdx)<<" new"<<endl;
//        cout<<newVg.row(cve->vert)<<" vert"<<endl;

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
        addToDulicateList(cve,seamsList, minusOneSeams, newVertIdx );
        toPattern_boundaryVerticesSet.insert(newVertIdx);// the duplicate is also on the boundary, hence insert it
        cve->leftdirection = toLeft;
        cve->rightdirection = toRight;
        cve-> leftCorner =  cve-> vert;
        cve-> rightCorner = newVertIdx;
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
void splitCounterPart(vector<cutVertEntry*>& cutPositions, int idxOfCVE,  cutVertEntry*& cve, MatrixXd& Vg, MatrixXi& Fg, vector<vector<int> >& vfAdj,
                      std::vector<std::vector<int> >& boundaryL,  vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams,
                      map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  MatrixXd& lengthsCurr,
                      MatrixXi& Fg_pattern, set<int> & cornerSet,  set<int>& handledVerticesSet){
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
                            boundaryL, seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet, lengthsCurr, cornerSet, handledVerticesSet);
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
                        boundaryL, seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet, lengthsCurr, cornerSet, handledVerticesSet);
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
//     cout<<"searching for index to open parallel of "<<cornerInitial<<endl;
    if(seamType<0){
//        cout<<" seam Type -1"<<endl;
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
//    cout<<searchedVert<<" the searched vert "<<endl;

    int size =  cornerToSeams[searchedVert].size();
    if( size != 2) {
        cout<<searchedVert<<" the size is not 2, it's "<<size<<endl;
    }
    int otherSeamId;
    if(cornerToSeams[searchedVert][0] == seamType){

        otherSeamId = cornerToSeams[searchedVert][1];
//        cout<<" other seam id is "<<otherSeamId<<endl;

    }else{
        otherSeamId = cornerToSeams[searchedVert][0];
    }

    for(int i=0; i<cutPositions.size(); i++){
        if(cutPositions[i]->vert == searchedVert || cutPositions[i]->cornerInitial == searchedVert){
            // then we should split herre!!

//            cout<<"we should split cut position "<<cutPositions[i]->vert<<" next, even if it has stress "<<cutPositions[i]->stress<<endl;
            return i;

        }
    }
//    for(int i=0; i< cornerToSeams[searchedVert].size(); i++){
////        cout<<cornerToSeams[searchedVert][i]<<" ";
//    }
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
//        cout<<cve->vert<<" "<< cutPositions[idx]->vert<<" matched"<<endl;

//        cout<<"found counter of middle, it's "<< cve-> vert <<" and "<< cutPositions[cve->counterPartIdx]->vert<<endl;
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

//    cout<<cve->vert<<" "<< cutPositions[idx]->vert<<" matched"<<endl;
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
        MatrixXd& lengthsOrig, MatrixXd& lengthsCurr,const std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary,
        map<int,
        vector<pair<int, int>>>& seamIdPerCorner, vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,
        double tailor_lazyness, double minConstrained, vector <cutVertEntry*>& cutPositions, VectorXd& cornerVert, MatrixXd& Vg){

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
    map<int, int> trackCornerIds; // keeps track of the var Id per corner - to understand if it exists already or not
    map <int, cutVertEntry*>  mapVarIdToVertId;
    map <pair<int, int>, int> mapVertAndSeamToVar;

    // a map to track the gurobi id of a corner.  whenever we set a corner we add the corner id, so we know once we set the second and connect them
    GRBVar* cutVar = model.addVars(numVar, GRB_BINARY);

    for(int i = 0; i < edgesPerBoundary.size(); i++) {
        int boundSize = boundaryL[i].size();
        int startVarPatch = varCount;
        for (int j = 0; j < edgesPerBoundary[i].size(); j++) {

            // first is the absolute index, second the index wrt the boundary loop
            auto cornerPair = edgesPerBoundary[i][j];
            if(seamIdPerCorner.find(cornerPair.first) == seamIdPerCorner.end()) continue;

            vector<pair<int, int>> seamId = seamIdPerCorner[cornerPair.first];// all seams that start at this corner, this can be max 2

            if(seamId.size()>2) cout<<" something is veryy odd!! we have more than two seams for a corner. impossible."<<endl;

            for(int si = 0; si < seamId.size(); si++) {
                cout<<" info "<<seamId[si].first<<" "<<seamId[si].second<<endl;

                GRBLinExpr innerSumConstr = 0; // interior sum
                GRBLinExpr lSumConstr = 0; // left sum
                GRBLinExpr rSumConstr = 0; // right sum
                int startVarOfThis = varCount;
                int count = 0;
                double distStartToEnd = computeSeamLength(seamId[si], seamsList, minusOneSeams, boundaryL, Vg);
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
            cve-> leftCorner =  -1;
            cve-> rightCorner = -1;
            int nextVert, prevVert;
            double nextStress, prevStress ;

            bool inverted = false;

            if(seamId<0) inverted = seamsList[(-1)*(seamId+1)]->inverted;
             getPrevAndNextVertAndStress( seamType, seamId, vert, prevVert, nextVert, prevStress, nextStress,
                                              seamsList, minusOneSeams, boundaryL, Fg_pattern, lengthsOrig, lengthsCurr, vfAdj, inverted );

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

                    cve->continuedDirection = Vg.row(nextVert)- Vg.row(vert);
                }else{
                    cve-> endCorner = true;
                    cve->stress = prevStress;
                    //experiment
                    if(seamId<0 &&  !inverted) cve->stress = nextStress;
                    cve->continuedDirection = Vg.row(prevVert)- Vg.row(vert);
                }

                    cout<<endl<<(firstInSeam==vert) <<" first ? "<<inverted<<endl;
                    cout<<nextStress<<" nextStress and Vert "<<nextVert<<endl;
                    cout<<prevStress<<" prev stress and vert "<<prevVert<<endl;

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
                  std::vector<std::vector<int> >& boundaryL, MatrixXi& Fg_pattern, vector<vector<int>> & vfAdj, MatrixXd& lengthsCurr
                  ){
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
    }

}

void tearFurther(vector<cutVertEntry*>& cutPositions, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams,
                 map<int, pair<int, int>> & releasedVert, set<int>& toPattern_boundaryVerticesSet,  std::vector<std::vector<int> >& boundaryL,
                 set<int> & cornerSet, set<int>& handledVerticesSet,  bool& prevFinished, const bool & preferManySmallCuts ){
    cout<<endl<<endl<<"-----------------------"<<endl<<endl;
    //when releasing the boundary it can turn into a non manifold mesh. not sure if this causes further problems
    Eigen::MatrixXi B;
    bool isManifold = igl::is_vertex_manifold( Fg_pattern, B);

    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);

    updateStress( cutPositions, seamsList, minusOneSeams, boundaryL,  Fg_pattern, vfAdj, lengthsCurr);

    if(prevFinished || preferManySmallCuts){
        // if we want many small cuts we sort always and there is no need to finish a seam before handling the next one!
        cout<<"It's time to sort again"<<endl;
        sort(cutPositions.begin(), cutPositions.end(), []( cutVertEntry* &a,  cutVertEntry* &b) { return a->stress > b-> stress; });

    }
//
    int  count = 0 ;
    for(int i = 0; i < 1; i++){
//    cout<<"current count "<<count<<endl;
        int currVert = cutPositions[count]->vert;
        bool parallelFinFlag = true;
        int parallel;
        //it cannot find the first one if it has not been released yet. need tto rewrite it

        if(cutPositions[count]->finFlag ){

//            cout<<"fin initial fin flag, check other "<<endl;
            // if it is a corner and it has been released
            if((cutPositions[count]->startCorner || cutPositions[count]->endCorner) &&
            releasedVertNew.find( cutPositions[count]->cornerInitial) != releasedVertNew.end()){
//                cout<<count<<" it was released, other exists  cornerInitial "<<cutPositions[count]-> cornerInitial<<endl;
                parallel = openParallelPosition(cutPositions[count]-> cornerInitial, releasedVertNew[cutPositions[count]->cornerInitial], seamsList, cutPositions);
//                cout<<"parallel idx is "<<parallel<<endl;
                if(parallel >= 0) parallelFinFlag = cutPositions[parallel]->finFlag;
                if(parallelFinFlag){
                    i--;
                    count ++;
//                    cout<<"other finished too"<<endl;
                }else{
//                    cout<<"the other is not finished yet. Continue cutting there"<<endl;

                    splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                                       minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, cornerSet, handledVerticesSet);
                    prevFinished = cutPositions[parallel]->finFlag;
//                    cout<<"finished  p ? "<<cutPositions[parallel]->finFlag<<endl;
                }
            }else{
//                cout<<"it has not been released, thus no parallel"<<endl;
                i--;
                count ++;
            }

        }else{

            cout<<endl<< cutPositions[count]->vert<<" vertex up next handling with i= "<<count<<" /"<<cutPositions.size()<<endl;
            if(!cutPositions[count]->finFlag){
                cout<<"split origninal "<<endl;
                splitVertexFromCVE(cutPositions[count], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList, minusOneSeams, releasedVert,
                                   toPattern_boundaryVerticesSet,lengthsCurr, cornerSet, handledVerticesSet );
                cout<<"finished ? "<<cutPositions[count]->finFlag<<endl;
                prevFinished = cutPositions[count] -> finFlag;
            }

//         int counterPart = cutPositions[count]->counterPartIdx;
//        if(counterPart >0){
//            cout<<endl<<" and corresponding "<<cutPositions[count]->counterPartIdx<<" "<<cutPositions[cutPositions[count]->counterPartIdx]->vert<<" with stress "<<cutPositions[cutPositions[count]->counterPartIdx]->stress<<endl;
//            splitVertexFromCVE(cutPositions[counterPart], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
//                             minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, Fg_pattern, cornerSet, handledVerticesSet);
//
//        }else {
//            cout<<"-1 seam, no counter to split"<<endl;
//        }
            parallel = -1;
            if(cutPositions[count]->startCorner || cutPositions[count]->endCorner){
                // once we finished cutting one side, check if it was a side opening. If so we can go on with the other side
                parallel = openParallelPosition(cutPositions[count]-> cornerInitial, releasedVertNew[cutPositions[count]-> cornerInitial], seamsList, cutPositions);
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
                                   minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, cornerSet, handledVerticesSet);
                cout<<"finished  p ? "<<cutPositions[parallel]->finFlag<<endl;
                prevFinished = (prevFinished && cutPositions[parallel] -> finFlag);

            }else{
                cout<<"parallel already finished"<<endl;
            }

        }

    }
    cout<<prevFinished<<" --------------------"<<endl;
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
void computeTear(Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern, MatrixXi& Fg_pattern_orig,
                 vector<seam*>& seamsList, vector<minusOneSeam*>& minusOneSeams, std::vector<std::vector<int> >& boundaryL, bool & finished,
                 const std::vector<std::vector<std::pair<int, int>>>& cornersPerBoundary, map<int, vector<pair<int, int>>>& seamIdPerCorner,
                 VectorXd& cornerVert, vector<cutVertEntry*>& cutPositions, map<int, pair<int, int>> & releasedVert,
                 set<int>& toPattern_boundaryVerticesSet, set<int> & cornerSet, set<int>& handledVerticesSet , MatrixXd& Vg, bool& prevFinished )
                 {

    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    MatrixXd lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    igl::edge_lengths(fromPattern, Fg_pattern_orig, lengthsOrig);

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
    setLP(boundaryL, vfAdj, Fg_pattern, lengthsOrig, lengthsCurr, cornersPerBoundary, seamIdPerCorner,
          seamsList, minusOneSeams, tailor_lazyness, minConstrained, cutPositions,
          cornerVert, Vg);

    //  here we need to sort and check if handled already
    sort(cutPositions.begin(), cutPositions.end(), []( cutVertEntry* &a,  cutVertEntry* &b) { return a->stress > b-> stress; });
    for(int i =0; i < cutPositions.size(); i++){
        findCorrespondingCounterCutPosition(cutPositions, i, cutPositions[i], currPattern, Fg_pattern, vfAdj, boundaryL,
                                            seamsList, minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr,Fg_pattern, cornerSet, handledVerticesSet );

    }

    int count=0;
    prevFinished = false;
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
        splitVertexFromCVE(cutPositions[count], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                           minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, cornerSet, handledVerticesSet);
        prevFinished = cutPositions[count] -> finFlag;

        if( releasedVertNew.find(currVert) != releasedVertNew.end()){
            // dann können müssen wir ja auch die passende andere seite des cuts öffnen
            // achtung, das ist nicht die corresponding seam
            int parallel = openParallelPosition(cutPositions[count]-> cornerInitial, releasedVertNew[currVert], seamsList, cutPositions);
            if(parallel<0) {cout<<"no proper parallel found"<<endl;continue;}
            splitVertexFromCVE(cutPositions[parallel], currPattern, Fg_pattern, vfAdj, boundaryL, seamsList,
                               minusOneSeams, releasedVert, toPattern_boundaryVerticesSet,lengthsCurr, cornerSet, handledVerticesSet);
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
}

void updatePositionToIntersection(MatrixXd& p,int next, const MatrixXd& Vg_bound){

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
//    if(next==62|| next == 63){
//        cout<<next<<": "<<Vg_bound.row(0)<<"    , "<<Vg_bound.row(1)<<"    ,"<<Vg_bound.row(2)<<endl;
//        cout<<mint<<" [pos] "<<p.row(next)<<" , i="<<mini<<endl<<" min dist target"<<minDistTarget.transpose() <<endl;
//    }
    p.row(next) += stiffness * (minDistTarget.transpose()-p.row(next));
//    if(next==62){cout<<"updated to "<<p.row(next)<<endl<<endl; }

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
void projectBackOnBoundary(const MatrixXd & Vg_to, MatrixXd& p, const vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,  const MatrixXi& Fg_pattern,
                           const MatrixXi& Fg_pattern_orig, const std::vector<std::vector<int> >& boundaryL_toPattern, map<int, pair<int, int>> & releasedVert ,bool visFlag){

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
            int v1_otherSide_idx = (stP2.first-i);// % boundLen2;
            if(v1_otherSide_idx < 0) {
                v1_otherSide_idx +=boundLen2;
            }
            if(seamsList[j]->inverted) v1_otherSide_idx = (stP2.first + i) % boundLen2;
            int v1_otherSide = boundaryL_toPattern[stP2.second][v1_otherSide_idx];

            Vg_seam1to.row(i) = Vg_to.row(v1_oneSide);
            Vg_seam2to.row(i) = Vg_to.row(v1_otherSide);

        }


        // for each interior (=not corner) vertex of the new boundary we need to find the closest position on the polyline and map it there
        // todo never ever cut the corner or change the corner index
        pair<int, int> ends = currSeam->getEndCornerIds();
        int i1=0;
        int bsize = boundaryL_toPattern[stP1.second].size();
        int next = boundaryL_toPattern[stP1.second][stP1.first];
        pair<int, int> compPair = make_pair(1,j );

        while( next!= ends.first ){
            // it is not released, project on boundary
            if(releasedVert.find(next) == releasedVert.end()){
                updatePositionToIntersection( p, next,Vg_seam1to);
            }
            // else it is released somehow. But from which seam? If it is released from another seam then pull it to this boundary still
            else if( releasedVertNew[next] != j){
               updatePositionToIntersection( p, next,Vg_seam1to);

            }
            i1++;
            next = boundaryL_toPattern[stP1.second][(stP1.first + i1) % bsize];
        }
        // the last corner. Again if it is constrained from another side pull it to boundary, else ignore since handled by corner
        if(releasedVert.find(next) != releasedVert.end() && releasedVertNew[next] != j){
            updatePositionToIntersection( p, next,Vg_seam1to);
        }

        int i2 = 0;
        bsize = boundaryL_toPattern[stP2.second].size();
        int nextidx = stP2.first;
        next = boundaryL_toPattern[stP2.second][nextidx];
        compPair = make_pair(1,-j-1 );

        while( next!= ends.second ){
            // general case an interior vertex , if it is not constrained pull it to boundary
            if(releasedVert.find(next) == releasedVert.end() ){// && i2!=0
                updatePositionToIntersection( p, next,Vg_seam2to);
            } else if( releasedVertNew[next] != j){
                updatePositionToIntersection( p, next,Vg_seam2to);

            }

            i2++;
            nextidx = (stP2.first - i2) % (bsize);

            if(nextidx < 0) {nextidx += bsize;}
            if(seamsList[j]->inverted) nextidx = (stP2.first + i2) % bsize;
            next = boundaryL_toPattern[stP2.second][nextidx];
        }

        if(releasedVert.find(next) != releasedVert.end() && releasedVertNew[next] != j){
            updatePositionToIntersection( p, next,Vg_seam2to);
        }


        // also project all duplicates of interior cut vertices
        for(const auto & addedVert : currSeam->duplicates){
            updatePositionToIntersection(p, addedVert.second, Vg_seam1to);
        }
        for(const auto & addedVert : currSeam->duplicates2){
            updatePositionToIntersection(p, addedVert.second, Vg_seam2to);
        }

    }
    for(int j = 0; j < minusOneSeams.size(); j++){

        minusOneSeam* currSeam  = minusOneSeams[j];
        int patch = currSeam -> getPatch();
        int startVert = currSeam -> getStartVert();
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

        int i1=0;
        int next = boundaryL_toPattern[patch][(startidx+i1) % boundLen];
        pair<int, int> compPair = make_pair(-1,j );
        while(next != endVert){
            // general case, it is not released hence pull it to the boundary
            if(releasedVert.find(next) == releasedVert.end()){
                updatePositionToIntersection( p, next,Vg_seamto);
            }else if(releasedVertNew[next] != (-1)*(j+1)){
                // it is released but not from this seam,thus it has to stay on the projection
                updatePositionToIntersection( p, next,Vg_seamto);
            }
            i1++;
            next = boundaryL_toPattern[patch][( startidx + i1) % boundLen];
        }
        // it is released for another side hence we have to pull it to our side
        if(releasedVert.find(next) != releasedVert.end() && releasedVertNew[next] != (-1)*(j+1)){
            updatePositionToIntersection( p, next,Vg_seamto);

        }
        // also map all projections
        for(const auto & addedVert : currSeam -> duplicates){
            updatePositionToIntersection(p, addedVert.second, Vg_seamto);

        }

    }

}

void updatePatchId(vector<cutVertEntry*>& cutPositions, const std::vector<std::vector<int> >& boundaryLnew, vector<seam*>& seamsList, vector<minusOneSeam*> & minusOneSeams){
    cout<<" The number of patches has changed. Create a mapping from one to the other and update all cve "<<endl;
    map<int, int> mapVertToNewPatch;
    for(int i = 0; i < boundaryLnew.size(); i++){
//        cout<<"Patch "<<i<<endl;
        for(int j = 0; j < boundaryLnew[i].size(); j++){
//            cout<<boundaryLnew[i][j]<<" ";
            mapVertToNewPatch[boundaryLnew[i][j]] = i;
        }
    }
    for(int i = 0; i < cutPositions.size(); i++){
//        cout<<"Vert "<<cutPositions[i] -> vert<<" was on patch "<<cutPositions[i] -> patch;
        cutPositions[i] -> patch = mapVertToNewPatch[cutPositions[i] -> vert];
//        cout<<" now it's on "<<cutPositions[i] -> patch<<endl;
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
//    vec += b;


}
