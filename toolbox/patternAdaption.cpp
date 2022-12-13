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
#include <igl/signed_distance.h>
#include "adjacency.h"

using namespace std;
using namespace Eigen;

class Node{
public:
    int vertId;
    Node *next;
    Node *prev;
    int faceId;
    int edgeOfFaceId;
//    int faceIdRight;
//    int edgeOfFaceIdRight;
    int duplicate;
    double stretchLeft;
//    double stretchRight;
};

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

void splitVertex(Node** head, int & listLength, int  whichTear, MatrixXd& Vg, MatrixXi& Fg, vector<vector<int> >& vfAdj, MatrixXd& lengthsOrig, MatrixXd& lengthsCurr){

//    Node splitNodetemp  = getDataOfIdx((*head), whichTear);
//    Node* splitNode = &splitNodetemp;
    cout<<"splitting "<<endl;
//    cout<<splitNode->vertId<<" "<< splitNode->stretchLeft + splitNode->prev->stretchLeft<<endl;
//    cout<<splitNode -> faceId<<" "<< splitNode -> edgeOfFaceId<<" "<< splitNode -> duplicate<<" the left side"<<endl;
//    cout<<  splitNode ->prev->faceId  <<" "<< splitNode-> prev->edgeOfFaceId  <<" the right side"<<endl;


    Node* node = (*head);
    int currIdx = 0;
    while (currIdx<whichTear){
        node = node-> next;
        currIdx++;
    }
    whichTear = node->vertId;

    // to which vert it the vertex adjacent? from there we get all the adjacent edges
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
    // we need to know which is closest to normal to insert this and new vertex and in right order
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
//        cout<<dist <<" "<<(midVec - edgeVec).norm()<<" for vert"<<adjVert<<endl;
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
    // if it coudl then it is a) not stressed or b) we found no splitnormal

    // right gets  newVertIdx
    double dRight = (Vg(rightId, 0) - x1) * (y2-y1) - (Vg(rightId, 1) - y1) * (x2 - x1);
    // defines the side of right, left is the opposite
    bool rightDSmaller = (dRight<0);

//    if(whichTear== 2898)
//        cout<<dRight<<" right of 2898  smaller 0 ? "<<rightDSmaller<<endl;

    if(Fg(node->prev->faceId , 1) == whichTear){
        Fg(node->prev->faceId, 1) = newVertIdx;
    }
    else if(Fg(node-> prev->faceId , 0)== whichTear){
        Fg(node->prev->faceId, 0) = newVertIdx;
    }
    else {
        Fg(node->prev->faceId, 2) = newVertIdx;
    }
    cout<<Fg.row(5319)<<" cout TEST FACE CHANGED 5319"<<endl;
    cout<<Fg.row(5328 )<<" cout TEST FACE CHANGED 5328"<<endl;


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

//    int nextToTear = printList((*head), listLength);
// cout<<nextToTear<<" fin split function "<<endl;

}
void computeTear(Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern,MatrixXi& Fg_pattern_orig,
                 vector<seam*>& seamsList, std::vector<std::vector<int> >& boundaryL, bool & finished){
    cout<<" in tear function"<<endl;

    // idea: we iterate over the first patch and see where it breaks apart
// for practical reasons we start with a small patch instead of the first on e
    vector<int> boundary = boundaryL[4];
    cout<<boundary.size()<<" the size "<<endl;
    Node* head = NULL;
    Node* tail = NULL;
    int listLength = 0;
    MatrixXd lengthsOrig, lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    igl::edge_lengths(fromPattern, Fg_pattern_orig, lengthsOrig);
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    cout<<"creating list"<<endl;
    for(int i=0; i< boundary.size(); i++){
        // for each edge we check how stressed it is
        // get old length
        int idx1 = boundary[i];
        int idx2 = boundary[(i+1) % boundary.size()];
        int faceIdx =  adjacentFaceToEdge( idx1, idx2, -1, vfAdj );

        int whichEdge =  findWhichEdgeOfFace(faceIdx, idx1, idx2, Fg_pattern);
        // assuming the new one is stretched it is certainly longer
        double stretch = lengthsCurr ( faceIdx, whichEdge)/lengthsOrig(faceIdx, whichEdge);
        cout<<stretch<<" i="<<i<<" "<<idx1<<endl;

        push(&tail, idx1, -1, faceIdx, whichEdge, stretch, listLength);

        if(head == NULL){
            head = tail; // now head points to the first element
        }
    }
    // need to add the right stretch of the last element
    head->prev = tail;
    tail-> next = head;
cout<<" stitched list"<<endl;
    int whichTear = printList(head, listLength, vfAdj);
    cout<<endl;
    cout<<whichTear<<" Tearing index"<<endl;
    printDataOfIdx(head, whichTear);
    cout<<currPattern.rows()<<"initial rows"<<endl;
    splitVertex(&head,listLength, whichTear,currPattern, Fg_pattern,vfAdj, lengthsOrig,lengthsCurr);
    cout<<listLength<<" list length"<<endl;
    cout<<currPattern.rows()<<" after rows"<<endl;
    std::vector<std::vector<int> > boundaryLnew;
    igl::boundary_loop(Fg_pattern, boundaryLnew);
    cout<<boundaryLnew[4].size()<<" bound size"<<endl;
    if(boundaryLnew[4].size()==24){
        for(int i=0; i<24; i++){
            cout<<boundaryLnew[4][i]<<" ";
        }
    }cout<<endl;
    if(printList(head, listLength, vfAdj)== -1){
        cout<<"finished tearing "<<endl;
        finished = true;
    }

}

void updatePositionToIntersection(MatrixXd& p,int next,  const MatrixXd& Vg_bound, const MatrixXi &Fg_bound){
    VectorXd S, I; MatrixXi C; MatrixXd N;
    igl::signed_distance(p.row(next), Vg_bound, Fg_bound, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N );
    // we know that the first two indices of the face define the egdge we want to intersect with
    Vector3d R = Vg_bound.row(Fg_bound(I(0), 0));
    Vector3d Q =  Vg_bound.row(Fg_bound(I(0), 1));
    // derive where QR and P meet = t, https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line
    double t = (R-Q).dot(Q-p.row(next).transpose())/((R-Q).dot(R-Q));
    Vector3d targetPos = Q-t*(R-Q);
    double stiffness = 0.8; //todo
    p.row(next)+= stiffness * (targetPos.transpose()-p.row(next));
}

/*we build a new structure to find the closest position on the original boundary
        *     v0*-------------------*v1*---------------------*v2
        *           - v_(len+1)-          - v_(len+1+1) -
        *
        *           where v_(len+i) are defined as (vi + v(i+1))/2 + eps* n w
        *           where n is the normal
        *
*/
void projectBackOnBoundary(const MatrixXd & Vg_to, MatrixXd& p, const vector<seam*>& seamsList, const MatrixXi& Fg_pattern,
                           const MatrixXi& Fg_pattern_orig, const std::vector<std::vector<int> >& boundaryL_toPattern, const std::vector<std::vector<int> >& boundaryL ){
    // Idea:
    // iterate over every seam, create a triangle mesh of the original boundary of that seam
    // then for each boundary vertex, project it's current position to the boundary
    // update the current position towards the projected
    int numSeams = seamsList.size();
    for (int j = 0; j<numSeams; j++){
        cout<<" handling seam "<<j<<endl;

        seam* currSeam  = seamsList[j];
        auto stP1= currSeam-> getStartAndPatch1();
        auto stP2 =  currSeam -> getStartAndPatch2ForCorres(); // attention this is with respect to the original pattern
        int len  = currSeam -> seamLength();
        int boundLen1 = boundaryL_toPattern[stP1.second].size();
        int boundLen2 = boundaryL_toPattern[stP2.second].size();

        // build the structure for closest search
        MatrixXd Vg_seam1to(len+1 + (len), 3);
        MatrixXd Vg_seam2to(len+1 + len, 3);
        MatrixXi Fg_seamto(len, 3);
//        MatrixXi Fg_seam2to(len, 3);

        for(int i = 0; i<= len ; i++){
            int v1_oneSide = boundaryL_toPattern[stP1.second][(stP1.first+i)% boundLen1];
            int v1_otherSide_idx = (stP2.first-i)% boundLen2;
            if(v1_otherSide_idx < 0) {
                v1_otherSide_idx +=boundLen2;
            }
            if(seamsList[j]->inverted) v1_otherSide_idx = (stP2.first + i) % boundLen2;
            int v1_otherSide = boundaryL[stP2.second][v1_otherSide_idx];

            Vg_seam1to.row(i) = Vg_to.row(v1_oneSide);
            Vg_seam2to.row(i) = Vg_to.row(v1_otherSide);

            if(i>0){
                Fg_seamto(i-1, 0)= i;
                Fg_seamto(i-1, 1)= i-1;
                Fg_seamto(i-1, 2)= len + i;// from len+1 + (i-1)

                Vg_seam1to.row(len+i ) = (Vg_seam1to.row(i)+  Vg_seam1to.row(i-1)) / 2; // still need the normal
                Vg_seam2to.row(len+i ) = (Vg_seam2to.row(i)+  Vg_seam2to.row(i-1)) / 2; // still need the normal
                // todo very ugly, just pick one side, don't care which ... this is not 100% correct
                Vector3d d =  (Vg_seam1to.row(i)-  Vg_seam1to.row(i-1));
                Vector3d normal = d; normal(0)= -d(1); normal(1)= d(0);
                Vg_seam1to.row(len+i )+= 0.001*normal;

                d =  (Vg_seam2to.row(i)-  Vg_seam2to.row(i-1));
                normal = d; normal(0)= -d(1); normal(1)= d(0);
                Vg_seam2to.row(len+i )+= 0.001*normal;
            }
        }
        cout<<" set up both helper meshes"<<endl;

        // for each vertex of the new boundary we need to check the closest position
        // todo never ever cut the corner or change the corner index
        auto ends = currSeam->getEndCornerIds();
        int i1=0;
        int bsize = boundaryL[stP1.second].size();
        int next = boundaryL[stP1.second][(stP1.first+i1)% bsize];

        while( next!= ends.first ){
            updatePositionToIntersection( p, next,Vg_seam1to,Fg_seamto);
            i1++;
            next = boundaryL[stP1.second][(stP1.first+i1)% bsize];
        }// and one more!!
        updatePositionToIntersection(p, next,Vg_seam1to,Fg_seamto);
        int sizeOneSide = i1+1; // account for  0and last


        int i2 = 0;
        bsize = boundaryL[stP2.second].size();
        int nextidx = (stP2.first) % bsize;
        if(nextidx < 0) nextidx += bsize;
        next = boundaryL[stP2.second][nextidx];

        while( next!= ends.second ){
            updatePositionToIntersection(p, next,Vg_seam2to,Fg_seamto);

            i2++;
            nextidx = (stP2.first - i2) % (bsize);

            if(nextidx < 0) {nextidx += bsize;}
            if(seamsList[j]->inverted) nextidx = (stP2.first + i2) % bsize;
            next = boundaryL[stP2.second][nextidx];
        }// and one more!!
        updatePositionToIntersection(p, next,Vg_seam2to,Fg_seamto);
        cout<<endl<<len+1<<" len "<<sizeOneSide<<" size of one, other "<<i2+1<<" size of other side "<<endl;

    }




}
