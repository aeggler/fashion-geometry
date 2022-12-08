//
// Created by Anna Maria Eggler on 08.12.22.
//
#include <Eigen/Dense>
#include "patternAdaption.h"
#include <iostream>
#include "seam.h"
#include <igl/edge_lengths.h>
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
    int faceIdRight;
    int edgeOfFaceIdRight;
    int duplicate;
    double stretchLeft;
    double stretchRight;
};

// push item to the back of the list
void push(Node** tail_ref, int vert_Id, int dupl, int face_id ,int whichEdge,  double stretch_Left){

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
        new_node -> stretchRight  = new_node-> prev->stretchLeft;
        new_node -> faceIdRight = new_node -> prev -> faceId;
        new_node-> edgeOfFaceIdRight = new_node -> prev -> edgeOfFaceId;
    }
}

int printList(Node * node){
    // this is the thereshold
    double maxStretch = 0.0001;
    int idxMax = -1;
    int currIdx = 0;
    while (node != NULL){
//        cout<<node->vertId<<" "<< node->stretchLeft + node->stretchRight<<endl;
        if(node->stretchLeft + node->stretchRight > maxStretch){
            maxStretch = node->stretchLeft + node->stretchRight;
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
    cout<<node->vertId<<" "<< node->stretchLeft + node->stretchRight<<endl;
    cout<<node -> faceId<<" "<< node -> edgeOfFaceId<<" "<< node -> duplicate<<" the metadata"<<endl;
    cout<<  node -> faceIdRight <<" "<< node-> edgeOfFaceIdRight <<" the right side"<<endl;
}

void computeTear(Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, std::vector<std::vector<int> >& boundaryL){
    cout<<" in tear function"<<endl;
    // idea: we iterate over the first patch and see where it breaks apart
// for practical reasons we start with a small patch instead of the first on e
    vector<int> boundary = boundaryL[4];
    cout<<boundary.size()<<" the size "<<endl;
    Node* head = NULL;
    Node* tail = NULL;

    MatrixXd lengthsOrig, lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    igl::edge_lengths(fromPattern, Fg_pattern, lengthsOrig);
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);

    VectorXd relativeStretch(boundary.size());
//    double lastSum = -1;
//    int tearIdx = -1;
//    int tearBoundIdxl = -1;
//    int tearBoundIdxr = -1;
    for(int i=0; i< boundary.size(); i++){
        // for each edge we check how stressed it is
        // get old length
        int idx1 = boundary[i];
        int idx2 = boundary[(i+1) % boundary.size()];
        int faceIdx =  adjacentFaceToEdge( idx1, idx2, -1, vfAdj );
//        cout<<idx1<<" "<<faceIdx<<" face index and row "<<Fg_pattern.row(faceIdx)<<endl;

        int faceidxv1, faceidxv2;
        for(int j=0; j<3; j++){
            if(Fg_pattern(faceIdx, j)== idx1) faceidxv1 = j;
            if(Fg_pattern(faceIdx, j)== idx2) faceidxv2 = j;

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
        double origLength = lengthsOrig(faceIdx, whichEdge);
        double newLength = lengthsCurr ( faceIdx, whichEdge);

        relativeStretch(i) = newLength/origLength;// assuming the new one is stretched it is certainly longer

        push(&tail, idx1, -1, faceIdx, whichEdge, relativeStretch(i));
        if(head == NULL){
            head = tail; // now head points to the first element
        }
//        if(i>0){
//
//            if(lastSum< relativeStretch(i)+ relativeStretch(i-1)){
//                lastSum = relativeStretch(i)+ relativeStretch(i-1);
//                tearIdx = idx1;
//                tearBoundIdxl = i-1;
//                tearBoundIdxr = i;
//            }
//        }

    }

    // need to add the right stretch of the last element
    head->stretchRight = tail->stretchLeft;
    head -> faceIdRight = tail -> prev -> faceId;
    head-> edgeOfFaceIdRight = tail -> prev -> edgeOfFaceId;

//    if(lastSum< relativeStretch(0)+ relativeStretch(boundary.size()-1)){
//        lastSum = relativeStretch(0)+ relativeStretch(boundary.size()-1);
//        tearIdx = boundary[0];
//        tearBoundIdxl = 0;
//        tearBoundIdxr = boundary.size()-1;
//
//    }
//    cout<<relativeStretch<<" the relative stretches of pattern 4 "<<endl;
//    cout<<lastSum <<" highest sum "<<endl;
//    cout<< tearIdx <<" tearIdx sum "<<endl;
//    cout<<tearBoundIdxl <<" tearBoundIdxl  "<<endl;
//    cout<<tearBoundIdxr <<" tearBoundIdxr  "<<endl;
    int whichTear = printList(head);
    cout<<endl;
    cout<<whichTear<<" Tearing index"<<endl;
    printDataOfIdx(head, whichTear);

}
