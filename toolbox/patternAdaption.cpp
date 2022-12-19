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

void addVarToModel (int vert, int nextVert, vector<vector<int>> & vfAdj, bool isConstrained, int& varCount, GRBVar* & cutVar,
                      MatrixXi& Fg_pattern,MatrixXd& lengthsOrig, MatrixXd& lengthsCurr, map <int, int> & mapVarIdToVertId ){

    int faceIdx = adjacentFaceToEdge(vert, nextVert, -1, vfAdj );
    int whichEdge = findWhichEdgeOfFace(faceIdx, vert, nextVert, Fg_pattern);
    double w_init = lengthsCurr(faceIdx, whichEdge)/lengthsOrig(faceIdx, whichEdge);
//    cout<<varCount<<"is constrained ? "<<isConstrained<<" ";
    if(isConstrained){
        // we add it with coefficient 0, it should be ignored
        cutVar[varCount].set(GRB_DoubleAttr_Obj, 0);
    }else{
        cutVar[varCount].set(GRB_DoubleAttr_Obj, w_init);
        //cout << w_init << endl;
    }
    mapVarIdToVertId[varCount]= vert;
    varCount++;


}

void setLP(std::vector<std::vector<int> >& boundaryL , vector<vector<int>> & vfAdj, MatrixXi& Fg_pattern,
MatrixXd& lengthsOrig, MatrixXd& lengthsCurr,const std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary, map<int, vector<pair<int, int>>>& seamIdPerCorner,
           vector<seam*>& seamsList, const vector<minusOneSeam*> & minusOneSeams,
            int tearPatch, int tearVert, int tearVertInBoundaryIndex ){

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();
    // Create an empty model
    GRBModel model = GRBModel(env);
    GRBLinExpr obj = 0.0;
    double minConstrained = 0.25;
    int varCount =0;
    int numVar=0; // whole boundary and all corners duplicate
    for(int i=0; i<edgesPerBoundary.size(); i++){
        numVar += edgesPerBoundary[i].size();
        numVar+= boundaryL[i].size();
    }
    map<int, int> trackCornerIds;
    map <int, int>  mapVarIdToVertId;
    // a map to track the gurobi id of a corner.  whenever we set a corner we add the corner id, so we know once we set the second and connect them

    GRBVar* cutVar = model.addVars(numVar, GRB_BINARY);
    //edgesPerBoundary.size();
    for(int i=0; i< edgesPerBoundary.size(); i++) {
        int boundSize = boundaryL[i].size();
        int startVarPatch = varCount;
        cout<<"Processing patch "<<i<<endl ;
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
                    if (seamId[si].first >= 0) {
                        seam *seam = seamsList[seamId[si].second];
                        auto startAndPatch = seam->getStartAndPatch1();
                        if (startAndPatch.second != i)
                            cout << " now in th e+1 seams the patch does not match where we are in the loop "<< startAndPatch.second << endl;
                        idx = startAndPatch.first;
                        vert = boundaryL[startAndPatch.second][idx];
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
                        if (relId == 0 || relId > minConstrained && relId < (1 - minConstrained)) {
                            isConstrained = false;
                        }else{
                            model.addConstr(cutVar[varCount] == 0);
                        }
                        addVarToModel(vert, boundaryL[i][(idx + 1) % boundSize], vfAdj, isConstrained, varCount, cutVar,
                                      Fg_pattern, lengthsOrig, lengthsCurr, mapVarIdToVertId );
                        if(varCount-1 != startVarOfThis+count){
                            cout<<varCount<<"---------------------counting is off -----------------"<< startVarOfThis+count<<" "<<count<<endl;
                        }
                        if (!isConstrained) {
                            lSumConstr += cutVar[startVarOfThis + count];
                            if (relId != 0) {
                                rSumConstr += cutVar[startVarOfThis + count];
                                innerSumConstr += cutVar[startVarOfThis + count];
                            }
                        }
                        idx = (idx + 1) % boundSize;
                        vert = boundaryL[i][idx];
                        count++;
                    }// handle the last, add it to the right sum

                    //  make sure duplicate corner is chosen only once, and if we found the last consider its duplicate with start of this patch

                    if(trackCornerIds.find(end) != trackCornerIds.end()){
                        model.addConstr(cutVar[ trackCornerIds[end]] + cutVar[varCount] <= 1);
                    }else{
                        trackCornerIds[end] = varCount;
                    }
                    rSumConstr += cutVar[varCount];
                    addVarToModel(vert, boundaryL[i][(idx + 1) % boundSize], vfAdj, false, varCount, cutVar, Fg_pattern,
                                  lengthsOrig, lengthsCurr, mapVarIdToVertId);



                } else {
                    // iterate in inverse direction
                    seam *seam = seamsList[(seamId[si].second ) * -1 -1];
                    cout<<" seam id "<<(seamId[si].second ) * -1 -1<<endl;
                    auto startAndPatch = seam -> getStartAndPatch2ForCorres();
                    int idx = startAndPatch.first;

                    if (startAndPatch.second != i) cout << " something is wrong, it should be in patch " << i << " but the seam is from patch  " << startAndPatch.second << endl;
                    int vert = boundaryL[startAndPatch.second][idx];
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
                        if (relId == 0 || relId > minConstrained && relId < (1 - minConstrained)) {
                            isConstrained = false;
                        }else{
                            model.addConstr(cutVar[varCount] == 0);
                        }
                        addVarToModel(vert, nextvert, vfAdj, isConstrained, varCount, cutVar, Fg_pattern, lengthsOrig,
                                      lengthsCurr, mapVarIdToVertId);
                        if (!isConstrained) {
                            lSumConstr += cutVar[startVarOfThis + count];
                            if (relId != 0) {
                                rSumConstr += cutVar[startVarOfThis + count];
                                innerSumConstr += cutVar[startVarOfThis + count];
                            }
                        }

                        vert = nextvert;
                        count++;
                        nextidx = (startAndPatch.first - (1 + count)) % boundSize;
                        if (nextidx < 0) nextidx += boundSize;
                        if (seam->inverted) nextidx = (startAndPatch.first + (1 + count)) % boundSize;
                        nextvert = boundaryL[startAndPatch.second][nextidx];

                    }// handle least element and set constraints
                    if(trackCornerIds.find(end) != trackCornerIds.end()){
                        model.addConstr(cutVar[ trackCornerIds[end]] + cutVar[varCount] <= 1);
                    }else{
                        trackCornerIds[end] = varCount;
                    }

                    addVarToModel(vert, boundaryL[startAndPatch.second][(idx + 1) % boundSize], vfAdj, false, varCount,
                                  cutVar, Fg_pattern, lengthsOrig, lengthsCurr, mapVarIdToVertId);
                    rSumConstr += cutVar[startVarOfThis + count];

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
//    cout<<"read them all "<<varCount<<endl;
//    cout<<mapVarIdToVertId.size()<<" the size of the map, should be num var which is "<<numVar<<endl;
    for(int i =0; i< numVar; i++){
//        cout<<(cutVar[i].get(GRB_DoubleAttr_X ) )<<" i "<<i<<endl;
        if( cutVar[i].get(GRB_DoubleAttr_X ) >0.99){
            cout<<"bigger "<< mapVarIdToVertId[i] <<endl;// ALRIGHT BUT NOW WE NEED A SEAM ID 
        }

    }


}
void computeTear(Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern, MatrixXi& Fg_pattern_orig,
                 vector<seam*>& seamsList,  const vector<minusOneSeam*> & minusOneSeams, std::vector<std::vector<int> >& boundaryL, bool & finished,
                 const std::vector<std::vector<std::pair<int, int>>>& edgesPerBoundary,  map<int, vector<pair<int, int>>>& seamIdPerCorner
                 ){
    cout<<" in tear function"<<endl;
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    MatrixXd lengthsOrig, lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    igl::edge_lengths(fromPattern, Fg_pattern_orig, lengthsOrig);

    // first we need to know where to tear, set up LP for this
    // information we need: stress. For stress we need lengths old and ned
    // for lengths old an new we need which edge of face
    // for which edge of face we need face
    // for face we need vertices
    // for vertices we need boundary loop
    int tearPatch, tearVert, tearVertInBoundaryIndex;
    setLP( boundaryL, vfAdj, Fg_pattern, lengthsOrig, lengthsCurr, edgesPerBoundary, seamIdPerCorner,seamsList, minusOneSeams, tearPatch, tearVert, tearVertInBoundaryIndex);


    // idea: we iterate over the first patch and see where it breaks apart
// for practical reasons we start with a small patch instead of the first on e
//    vector<int> boundary = boundaryL[4];
//    cout<<boundary.size()<<" the size "<<endl;
//    Node* head = NULL;
//    Node* tail = NULL;
//    int listLength = 0;
//
//    cout<<"creating list"<<endl;
//    for(int i=0; i< boundary.size(); i++){
//        // for each edge we check how stressed it is
//        // get old length
//        int idx1 = boundary[i];
//        int idx2 = boundary[(i+1) % boundary.size()];
//        int faceIdx =  adjacentFaceToEdge( idx1, idx2, -1, vfAdj );
//
//        int whichEdge =  findWhichEdgeOfFace(faceIdx, idx1, idx2, Fg_pattern);
//        // assuming the new one is stretched it is certainly longer
//        double stretch = lengthsCurr ( faceIdx, whichEdge)/lengthsOrig(faceIdx, whichEdge);
//        cout<<stretch<<" i="<<i<<" "<<idx1<<endl;
//
//        push(&tail, idx1, -1, faceIdx, whichEdge, stretch, listLength);
//
//        if(head == NULL){
//            head = tail; // now head points to the first element
//        }
//    }
//    // need to add the right stretch of the last element
//    head->prev = tail;
//    tail-> next = head;
//    cout<<" stitched list"<<endl;
//
//    int whichTear = printList(head, listLength, vfAdj);
//    cout<<endl;
//    cout<<whichTear<<" Tearing index"<<endl;
//    printDataOfIdx(head, whichTear);
//    cout<<currPattern.rows()<<"initial rows"<<endl;
//    splitVertex(&head,listLength, whichTear,currPattern, Fg_pattern,vfAdj, lengthsOrig,lengthsCurr);
//    cout<<listLength<<" list length"<<endl;
//    cout<<currPattern.rows()<<" after rows"<<endl;
//    std::vector<std::vector<int> > boundaryLnew;
//    igl::boundary_loop(Fg_pattern, boundaryLnew);
//    cout<<boundaryLnew[4].size()<<" bound size"<<endl;
//    if(boundaryLnew[4].size()==24){
//        for(int i=0; i<24; i++){
//            cout<<boundaryLnew[4][i]<<" ";
//        }
//    }cout<<endl;
//    if(printList(head, listLength, vfAdj)== -1){
//        cout<<"finished tearing "<<endl;
//        finished = true;
//    }
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
                           const MatrixXi& Fg_pattern_orig, const std::vector<std::vector<int> >& boundaryL_toPattern, const std::vector<std::vector<int> >& boundaryL ){
    // Idea:
    // iterate over every seam, create a triangle mesh of the original boundary of that seam
    // then for each boundary vertex, project it's current position to the boundary
    // update the current position towards the projected
    int numSeams = seamsList.size();
    for (int j = 0; j<numSeams; j++){
//        cout<<" handling seam "<<j<<endl;

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
            int v1_otherSide = boundaryL[stP2.second][v1_otherSide_idx];

            Vg_seam1to.row(i) = Vg_to.row(v1_oneSide);
            Vg_seam2to.row(i) = Vg_to.row(v1_otherSide);

        }
        // for each interior (=not corner) vertex of the new boundary we need to check the closest position
        // todo never ever cut the corner or change the corner index
        auto ends = currSeam->getEndCornerIds();
        int i1=1;
        int bsize = boundaryL[stP1.second].size();
        int next = boundaryL[stP1.second][(stP1.first+i1)% bsize];

        while( next!= ends.first ){
            updatePositionToIntersection( p, next,Vg_seam1to);
            i1++;
            next = boundaryL[stP1.second][(stP1.first+i1)% bsize];
        }
        int sizeOneSide = i1+1; // account for  0and last

        int i2 = 1;
        bsize = boundaryL[stP2.second].size();
        int nextidx = (stP2.first- i2) % bsize;
        if(nextidx < 0) nextidx += bsize;
        next = boundaryL[stP2.second][nextidx];

        while( next!= ends.second ){
            updatePositionToIntersection(p, next,Vg_seam2to);

            i2++;
            nextidx = (stP2.first - i2) % (bsize);

            if(nextidx < 0) {nextidx += bsize;}
            if(seamsList[j]->inverted) nextidx = (stP2.first + i2) % bsize;
            next = boundaryL[stP2.second][nextidx];
        }

    }
    for(int j=0; j<minusOneSeams.size(); j++){
//        cout<<"minus one seam "<<j<<endl;
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
        int next = boundaryL[patch][(startidx+i1) % boundLen];
        while(next != endVert){
            updatePositionToIntersection( p, next,Vg_seamto);
            i1++;
            next = boundaryL[patch][( startidx+i1 ) % boundLen];
        } //updatePositionToIntersection( p, next,Vg_seamto);

    }


}
