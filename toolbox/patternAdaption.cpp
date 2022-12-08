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



void computeTear(Eigen::MatrixXd & fromPattern, MatrixXd&  currPattern, MatrixXi& Fg_pattern,vector<seam*>& seamsList, std::vector<std::vector<int> >& boundaryL){
    cout<<" in tear function"<<endl;
    // idea: we iterate over the first patch and see where it breaks apart
// for practical reasons we start with a small patch instead of the first on e
    vector<int> boundary = boundaryL[4];
    cout<<boundary.size()<<" the size "<<endl;

    MatrixXd lengthsOrig, lengthsCurr;
    igl::edge_lengths(currPattern, Fg_pattern, lengthsCurr);
    igl::edge_lengths(fromPattern, Fg_pattern, lengthsOrig);
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);

    VectorXd relativeStretch(boundary.size());
    for(int i=0; i< boundary.size(); i++){
        // for each edge we check how stressed it is
        // get old length
        int idx1 = boundary[i];
        int idx2 = boundary[(i+1) % boundary.size()];
        int faceIdx =  adjacentFaceToEdge( idx1, idx2, -1, vfAdj );
        cout<<idx1<<" "<<faceIdx<<" face index and row "<<Fg_pattern.row(faceIdx)<<endl;
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

    }

    cout<<relativeStretch<<" the reltaive stretches of pattern 4 "<<endl;


}
