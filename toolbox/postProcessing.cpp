//
// Created by Anna Maria Eggler on 17.01.23.
//

#include "postProcessing.h"
#include <Eigen/Dense>
#include <iostream>
#include "seam.h"
#include <igl/edge_lengths.h>
#include <igl/adjacency_list.h>
#include "igl/boundary_loop.h"
#include "adjacency.h"
#include <cmath>
//#include <igl/triangle/triangulate.h>


using namespace std;
using namespace Eigen;


void smoothBetweenVertices(MatrixXd& currPattern, MatrixXi& Fg_pattern, vector<int>& startAndEnd){
    vector<vector<int>> boundaryL;
    int patch ; int startIdx, endIdx;

    igl::boundary_loop(Fg_pattern, boundaryL);
    for (int i = 0; i < boundaryL.size(); ++i) {
        for (int j = 0; j < boundaryL[i].size(); ++j) {
            if(boundaryL[i][j] == startAndEnd[0]){
                patch = i;
                startIdx = j;
            }else if(boundaryL[i][j] == startAndEnd[1]){
                endIdx = j;
            }
        }
    }
    if(startIdx>endIdx){
        auto temp = startIdx;
        startIdx = endIdx;
        endIdx = temp;
    }
    auto dist = (endIdx - startIdx);
    cout<< "dist= "<<dist<<endl;

    if(fabs(startIdx- endIdx)<dist){
        cout<<"changed"<<endl;
        auto tem= startIdx;
        startIdx = endIdx;
        endIdx = tem;
    }

    cout<<"all the ones in between "<<endl;
    int idx= startIdx;
    Vector3d R = currPattern.row(boundaryL[patch][startIdx]);
    Vector3d Q =  currPattern.row(boundaryL[patch][endIdx]);

    while(idx!= endIdx){
        cout<<" "<<boundaryL[patch][idx]<<endl;
        idx++;
        idx%= boundaryL[patch].size();
        double t = (R-Q).dot(currPattern.row(boundaryL[patch][idx]).transpose()-Q)/((R-Q).dot(R-Q));
        Vector3d targetPos = Q+t*(R-Q);

        currPattern.row(boundaryL[patch][idx]) = targetPos;
    }


}
void startRetriangulation(vector<VectorXd>& polylineSelected){
    int n = polylineSelected.size();
    MatrixXd V (n, 2);
    for(int i=0; i<n; i++){
        V.row(i) = polylineSelected[i].leftCols(2);
    }
    MatrixXi E(n,2);
    for(int i=0; i < n; i++){
        E(i, 0) = i;
        E(i, 1) = (i+1) % n;
    }
    MatrixXd H;



}

