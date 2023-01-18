//
// Created by Anna Maria Eggler on 17.01.23.
//

#include "postProcessing.h"
#include <Eigen/Dense>
#include <iostream>
#include "seam.h"
#include <igl/edge_lengths.h>
#include <igl/signed_distance.h>
#include "igl/boundary_loop.h"
#include "adjacency.h"
#include <cmath>
#include "igl/barycentric_interpolation.h"
//#include <igl/triangle/triangulate.h>
//#include <>


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

void backTo3Dmapping(MatrixXd& adaptedPattern, MatrixXi& adaptedPattern_faces, MatrixXd& perfectPattern, MatrixXi& perfectPattern_faces ,
                     MatrixXd& perfectPatternIn3d, MatrixXi& perfectPatternIn3d_faces, MatrixXd& adaptedPatternIn3d){
    //idea: we have with perfectPatternForThisShape the perfect pattern and also in 3d
    // since the adapted pattern is a subset of the perfect pattern, we can locate every vertex of adapted pattern in perfectPattern and apply it using barycentric coordinates in 3d
    // maybe we have to do some manual stitching later but that should be ok.

    VectorXd S; VectorXi I;//face index of smallest distance
    MatrixXd C,N;
    cout<<"Starting with signed distance "<<endl;
    igl::signed_distance(adaptedPattern, perfectPattern, perfectPattern_faces, igl::SIGNED_DISTANCE_TYPE_UNSIGNED, S, I, C, N);
    cout<<"Finished signed distance "<<endl;
    for(int i=3149; i<=3158; i++){
        cout<<I(i)<<" face of"<<i<<endl;
    }
    MatrixXd B(adaptedPattern.rows(), 3); // contains all barycentric coordinates
    for(int i=0; i< adaptedPattern.rows(); i++){
        VectorXd bary;
        auto face = perfectPattern_faces.row(I(i));
        igl::barycentric_coordinates(adaptedPattern.row(i), adaptedPattern.row(face(0)), adaptedPattern.row(face(1)),
                                     adaptedPattern.row(face(2)), bary);
        B.row(i) = bary;
    }
    cout<<"Got all barycentric coords"<<endl;
    igl::barycentric_interpolation(perfectPatternIn3d, perfectPatternIn3d_faces, B, I, adaptedPatternIn3d);
    cout<<"Got interpolation"<<endl;

}
