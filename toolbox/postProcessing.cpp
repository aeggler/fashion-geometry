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
#include "igl/adjacency_list.h"
//#define ANSI_DECLARATORS
//#define TRILIBRARY
//extern "C"{
#include "../triangle/triangle.h"

//};

//#include "igl/triangle/triangulate.h"

using namespace std;
using namespace Eigen;

//tood attention this can create degenerate triangles
void smoothBetweenVertices(MatrixXd& currPattern, MatrixXi& Fg_pattern, vector<int>& startAndEnd){

    vector<vector<int>> boundaryL;
    int patch =-1; int startIdx, endIdx;

    igl::boundary_loop(Fg_pattern, boundaryL);
    vector<vector<int> > vvAdj;
    igl::adjacency_list(Fg_pattern,vvAdj);
    for (int i = 0; i < boundaryL.size(); ++i) {
        for (int j = 0; j < boundaryL[i].size(); ++j) {
            if(boundaryL[i][j] == startAndEnd[0] || boundaryL[i][j] == startAndEnd[1] ) {
                if (patch == -1) {
                    patch = i;
                    startIdx = j;
                } else {
                    endIdx = j;
                }
            }
        }
    }

    auto dist = (endIdx - startIdx);
    cout<< "dist= "<<dist<<endl;

    if((boundaryL[patch].size() - dist) < dist){
        cout<<"changed"<<endl;
        auto tem= startIdx;
        startIdx = endIdx;
        endIdx = tem;
    }

    cout<<"all the ones in between "<<endl;
    int idx= startIdx;
    Vector3d R = currPattern.row(boundaryL[patch][startIdx]);
    Vector3d Q =  currPattern.row(boundaryL[patch][endIdx]);
    vector<int> furtherChecks;
    while(idx!= endIdx){

        idx++;
        idx %= boundaryL[patch].size();
        int vert = boundaryL[patch][idx];

        // we need to find out if one neighbor of the current vertex is also not on the boundary
        double t = (R-Q).dot(currPattern.row(vert).transpose()-Q)/((R-Q).dot(R-Q));
        Vector3d targetPos = Q+t*(R-Q);

        currPattern.row(vert) = targetPos;
        for(int j =0; j<vvAdj[vert].size(); j++){
            int neigh = vvAdj[vert][j];
//            cout<<" "<<neigh<<endl;

            Vector3d C = currPattern.row(neigh);
            auto sign = ((R(0)-Q(0) )*(Q(1)-C(1))-(R(1)-Q(1))*(Q(0)-C(0)));

            if(sign < -1){
                double tt = (R-Q).dot(C-Q)/((R-Q).dot(R-Q));
                Vector3d targetPosNeigh = Q+tt*(R-Q);
                cout<<neigh<<" Neigh is on on the other side of the boundary "<<sign<<endl;
                furtherChecks.push_back(neigh);
                currPattern.row(neigh) = targetPosNeigh;

            }
        }
    }
    while(furtherChecks.size() > 0){
        int next = furtherChecks[furtherChecks.size()-1];
        furtherChecks.pop_back();
        for(int j =0; j < vvAdj[next].size(); j++){
            int neigh = vvAdj[next][j];
            Vector3d C = currPattern.row(neigh);
            auto sign = ((R(0)-Q(0) )*(Q(1)-C(1))-(R(1)-Q(1))*(Q(0)-C(0)));

            if(sign < -1){
                double tt = (R-Q).dot(C-Q)/((R-Q).dot(R-Q));
                Vector3d targetPosNeigh = Q+tt*(R-Q);
                cout<<neigh<<" iterative neigh is on on the other side of the boundary, maybe creating degenerate triangles! "<<sign<<endl;
                furtherChecks.push_back(neigh);
                currPattern.row(neigh) = targetPosNeigh;

            }
        }
    }

}

// see https://github.com/libigl/libigl/blob/main/include/igl/triangle/triangulate.cpp
void triangulateFAKE(MatrixXd& V, MatrixXi& E, MatrixXd& H, string flags, MatrixXd& V2, MatrixXi& F2){
    Eigen::VectorXi VM,EM,VM2,EM2;
    // "Vertex markers must be empty or same size as V");
//    "Segment markers must be empty or same size as E");
    Eigen::MatrixXi E2;
    // Prepare the flags
    string full_flags = flags + "pz" + (EM.size() || VM.size() ? "" : "B");
    typedef Map< Matrix<double,Dynamic,Dynamic,RowMajor> > MapXdr;
    typedef Map< Matrix<int,Dynamic,Dynamic,RowMajor> > MapXir;
    triangulateio in;
    in.numberofpoints = V.rows();
    in.pointlist = (double*)calloc(V.size(),sizeof(double));
    {
        MapXdr inpl(in.pointlist,V.rows(),V.cols());
        inpl = V.template cast<double>();
    }

    in.numberofpointattributes = 0;
    in.pointmarkerlist = (int*)calloc(V.size(),sizeof(int)) ;
    for(unsigned i=0;i<V.rows();++i) in.pointmarkerlist[i] = VM.size()?VM(i):1;

    in.trianglelist = NULL;
    in.numberoftriangles = 0;
    in.numberofcorners = 0;
    in.numberoftriangleattributes = 0;
    in.triangleattributelist = NULL;

    in.numberofsegments = E.size()?E.rows():0;
    in.segmentlist = (int*)calloc(E.size(),sizeof(int));
    {
        MapXir insl(in.segmentlist,E.rows(),E.cols());
        insl = E.template cast<int>();
    }
    in.segmentmarkerlist = (int*)calloc(E.rows(),sizeof(int));
    for (unsigned i=0;i<E.rows();++i) in.segmentmarkerlist[i] = EM.size()?EM(i):1;

    in.numberofholes = H.size()?H.rows():0;
    in.holelist = (double*)calloc(H.size(),sizeof(double));
    {
        MapXdr inhl(in.holelist,H.rows(),H.cols());
        inhl = H.template cast<double>();
    }
    in.numberofregions = 0;

    // Prepare the output struct
    triangulateio out;
    out.pointlist = NULL;
    out.trianglelist = NULL;
    out.segmentlist = NULL;
    out.segmentmarkerlist = NULL;
    out.pointmarkerlist = NULL;

//     Call triangle/
     triangulate((char *)full_flags.c_str(), &in, &out, (triangulateio *) NULL);
//    triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);
//
//    // Return the mesh
//    V2 = MapXdr(out.pointlist,out.numberofpoints,2).cast< double>();
//    F2 = MapXir(out.trianglelist,out.numberoftriangles,3).cast< int>();
//    E2 = MapXir(out.segmentlist,out.numberofsegments,2).cast< int>();
//    if(VM.size())
//    {
//        VM2 = MapXir(out.pointmarkerlist,out.numberofpoints,1).cast<int>();
//    }
//    if(EM.size())
//    {
//        EM2 = MapXir(out.segmentmarkerlist,out.numberofsegments,1).cast< int>();
//    }



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
    string flags = "";
    MatrixXd V2;
    MatrixXi F2;
    triangulateFAKE(V, E, H, flags, V2, F2 );

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

    MatrixXd B(adaptedPattern.rows(), 3); // contains all barycentric coordinates
    for(int i=0; i< adaptedPattern.rows(); i++){
        VectorXd bary;
        auto face = perfectPattern_faces.row(I(i));
        igl::barycentric_coordinates(adaptedPattern.row(i), perfectPattern.row(face(0)), perfectPattern.row(face(1)),
                                     perfectPattern.row(face(2)), bary);
        B.row(i) = bary;
    }
    cout<<"Got all barycentric coords"<<endl;
    igl::barycentric_interpolation(perfectPatternIn3d, perfectPatternIn3d_faces, B, I, adaptedPatternIn3d);
    cout<<"Got interpolation"<<endl;

}
