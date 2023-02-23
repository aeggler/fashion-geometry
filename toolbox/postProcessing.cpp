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
#include "../triangle/triangle.h"
#include <fstream>
#include <iterator>
#include <set>
#include "MathFunctions.h"
#include "seam.h"
#include <igl/vertex_components.h>
#include <igl/internal_angles.h>
#include <igl/writeOBJ.h>



using namespace std;
using namespace Eigen;

//Taubin smoothing according to
//https://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/06_smoothing.pdf
// from start to end in direction of the middle vertex

void smoothBetweenVertices(MatrixXd& currPattern, MatrixXi& Fg_pattern, vector<int>& startAndEnd){
    cout<<startAndEnd.size()<<" start and end size"<<endl;
    vector<vector<int>> boundaryL;
    int patch =-1; int startIdx = -1; int midIdx = -1; int endIdx = -1;

    igl::boundary_loop(Fg_pattern, boundaryL);
    vector<vector<int> > vvAdj;
    igl::adjacency_list(Fg_pattern,vvAdj);
    for (int i = 0; i < boundaryL.size(); ++i) {
        for (int j = 0; j < boundaryL[i].size(); ++j) {
            if(boundaryL[i][j] == startAndEnd[0] ) {
                patch = i;
                startIdx = j;
            }
        }
    }
    if(patch ==-1) cout<<"error vertex not found"<<endl ;
    bool otherDir= false;
    vector<int> boundary  = boundaryL[patch];
    int loopSize = boundary.size();

    for (int j = 0; j < loopSize ; ++j) {
        if( boundary[j % loopSize] == startAndEnd[2] ) {
            endIdx = j;
        }else if(boundary[j % loopSize] == startAndEnd[1] ){
            midIdx = j;
        }
    }

    if(startIdx < endIdx){
        if(midIdx<endIdx && midIdx>startIdx){
            otherDir= false;
        }else{
            otherDir = true;
        }
    }else{
        //starrt>mid
        if(midIdx<startIdx && midIdx>endIdx){
            otherDir = true;
        }
    }

    if(startIdx==-1 || midIdx ==-1 || endIdx ==-1){
        cout<<"Something is -1. stopping here."  <<endl; return ;
    }
    double lamda = 0.1;
    double mu = -0.1;
    int iterations = 100;
    for (int i = 0; i < iterations; i++){

        int curr, next, prev;
        // shirnking wiht lamnda
        curr = startIdx;
        next = (otherDir)? (curr-1) : ((curr +1)%loopSize);
        if(next<0)next+= loopSize;
        while(next!= endIdx){
            prev= curr;
            curr= next;
            next = (!otherDir)? (curr+1 )%loopSize : (curr -1);
            if(next<0)next+= loopSize;

            VectorXd deltaP= (0.5 * (currPattern.row(boundary[prev]) + currPattern.row(boundary[next])) - currPattern.row(boundary[curr])).transpose();
            currPattern.row(boundary[curr]) += lamda * deltaP.transpose();

        }

        //enlarging with mu
        curr = startIdx;
        next = (!otherDir)? ((curr+1 )%loopSize) : (curr -1);
        if(next<0)next+= loopSize;
        while(next!= endIdx){
            prev= curr;
            curr= next;
            next = (!otherDir)? (curr+1 )%loopSize : (curr -1);
            if(next<0)next+= loopSize;

            VectorXd deltaP= (0.5*(currPattern.row(boundary[prev])+ currPattern.row(boundary[next])) - currPattern.row(boundary[curr])).transpose();
            currPattern.row(boundary[curr]) += mu * deltaP.transpose();

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
     cout<<" calling triangle"<<endl;
     triangulate((char *)full_flags.c_str(), &in, &out, (triangulateio *) NULL);
//    // Return the mesh
    V2 = MapXdr(out.pointlist,out.numberofpoints,2).cast< double>();

    F2 = MapXir(out.trianglelist,out.numberoftriangles,3).cast< int>();
    E2 = MapXir(out.segmentlist,out.numberofsegments,2).cast< int>();
    if(VM.size())
    {
        VM2 = MapXir(out.pointmarkerlist,out.numberofpoints,1).cast<int>();
    }
    if(EM.size())
    {
        EM2 = MapXir(out.segmentmarkerlist,out.numberofsegments,1).cast< int>();
    }
    MatrixXd V2new = MatrixXd::Constant(V2.rows(), 3, 200);
    V2new.block(0,0,V2.rows(), 2) = V2;
    V2.resize(V2.rows(), V2new.cols());
    V2 = V2new;
}
void startRetriangulation(vector<VectorXd>& polylineSelected, MatrixXd& V2, MatrixXi& F2 ){
    cout<<endl<<endl;
    int n = polylineSelected.size();
    MatrixXd V (n, 2);
    for(int i=0; i<n; i++){
        V.row(i) = polylineSelected[i].transpose().leftCols(2);
    }
    MatrixXi E(n,2);
    for(int i=0; i < n; i++){
        E(i, 0) = i;
        E(i, 1) = (i+1) % n;
    }
    MatrixXd H;
    string flags = "qa234.64";
    cout<<" starting triangulation with "<<E.rows()<<" edges and "<<V.rows()<<" points"<<endl;

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

bool vertOnEdge(const VectorXd& R, const VectorXd& Q, VectorXd& p,int v, int v1){
    double eps = 0.5;
    auto QR = Q-R;
    auto Qp = Q-p;
//    cout<<"v= "<<v<<endl;
    VectorXd diff = (Qp).normalized() - (QR).normalized();
//    if((v==2967 && v1 ==2968) ||(v==2968 && v1 ==2967) ){
//        cout<<diff.transpose()<<" diff"<<endl;
//        double tt = (p-R)(0)/(QR)(0);
//        double t2 =  (p-R)(1)/(QR)(1);
//        cout<<tt<<" = t=  "<<t2<<", and makes "<<endl<<(R+tt*(Q-R)).transpose()<<endl<<(R+t2*(Q-R)).transpose()<<endl<<p.transpose()<<" =? "<<endl;
//    }
    if (abs(diff(0))+abs(diff(1)) >eps) return false;

    double t = (p-R)(0)/(QR)(0);
    double tt=  (p-R)(1)/(QR)(1);
// numerically instable as quite often the x or y corrdinate is the same... compute both and take better choice
    double finalT;
    if (0 > t || t > 1 ) {// t does not work anyways ,try with tt
        finalT =tt;
    }
    if (0> finalT || finalT>1 ) return false;
    // both are possible, choose the better
    if( ((R+t*(QR))-p).norm()< ((R+tt*(QR))-p).norm())
    {
        finalT = t;
    }else{
        finalT = tt;
    }
    VectorXd posdiff = R+finalT*(Q-R) - p;
    cout<<"diff from real pos " << abs(posdiff(0)) + abs(posdiff(1)) <<endl;
//    if(v==1556|| v==1557){
//        cout<<posdiff.transpose()<<" posdiff"<<endl;
//    }
    return (abs(posdiff(0)) + abs(posdiff(1)) < 1);


}
/*
 * V0 _________________ v9, v10, v11
 * |                 |V8
 * |                 |
 * |V1               |
 * |                 |V7
 * |                 |
 * |V2___v3___v4__v5_|V6
 * we assume V0-2 to be on the current pattern (patch 1), and V6-8 patch2
 * further v2-4 & 9-11 are on patch c on the vg_to pattern
 * if for some the vertices are the same, there is just not enough between and we're done, no need to traverse
 * */
bool isAsc(int s, int m, int t){
    if(m==t || m==s){
        // there is only tow
        if(!(s==0 || t==0)){
            return (s < t);
        }else{
            if(s==0){
                if ( t==1) {return true;}else {return false;}
            }else{
                // t==0
                if(s==1){
                    return false;
                }else{
                    return true;
                }
            }
        }
    }
    if(s< t){
        if (s < m && m < t){
            return true;
        }else{
            return false;
        }
    } else{
        if( s > m && m > t){
            return false;
        }else{
            return true;
        }
    }
}
void computeAllBetweensConnectPatches(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                           vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                           MatrixXd& currPattern, MatrixXd& Vg_to, vector<VectorXd>& polyLineInput, vector<vector<int>>& connectedVertVec, vector<int>& patchId, vector<bool>& isAscVec) {
    VectorXi idx(12);
    patchId.clear();
    isAscVec.clear();
    connectedVertVec.clear();

    int patchL, patchR, patchTo;
    for(int i = 0; i < boundaryL_adaptedFromPattern.size(); i++){
        vector<int> boundary = boundaryL_adaptedFromPattern[i];
        for(int j = 0; j < boundary.size(); j++){
            if(boundary[j]== polylineIndex[0]){
                idx(0) = j;
                patchL= i;
            }
            if(boundary[j]== polylineIndex[1]){
                idx(1) = j;
            }
            if(boundary[j]== polylineIndex[2]){
                idx(2) = j;
            }

            if(boundary[j]== polylineIndex[6]){
                idx(6) = j;
                patchR= i;
            }
            if(boundary[j]== polylineIndex[7]){
                idx(7) = j;
            }
            if(boundary[j]== polylineIndex[8]){
                idx(8) = j;
            }
        }
    }
    bool leftAsc = isAsc(idx(0), idx(1), idx(2));
    bool rightAsc = isAsc(idx(6), idx(7), idx(8));
    patchId.push_back(patchL); patchId.push_back(patchR);
    isAscVec.push_back(leftAsc);
    isAscVec.push_back(rightAsc);
    for(int i = 0; i < boundaryL_toPattern.size(); i++){
        vector<int> boundary = boundaryL_toPattern[i];
        for(int j = 0; j < boundary.size(); j++){
            if(boundary[j]== polylineIndex[3]){
                idx(3) = j;
                patchTo= i;
            }
            if(boundary[j]== polylineIndex[4]){
                idx(4) = j;
            }
            if(boundary[j]== polylineIndex[5]){
                idx(5) = j;
            }

            if(boundary[j]== polylineIndex[9]){
                idx(9) = j;
            }
            if(boundary[j]== polylineIndex[10]){
                idx(10) = j;
            }
            if(boundary[j]== polylineIndex[11]){
                idx(11) = j;
            }
        }
    }
    bool toAsc = isAsc(idx(3), idx(4), idx(5));
    cout<<idx.transpose()<<" the indices "<<endl;

    int currIdx = idx(0);
    vector<int> boundary = boundaryL_adaptedFromPattern[patchL];
    vector<int> connVert;
    while(boundary[currIdx]!= boundary[idx(2)]){
        //add them
        polyLineInput.push_back(currPattern.row(boundary[currIdx]).transpose());
        connVert.push_back(boundary[currIdx]);
        cout<<boundary[currIdx]<<" first"<<endl;

        currIdx = (leftAsc)? (currIdx+1) % boundary.size() : currIdx-1;
        if(currIdx<0) currIdx+= boundary.size();

    }// add last
    polyLineInput.push_back(currPattern.row(boundary[currIdx]).transpose());
    connVert.push_back(boundary[currIdx]);
    connectedVertVec.push_back(connVert);
    cout<<boundary[currIdx]<<" first"<<endl;

    // go on on boundary of to pattern
    currIdx = idx(3);
    boundary.clear();
    boundary = boundaryL_toPattern[patchTo];
    while(boundary[currIdx]!= boundary[idx(5)]){
        //add them
        polyLineInput.push_back(Vg_to.row(boundary[currIdx]).transpose());
        cout<<boundary[currIdx]<<" sec"<<endl;

        currIdx = (toAsc)? (currIdx+1) % boundary.size() : currIdx-1;
        if(currIdx<0) currIdx+= boundary.size();

    }// add last
    polyLineInput.push_back(Vg_to.row(boundary[currIdx]).transpose());
    cout<<boundary[currIdx]<<" sec"<<endl;

    currIdx = idx(6);
    boundary.clear();
    boundary = boundaryL_adaptedFromPattern[patchR];
    connVert.clear();
    while(boundary[currIdx]!= boundary[idx(8)]){
        //add them
        polyLineInput.push_back(currPattern.row(boundary[currIdx]).transpose());
        cout<<boundary[currIdx]<<" third"<<endl;
        connVert.push_back(boundary[currIdx]);
        currIdx = (rightAsc)? (currIdx+1) % boundary.size() : currIdx-1;
        if(currIdx<0) currIdx+= boundary.size();

    }// add last
    polyLineInput.push_back(currPattern.row(boundary[currIdx]).transpose());
    connVert.push_back(boundary[currIdx]);
    connectedVertVec.push_back(connVert);
    cout<<boundary[currIdx]<<" third"<<endl;


    // go on on boundary of to pattern
    currIdx = idx(9);
    boundary.clear();
    boundary = boundaryL_toPattern[patchTo];
    while(boundary[currIdx]!= boundary[idx(11)]){
        //add them
        polyLineInput.push_back(Vg_to.row(boundary[currIdx]).transpose());
        cout<<boundary[currIdx]<<" four"<<endl;

        currIdx = (toAsc)? (currIdx+1) % boundary.size() : currIdx-1;
        if(currIdx<0) currIdx+= boundary.size();

    }// add last
    polyLineInput.push_back(Vg_to.row(boundary[currIdx]).transpose());
    cout<<boundary[currIdx]<<" four"<<endl;

}
void computeAllBetweensNew(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                           vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                           MatrixXd& currPattern, MatrixXd& Vg_to, vector<VectorXd>& polyLineInput, vector<vector<int>>& connectedVertVec, vector<int>& patchId, vector<bool>& isAscVec) {
    polyLineInput.clear();
    connectedVertVec.clear();
    cout<<endl<<"Seam Size "<<polylineSelected.size()<<endl<<endl;
    if(polylineSelected.size() == 12){
        computeAllBetweensConnectPatches(polylineSelected, polylineIndex, polyLineMeshIndicator,
                               boundaryL_adaptedFromPattern, boundaryL_toPattern, currPattern, Vg_to, polyLineInput, connectedVertVec, patchId, isAscVec);
        return;
    }
    cout<<"--------------------- -------------------------"<<endl;
    cout<<"--------------------- -------------------------"<<endl;
    cout<<"--------------------- -------------------------"<<endl;
    cout<<"--------------------- -------------------------"<<endl;
    cout<<" WE ARE IN THE CASE OF NOT HAVING not 12 nor 2. WHAT TO DO NOW IS UNCLEAR"<<endl;
    cout<<"--------------------- -------------------------"<<endl;
    cout<<"--------------------- -------------------------"<<endl;
    cout<<"--------------------- -------------------------"<<endl;
    cout<<"--------------------- -------------------------"<<endl;

    /* given 6 positions in total
     *   we assume v0 is on the from mesh ,adapted pattern
     *   v1 is a corner that should intersect the to pattern
     *   v2 is on the to pattern (but need not be a vertex, but shouuld be far away in another face)
     *      |           |
     *      v2          v3
     *      |           |
     * --v0--v1         v4 -- v5
     *      |           |
     *
     * */
    //step 1 locate v1 on the toPattern by checking on which edge it is
    int closer1, closer2, far1, far2; // closer is the one in direction of v2
    int patch;
    vector<vector<int>> boundaryToSearch = boundaryL_toPattern;
    for (int j = 0; j < boundaryToSearch.size(); j++) {
        bool found = false;
        cout<<"searching patch j="<<j<<endl;
        for (int k = 0; k < boundaryToSearch[j].size(); k++) {
            int v = boundaryToSearch[j][k];
            int v1 = boundaryToSearch[j][(k+1) % boundaryToSearch[j].size() ];
            if (vertOnEdge(Vg_to.row(v).transpose(), Vg_to.row(v1).transpose(), polylineSelected[1], v, v1)) {
                cout << "found vertex 1 on patch " << j<<" between vertices "<<v<<" and "<<v1 << endl;
                cout<<Vg_to.row(v)<<", " << Vg_to.row(v1)<<", "<< polylineSelected[1].transpose()<<endl;
                closer1 = k;
                far1 = (k+1) % boundaryToSearch[j].size();
                VectorXd v21 = (polylineSelected[2]-polylineSelected[1]).normalized();
                VectorXd vv1 = ( Vg_to.row(v)-polylineSelected[1]).normalized()- v21;
                VectorXd vv11 =  (Vg_to.row(v1)-polylineSelected[1]).normalized()- v21;
                if( (abs(vv11(0))+abs(vv11(1))) > (abs(vv1(0))+abs(vv1(1))) ){
//                if ((Vg_to.row(v)- polylineSelected[2]).norm() < ( Vg_to.row(v1)-polylineSelected[2]).norm() ) {
                    // swap
                    swap(closer1, far1);
                }
                patch = j;
                found = true;
                break;
            }
        }
        if (!found) continue;
        cout<<"***************** second searched"<<endl;
        // else we look for the other one
        for (int k = 0; k < boundaryToSearch[j].size(); k++) {
            int v = boundaryToSearch[j][k];
            int v1 = boundaryToSearch[j][k+1 % boundaryToSearch[j].size() ];
            if (vertOnEdge(Vg_to.row(v).transpose(), Vg_to.row(v1).transpose(), polylineSelected[4], v ,v1)) {
                cout << "found vertex  4 on patch " << j<<" between vertices "<<v<<" and "<<v1 << endl;
                cout<<Vg_to.row(v)<<" " << Vg_to.row(v1)<<" "<< polylineSelected[4]<<endl;
                closer2 = k;
                far2 = (k+1) % boundaryToSearch[j].size();
//                if ((Vg_to.row(v)- polylineSelected[3]).norm() < ( Vg_to.row(v1)-polylineSelected[3]).norm() ) {
                VectorXd v21 = (polylineSelected[4]-polylineSelected[3]).normalized();
                VectorXd vv1 = ( Vg_to.row(v)-polylineSelected[3]).normalized()- v21;
                VectorXd vv11 =  (Vg_to.row(v1)-polylineSelected[3]).normalized()- v21;
                if( (abs(vv11(0))+abs(vv11(1))) < (abs(vv1(0))+abs(vv1(1))) ){
                    // swap
                    swap(closer2, far2);
                }
                break;
            }
        }
    }
    // hopefully we found both by now!
    cout<<boundaryToSearch[patch][closer1]<<" "<<boundaryToSearch[patch][far1]<<" and other side "<<boundaryToSearch[patch][far2]<<" "<<boundaryToSearch[patch][closer2]<<endl;
    bool asc1;
    if(closer1==0){
        asc1 = closer1 + boundaryToSearch[patch].size() < far1;
    }else if(far1==0){
        asc1 = closer1 < boundaryToSearch[patch].size() +far1;
    }else{
        asc1 = closer1<far1;
    }

    if(asc1){// ascending order
        int i= far1;// for safety of tri better not?
//        i++;
//        i= i % boundaryToSearch[patch].size();
        int count = 0;
        while( i != far2 && i!=closer2 && count < boundaryToSearch[patch].size()){
            cout<<boundaryToSearch[patch][i]<<" 1, asc "<<Vg_to.row(boundaryToSearch[patch][i])<<endl;
            polyLineInput.push_back(Vg_to.row(boundaryToSearch[patch][i]).transpose());
            i++;
            i= i % boundaryToSearch[patch].size();
        }
        cout<<boundaryToSearch[patch][i]<<" "<<Vg_to.row(boundaryToSearch[patch][i])<<endl;
        polyLineInput.push_back(Vg_to.row(boundaryToSearch[patch][i]).transpose());


    }else{// descending
        int i= far1;
//        i--;
//        if(i<0) i+= boundaryToSearch[patch].size();
        int count = 0;

        while( i != far2 && i!=closer2 && count < boundaryToSearch[patch].size()){
            count++;
            cout<<boundaryToSearch[patch][i]<<" 1, desc "<<Vg_to.row(boundaryToSearch[patch][i])<<endl;
            polyLineInput.push_back(Vg_to.row(boundaryToSearch[patch][i]).transpose());
            i--;
            if(i<0) i+= boundaryToSearch[patch].size();
        }
        cout<<boundaryToSearch[patch][i]<<" "<<Vg_to.row(boundaryToSearch[patch][i])<<endl;

        cout<<"end "<<endl;
        polyLineInput.push_back(Vg_to.row(boundaryToSearch[patch][i]).transpose());
    }

    // search 0,1,4,5 on the other smaller pattern
    int patchFrom;
    int idx0, idx1, idx4, idx5;
    boundaryToSearch.clear();
    boundaryToSearch = boundaryL_adaptedFromPattern;
    for(int j=0; j< boundaryToSearch.size(); j++){
        for(int i=0; i< boundaryToSearch[j].size(); i++){
            if(boundaryToSearch[j][i] == polylineIndex[0]){
                cout<<i<<"found vertex 0 "<<polylineIndex[0]<<" on patch "<<j<<endl;
                idx0 = i;
                patchFrom = j;

            }else if(boundaryToSearch[j][i] == polylineIndex[1]){
                cout<<i<<"found vertex 1 "<<polylineIndex[1]<<" on patch "<<j<<endl;
                idx1 = i;
            }else if(boundaryToSearch[j][i] == polylineIndex[4]){
                cout<<i<<"found vertex 4 "<<polylineIndex[4]<<" on patch "<<j<<endl;
                idx4 = i;
            }else if(boundaryToSearch[j][i] == polylineIndex[5]){
                cout<<i<<"found vertex 5 "<<polylineIndex[5]<<" on patch "<<j<<endl;
                idx5 = i;
            }
        }
    }
    // hopefully now all indices are found
    bool asc = (idx5>idx4);
    if(idx4==0) asc = (idx5>(idx4+boundaryToSearch[patchFrom].size() ));
    cout<<"asc? "<<asc<<endl;

    connectedVertVec.clear();
    vector<int> connVec;
    patchId.clear(); patchId.push_back(patchFrom);
    isAscVec.clear(); isAscVec.push_back(asc);
    if(asc){// ascending
        int i = idx4;
        int count = 0;
        while(i !=idx1 && count < boundaryToSearch[patchFrom].size()){
            polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
            connVec.push_back(boundaryToSearch[patchFrom][i]);
            cout<<boundaryToSearch[patchFrom][i]<<" 2, asc"<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;
            i++; count++;
            i = i %boundaryToSearch[patchFrom].size();
        }
        polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
        cout<<boundaryToSearch[patchFrom][i]<<" "<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;
        connVec.push_back(boundaryToSearch[patchFrom][i]);


    }else{
        int i = idx4;
        int count = 0;
        while(i !=idx1 && count < boundaryToSearch[patchFrom].size()){
            polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
            connVec.push_back(boundaryToSearch[patchFrom][i]);
            cout<<boundaryToSearch[patchFrom][i]<<" 2, desc "<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;
            i--; count++;
            if(i <0) i+= boundaryToSearch[patchFrom].size();
        }
        polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
        cout<<boundaryToSearch[patchFrom][i]<<" "<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;
        connVec.push_back(boundaryToSearch[patchFrom][i]);


    }
    connectedVertVec.push_back(connVec);


}
void computeAllBetweens(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                   vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                   MatrixXd& currPattern, MatrixXd& Vg_pattern_orig, vector<VectorXd>& polyLineInput, vector<vector<int>>& connectedVertVec, vector<int>& patchId, vector<bool>& isAscVec){
    if(polylineSelected.size()>2){
        computeAllBetweensNew(polylineSelected, polylineIndex, polyLineMeshIndicator,
                              boundaryL_adaptedFromPattern, boundaryL_toPattern,
                              currPattern,  Vg_pattern_orig, polyLineInput, connectedVertVec, patchId, isAscVec);
        return;
    }


    polyLineInput.clear();
    connectedVertVec.clear();
    vector<int> connectedVert; int p; bool inverted;

    for(int i=0; i< polylineIndex.size(); i++){
        //add the current one
        polyLineInput.push_back(polylineSelected[i]);

        if(i+1 == polylineIndex.size()) continue; // there is no after
        if(polyLineMeshIndicator[i] ==2 || polyLineMeshIndicator[i+1]==2 ){
           // check if there are ot be inserted vertices on the other boundary line
            continue;
        }else{
            // add all vertices in between to the vector , but attention,exclude start and end
            int start = polylineIndex[i]; int startIdx, endIdx, patch;
            int end = polylineIndex[i+1];
            cout<<"both from same patch"<<endl;
            vector<vector<int>> boundaryToSearch = boundaryL_adaptedFromPattern;
            MatrixXd v_used = currPattern ;

            for(int j=0; j<boundaryToSearch.size(); j++){
                bool found = false;
                for (int k =0; k<boundaryToSearch[j].size(); k++){
                    if(boundaryToSearch[j][k] == start){
                        cout<<"found vertex "<<start<<" on patch "<<j<<endl;
                        startIdx = k;
                        patch = j;
                        found = true;
                    }
                }
                if(!found) continue;

                for (int k =0; k<boundaryToSearch[j].size(); k++){
                    if(boundaryToSearch[j][k] == end){
                        endIdx = k;
                        cout<<"found vertex "<<end<<endl;
                    }
                }
                // we have both start and end , their absolute distance should be
                int smaller = (endIdx > startIdx) ? startIdx : endIdx;
                inverted = false;
                if(smaller == endIdx) {
                    inverted = true;
                }
                int greater = (endIdx < startIdx) ? startIdx : endIdx;

                int dist = (greater - smaller);
                int otherdist = boundaryToSearch[j].size()-greater + smaller;
                cout<<otherdist<<" betweens "<<dist<<endl;
                p = j;

                if(dist<otherdist){
                    if(!inverted){
                        connectedVert.push_back(boundaryToSearch[j][smaller]);
                        for(int k = smaller+1; k < greater; k++){
                            polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
                            cout<<v_used.row(boundaryToSearch[j][k])<<endl;
                            connectedVert.push_back(boundaryToSearch[j][k]);

                        }
                        connectedVert.push_back(boundaryToSearch[j][greater]);

                    }else{
                        connectedVert.push_back(boundaryToSearch[j][greater]);
                        for(int k= greater-1; k> smaller ; k--){
                            polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
                            cout<<v_used.row(boundaryToSearch[j][k])<<endl;
                            connectedVert.push_back(boundaryToSearch[j][k]);

                        }
                        connectedVert.push_back(boundaryToSearch[j][smaller]);

                    }

                }else{
                    if(!inverted){
                        connectedVert.push_back(boundaryToSearch[j][smaller]);
                        int k= smaller; k--;
                        if(k<0) k += boundaryToSearch[j].size();
                        while( k != greater ){
                            polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
                            cout<<v_used.row(boundaryToSearch[j][k])<<endl;
                            connectedVert.push_back(boundaryToSearch[j][k]);

                            k--;
                            if(k<0) k += boundaryToSearch[j].size();

//                        k = k % boundaryToSearch[j].size();
                        }// and one more
//                        polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
//                        cout<<v_used.row(boundaryToSearch[j][k])<<"w-i"<<endl;
                        connectedVert.push_back(boundaryToSearch[j][k]);
                    }else{
                        connectedVert.push_back(boundaryToSearch[j][greater]);
                        int k = greater; k++;
                        k = k % boundaryToSearch[j].size();
                        while (k!= smaller){
                            polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
                            cout<<v_used.row(boundaryToSearch[j][k])<<endl;
                            connectedVert.push_back(boundaryToSearch[j][k]);
                            k++;
                            k = k % boundaryToSearch[j].size();
                        }// last one k==greater
//                        polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
//                        cout<<v_used.row(boundaryToSearch[j][k])<<"wi"<<endl;
                        connectedVert.push_back(boundaryToSearch[j][k]);
                    }
                }
            }
        }
    }
    connectedVertVec.push_back(connectedVert);
    isAscVec.push_back(!inverted);
    patchId.push_back(p);
    for(int i=0; i<isAscVec.size(); i++){
        cout<<isAscVec[i]<<" asc ? "<<endl ;
    }
    for(int i=0; i<patchId.size(); i++){
        cout<<patchId[i]<<" patch ID "<<endl ;
    }
    cout<<polyLineInput.size()<<" =?= "<<connectedVert.size()<<endl;
    for(int i=0; i<connectedVert.size(); i++){
        cout<< connectedVert[i] <<" is index i="<<i<<endl;
    }
}

int checkIfMatchesOne(VectorXd v,const vector<int>& connectedVert, const MatrixXd& currPattern ){
    for(int i=0; i<connectedVert.size(); i++){
        if(v == currPattern.row(connectedVert[i])){
            return connectedVert[i];
        }
    }
    return -1;
}
void replaceInFaces(int id, int newId, MatrixXi& Fg){
    for(int i = 0; i< Fg.rows(); i++){
        for(int j = 0; j < 3; j++){
            if(Fg(i, j) == id){
                Fg(i, j) = newId;
            }
        }
    }
}
void mergeTriagulatedAndPattern(const vector<vector<int>>& connectedVertVec, const vector<int>& patchId,const vector<bool>& isAscVec,
                           MatrixXd& Vg_retri, MatrixXi& Fg_retri, MatrixXd& currPattern, MatrixXi& Fg_pattern, vector<int> & newFaces){
    currPattern.col(2).setConstant(200);
    vector<vector<int>> boundaryLinsert, boundaryLGar;
    igl::boundary_loop(Fg_retri, boundaryLinsert);
    igl::boundary_loop(Fg_pattern, boundaryLGar);
    VectorXi mergedIds(Vg_retri.rows());
    mergedIds.setConstant(-1);
    igl::writeOBJ("retriPatch.obj", Vg_retri, Fg_retri);
    // insertion on each pattern individually
    for(int i=0; i<isAscVec.size(); i++){
        int vertId = connectedVertVec[i][0];
        vector<int> boundary = boundaryLGar[patchId[i]];
        vector<int> boundaryIns = boundaryLinsert[0];
        int patchSize = boundaryLGar[patchId[i]].size();
        int patchSizeIns = boundaryLinsert[0].size();

        int garIdx =0;
        while(boundary[garIdx] != vertId && garIdx < patchSize){
            garIdx++;
        }
        if(boundary[garIdx] != vertId){cout<<"Something went wrong, index not found. "<< endl; return; }
        int insertIdx=0;
        while(Vg_retri.row(boundaryIns[insertIdx]) != currPattern.row(vertId) && insertIdx < patchSizeIns){
            insertIdx++;
        }
        if(Vg_retri.row(boundaryIns[insertIdx]) != currPattern.row(vertId)){
            cout<<"Something went wrong. Not found in insert patch. "<<endl; return;
        }

        mergedIds(boundaryIns[insertIdx])= boundary[garIdx];

        int currVertIdGar = boundary[garIdx];
        int currVertIdInsert = boundaryIns[insertIdx];
        int jcounter = 1; int j=1;
        for(int jj=1; jj < connectedVertVec[i].size(); jj++){//
            int nextVertId, nextVertIdInserted;
            if(isAscVec[i]){
                nextVertId = boundary[ (garIdx + jcounter) % patchSize];
                nextVertIdInserted = (insertIdx - j) ;
                if(nextVertIdInserted < 0) nextVertIdInserted += patchSizeIns;
                nextVertIdInserted = boundaryIns[nextVertIdInserted ];

            }else{
                nextVertId = (garIdx - jcounter) ;
                if(nextVertId < 0) nextVertId += patchSize;
                nextVertId = boundary[nextVertId];
                nextVertIdInserted = boundaryIns[(insertIdx + j) % patchSizeIns];

            }
            if((Vg_retri.row(nextVertIdInserted) - currPattern.row(nextVertId )).norm() >0.0001){
                cout<<" We have to insert a vertex "<<endl ;
                vector<vector<int>> vfAdj;
                createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
                int faceToDupl = adjacentFaceToEdge(nextVertId,currVertIdGar, -1, vfAdj);
                MatrixXi Fg_new (Fg_pattern.rows()+1, 3);
                MatrixXd currPattern_new (currPattern.rows()+1, 3);
                Fg_new.block(0,0,Fg_pattern.rows(), 3 ) = Fg_pattern;
                Fg_new.row(Fg_pattern.rows()) = Fg_pattern.row(faceToDupl);
                currPattern_new.block(0,0,currPattern.rows(), 3 ) = currPattern;
                currPattern_new.row(currPattern.rows()) = Vg_retri.row(nextVertIdInserted);

                if(Fg_new(Fg_pattern.rows(), 0) == nextVertId){
                    Fg_new(Fg_pattern.rows(), 0) = currPattern.rows();
                }else if(Fg_new(Fg_pattern.rows(), 1) == nextVertId ){
                    Fg_new(Fg_pattern.rows(), 1) = currPattern.rows();
                }else{
                    Fg_new(Fg_pattern.rows(), 2) = currPattern.rows();
                }
                if(Fg_new(faceToDupl, 0) == currVertIdGar ){
                    Fg_new(faceToDupl, 0) = currPattern.rows();
                }else if(Fg_new(Fg_pattern.rows(), 1) == currVertIdGar ){
                    Fg_new(faceToDupl, 1) = currPattern.rows();
                }else{
                    Fg_new(faceToDupl, 2) = currPattern.rows();
                }

                mergedIds(nextVertIdInserted) = currPattern.rows();
                Fg_pattern.resize(Fg_new.rows(), 3); Fg_pattern = Fg_new;
                currPattern.resize(currPattern_new.rows(), 3); currPattern = currPattern_new;
                jj--;
                currVertIdGar = currPattern.rows()-1;//  nextVertId;
            }else{
                mergedIds(nextVertIdInserted)= nextVertId;
                jcounter++;
                currVertIdGar = nextVertId;
            }
            j++;

            currVertIdInsert = nextVertIdInserted;

        }

    }
    cout<<"starting the merger"<<endl;
    //merging part
    int count = 0;
    VectorXi newId(Vg_retri.rows());
    MatrixXd Vg_help (Vg_retri.rows(), 3);
    for (int i=0 ; i<Vg_retri.rows(); i++){
        if(mergedIds(i) == -1){
            newId(i) = count;
            Vg_help.row(count) = Vg_retri.row(i);
            count++;
        }
    }
    int offset = currPattern.rows();
    for (int i = 0; i < Fg_retri.rows(); i++){
        for(int j=0; j<3; j++){
            if(mergedIds(Fg_retri(i, j)) != -1){
                Fg_retri(i, j) = mergedIds(Fg_retri(i, j));
            }else {
                Fg_retri(i, j) = newId(Fg_retri(i,j)) + offset;
            }
        }
    }

    int faceOffset = Fg_pattern.rows();
    MatrixXi Fg_new (faceOffset+ Fg_retri.rows(), 3);
    MatrixXd currPattern_new (offset + count, 3);
    Fg_new.block(0,0,Fg_pattern.rows(), 3 ) = Fg_pattern;
    Fg_new.block(faceOffset, 0, Fg_retri.rows(), 3) = Fg_retri;
    currPattern_new.block(0,0,currPattern.rows(), 3 ) = currPattern;
    currPattern_new.block(offset, 0,count, 3 ) = Vg_help.block(0,0, count, 3);

    Fg_pattern.resize(Fg_new.rows(), 3); Fg_pattern = Fg_new;
    currPattern.resize(currPattern_new.rows(), 3); currPattern = currPattern_new;
    for(int i= faceOffset; i<Fg_pattern.rows(); i++){
        newFaces.push_back(i);
    }
}

vector<int> toVecInt(VectorXi& v){
    vector<int> vec;
    for(int i=0; i<v.rows(); i++){
        vec.push_back(v(i));
    }
    return vec;
}

vector<double> toVecDouble(VectorXd& v){
    vector<double> vec;
    for(int i=0; i<v.rows(); i++){
        vec.push_back(v(i));
    }
    return vec;
}
void createHalfAvatarMap(MatrixXd& testMorph_V1, MatrixXi& testMorph_F1,
    MatrixXd& testMorph_V1left, MatrixXi& testMorph_F1left,
    MatrixXd& testMorph_V1right, MatrixXi& testMorph_F1right,
    map<int, int>& leftHalfToFullFaceMap,  map<int, int>& rightHalfToFullFaceMap){

    map<vector<double>, int> posToVertIdFull;
    map<vector<int>, int> vertToFaceIdFull;

//    vector<int> indicesOfLeft;
//    vector<int> indicesOfRight;
//    set<vector<double>> setOfLeftPos, setOfRightPos;
//
//    for(int i=0; i< testMorph_V1left.rows(); i++){
//        VectorXd v = testMorph_V1left.row(i).transpose();
//        vector<double> vert = toVecDouble(v);
//        setOfLeftPos.insert(vert);
//    }
//    for(int i=0; i< testMorph_V1right.rows(); i++){
//        VectorXd v = testMorph_V1right.row(i).transpose();
//        vector<double> vert = toVecDouble(v);
//        setOfRightPos.insert(vert);
//    }

//    string fileName = "./rightHalfFaces.txt";
//    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
//    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
//
//    ofstream file(fileName);
//    if (file.is_open())
//    {
//        file << testMorph_F1right.format(CSVFormat);
//        file.close();
//    }

    for(int i=0; i< testMorph_V1.rows(); i++){
        VectorXd v = testMorph_V1.row(i).transpose();
        vector<double> vert = toVecDouble(v);
        posToVertIdFull[vert] = i;
//        if(setOfRightPos.contains(vert)){
//        if(auto search = setOfRightPos.find(vert); search != setOfRightPos.end()){
//            indicesOfRight.push_back(i);
//        }else if (auto search = setOfLeftPos.find(vert); search != setOfLeftPos.end()){
//            indicesOfLeft.push_back(i);
//        }else{
//            cout<<"NONE CONTAINS VERTEX "<<i<<endl;
//        }
    }
//    std::ofstream output_file("./rightHalfIndices.txt");
//    std::ostream_iterator<int> output_iterator(output_file, "\n");
//    std::copy(std::begin(indicesOfRight), std::end(indicesOfRight), output_iterator);

//    std::ostream_iterator<std::string> output_iterator2(output_file2, "\n");
//    std::copy(std::begin(indicesOfRight), std::end(indicesOfRight), output_iterator2);

    for(int i=0; i< testMorph_F1.rows(); i++){
        VectorXi f = testMorph_F1.row(i).transpose();
        vector<int> face =  toVecInt(f);
        vertToFaceIdFull[face]= i;
    }

    for(int i=0; i< testMorph_F1left.rows(); i++){
        Vector3i idLeft = testMorph_F1left.row(i);
        vector<int> idTotal;
        VectorXd v = testMorph_V1left.row(idLeft(0)).transpose();
        idTotal.push_back( posToVertIdFull[toVecDouble(v)]);
        v = testMorph_V1left.row(idLeft(1)).transpose();
        idTotal.push_back( posToVertIdFull[toVecDouble(v)]);
        v = testMorph_V1left.row(idLeft(2)).transpose();
        idTotal.push_back( posToVertIdFull[toVecDouble(v)]);

        leftHalfToFullFaceMap[i] = vertToFaceIdFull[idTotal] ;
    }

    for(int i=0; i< testMorph_F1right.rows(); i++){
        Vector3i idRight = testMorph_F1right.row(i);
        vector<int> idTotal;
        VectorXd v = testMorph_V1right.row(idRight(0)).transpose();
        idTotal.push_back( posToVertIdFull[toVecDouble(v)]);

        v = testMorph_V1right.row(idRight(1)).transpose();
        idTotal.push_back( posToVertIdFull[toVecDouble(v)]);

        v = testMorph_V1right.row(idRight(2)).transpose();
        idTotal.push_back( posToVertIdFull[toVecDouble(v)]);

        rightHalfToFullFaceMap[i] = vertToFaceIdFull[idTotal] ;

    }

}

void initialGuessAdaption(MatrixXd& currPattern_nt, MatrixXd& mapToVg_nt, MatrixXd& perfectPattern_nt,  MatrixXi& Fg_pattern_curr, MatrixXi& mapToFg, bool symetry,
                          set<int> & cornerSet, map<int, int >& mapCornerToCorner, int origHalfSize, map<int, int>& halfPatternVertToFullPatternVertT){
//perfect pattern is the initial perfect pattern for the new shape, and curr is the existing adapted to the perfect, thus is has more faces and its own Fg
// note that mapToVg was the one we started with, and it has face correspondence with perfect pattern, therefore they both use mapToFg
    cout<<"start init guess"<<endl ;
    Eigen::VectorXi componentIdPerVert_curr, componentIdPerVert_other;
    igl::vertex_components(mapToFg, componentIdPerVert_other);
    igl::vertex_components(Fg_pattern_curr,componentIdPerVert_curr );
    MatrixXd currPattern = currPattern_nt;

    vector<int> right, left, rightFull, leftFull;
    rightFull.push_back(1); rightFull.push_back(5);
    leftFull.push_back(3); leftFull.push_back(6);
    if(!symetry){
        right = rightFull;
        left = leftFull;
    }else{
        right.push_back(1); right.push_back(3); // maybe switched!!
        left.push_back(100); left.push_back(3000);
    }

    for(int i=0; i < currPattern_nt.rows(); i++){
        if(componentIdPerVert_curr(i)== right[0] ||componentIdPerVert_curr(i)== right[1]){
            currPattern(i, 0) += 100;
        }
        else if(componentIdPerVert_curr(i)== left[0] ||componentIdPerVert_curr(i)== left[1] ){
            currPattern(i, 0) -= 100;
        }
    }

    MatrixXd mapToVg = mapToVg_nt;
    MatrixXd perfectPattern = perfectPattern_nt;
    for(int i=0; i < mapToVg.rows(); i++){
        if(componentIdPerVert_other(i)== rightFull[0] ||componentIdPerVert_other(i)== rightFull[1]){
            mapToVg(i, 0) += 100;
            perfectPattern(i, 0) += 100;
        }
        else if(componentIdPerVert_other(i)== leftFull[0] ||componentIdPerVert_other(i)== leftFull[1] ){
            mapToVg(i, 0) -= 100;
            perfectPattern(i, 0) -= 100;

        }
    }

    VectorXd S;
    VectorXi I;//face index of smallest distance
    MatrixXd C,N, initGuess;
    cout<<"Starting with signed distance "<<endl;

    igl::signed_distance(currPattern, perfectPattern, mapToFg, igl::SIGNED_DISTANCE_TYPE_UNSIGNED, S, I, C, N);
//    cout<<"Finished signed distance "<<endl;
    vector<set<pair<int, double>>> newCornerCandidate(perfectPattern.rows());
    MatrixXd B(currPattern.rows(), 3); // contains all barycentric coordinates
    for(int i = 0; i < currPattern.rows(); i++){
        VectorXd bary;

        auto face = mapToFg.row(I(i));
        igl::barycentric_coordinates(currPattern.row(i), perfectPattern.row(face(0)), perfectPattern.row(face(1)),
                                     perfectPattern.row(face(2)), bary);
        for(int c=0; c<3; c++){
            if(cornerSet.find(mapToFg(I(i),c)) != cornerSet.end()){
                newCornerCandidate[mapToFg(I(i),c) ].insert(make_pair(i, (currPattern.row(i) - perfectPattern.row( mapToFg(I(i),c) )).norm() ));
            }
        }

        B.row(i) = bary;
    }
    cout<<"Got all barycentric coords "<< newCornerCandidate.size()<<endl;
    igl::barycentric_interpolation(mapToVg, mapToFg, B, I, initGuess);
    cout<<"Got interpolation "<<initGuess.rows()<<" and before we had "<<currPattern.rows()<<endl;
    currPattern  = initGuess;
    //undo the transposition
    currPattern_nt = currPattern;
    for(int i=0; i < currPattern_nt.rows(); i++){
        if(componentIdPerVert_curr(i)== right[0] ||componentIdPerVert_curr(i)== right[1]){
            currPattern_nt(i, 0) -= 100;
        }
        else if(componentIdPerVert_curr(i)== left[0] ||componentIdPerVert_curr(i)== left[1] ){
            currPattern_nt(i, 0) += 100;
        }
    }
    for(int i = 0; i<newCornerCandidate.size(); i++){
        double min = 10000;
        int newCorner = -1;
        if(newCornerCandidate[i].size()<1)continue;
        for(auto sugg : newCornerCandidate[i]){
            if(sugg.second<min){
               newCorner = sugg.first;
               min = sugg.second;
            }
        }
        if(newCorner ==-1){
            cout<<"no new corner found "<<i<<endl;
        }else{
            cout<<i<<" new corner is "<<newCorner<<endl;
            if(newCorner>origHalfSize){
                cout<<"added neg"<<endl;
                newCorner *= -1;
            }else{
                newCorner = halfPatternVertToFullPatternVertT[newCorner];
                cout<<"updated to "<<newCorner<<endl;
            }
            mapCornerToCorner[i]= newCorner;
        }
    }

cout<<"fin init guess "<<currPattern_nt.rows()<<endl;
}
void initialGuessAdaptionWithoutT(MatrixXd& currPattern, MatrixXd& mapToVg, MatrixXd& perfectPattern,  MatrixXi& Fg_pattern_curr, MatrixXi& mapToFg){
//perfect pattern is the initial perfect pattern for the new shape, and curr is the existing adapted to the perfect, thus is has more faces and its own Fg
// note that mapToVg was the one we started with, and it has face correspondance with perfect pattern, therefore they both use mapToFg
    VectorXd S;
    VectorXi I;//face index of smallest distance
    MatrixXd C,N, initGuess;
    cout<<"Starting with signed distance "<<endl;
    igl::signed_distance(currPattern, perfectPattern, mapToFg, igl::SIGNED_DISTANCE_TYPE_UNSIGNED, S, I, C, N);

    MatrixXd B(currPattern.rows(), 3); // contains all barycentric coordinates
    for(int i = 0; i < currPattern.rows(); i++){
        VectorXd bary;
        auto face = mapToFg.row(I(i));
        igl::barycentric_coordinates(currPattern.row(i), perfectPattern.row(face(0)), perfectPattern.row(face(1)),
                                     perfectPattern.row(face(2)), bary);
        B.row(i) = bary;
    }
    cout<<"Got all barycentric coords"<<endl;
    igl::barycentric_interpolation(mapToVg, mapToFg, B, I, initGuess);
    cout<<"Got interpolation "<<initGuess.rows()<<" and before we had "<<currPattern.rows()<<endl;
    currPattern  = initGuess;

}

set<int> boundaryVerticesP;
MatrixXd internalAngles;
void ensureAngle(MatrixXd& p, MatrixXd& fromPattern, MatrixXi& Fg_pattern, MatrixXi& fromPatternFg){
    if(boundaryVerticesP.empty()){
        vector<vector<int>> boundaryVP;
        igl::boundary_loop(Fg_pattern, boundaryVP);
        for(int i=0; i<boundaryVP.size(); i++){
            for(int j=0; j<boundaryVP[i].size(); j++){
                boundaryVerticesP.insert(boundaryVP[i][j]);
            }
        }
    }
    if(internalAngles.rows() != fromPatternFg.rows()){
        igl::internal_angles(fromPattern, fromPatternFg, internalAngles);

    }
    MatrixXd internalAngleCurr;
    igl::internal_angles(p, Fg_pattern, internalAngleCurr);
    for(int i = 0; i < Fg_pattern.rows(); i++){
        for(int j=0; j<3; j++){

            double newAngle = internalAngleCurr(i,j);
            double newdegree = newAngle*180/M_PI;

            if(newdegree<10){//only an issue if the original angle was better

                double oldAngle = internalAngles(i,j);
                double olddegree = oldAngle*180/M_PI;
                if(newdegree/olddegree > 0.75) continue;

                MatrixXd rot = MatrixXd::Identity(3, 3);
                double damp = (15-newdegree)/15;
                double rad =  ((olddegree-newdegree) / 2) / 180 *M_PI; //converting to radian value
                rot(0,0)= cos(rad);
                rot(1, 1) = rot(0,0) ;
                rot(0,1) =  - sin(rad);
                rot(1,0) = sin(rad);
                double stiffness = 0.5;

                Vector3d e2 = p.row(Fg_pattern(i, (j + 1) % 3))-p.row(Fg_pattern(i, (j) % 3));
                Vector3d e1 = p.row(Fg_pattern(i, (j + 2) % 3 ))-p.row(Fg_pattern(i, (j) % 3));
                VectorXd e1rot= rot * e1;
                VectorXd dir = e1rot - e1;
                p.row(Fg_pattern(i, j)) += damp*stiffness * dir;

                VectorXd e2rot= rot.transpose() * e2;
                dir = e2rot - e2;
                p.row(Fg_pattern(i, (j+2) % 3)) += damp* stiffness * dir;

//                e1 = p.row(Fg_pattern(i, j))-p.row(Fg_pattern(i, (j+1) % 3));
//                e2 = p.row(Fg_pattern(i, (j+2) % 3 ))-p.row(Fg_pattern(i, (j+1) % 3));
//                newAngle = acos(min(max((e1.normalized()).dot(e2.normalized()), -1.), 1.));
//                newdegree = newAngle*180/M_PI;
////                cout<<"new degree: "<<newdegree<<endl;

            }


        }


    }
}

void ensurePairwiseDist(MatrixXd& p, MatrixXd& toPattern, MatrixXi& Fg_pattern){
    vector<vector<int>> faceFaceAdjecencyList;

    createFaceFaceAdjacencyList(Fg_pattern,faceFaceAdjecencyList);
    for(int i = 0; i < Fg_pattern.rows(); i++){
        vector<int> neigh = faceFaceAdjecencyList[i];
        VectorXi face = Fg_pattern.row(i);
        Vector3d e1 = toPattern.row(Fg_pattern(i, 0))-toPattern.row(Fg_pattern(i, 1));
        Vector3d e2 = toPattern.row(Fg_pattern(i, 2))-toPattern.row(Fg_pattern(i, 1));
        auto crossp = e1.cross(e2);

        Vector3d e1old = p.row(Fg_pattern(i, 0))-p.row(Fg_pattern(i, 1));
        Vector3d e2old = p.row(Fg_pattern(i, 2))-p.row(Fg_pattern(i, 1));
        auto crosspold = e1old.cross(e2old);

        if(crossp.dot(crosspold)<0){
//            cout<<"face "<<i<<" flipped"<<endl;
        }

//        for(int j = 0; j < neigh.size(); j++){
//            VectorXi other = Fg_pattern.row(neigh[j]);
//            // find the two different vertices
//            int first=0; int second =0;
//            while(other(0)== face(first) || other(1)== face(first) || other(2) == face(first) ){
//                first++; // first is not in the other
//            }
//            while(face(0)== other(second) || face(1)== other(second) || face(2) == other(second) ){
//                second++; // second is not in face
//            }
//            VectorXd distOrig = (toPattern.row(face(first)) - toPattern.row(other(second)));
//            VectorXd distNew = (p.row(face(first)) - p.row(other(second)));
//            if(i== 2957 && neigh[j]== 2904){
////                cout<<face(first)<<" vertices " <<other(second)<<endl;
////                cout<< (p.row(face(first)))<<", "<<endl<<( p.row(other(second)))<<", "<<endl<<(p.row(face(first)) - p.row(other(second)))<<endl;
//                cout<<"face "<<neigh[j]<<" :"<<distNew.norm()/ distOrig.norm()<<", "<< distNew.norm()<< ", "<<distOrig.norm()<<endl;
//            }
//            if(distOrig.dot(distNew) < 0){//|| abs(distNew.norm()/ distOrig.norm()-1 ) < 0.5
//                cout<<"face "<<i<<" generates a flip with face "<<neigh[j]<<endl;
//            }
//        }
    }
}
//map<int, int> patchMapToHalf;
map<int, int> htFFace ,pMapToHalf,
//        fullPatternFaceToHalfPatternFace,
        fTHVert;
//        halfPatternVertToFullPatternVert;

void createMapCornersToNewCorner(MatrixXd& currPattern,MatrixXd& mapToVg, vector<vector<pair<int, int>>>& cornerPerBoundary,// first is vert id, second ins loop id, but thats bullshit
                                 map<int, int>& mapCornerToCorner, vector<vector<int>>& boundaryL, map<int, int>& halfPatternVertToFullPatternVertT,
                                 map<int, int>& fullPatternVertToHalfPatternVertT, bool symetry, map<int, int>&  halfPatternFaceToFullPatternFaceT, map<int, int>&
                                 fullPatternFaceToHalfPatternFaceT ){
    currPattern.col(2).setConstant(200);
    mapToVg.col(2).setConstant(200);
    fTHVert = fullPatternVertToHalfPatternVertT;
    double eps = 0.1; int count = 0;
//    for(int i = 0; i<cornerPerBoundary.size(); i++){
//        vector<pair<int, int>> cornersOfB = cornerPerBoundary[i];
//        for(int j = 0; j< cornersOfB.size(); j++){
//            int newId;
//
//            if(fullPatternVertToHalfPatternVertT.find(cornersOfB[j].first ) == fullPatternVertToHalfPatternVertT.end()){continue; }
//
//            if(pMapToHalf.find(i)== pMapToHalf.end()){
//                pMapToHalf[i]=count;
//                cout<<i<<" is now patch "<<count<<endl;
//                count++;
//            }
//            cout<<j<<" / "<<endl;
//            int halfId= fullPatternVertToHalfPatternVertT[cornersOfB[j].first];
//            if((currPattern.row(halfId )- mapToVg.row(halfId)).norm()> eps){
//                // we need to find it
//                bool found = false;
//                double minDist = (currPattern.row(boundaryL[pMapToHalf[i]][0])- mapToVg.row(halfId)).norm();
//                for(int ii=0; ii< boundaryL[pMapToHalf[i]].size(); ii++){
//
//                    if((currPattern.row(boundaryL[pMapToHalf[i]][ii])- mapToVg.row(halfId)).norm()< minDist){
//                        found = true;
//                        newId= boundaryL[i][ii];
//                        minDist = (currPattern.row(boundaryL[pMapToHalf[i]][ii])- mapToVg.row(halfId)).norm();
//                    }
//                }
//                if(!found){
//                    cout<<" no alternative corner found!"<<endl ;
//                }else{cout<<"taking "<<newId<<" with dist "<<minDist<<" as heuristic"<<endl;
//                    if(halfPatternVertToFullPatternVertT.find(newId) == halfPatternVertToFullPatternVertT.end() ){
//                        newId *= (-1);
//                        fullPatternVertToHalfPatternVertT[newId] = (-1) * newId;
//                        cout<<newId<<" a new introduced vertex iis corner, map it negatively"<<endl;
//                    }else{
//                        cout<<"Full Id compute"<<endl;
//                        newId = halfPatternVertToFullPatternVertT[newId];
//                        cout<<"Full Id computed * "<<endl;
//
//                    }
//                }
//            }else{
//                newId = cornersOfB[j].first;
//            }
//
//            mapCornerToCorner[cornersOfB[j].first] = newId;
//        }
//    }
//    if(symetry){
//        halfPatternFaceToFullPatternFace= halfPatternFaceToFullPatternFaceT;
//        fullPatternFaceToHalfPatternFace =fullPatternFaceToHalfPatternFaceT;
//        fullPatternVertToHalfPatternVert= fullPatternVertToHalfPatternVertT;
//        halfPatternVertToFullPatternVert=halfPatternVertToFullPatternVertT;
//
//        for(int i=0; i<cornerPerBoundary.size(); i++){
//            for(int j=0; j<cornerPerBoundary[i].size(); j++){
//                int ver = cornerPerBoundary[i][j].first;
//                if(fullPatternVertToHalfPatternVert.find(ver) == fullPatternVertToHalfPatternVert.end()){
//                    continue;
//                }
//                int newVer = fullPatternVertToHalfPatternVert[ver];
//                for(int ii=0; ii < boundaryL.size(); ii++){
//                    for(int jj= 0; jj<boundaryL[ii].size(); jj++){
//                        if(boundaryL[ii][jj] == newVer){
//                            patchMapToHalf[i]= ii;
//                        }
//                    }
//                }
//            }
//        }
//    }
}

void updateCornerUtils(set<int>& cornerSet , // a set containing all corner vertices
                       vector<vector<pair<int, int>>>& cornerPerBoundary,
                       map<int, vector<pair<int, int>>>& seamIdPerCorner,    // contains corner id and a list of which seams start here (max 2),
                        map<int, int>& mapCornerToCorner,  VectorXd&   cornerVertices // 1 for each corner, 0 for all other vertices
){
todo this is wrong. update with old and new for corner mapping

    cornerSet.clear();

    map<int, vector<pair<int, int>>> newSeamIdPerCorner;

    for(int i=0; i<cornerPerBoundary.size(); i++){
        for(int j=0; j<cornerPerBoundary[i].size(); j++){
            newSeamIdPerCorner[mapCornerToCorner[ cornerPerBoundary[i][j].first]] = seamIdPerCorner[cornerPerBoundary[i][j].first];
            cornerPerBoundary[i][j].first = mapCornerToCorner[ cornerPerBoundary[i][j].first];
            cornerSet.insert(cornerPerBoundary[i][j].first);
            cornerVertices(cornerPerBoundary[i][j].first) = 1;
        }
    }
    seamIdPerCorner.clear();
    seamIdPerCorner = newSeamIdPerCorner;
}

void updateSeamCorner( vector<seam*>& seamsList,  vector<minusOneSeam*> & minusOneSeams, map<int, int>& mapCornerToCorner,
                       vector<vector<int>>& boundaryL){
    for(int i=0; i<seamsList.size(); i++){
        int start1 =  seamsList[i]->getStart1();
        int start2 =  seamsList[i]->getStart2();


        //issue: some verts are in half-pattern
        // half-pattern gets some new verts
        // but these new verts of half pattern are not in halfMap
        // which is ok bc they are not in full, but maybe it is needed in some mappings?


        auto ends =  seamsList[i]->getEndCornerIds();
        seamsList[i]->updateStartEnd( mapCornerToCorner[start1], mapCornerToCorner[start2],
                                      mapCornerToCorner[ends.first],  mapCornerToCorner[ends.second]) ;
        cout<<"seam "<<i<<" "<<mapCornerToCorner[start1]<<" "<< mapCornerToCorner[start2]<<" "<<
                mapCornerToCorner[ends.first]<<" "<<  mapCornerToCorner[ends.second]<<endl;
        int patch1 = seamsList[i]->getPatch1();
        int patch2 = seamsList[i]->getPatch2();

//        int start1idx=0;
//        while(boundaryL[patchMapToHalf[patch1]][start1idx] != mapCornerToCorner[start1] && start1idx <= boundaryL[patch1].size()){
//            start1idx ++;
//        } if (boundaryL[patch1][start1idx] != mapCornerToCorner[start1]) { cout<<"start1 not found"<<endl;}
//
//        int start2idx=0;
//        while(boundaryL[patch2][start2idx] != mapCornerToCorner[start2] && start2idx <= boundaryL[patch2].size() ){
//            start2idx ++;
//        } if(boundaryL[patch2][start2idx] != mapCornerToCorner[start2]){ cout<<"start2 not found"<<endl;}
//
//        int end1idx=0;
//        while(boundaryL[patch1][end1idx] != mapCornerToCorner[ends.first] && end1idx <= boundaryL[patch1].size()){
//            end1idx ++;
//        } if(boundaryL[patch1][end1idx] != mapCornerToCorner[ends.first]){cout<<"end1 not found"<<endl; }
//
//        int end2idx=0;
//        while(boundaryL[patch2][end2idx] != mapCornerToCorner[ends.second] && end2idx <= boundaryL[patch2].size()){
//            end2idx ++;
//        } if(boundaryL[patch2][end2idx] != mapCornerToCorner[ends.second]){ cout<<"end2 not found"<<endl;}


//        seamsList[i]->updateStartEndIdx( start1idx, start2idx, end1idx, end2idx);

    }
    for(int i=0; i<minusOneSeams.size(); i++){
        int start =  minusOneSeams[i]->getStartVert();
        int end =  minusOneSeams[i]->getEndVert();
        if(fTHVert.find(start)== fTHVert.end())continue; // no need to update, it is not in the half pattern

        minusOneSeams[i]->updateStartEnd( mapCornerToCorner[start], mapCornerToCorner[end]) ;
        cout<<"From "<<start<<" "<<end<<" to "<<mapCornerToCorner[start]<<" "<<mapCornerToCorner[end]<<endl;
        int patch = minusOneSeams[i]->getPatch();
//        int startidx=0;
//        while(boundaryL[patch][startidx] != mapCornerToCorner[start] && startidx <= boundaryL[patch].size()){
//            startidx ++;
//        }if (boundaryL[patch][startidx] != mapCornerToCorner[start] ){ cout<<"start not found" <<endl; }
//        int endidx=0;
//        while(boundaryL[patch][endidx] != mapCornerToCorner[end] && endidx <= boundaryL[patch].size()){
//            endidx ++;
//        } if (boundaryL[patch][endidx] != mapCornerToCorner[end]){ cout<<"end not fond "<<endl ; }
//        minusOneSeams[i]->updateStartEndIdx( startidx, endidx) ;

    }
    cout<<"updated all corners :) "<<endl;
}

void stitchSeam(vector<int>& startAndEnd, MatrixXd& currPattern, MatrixXi& Fg_pattern_curr){
    Eigen::VectorXi componentIdPerVert;
    igl::vertex_components(Fg_pattern_curr, componentIdPerVert);
    vector<vector<int>> boundaryLCurr ;
    igl::boundary_loop(Fg_pattern_curr, boundaryLCurr);
    if(startAndEnd.size() != 6) {
        cout<<" Sorry you should select 6 vertices, 4 corners and 4 intermediate.  "<<endl;
    }
    VectorXi idxOf(6);
    int patch1 = componentIdPerVert(startAndEnd[0]);
    int patch2 = componentIdPerVert(startAndEnd[3]);
    // search for the indices of start and end elements

    for(int l = 0; l < boundaryLCurr.size(); l++){
        for (int m = 0; m < boundaryLCurr[l].size(); m++){
            for(int i = 0; i < startAndEnd.size(); i++) {
                if (startAndEnd[i] == boundaryLCurr[l][m] ){
                    idxOf(i) = m;

                }
            }
        }
    }

    bool patch2Asc = isAsc(idxOf(3), idxOf(4), idxOf(5));
    bool patch1Asc = isAsc(idxOf(0), idxOf(1), idxOf(2));

    // we need to translate the whole component
    VectorXd offset = currPattern.row(startAndEnd[0])-currPattern.row(startAndEnd[3]);
    int componentOfOne = componentIdPerVert(startAndEnd[3]);
    for(int i =0; i<currPattern.rows(); i++){
        if(componentIdPerVert(i) == componentOfOne){
            currPattern.row(i)+= offset.transpose();
        }
    }
    cout<<"For now we assume the number of vertices on both sides is the same. Maybe this has to be adapted later"<<endl;

    vector<vector<int>> vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern_curr, vfAdj);

    int i=0;
    int currIdx1 = idxOf(0);
    int currIdx2 = idxOf(3);

    while (boundaryLCurr[patch1][currIdx1] != boundaryLCurr[patch1][idxOf(2)]){
        // do I need to remove it fully or unreference it or just move the posititions
        currPattern.row(boundaryLCurr[patch2][currIdx2]) = currPattern.row(boundaryLCurr[patch1][currIdx1]);
        i++;
        if(patch1Asc){
            currIdx1 = (idxOf(0)+ i) % boundaryLCurr[patch1].size();
        }else{
            currIdx1 = ((idxOf(0)- i)>= 0) ?(idxOf(0)- i) : (idxOf(0)- i)+boundaryLCurr[patch1].size() ;
        }

        if(patch2Asc){
            currIdx2 = (idxOf(3)+ i) % boundaryLCurr[patch2].size();
        }else{
            currIdx2 = ((idxOf(3)- i)>= 0) ?(idxOf(3)- i) : (idxOf(3)- i)+boundaryLCurr[patch2].size() ;
        }
    }
    currPattern.row(boundaryLCurr[patch2][currIdx2]) = currPattern.row(boundaryLCurr[patch1][currIdx1]);

}
