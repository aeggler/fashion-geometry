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



using namespace std;
using namespace Eigen;

//Taubin smoothing according to
//https://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/06_smoothing.pdf
// from start to end in direction of the middle vertex

void smoothBetweenVertices(MatrixXd& currPattern, MatrixXi& Fg_pattern, vector<int>& startAndEnd){

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

    for (int j = startIdx; j < loopSize + startIdx; ++j) {
        if( boundary[j % loopSize] == startAndEnd[2] ) {
            otherDir= true;
            break;
        }else if(boundary[j % loopSize] == startAndEnd[1] ){
            midIdx = j;
        }
    }
    if(!otherDir){
        for (int j = midIdx; j < loopSize + midIdx; ++j) {
            if( boundary[j % loopSize] == startAndEnd[2] ) {
               endIdx = j;
            }
        }
    }else{
        // we go in counter direction
         int j = startIdx;
        while ( boundary[j] != startAndEnd[2]){
            if(boundary[j] == startAndEnd[1] ){
                midIdx = j;
            }
            j--;
            if(j<0) j+= loopSize;
        }
        endIdx = j;

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
        next = (otherDir)? (curr-1) : (curr +1)%loopSize;
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
        next = (!otherDir)? (curr+1 )%loopSize : (curr -1);
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

bool vertOnEdge(const VectorXd& R, const VectorXd& Q, VectorXd& p,int v){
    double eps = 0.5;
    auto QR = Q-R;
    auto Qp = Q-p;
cout<<"v= "<<v<<endl;
    VectorXd diff = (Qp).normalized() - (QR).normalized();
    if(v==2918|| v==2938){
        cout<<diff.transpose()<<" diff"<<endl;
        double tt = (p-R)(0)/(QR)(0);
        double t2 =  (p-R)(1)/(QR)(1);
        cout<<tt<<" = t=  "<<t2<<", and makes "<<endl<<(R+tt*(Q-R)).transpose()<<endl<<p.transpose()<<" =? "<<endl;
    }
    if (abs(diff(0))+abs(diff(1)) >eps) return false;
    // p is on the line, but also on the segment?
    // R + t* (Q-R)= p
//    double t = (R-Q).dot(p-Q)/((R-Q).dot(R-Q));
    double t = (p-R)(0)/(QR)(0);
    double tt=  (p-R)(1)/(QR)(1);
//    if((p-R)(0)-(QR)(0)< 0.01){// numerically instable
    if (0 > t || t > 1 ) {// t does not work anyways ,try with tt
        t =tt;
    }
    if (0<=t && t<=1 ) cout<<v<<" : t= "<<t<<", and makes "<<endl<<(R+t*(QR)).transpose()<<endl<<p.transpose()<<" =? "<<endl;
    if(v==2918|| v==2938){
        cout<<"  R "<< R.transpose() <<endl;
        cout<<"Q "<<Q.transpose()<<endl;
        cout<<"  QR "<< QR.transpose() <<endl;


    }
    if (0>t || t>1 ) return false;
    VectorXd posdiff = R+t*(Q-R) - p;
    cout<<"diff from real pos " << abs(posdiff(0)) + abs(posdiff(1)) <<endl;
//    if(v==1556|| v==1557){
//        cout<<posdiff.transpose()<<" posdiff"<<endl;
//    }
    return (abs(posdiff(0)) + abs(posdiff(1)) < 1);


}
void computeAllBetweensNew(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                           vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                           MatrixXd& currPattern, MatrixXd& Vg_to, vector<VectorXd>& polyLineInput, vector<int>& connectedVert) {
    polyLineInput.clear();
    connectedVert.clear();
    cout<<endl<<"Seam Size "<<polylineSelected.size()<<endl<<endl;
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
            if (vertOnEdge(Vg_to.row(v).transpose(), Vg_to.row(v1).transpose(), polylineSelected[1], v)) {
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
            if (vertOnEdge(Vg_to.row(v).transpose(), Vg_to.row(v1).transpose(), polylineSelected[4], v)) {
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
    if(asc){// ascending
        int i = idx4;
        int count = 0;
        while(i !=idx1 && count < boundaryToSearch[patchFrom].size()){
            polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
            cout<<boundaryToSearch[patchFrom][i]<<" 2, asc"<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;
            i++; count++;
            i = i %boundaryToSearch[patchFrom].size();
        }
        polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
        cout<<boundaryToSearch[patchFrom][i]<<" "<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;

    }else{
        int i = idx4;
        int count = 0;
        while(i !=idx1 && count < boundaryToSearch[patchFrom].size()){
            polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
            cout<<boundaryToSearch[patchFrom][i]<<" 2, desc "<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;
            i--; count++;
            if(i <0) i+= boundaryToSearch[patchFrom].size();
        }
        polyLineInput.push_back(currPattern.row(boundaryToSearch[patchFrom][i]).transpose());
        cout<<boundaryToSearch[patchFrom][i]<<" "<<currPattern.row(boundaryToSearch[patchFrom][i])<<endl;

    }



}
void computeAllBetweens(vector<VectorXd>& polylineSelected,vector<int>& polylineIndex, vector<int>& polyLineMeshIndicator,
                   vector<vector<int>>& boundaryL_adaptedFromPattern, vector<vector<int>>& boundaryL_toPattern,
                   MatrixXd& currPattern, MatrixXd& Vg_pattern_orig, vector<VectorXd>& polyLineInput, vector<int>& connectedVert){
    if(polylineSelected.size()>2){
        computeAllBetweensNew(polylineSelected, polylineIndex, polyLineMeshIndicator,
                              boundaryL_adaptedFromPattern, boundaryL_toPattern,
                              currPattern,  Vg_pattern_orig, polyLineInput, connectedVert);
        return;
    }


    polyLineInput.clear();
    connectedVert.clear();

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
                bool inverted = false;
                if(smaller == endIdx) {
                    inverted = true;
                }
                int greater = (endIdx < startIdx) ? startIdx : endIdx;

                int dist = (greater - smaller);
                int otherdist = boundaryToSearch[j].size()-greater + smaller;
                cout<<otherdist<<" betweens "<<dist<<endl;


                if(dist<otherdist){
                    if(!inverted){
                        for(int k = smaller+1; k < greater; k++){
                            polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
                            cout<<v_used.row(boundaryToSearch[j][k])<<endl;
//                            connectedVert.push_back(boundaryToSearch[j][k]);

                        }
                    }else{
                        for(int k= greater-1; k> smaller ; k--){
                            polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
                            cout<<v_used.row(boundaryToSearch[j][k])<<endl;
//                            connectedVert.push_back(boundaryToSearch[j][k]);

                        }
                    }

                }else{
                    if(!inverted){
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
//                        connectedVert.push_back(boundaryToSearch[j][k]);
                    }else{
                        int k = greater; k++;
                        k = k % boundaryToSearch[j].size();
                        while (k!= smaller){
                            polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
                            cout<<v_used.row(boundaryToSearch[j][k])<<endl;
//                            connectedVert.push_back(boundaryToSearch[j][k]);
                            k++;
                            k = k % boundaryToSearch[j].size();
                        }// last one k==greater
//                        polyLineInput.push_back(v_used.row(boundaryToSearch[j][k]));
//                        cout<<v_used.row(boundaryToSearch[j][k])<<"wi"<<endl;
//                        connectedVert.push_back(boundaryToSearch[j][k]);
                    }

                }
            }

        }


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

void mergeTriagulatedAndPattern(const vector<int> &connectedVert, MatrixXd& Vg_retri, MatrixXi& Fg_retri, MatrixXd& currPattern, MatrixXi& Fg_pattern){
    int offset = currPattern.rows();
//    int count = 0;
    MatrixXd newVg (offset+Vg_retri.rows(), 3);
    newVg.block(0,0,offset, 3) = currPattern;

    for(int i =0; i< Vg_retri.rows(); i++){
//        mathcesOne = checkIfMatchesOne(Vg_retri.row(i), connectedVert, currPattern );
        newVg.row(offset+i) = Vg_retri.row(i);
//        cout<<newVg.row(offset+i)<<endl;
        replaceInFaces(i, offset+i, Fg_retri);


    }

    MatrixXi newFg (Fg_pattern.rows()+ Fg_retri.rows(), 3);
    newFg.block(0,0,Fg_pattern.rows(), 3) = Fg_pattern;
    newFg.block(Fg_pattern.rows(), 0, Fg_retri.rows(), 3) = Fg_retri;

    currPattern.resize(newVg.rows(), 3);
    currPattern = newVg;
    Fg_pattern.resize(newFg.rows(), 3);
    Fg_pattern = newFg;

    vector<vector<int>> newBound;
    igl::boundary_loop(Fg_pattern, newBound);
    cout<<"We now have "<<newBound.size()<<" patches. "<<endl;

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

void initialGuessAdaption(MatrixXd& currPattern, MatrixXd& toPattern, MatrixXi& Fg_pattern){

    currPattern.block(0, 0,toPattern.rows(), toPattern.cols())= toPattern;
    cout<<currPattern.row(66)<<" v 66"<<endl ;
    cout<<currPattern.row(368)<<" v 368"<<endl ;

    if(currPattern.rows()>toPattern.rows()){
        // we added some vertices already, map the vertices to the duplicates
        cout<<"add some more. TODO "<<endl;
    }
}
set<int> boundaryVerticesP;
void ensureAngle(MatrixXd& p, MatrixXd& toPattern, MatrixXi& Fg_pattern){
    if(boundaryVerticesP.empty()){
        vector<vector<int>> boundaryVP;
        igl::boundary_loop(Fg_pattern, boundaryVP);
        for(int i=0; i<boundaryVP.size(); i++){
            for(int j=0; j<boundaryVP[i].size(); j++){
                boundaryVerticesP.insert(boundaryVP[i][j]);
            }
        }
    }
    for(int i = 0; i < Fg_pattern.rows(); i++){
        for(int j=0; j<3; j++){
            Vector3d e1 = p.row(Fg_pattern(i, j))-p.row(Fg_pattern(i, (j+1) % 3));
            Vector3d e2 = p.row(Fg_pattern(i, (j+2) % 3 ))-p.row(Fg_pattern(i, (j+1) % 3));
            double newAngle = acos(min(max((e1.normalized()).dot(e2.normalized()), -1.), 1.));

            //https://www.euclideanspace.com/maths/algebra/vectors/angleBetween/ and https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
//            // praise stackoverflow
//
//            auto crossVec = alignFrom.cross(alignTo); // or other way round?
//            if(oldNormalVec.dot(crossVec)<0){
//                newAngle = -newAngle;
//            }

            double newdegree = newAngle*180/M_PI;
            if(newdegree<10){
//                cout<<"old degree: "<<newdegree<<endl;

                MatrixXd rot = MatrixXd::Identity(3, 3);
                double rad =  newAngle; //converting to radian value
                rot(0,0)= cos(rad);
                rot(1, 1) = rot(0,0) ;
                rot(0,1) =  - sin(rad);
                rot(1,0) = sin(rad);
                double stiffness = 0.5;
//                if(boundaryVerticesP.find(Fg_pattern(i, j))!= boundaryVerticesP.end()){
                    VectorXd e1rot= rot * e1;
                    VectorXd dir = e1rot - e1;
                    p.row(Fg_pattern(i, j)) += stiffness * dir;

//                }
//                if(boundaryVerticesP.find(Fg_pattern(i, (j+2) % 3))!= boundaryVerticesP.end()){
                    VectorXd e2rot= rot.transpose() * e2;
                     dir = e2rot - e2;
                    p.row(Fg_pattern(i, (j+2) % 3)) += stiffness * dir;

//                }


                ///ttest
                e1 = p.row(Fg_pattern(i, j))-p.row(Fg_pattern(i, (j+1) % 3));
                e2 = p.row(Fg_pattern(i, (j+2) % 3 ))-p.row(Fg_pattern(i, (j+1) % 3));
                newAngle = acos(min(max((e1.normalized()).dot(e2.normalized()), -1.), 1.));
                newdegree = newAngle*180/M_PI;
//                cout<<"new degree: "<<newdegree<<endl;

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
            cout<<"face "<<i<<" flipped"<<endl;
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
