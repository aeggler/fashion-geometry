//
// Created by Anna Maria Eggler on 30.01.23.
//

#include "preProcessing.h"

#include <Eigen/Dense>
#include <iostream>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/signed_distance.h>
#include <igl/AABB.h>
#include <igl/exact_geodesic.h>
#include <map>


using namespace std;
using namespace Eigen;

void setupCollisionConstraintsCall(Eigen::MatrixXi& collisionVert, vector<int> & pureCollVert, MatrixXd& testMorph_V1left, MatrixXi& testMorph_F1left, MatrixXd& p, int& numVert, double coll_EPS,
                               std::map<int, int> & leftHalfToFullFaceMap, vector<VectorXd> & CleftRight , vector<VectorXd>& NleftRight, VectorXi& closestFaceId,
                               MatrixXd& Vm, MatrixXi& Fm, MatrixXi& Fg){
    CleftRight.clear();
    NleftRight.clear();
    VectorXd S;
    MatrixXd Cleft, Nleft;
    collisionVert = Eigen::MatrixXi::Zero(numVert, 1);
    pureCollVert.clear();
    //new part
    VectorXi closestFaceIdCollision;
    igl::AABB<Eigen::MatrixXd, 3> col_treeLeft;
    col_treeLeft.init(testMorph_V1left, testMorph_F1left);
    MatrixXd Vmleft, FN_mleft,FN_gar, EN_mleft, VN_mleft, VN_gar;
    MatrixXi Fmleft, EMAP_mleft, E_mleft;


    igl::per_face_normals(p, Fg, FN_gar);
    igl::per_vertex_normals(p, Fg, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_gar, VN_gar);


    Vmleft = testMorph_V1left;
    Fmleft = testMorph_F1left;
    igl::per_face_normals(Vmleft, Fmleft, FN_mleft);
    igl::per_vertex_normals(Vmleft, Fmleft, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_mleft, VN_mleft);
    igl::per_edge_normals(Vmleft, Fmleft, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN_mleft,
                          EN_mleft, E_mleft, EMAP_mleft);
    igl::signed_distance_pseudonormal(p, Vmleft, Fmleft, col_treeLeft, FN_mleft, VN_mleft, EN_mleft,
                                      EMAP_mleft, S, closestFaceIdCollision, Cleft, Nleft);

    cout<<numVert<<" num vert "<<Nleft.rows()<<" "<<Nleft.cols()<<" and c "<<Cleft.rows()<<" "<<Cleft.cols()<<endl;
    int collCount = 0;
    for(int i=0; i<numVert; i++){
        // to make sure negative does not count
        if(i==27|| i==28||i ==220){
            cout<<i<<" "<<closestFaceId(i)<<" expecting 51408 "<<S(i)<<endl;
        }
        if(S(i) < coll_EPS && abs(S(i))<50){

            //TODO it might well be that closest Face Id is not always correct?
            //todo or its neighbors!!!
            if(leftHalfToFullFaceMap[closestFaceIdCollision(i)]== closestFaceId(i)){
                collCount++;
                collisionVert(i) = 1;
                CleftRight.push_back(Cleft.row(i));
                NleftRight.push_back(Nleft.row(i));
                pureCollVert.push_back(i);
            }else{
                // wenn die vertex normal von i und die face normal von closestFaceIdCollision
                cout<<"Vert "<<i<<" intersecting face "<<leftHalfToFullFaceMap[closestFaceIdCollision(i)]<<endl;
                double cosVal = VN_gar.row(i).dot(FN_mleft.row(closestFaceIdCollision(i)));
                cout<<S(i)<<" dist and cos val "<<cosVal<<endl;

                // new idea use cosine to check
//                igl::exact_geodesic(Vm, Fm, FS, VS, FT,  VT, D);
//                cout<<"doubt: "<<i<<"  would intersect face "<<closestFaceIdCollision(i)<<" on half model "<<endl;
                if(cosVal<0){//D(0)>50
                    cout<<"don't consider, it's a intersection we dont want to handle, let other side do this "<<endl;

                }else{
                    collCount++;
                    collisionVert(i) = 1;
                    pureCollVert.push_back(i);
                    CleftRight.push_back(Cleft.row(i));
                    NleftRight.push_back(Nleft.row(i));
                }
                cout<<endl;
            }
        }
    }

    //trial adding right side
//    igl::AABB<Eigen::MatrixXd, 3> col_treeRight;
//    col_treeRight.init(testMorph_V1right, testMorph_F1right);
//    MatrixXd Vmleftr, FN_mleftr, EN_mleftr, VN_mleftr;
//    MatrixXi Fmleftr, EMAP_mleftr, E_mleftr;
//
//    Vmleftr = testMorph_V1right;
//    Fmleftr = testMorph_F1right;
//    igl::per_face_normals(Vmleftr, Fmleftr, FN_mleftr);
//    igl::per_vertex_normals(Vmleftr, Fmleftr, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_mleftr, VN_mleftr);
//    igl::per_edge_normals(Vmleftr, Fmleftr, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN_mleftr,EN_mleftr, E_mleftr, EMAP_mleftr);
//    igl::signed_distance_pseudonormal(p, Vmleftr, Fmleftr, col_treeRight, FN_mleftr, VN_mleftr, EN_mleftr,EMAP_mleftr, S, closestFaceIdCollision, Cright, Nright);
//    igl::signed_distance_pseudonormal(p, Vmleftr, Fmleftr, S, closestFaceIdCollision, Cright, Nright);
//
//
//    for(int i=0; i<numVert; i++){
//        if(S(i) < coll_EPS && rightHalfToFullFaceMap[closestFaceIdCollision(i)]== closestFaceId(i)){
//            if(collisionVert(i) == 1){
//                cout<<"it is in both. What to do? "<<endl;
//                continue;
//            }
//            collisionVert(i) = 1;
//            pureCollVert.push_back(i);
//            lorr.push_back(2);
//
//            cout<<"Vertex "<<i<<" collides right"<<endl;
//        }
//    }
//


//    igl::signed_distance_pseudonormal(p, Vm, Fm, col_tree, FN_m, VN_m, EN_m, EMAP_m, S, closestFaceId, C, N);
//    int collCount=0;
//    for(int i=0; i<numVert; i++){
//        if(S(i)<coll_EPS){
//            collCount++;
//            collisionVert(i)=1;
//            pureCollVert.push_back(i);
//
//        }
//    }
    if(pureCollVert.size()!= collCount){
        cout<<" size problem"<<endl;
    }
}