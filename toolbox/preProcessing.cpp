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

MatrixXd FN_mleft,FN_mright, FN_gar, EN_mleft, EN_mright, VN_mleft,VN_mright, VN_gar;
MatrixXi EMAP_mleft, E_mleft, E_mright, EMAP_mright;
MatrixXd Cleft, Cright, Nleft, Nright;
VectorXd Sleft, Sright;
igl::AABB<Eigen::MatrixXd, 3> col_treeLeft, col_treeRight;
VectorXi closestFaceIdLeft, closestFaceIdright;


void initCollMeshCall( MatrixXd& Vm_left, MatrixXi& Fm_left,
                      MatrixXd& Vm_right, MatrixXi& Fm_right){

    col_treeLeft.init(Vm_left, Fm_left);
    col_treeRight.init(Vm_right, Fm_right);

    igl::per_face_normals(Vm_left, Fm_left, FN_mleft);
    igl::per_vertex_normals(Vm_left, Fm_left, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_mleft, VN_mleft);
    igl::per_edge_normals(Vm_left, Fm_left, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN_mleft,
                          EN_mleft, E_mleft, EMAP_mleft);
    // right half mesh
    igl::per_face_normals(Vm_right, Fm_right, FN_mright);
    igl::per_vertex_normals(Vm_right, Fm_right, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_mright, VN_mright);
    igl::per_edge_normals(Vm_right, Fm_right, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN_mright, EN_mright, E_mright, EMAP_mright);

}

void setupCollisionConstraintsCall(Eigen::MatrixXi& collisionVert, vector<int> & pureCollVert,
                                   MatrixXd& Vm_left, MatrixXi& Fm_left,
                                   MatrixXd& Vm_right, MatrixXi& Fm_right, MatrixXd& p, int& numVert, double coll_EPS,
                               std::map<int, int> & leftHalfToFullFaceMap,   std::map<int, int> & rightHalfToFullFaceMap, vector<VectorXd> & CleftRight , vector<VectorXd>& NleftRight, VectorXi& closestFaceId,
                               MatrixXd& Vm, MatrixXi& Fm, MatrixXi& Fg){
    CleftRight.clear();
    NleftRight.clear();
    double absThereshold = 10;
    collisionVert = Eigen::MatrixXi::Zero(numVert, 1);
    pureCollVert.clear();
    //new part
    igl::per_face_normals(p, Fg, FN_gar);
    igl::per_vertex_normals(p, Fg, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_gar, VN_gar);


    igl::signed_distance_pseudonormal(p, Vm_left, Fm_left, col_treeLeft, FN_mleft, VN_mleft, EN_mleft,
                                      EMAP_mleft, Sleft, closestFaceIdLeft, Cleft, Nleft);
    int collCount = 0;
    for(int i=0; i<numVert; i++){
        // to make sure negative does not count

        if(Sleft(i) < coll_EPS && abs(Sleft(i))<absThereshold){// assuming no collision is bad enough to cause -10 problem!

            //TODO it might well be that closest Face Id is not always correct?
            //todo or its neighbors!!!
            if(leftHalfToFullFaceMap[closestFaceIdLeft(i)]== closestFaceId(i)){
                collCount++;
                collisionVert(i) = 1;
                CleftRight.push_back(Cleft.row(i));
                NleftRight.push_back(Nleft.row(i));
                pureCollVert.push_back(i);
            }else{
                // wenn die vertex normal von i und die face normal von closestFaceIdCollision
//                cout<<"Vert "<<i<<" intersecting face "<<leftHalfToFullFaceMap[closestFaceIdLeft(i)]<<endl;
                double cosVal = VN_gar.row(i).dot(FN_mleft.row(closestFaceIdLeft(i)));

                if(cosVal<0){//D(0)>50
//                    cout<<"don't consider, it's a intersection we dont want to handle, let other side do this "<<endl;

                }else{
                    collCount++;
                    collisionVert(i) = 1;
                    pureCollVert.push_back(i);
                    CleftRight.push_back(Cleft.row(i));
                    NleftRight.push_back(Nleft.row(i));
                }
            }
        }
    }

    //trial adding right side
    igl::signed_distance_pseudonormal(p, Vm_right, Fm_right, col_treeRight, FN_mright, VN_mright, EN_mright,EMAP_mright, Sright, closestFaceIdright, Cright, Nright);
    for(int i=0; i<numVert; i++){
        if(Sright(i) < coll_EPS && abs(Sright(i))<absThereshold) {

            if(rightHalfToFullFaceMap[closestFaceIdright(i)]== closestFaceId(i)){
                if (collisionVert(i) == 1) {
                    cout << i<<" with same index it is in both. What to do? " << endl;// ignoore it. most likely it is because left and right are not properly seperated
                    continue;
                }
                collCount++;
                collisionVert(i) = 1;
                CleftRight.push_back(Cright.row(i));
                NleftRight.push_back(Nright.row(i));
                pureCollVert.push_back(i);
            }else{
                // wenn die vertex normal von i und die face normal von closestFaceIdCollision
                double cosVal = VN_gar.row(i).dot(FN_mright.row(closestFaceIdright(i)));
                if(cosVal>0){//D(0)>50
                    if (collisionVert(i) == 1) {
                        cout << i<<" with cos it is in both. What to do? " <<rightHalfToFullFaceMap[closestFaceIdright(i)]<<" "<<leftHalfToFullFaceMap[closestFaceIdLeft(i)]<< endl;
                        continue;
                    }
                    collCount++;
                    collisionVert(i) = 1;
                    pureCollVert.push_back(i);
                    CleftRight.push_back(Cright.row(i));
                    NleftRight.push_back(Nright.row(i));
                }
            }
        }
    }

    if(pureCollVert.size()!= collCount){
        cout<<" size problem"<<endl;
    }
}