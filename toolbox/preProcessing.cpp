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
#include<Eigen/SparseCholesky>
#include "constraint_utils.h"
#include <igl/writeOBJ.h>

using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double, RowMajor> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
Eigen::SparseMatrix<double, RowMajor> A;
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholSolver;
int numVertPattern =0;

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
// attention we lose face to face correspondance between the original face and the new -> create a map from halfPatternFaceToFullPatternFace
// decision if it is in the half pattern depends on x- coordinate on 3D garment
void createHalfSewingPattern(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, MatrixXd& Vg_pattern_half, MatrixXi& Fg_pattern_half,
                             map<int, int>& halfPatternFaceToFullPatternFace, map<int, int>& fullPatternFaceToHalfPatternFace, map<int, int>& halfPatternVertToFullPatternVert ,
                             map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& insertedIdxToPatternVert, VectorXi& isLeftVertPattern,
                             MatrixXd& R_sym, VectorXd& T_sym, MatrixXd& rightVert){
   int n = Vg.rows();
   int m  = Fg.rows();
//   cout<<endl<<"in half sewing pattern"<<endl;
    VectorXi isLeftVert(Vg.rows());
    VectorXi isRightVert(Vg.rows());

    isLeftVert.setConstant(-1);
    isRightVert.setConstant(-1);

    int leftCount = 0; int rightCount = 0 ;
    for(int i=0; i<n; i++){
        if(Vg(i, 0) <= 0){
            isLeftVert(i) = 1;
            leftCount++;
        }
        if(Vg(i, 0) >= 0){
            isRightVert(i) = 1;
            rightCount++;
        }
    }
    cout<<leftCount<<" left and right in 3D "<<rightCount<<endl;
    MatrixXd Vg_half(leftCount, 3);
    map<int, int> halfVertToFullVert, fullVertToHalfVert;
    int idx=0;
    for(int i=0; i<n; i++){
        if(isLeftVert(i)){
            Vg_half.row(idx) = Vg.row(i);
            halfVertToFullVert[idx] = i;
            fullVertToHalfVert[i] = idx;
            idx++;

        }
    }
    cout<<leftCount<<" left and right in 3D after "<<rightCount<<endl;

    VectorXi isRightVertPattern(Vg_pattern.rows());
    isRightVertPattern.setConstant(-1);
    isLeftVertPattern.resize(Vg_pattern.rows());
    isLeftVertPattern.setConstant(-1);
    VectorXi isLeftFace(m);
    isLeftFace.setConstant(0);
    int faceCount = 0;
    cout<<"start iterating"<<endl;
    for(int i=0; i<m; i++) {
        int v0 = Fg(i, 0);
        int v1 = Fg(i, 1);
        int v2 = Fg(i, 2);
        isLeftVertPattern(Fg_pattern(i, 0)) = (isLeftVert(v0) == 1) ? 1 : 0;
        isLeftVertPattern(Fg_pattern(i, 1)) = (isLeftVert(v1) == 1) ? 1 : 0;
        isLeftVertPattern(Fg_pattern(i, 2)) = (isLeftVert(v2) == 1) ? 1 : 0;

        isRightVertPattern(Fg_pattern(i, 0)) = (isRightVert(v0) == 1) ? 1 : 0;
        isRightVertPattern(Fg_pattern(i, 1)) = (isRightVert(v1) == 1) ? 1 : 0;
        isRightVertPattern(Fg_pattern(i, 2)) = (isRightVert(v2) == 1) ? 1 : 0;

        if(isLeftVert(v0) == -1 && isLeftVert(v1) == -1 && isLeftVert(v2) == -1){
            // all on the right, ignore it
            isLeftFace(i) = -1;
        }else if(isLeftVert(v0) == 1 && isLeftVert(v1) == 1 && isLeftVert(v2) == 1){
            // all on the left -> just take it
            isLeftFace(i) = 1;
            faceCount ++;
        }else if (Vg(v0, 0)<0 ||Vg(v1, 0)<0 ||Vg(v2, 0)<0  ){
            faceCount++;
        }else{
            isLeftFace(i) = -1;
            isLeftVertPattern(Fg_pattern(i, 0))=0;
            isLeftVertPattern(Fg_pattern(i, 1))=0;
            isLeftVertPattern(Fg_pattern(i, 2))=0;
        }

        if(isRightVert(v0) == -1 && isRightVert(v1) == -1 && isRightVert(v2) == -1){

        }else if(isRightVert(v0) == 1 && isRightVert(v1) == 1 && isRightVert(v2) == 1){

        }else if (Vg(v0, 0)>0 ||Vg(v1, 0)>0 ||Vg(v2, 0)>0  ){

        }else{
            isRightVertPattern(Fg_pattern(i, 0))=0;
            isRightVertPattern(Fg_pattern(i, 1))=0;
            isRightVertPattern(Fg_pattern(i, 2))=0;
        }
    }
    idx=0; int rightIdx = 0;

    cout<<isLeftVertPattern.sum()<<" and right side sum "<<isRightVertPattern.sum()<<endl;
    int patternHalfVert = isLeftVertPattern.sum();
    Vg_pattern_half.resize(patternHalfVert, 3);
     rightVert.resize(isRightVertPattern.sum(), 3);

    for(int i=0; i< Vg_pattern.rows(); i++){
        if(isLeftVertPattern(i)){
            Vg_pattern_half.row(idx) = Vg_pattern.row(i);
            halfPatternVertToFullPatternVert[idx] = i;
            fullPatternVertToHalfPatternVert[i] = idx;
            idx++;

        }
        if(isRightVertPattern(i)){
            rightVert.row(rightIdx) = Vg_pattern.row(i);
            rightIdx ++;
        }
    }



    vector<VectorXd> addedVert;
    int newIdx = leftCount;
    Fg_pattern_half.resize(faceCount, 3);

    idx = 0;
    for(int i = 0; i<m; i++){
        if(isLeftFace(i) == 1){
            Fg_pattern_half(idx, 0)= fullPatternVertToHalfPatternVert[Fg_pattern(i, 0)];
            Fg_pattern_half(idx, 1)= fullPatternVertToHalfPatternVert[Fg_pattern(i, 1)];
            Fg_pattern_half(idx, 2)= fullPatternVertToHalfPatternVert[Fg_pattern(i, 2)];

            halfPatternFaceToFullPatternFace[idx] = i;
            fullPatternFaceToHalfPatternFace[i] = idx;
            idx++;

        }
    }


    // we hve rightVert adn Vg_pattern_half for right and left vertices.
    // now find the syymetry using procrustes with reflection
    cout<<rightVert.rows()<<" right and left 2D verts "<<Vg_pattern_half.rows()<<endl;
    MatrixXd R; VectorXd T;
    MatrixXd input = Vg_pattern_half.block(0,0,Vg_pattern_half.rows(), 2);
    MatrixXd input2 = rightVert.block(0,0,rightVert.rows(), 2);

    procrustes( input, input2,R, T);
//    cout<<R<<" rotation pro "<<endl<<T<<endl;
//    cout<<R<<" rotation wo ref "<<endl;
    R_sym.resize(3, 3);
    R_sym.block(0,0,2,2) = R;
    R_sym.row(2).setConstant(0);
    R_sym.col(2).setConstant(0);
    R_sym(2,2) = 1;

    T_sym.resize(3);T_sym(0)= T(0); T_sym(1) = T(1);  T_sym(2) = 0;
        cout<<R_sym<<" after"<<endl;


}

void preProcessGarment(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern){
//    map<double, std::pair<int, int>> yToFaceAndIdx;
//    int vgsize = Vg_pattern.rows();
//    for(int i=0; i<Fg.rows(); i++){
//        bool hasLeft = false;
//        bool hasRight = false;
//        Vector3i LR;
//        for (int j=0; j<3; j++){
//            double x = Vg(Fg(i, j), 0);
//            if(x==0) {
//                continue;
//            }
//            if(x<0){
//                hasLeft = true;
//                LR(j)= 1;
//            }else{
//                hasRight = true;
//                LR(j) = 0;
//            }
//        }
//        if(hasLeft && hasRight){
//            cout<<i<<" is in the middle"<<endl;
//            int otherSide=0;
//
//            if(LR.sum()==1){
//                // search for the single left
//                while(LR(otherSide) != 1){
//                    otherSide++;
//                }
//            }else{
//                while(LR(otherSide) != 0){
//                    otherSide++;
//                }
//            }
//            // otherSide is the index that has cuts on both sideds!
//            int v1 = Fg(i, otherSide);
//            int v2 = Fg(i, (otherSide+1) % 3 );
//            int v3 = Fg(i, (otherSide+2) % 3 );
//            VectorXd edge2 = Vg.row(v2) - Vg.row(v1);
//            VectorXd edge3 = Vg.row(v3) - Vg.row(v1);
//            double t2 = -Vg(v1,0) / edge2(0);
//            double t3 = -Vg(v1,0) / edge3(0);
//            VectorXd newPos1 = Vg.row(v1) + t2 * edge2.transpose();
//            newPos1(0)= 0;
//            VectorXd newPos2 = Vg.row(v1) + t3 * edge3.transpose();
//            newPos2(0) = 0;
//
//            int idx1 , idx2; int fac1, idfac1, fac2, idfac2;
//            bool new1= false; bool new2 = false;
//            if(i==1132){
//                cout<<newPos1(1)<<" "<<(yToFaceAndIdx.find(newPos1(1)) == yToFaceAndIdx.end())<<" ad vert "<<Vg(1804,1)<< " "<<(newPos1(1) == Vg(1804,1))<<endl;
//                fac1= 457;
//                idfac1 = 0;
//                idx1= 1804;
//            }else
//            if(yToFaceAndIdx.find(newPos1(1)) == yToFaceAndIdx.end()){
//                yToFaceAndIdx[newPos1(1)] = std::make_pair(i, (otherSide+1) % 3 );
//
//                int vgrow =  Vg.rows();
//                MatrixXd Vgnew( vgrow + 1, 3);
//                Vgnew.block(0,0,vgrow, 3) = Vg;
//                Vgnew.row(vgrow ) = newPos1;
//                Vg.resize(vgrow+1, 3);
//                Vg = Vgnew;
//                idx1 = vgrow; new1= true;
//            }else{
//                fac1 = yToFaceAndIdx[newPos1(1)].first;
//                idfac1 = yToFaceAndIdx[newPos1(1)].second;
//
//                idx1 =Fg (fac1, idfac1 );
////                cout<<"found at idx "<<idx1<<endl;
//
//            }
//            if(yToFaceAndIdx.find(newPos2(1)) == yToFaceAndIdx.end()){
//                yToFaceAndIdx[newPos2(1)] = std::make_pair(i, (otherSide+2) % 3 );
//
//                int vgrow =  Vg.rows();
//                MatrixXd Vgnew( vgrow + 1, 3);
//                Vgnew.block(0,0,vgrow, 3) = Vg;
//                Vgnew.row(vgrow ) = newPos2;
//                Vg.resize(vgrow+1, 3);
//                Vg = Vgnew;
//                idx2 = vgrow;
//                new2= true;
//            }else{
//                fac2 = yToFaceAndIdx[newPos2(1)].first;
//                idfac2 = yToFaceAndIdx[newPos2(1)].second;
//                idx2 =Fg ( fac2, idfac2 );
////                cout<<"found at idx "<<idx2<<endl;
//
//            }
//
//            int fgrow = Fg.rows();
//            MatrixXi Fgnew (fgrow+2, 3);
//            Fgnew.block(0,0,fgrow, 3) = Fg;
//
//            Fgnew(i, (otherSide+1) % 3 ) = idx1;
//            Fgnew(i, (otherSide+2) % 3 ) = idx2;
//            Fgnew(fgrow, 0) = v2;
//            Fgnew(fgrow, 1) = v3;
//            Fgnew(fgrow, 2) = idx1;
//            Fgnew(fgrow + 1, 0) = idx1;
//            Fgnew(fgrow + 1, 1) = v3;
//            Fgnew(fgrow + 1, 2) = idx2;
//            Fg.resize(Fgnew.rows(), 3);
//            Fg= Fgnew;
//
//            // now change the pattern too!
//            MatrixXd bary, bary2;
//            MatrixXd input (1, 3);
//            input.row(0) = newPos1;
//
//            igl::barycentric_coordinates(input, Vg.row(v1), Vg.row(v2), Vg.row(v3), bary);
//            newPos1 = bary(0,0) * Vg_pattern.row(Fg_pattern(i, otherSide))+
//                    bary(0,1) * Vg_pattern.row(Fg_pattern(i, (otherSide+1) % 3 ))+
//                    bary(0,2) * Vg_pattern.row(Fg_pattern(i, (otherSide+2) % 3 )) ;
//
//            input.row(0) = newPos2;
//            igl::barycentric_coordinates(input, Vg.row(v1), Vg.row(v2), Vg.row(v3), bary2);
//
//            newPos2 = bary2(0,0) * Vg_pattern.row(Fg_pattern(i, otherSide))+
//                      bary2(0,1) * Vg_pattern.row(Fg_pattern(i, (otherSide+1) % 3 ))+
//                      bary2(0,2) * Vg_pattern.row(Fg_pattern(i, (otherSide+2) % 3 )) ;
//
//            if(new1){
//                int vgrow =  Vg_pattern.rows();
//                MatrixXd Vgnew( vgrow + 2, 3);
//                Vgnew.block(0,0,vgrow, 3) = Vg_pattern;
//
//                Vgnew.row(vgrow ) = newPos1;
//                Vgnew.row(vgrow+1 ) = newPos1;
//                Vg_pattern.resize(vgrow+2, 3);
//                Vg_pattern = Vgnew;
//                idx1 = vgrow;
//            }
//
//            if(new2){
//                int vgrow =  Vg_pattern.rows();
//                MatrixXd Vgnew( vgrow + 2, 3);
//                Vgnew.block(0,0,vgrow, 3) = Vg_pattern;
//                Vgnew.row(vgrow ) = newPos2;
//                Vgnew.row(vgrow +1) = newPos2;
//
//                Vg_pattern.resize(vgrow + 2, 3);
//                Vg_pattern = Vgnew;
//                idx2 = vgrow;
//            }
//
//
//            fgrow = Fg_pattern.rows();
//            Fgnew.resize (fgrow+2, 3);
//            Fgnew.block(0,0,fgrow, 3) = Fg_pattern;
//            if (!new1)  idx1 = Fg_pattern(fac1, idfac1);
//            if (!new2)  idx2 = Fg_pattern(fac2, idfac2);
//
//            Fgnew(i, (otherSide+1) % 3 ) = idx1;   // if(LR.sum()==1)   Fgnew(i, (otherSide+1) % 3 )++;
//            Fgnew(i, (otherSide+2) % 3 ) = idx2; //   if(LR.sum()==1)   Fgnew(i, (otherSide+2) % 3 )++;
//
//            Fgnew(fgrow, 0) = Fg_pattern(i, (otherSide + 1) % 3 );
//            Fgnew(fgrow, 1) = Fg_pattern(i, (otherSide + 2) % 3 );
//            Fgnew(fgrow, 2) = idx1;           //      if(LR.sum()==2 && new1) Fgnew(fgrow, 2)= idx1+1;
//////
//            Fgnew(fgrow + 1, 0) = idx1;          //   if(LR.sum()==2) Fgnew(fgrow +1 , 0)++;
//            Fgnew(fgrow + 1, 1) = Fg_pattern(i, (otherSide + 2) % 3 );
//            Fgnew(fgrow + 1, 2) = idx2;           //  if(LR.sum()==2) Fgnew(fgrow +1 , 2)++;
//
//            Fg_pattern.resize(Fgnew.rows(), 3);
//            Fg_pattern= Fgnew;
//
//
//
//        }
//    }
//
//    igl::writeOBJ("dress_3d.obj", Vg, Fg);
//    igl::writeOBJ("dress_2d.obj", Vg_pattern, Fg_pattern);
//    int added = Vg_pattern.rows()- vgsize;
//    int vgnewsize = Vg_pattern.rows();
//    MatrixXd Vgp (Vg_pattern.rows()+ added, 3);
//    Vgp.block(0,0,vgnewsize, 3) = Vg_pattern;
//    Vgp.block(vgnewsize, 0, added, 3) = Vg_pattern.block(vgsize, 0, added, 3);
//    cout<<"continue"<<endl;
//    // we duplicate the new vertices to split them
//    for (int i=0; i<Fg_pattern.rows(); i++){
//        for(int j=0; j<3; j++){
//            cout<<i<<" j"<<j<<" " <<Fg_pattern(i,j)<<endl;
//            if(Fg_pattern(i,j)>= vgsize){
//
//                int other = (j+1)%3; if(i== 553){
//                    cout<< " other "<<other <<endl;
//                }
//                if(i== 553){
//                    cout<<Fg.rows()<< " rows "<<endl;
//                    cout<<Fg(i, other)<< " element "<<endl;
//                    cout<< Vg.row(Fg(i, other))<< endl;
//
//                }
//                while(Vg(Fg(i, other), 0) == 0){
//                    other++;
//                    other %= 3;
//                    cout<<"update "<<other<<endl;
//                }// find one that is not 0
//                if(i== 553){
//                    cout<< " other "<<endl;
//                }
//                bool isLeft = false;
//                if(Vg(Fg(i, other), 0)<0){
//                    isLeft= true;
//                    if(i== 553){
//                        cout<< " left "<<endl;
//                    }
//                }
//                if(i== 553){
//                    cout<< " before "<<endl;
//                }
//                if(isLeft){
//                    Fg_pattern(i,j) +=added;
//                }
//                if(i== 553){
//                    cout<< " fin "<<endl;
//                }
//            }
//        }
//    }
//    igl::writeOBJ("dress_2d_dupl.obj", Vgp, Fg_pattern);


}