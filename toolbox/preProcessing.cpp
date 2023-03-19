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
#include <igl/readOBJ.h>
#include <igl/signed_distance.h>
#include <igl/barycentric_coordinates.h>
#include <igl/barycentric_interpolation.h>
#include <igl/vertex_components.h>
#include "igl/adjacency_list.h"
#include "adjacency.h"
#include <igl/vertex_components.h>
#include <igl/facet_components.h>

using namespace std;
using namespace Eigen;
using namespace igl;

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
void vertex_componentsBasedOnFacet(MatrixXi& Fg_pattern, VectorXi& componentIdPerFace,VectorXi& componentIdPerVert, int n){
    componentIdPerVert.resize(n);
    for(int i=0; i<Fg_pattern.rows(); i++){
        for(int j=0; j<3; j++){
            componentIdPerVert(Fg_pattern(i,j)) = componentIdPerFace(i);
        }
    }

}

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
                             map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& insertedIdxToPatternVert, VectorXi& isLeftVertPattern,MatrixXd& rightVert){
   int n = Vg.rows();
   int m  = Fg.rows();
//   cout<<endl<<"in half sewing pattern"<<endl;
    VectorXi isLeftVert(n);
    VectorXi isRightVert(n);

    isLeftVert.setConstant(-1);
    isRightVert.setConstant(-1);

    int leftCount = 0; int rightCount = 0 ;
    for(int i=0; i<n; i++){
        if(Vg(i, 0) <= 0.2){
            isLeftVert(i) = 1;
            leftCount++;
        }
        if(Vg(i, 0) >= -0.2){
            isRightVert(i) = 1;
            rightCount++;
        }
    }
    cout<<leftCount<<" left and right in 3D "<<rightCount<<endl;
    MatrixXd Vg_half(leftCount, 3);
    map<int, int> halfVertToFullVert, fullVertToHalfVert;
    int idx=0;
    for(int i=0; i<n; i++){
        if(isLeftVert(i)==1){
            Vg_half.row(idx) = Vg.row(i);
            halfVertToFullVert[idx] = i;
            fullVertToHalfVert[i] = idx;
            idx++;

        }
    }
//    cout<<leftCount<<" left and right in 3D after "<<rightCount<<endl;

    VectorXi isRightVertPattern(Vg_pattern.rows());
    isRightVertPattern.setConstant(-1);
    isLeftVertPattern.resize(Vg_pattern.rows());
    isLeftVertPattern.setConstant(-1);
    VectorXi isLeftFace(m);
    isLeftFace.setConstant(0);
    int faceCount = 0;
//    cout<<"start iterating"<<endl;
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
//            faceCount++;
        }else{
            cout<<"EVER HERE??"<<endl; 
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
cout<<"left face pattern count "<< faceCount<<endl;
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
igl::writeOBJ("halfPattern.obj", Vg_pattern_half,  Fg_pattern_half);
    cout<<" after"<<endl;

}
void insertPlane(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, string garment){
    // we insert a plane along the x=0 axis, add vertices and introduce new faces
    // this ensures we can split the garment safely in the pre processing
    // we do the same for the pattern and duplicate x=0 vertices -> to make sure the patch is acutally disconnected and we can take only the half patch for symetry
    map<double, std::pair<int, int>> yToFaceAndIdx;
    int vgsize = Vg_pattern.rows();
    std::vector< std::vector<int> > vfAdj,vfAdjG, vvAdj;
    createVertexFaceAdjacencyList(Fg, vfAdjG);

    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    igl::adjacency_list(Fg_pattern, vvAdj);


    for(int i=0; i<Fg.rows(); i++){
        bool hasLeft = false;
        bool hasRight = false;
        Vector3i LR;
        for (int j=0; j<3; j++){
            double x = Vg(Fg(i, j), 0);
            if(x==0) {
                continue;
            }
            if(x<0){
                hasLeft = true;
                LR(j)= 1;
            }else{
                hasRight = true;
                LR(j) = 0;
            }
        }
        if(hasLeft && hasRight){
//            cout<<i<<" is in the middle"<<endl;
            int otherSide=0;

            if(LR.sum()==1){
                // search for the single left
                while(LR(otherSide) != 1){
                    otherSide++;
                }
            }else{
                while(LR(otherSide) != 0){
                    otherSide++;
                }
            }
            // otherSide is (one of the two) index that has cuts on both sides!
            int v1 = Fg(i, otherSide);
            int v2 = Fg(i, (otherSide+1) % 3 );
            int v3 = Fg(i, (otherSide+2) % 3 );
            VectorXd edge2 = Vg.row(v2) - Vg.row(v1);
            VectorXd edge3 = Vg.row(v3) - Vg.row(v1);
            double t2 = -Vg(v1,0) / edge2(0);
            double t3 = -Vg(v1,0) / edge3(0);
            VectorXd newPos1 = Vg.row(v1) + t2 * edge2.transpose();
            newPos1(0)= 0;
            VectorXd newPos2 = Vg.row(v1) + t3 * edge3.transpose();
            newPos2(0) = 0;

            int idx1 , idx2; int fac1, idfac1, fac2, idfac2;
            bool new1= false; bool new2 = false;
            bool extra1= false; bool extra2 = false;
            if(yToFaceAndIdx.find(newPos1(1)) == yToFaceAndIdx.end()){
                yToFaceAndIdx[newPos1(1)] = std::make_pair(i, (otherSide+1) % 3 );

                int vgrow =  Vg.rows();
                MatrixXd Vgnew( vgrow + 1, 3);
                Vgnew.block(0,0,vgrow, 3) = Vg;
                Vgnew.row(vgrow ) = newPos1;
                Vg.resize(vgrow+1, 3);
                Vg = Vgnew;
                idx1 = vgrow; new1= true;
            }else{
                fac1 = yToFaceAndIdx[newPos1(1)].first;
                idfac1 = yToFaceAndIdx[newPos1(1)].second;

                idx1 =Fg (fac1, idfac1 );
                if(i== 1128 && garment== "skirt_no2") extra1 = true;
            }
            if(yToFaceAndIdx.find(newPos2(1)) == yToFaceAndIdx.end()){
                yToFaceAndIdx[newPos2(1)] = std::make_pair(i, (otherSide+2) % 3 );

                int vgrow =  Vg.rows();
                MatrixXd Vgnew( vgrow + 1, 3);
                Vgnew.block(0,0,vgrow, 3) = Vg;
                Vgnew.row(vgrow ) = newPos2;
                Vg.resize(vgrow+1, 3);
                Vg = Vgnew;
                idx2 = vgrow;
                new2= true;
            }else{
                fac2 = yToFaceAndIdx[newPos2(1)].first;
                idfac2 = yToFaceAndIdx[newPos2(1)].second;
                idx2 =Fg ( fac2, idfac2 );
                if(i== 1256 && garment== "skirt_no2") extra2 = true;

//                cout<<"found at idx "<<idx2<<endl;

            }
            if(idx1 == 1256) cout<<extra1<<" "<<extra2<<endl;
            int fgrow = Fg.rows();
            MatrixXi Fgnew (fgrow+2, 3);
            Fgnew.block(0,0,fgrow, 3) = Fg;

            Fgnew(i, (otherSide+1) % 3 ) = idx1;
            Fgnew(i, (otherSide+2) % 3 ) = idx2;
            Fgnew(fgrow, 0) = v2;
            Fgnew(fgrow, 1) = v3;
            Fgnew(fgrow, 2) = idx1;
            Fgnew(fgrow + 1, 0) = idx1;
            Fgnew(fgrow + 1, 1) = v3;
            Fgnew(fgrow + 1, 2) = idx2;
            Fg.resize(Fgnew.rows(), 3);
            Fg= Fgnew;

            // now change the pattern too!
            MatrixXd bary, bary2;
            MatrixXd input (1, 3);
            input.row(0) = newPos1;

            igl::barycentric_coordinates(input, Vg.row(v1), Vg.row(v2), Vg.row(v3), bary);
            newPos1 = bary(0,0) * Vg_pattern.row(Fg_pattern(i, otherSide))+
                    bary(0,1) * Vg_pattern.row(Fg_pattern(i, (otherSide+1) % 3 ))+
                    bary(0,2) * Vg_pattern.row(Fg_pattern(i, (otherSide+2) % 3 )) ;

            input.row(0) = newPos2;
            igl::barycentric_coordinates(input, Vg.row(v1), Vg.row(v2), Vg.row(v3), bary2);

            newPos2 = bary2(0,0) * Vg_pattern.row(Fg_pattern(i, otherSide))+
                      bary2(0,1) * Vg_pattern.row(Fg_pattern(i, (otherSide+1) % 3 ))+
                      bary2(0,2) * Vg_pattern.row(Fg_pattern(i, (otherSide+2) % 3 )) ;

            if(new1 || extra1){
                int vgrow =  Vg_pattern.rows();
                MatrixXd Vgnew( vgrow + 2, 3);
                Vgnew.block(0,0,vgrow, 3) = Vg_pattern;

                Vgnew.row(vgrow ) = newPos1;
                Vgnew.row(vgrow+1 ) = newPos1;
                Vg_pattern.resize(vgrow+2, 3);
                Vg_pattern = Vgnew;
                idx1 = vgrow;
            }

            if(new2|| extra2 ){
                int vgrow =  Vg_pattern.rows();
                MatrixXd Vgnew( vgrow + 2, 3);
                Vgnew.block(0,0,vgrow, 3) = Vg_pattern;
                Vgnew.row(vgrow ) = newPos2;
                Vgnew.row(vgrow +1) = newPos2;

                Vg_pattern.resize(vgrow + 2, 3);
                Vg_pattern = Vgnew;
                idx2 = vgrow;
            }


            fgrow = Fg_pattern.rows();
            Fgnew.resize (fgrow+2, 3);
            Fgnew.block(0,0,fgrow, 3) = Fg_pattern;
            if (!new1  && !extra1)  idx1 = Fg_pattern(fac1, idfac1);
            VectorXi pattComp;
            igl::facet_components(Fg_pattern, pattComp);


            if (!new2 && !extra2)  idx2 = Fg_pattern(fac2, idfac2);

            Fgnew(i, (otherSide+1) % 3 ) = idx1;
            Fgnew(i, (otherSide+2) % 3 ) = idx2;

            Fgnew(fgrow, 0) = Fg_pattern(i, (otherSide + 1) % 3 );
            Fgnew(fgrow, 1) = Fg_pattern(i, (otherSide + 2) % 3 );
            Fgnew(fgrow, 2) = idx1;
////
            Fgnew(fgrow + 1, 0) = idx1;
            Fgnew(fgrow + 1, 1) = Fg_pattern(i, (otherSide + 2) % 3 );
            Fgnew(fgrow + 1, 2) = idx2;

            Fg_pattern.resize(Fgnew.rows(), 3);
            Fg_pattern= Fgnew;

        }
    }

    igl::writeOBJ("dress_3d.obj", Vg, Fg);
    igl::writeOBJ("dress_2d.obj", Vg_pattern, Fg_pattern);
    int added = Vg_pattern.rows()- vgsize;
    int vgnewsize = Vg_pattern.rows();
    MatrixXd Vgp (Vg_pattern.rows()+ added, 3);
    Vgp.block(0,0,vgnewsize, 3) = Vg_pattern;
    Vgp.block(vgnewsize, 0, added, 3) = Vg_pattern.block(vgsize, 0, added, 3);
    cout<<"continue"<<endl;
    // we duplicate the new vertices to split them
    for (int i=0; i<Fg_pattern.rows(); i++){
        for(int j=0; j<3; j++){
            if(Fg_pattern(i,j)>= vgsize){

                int other = (j+1)%3;

                while(Vg(Fg(i, other), 0) == 0){
                    other++;
                    other %= 3;
                }// find one that is not 0

                bool isLeft = false;
                if(Vg(Fg(i, other), 0)<0){
                    isLeft= true;
                }

                if(isLeft){
                    Fg_pattern(i,j) +=added;
                }

            }
        }
    }
    igl::writeOBJ("dress_2d_dupl.obj", Vgp, Fg_pattern);
    Vg_pattern.resize(Vgp.rows(), Vgp.cols()); Vg_pattern = Vgp;

}
void laplFilter(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern){
    MatrixXd newPattern = Vg_pattern;
    vector<vector<int>> vvAdj, vfAdj;
    igl::adjacency_list(Fg_pattern, vvAdj);
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);

    for(int i=0; i<Vg_pattern.rows(); i++){
        if(isBoundaryVertex(Vg_pattern,i, vvAdj,vfAdj))continue;

        VectorXd newPos = VectorXd::Zero(3);
        int count = 0;
        for(auto it: vvAdj[i]){
            newPos += Vg_pattern.row(it);
            count++;
        }
        newPos /= count;
        newPattern.row(i) = newPos;
    }
    igl::writeOBJ("laplFiltered.obj", newPattern, Fg_pattern);
    Vg_pattern= newPattern;
    MatrixXd newnewPattern = newPattern;
    for(int i=0; i<Vg_pattern.rows(); i++){
        if(isBoundaryVertex(newPattern,i, vvAdj,vfAdj))continue;
        VectorXd newPos = VectorXd::Zero(3);
        int count = 0;
        for(auto it: vvAdj[i]){
            newPos += newPattern.row(it);
            count++;
        }
        newPos /= count;
        newnewPattern.row(i) = newPos;
    }
    igl::writeOBJ("laplFiltered2.obj", newnewPattern, Fg_pattern);

}
void smoothLaplacian(MatrixXd& Vg, MatrixXi& Fg){
    laplFilter(Vg, Fg, Vg, Fg);
}
void splitAndSmooth(MatrixXd& Vg,MatrixXi& Fg,MatrixXd& Vg_pattern,MatrixXi& Fg_pattern,
                    MatrixXd& VgPatternRet,MatrixXi& FgPatternRet,
                    MatrixXd& VgRet, MatrixXi& FgRet, string garment ){

    insertPlane(Vg, Fg, Vg_pattern, Fg_pattern, garment);
    VectorXd leftFaces = VectorXd::Zero(Fg.rows());
    VectorXd leftVert = VectorXd::Zero(Vg.rows());
    VectorXd leftVert_pattern = VectorXd::Zero(Vg_pattern.rows());

    for(int i=0; i<Fg.rows(); i++){
        bool isLeft= false;
        for (int j=0; j<3; j++){
            if(Vg(Fg(i,j),0) < 0){
                isLeft = true;
            }
        }
        if(isLeft) {
            leftFaces(i) = 1;
            for (int j=0; j<3; j++) {
                leftVert(Fg(i, j)) = 1;
                leftVert_pattern(Fg_pattern(i, j)) = 1;
            }
        }
    }

    VectorXi mapVert(Vg.rows());
    VectorXi mapVert_pattern(Vg_pattern.rows());

    int newVert = leftVert.sum();
    int newVert_pattern = leftVert_pattern.sum();

    MatrixXd newVg (newVert, 3);
    MatrixXd newVg_pattern(newVert_pattern, 3);

    // for now no on seam for the pattern. just copy them all
    VectorXd onSeam = VectorXd::Zero(newVert);
    int count=0; int count_pattern= 0;
    for(int i=0; i<Vg.rows(); i++){
        if(leftVert(i) == 1){
            newVg.row(count) = Vg.row(i);
            mapVert(i) = count;
            count++;
        }
    }
    for(int i=0; i<Vg_pattern.rows(); i++){
        if(leftVert_pattern(i) ==1){
            newVg_pattern.row(count_pattern) = Vg_pattern.row(i);
            mapVert_pattern(i) = count_pattern;
            count_pattern++;
        }
    }


    int newFace = leftFaces.sum();
    MatrixXi newFg (newFace, 3);
    MatrixXi newFg_pattern(newFace, 3);
    count = 0;
    for(int i=0; i<Fg.rows(); i++){
        if(leftFaces(i) == 1){
            for(int j = 0; j < 3; j++){
                newFg(count, j) = mapVert(Fg(i,j));
                newFg_pattern(count, j) = mapVert_pattern(Fg_pattern(i,j));
                if(Vg(Fg(i,j),0) == 0){
                    onSeam( mapVert(Fg(i,j))) = 1;
                }
            }
            count++;
        }
    }
    if(newVg_pattern.rows() == 691){// risky tweak for shirt
        VectorXi componentIdPerVert;
        igl::vertex_components(newFg_pattern, componentIdPerVert);
        for(int i=0; i<newVg_pattern.rows(); i++){
            if(componentIdPerVert(i) == 1){
                newVg_pattern(i,0) -= 50;
                newVg_pattern(i,1) -= 50;
            }
        }
    }

    VgRet.resize(newVg.rows(), 3);
    VgRet= newVg;
    FgRet.resize(newFg.rows(), 3);
    FgRet= newFg;
    igl::writeOBJ("leftGarmentBeforeSmooth.obj", newVg, newFg );
    VgPatternRet.resize(newVg_pattern.rows(), 3);
    VgPatternRet= newVg_pattern;
    FgPatternRet.resize(newFg_pattern.rows(), 3);
    FgPatternRet= newFg_pattern;
    igl::writeOBJ("leftPatternBeforeSmooth.obj", newVg_pattern, newFg_pattern);
}

void preProcessGarment(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, bool insPlane, int symVert1, int symVert2 ,VectorXd& T_sym, string garment){
    MatrixXd newVg, origVg;
    MatrixXd newVg_pattern;
    MatrixXi newFg;
    MatrixXi newFg_pattern;
    VectorXd onSeam;
    if(!insPlane) {
       VectorXd leftFaces = VectorXd::Zero(Fg.rows());
       VectorXd leftVert = VectorXd::Zero(Vg.rows());
       VectorXd leftVert_pattern = VectorXd::Zero(Vg_pattern.rows());

       for (int i = 0; i < Fg.rows(); i++) {
           bool isLeft = false;
           for (int j = 0; j < 3; j++) {
               if (Vg(Fg(i, j), 0) < 0) {
                   isLeft = true;
               }
           }
           if (isLeft) {
               leftFaces(i) = 1;
               for (int j = 0; j < 3; j++) {
                   leftVert(Fg(i, j)) = 1;
                   leftVert_pattern(Fg_pattern(i, j)) = 1;
               }
           }
       }

       VectorXi mapVert(Vg.rows());
       VectorXi mapVert_pattern(Vg_pattern.rows());

       int newVert = leftVert.sum();
       int newVert_pattern = leftVert_pattern.sum();

       newVg.resize(newVert, 3);
       newVg_pattern.resize(newVert_pattern, 3);

       // for now no on seam for the pattern. just copy them all
       onSeam = VectorXd::Zero(newVert);
       int count = 0;
       int count_pattern = 0;
       for (int i = 0; i < Vg.rows(); i++) {
           if (leftVert(i) == 1) {
               newVg.row(count) = Vg.row(i);
               mapVert(i) = count;
               count++;
           }
       }
       for (int i = 0; i < Vg_pattern.rows(); i++) {
           if (leftVert_pattern(i) == 1) {
               newVg_pattern.row(count_pattern) = Vg_pattern.row(i);
               mapVert_pattern(i) = count_pattern;
               count_pattern++;
           }
       }


       int newFace = leftFaces.sum();
       newFg.resize(newFace, 3);
       newFg_pattern.resize(newFace, 3);
       count = 0;
       for (int i = 0; i < Fg.rows(); i++) {
           if (leftFaces(i) == 1) {
               for (int j = 0; j < 3; j++) {
                   newFg(count, j) = mapVert(Fg(i, j));
                   newFg_pattern(count, j) = mapVert_pattern(Fg_pattern(i, j));
                   if (Vg(Fg(i, j), 0) == 0) {
                       onSeam(mapVert(Fg(i, j))) = 1;
                   }
               }
               count++;
           }
       }
       if (newVg_pattern.rows() == 691) {// risky tweak for shirt
           VectorXi componentIdPerVert;
           igl::vertex_components(newFg_pattern, componentIdPerVert);
           for (int i = 0; i < newVg_pattern.rows(); i++) {
               if (componentIdPerVert(i) == 1) {
                   newVg_pattern(i, 0) -= 50;
                   newVg_pattern(i, 1) -= 50;
               }
           }
       }

       igl::writeOBJ("leftGarment.obj", newVg, newFg);

       igl::writeOBJ("leftPattern.obj", newVg_pattern, newFg_pattern);
   }else{
       igl::readOBJ("leftPattern_SmoothBound.obj", newVg_pattern, newFg_pattern);
       igl::readOBJ("leftGarmentBeforeSmooth.obj", newVg, newFg);
       cout<<newFg_pattern.rows()<<" pattern rows and garment rows"<<newFg.rows()<<endl;
        origVg = newVg;
       laplFilter(newVg, newFg , newVg_pattern, newFg_pattern);
        // addapt the garment similarly
        VectorXd S;
        VectorXi I;//face index of smallest distance
        MatrixXd C,N, initGuess;
        MatrixXd pattLeftV;
        MatrixXi pattLeftF;
        igl::readOBJ("leftPatternBeforeSmooth.obj", pattLeftV, pattLeftF);

        igl::signed_distance(newVg_pattern, pattLeftV, pattLeftF, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N);

        MatrixXd B(newVg_pattern.rows(), 3); // contains all barycentric coordinates
        for(int i = 0; i < newVg_pattern.rows(); i++){
            MatrixXd bary;
            auto face = pattLeftF.row(I(i));
            igl::barycentric_coordinates(newVg_pattern.row(i), pattLeftV.row(face(0)), pattLeftV.row(face(1)),
                                         pattLeftV.row(face(2)), bary);
            B.row(i) = bary.row(0);
        }
        initGuess.resize(newVg.rows(), 3);
        for(int i=0; i<newFg.rows(); i++){
            for(int j=0; j<3; j++){

                int vidx = pattLeftF(i,j);
                VectorXi inFace = newFg.row(I(vidx));
                initGuess.row(newFg(i,j)) = newVg.row(inFace( 0))* B(vidx, 0)+
                        newVg.row(inFace( 1))* B(vidx, 1)+ newVg.row(inFace( 2))* B(vidx, 2);

            }
        }

        writeOBJ("interpolLaplGar.obj", initGuess, newFg);
        newVg= initGuess;

   }

    // finished the split, now duplicate to make it one again

    //first duplicate all vertices and faces, then remove the ones that are on the seam

    VectorXi mapDupl(newVg.rows());
    int count = 0;
    VectorXd onSeam2 = VectorXd::Zero(newVg.rows());
    if(!insPlane){
        for(int i =0; i<newVg.rows(); i++){
            if( onSeam(i) == 0){
                mapDupl(i) = count;
                count++;
                onSeam2(i)= onSeam(i);
            }
        }
    }else{
        VectorXi checked= VectorXi::Zero(newVg.rows());

        for(int i =0; i< newFg.rows(); i++){
            for(int j =0; j<3; j++){

                if(checked(newFg(i,j))!=0)continue;
                if(newFg(i,j)==320)cout<<origVg.row(newFg(i,j))<<" the orig row"<<endl;
                 if(abs( origVg(newFg(i,j),0)) >= 0.05){
                    mapDupl(newFg(i,j)) = count;
                    count++;
                    onSeam2(newFg(i,j))= 0;
                     if(newFg(i,j)==320)cout<<onSeam2(newFg(i,j))<<" on seam? and map "<<  mapDupl(newFg(i,j))<<endl;
                 }else{// if( newVg_pattern(newFg_pattern(i,j),0) == 0){
                    onSeam2(newFg(i,j))= 1;
                }
                 checked(newFg(i,j))++;
            }
        }
    }

    MatrixXd VgDupl (count ,3);
    count = 0;
    if(!insPlane){
    for(int i=0; i<newVg.rows(); i++){
        if( onSeam2(i) == 0){
            VgDupl.row(count)= newVg.row(i);
            count++;
        }
    }
    }else{
        VectorXi checked= VectorXi::Zero(newVg.rows());
        for(int i =0; i< newFg.rows(); i++){
            for(int j =0; j<3; j++){
                if(checked(newFg(i,j))!=0)continue;
                if( onSeam2(newFg(i,j))== 0){
                    VgDupl.row(count)= newVg.row(newFg(i,j));
                    count++;
                }
                checked(newFg(i,j))++;
            }
        }
    }

    MatrixXd rot = MatrixXd::Identity(3,3); rot(0,0)= -1; // reflection

    VgDupl = (rot * VgDupl.transpose()).transpose();
    MatrixXi FgDupl( newFg.rows(), 3);
    int offset = newVg.rows();

    for(int i=0; i<newFg.rows(); i++){
        for(int j=0; j<3; j++){
            // if the vertex is on the seam, we take the same one!
            if( onSeam2(newFg(i,j )) != 0){
                FgDupl(i,j) = newFg(i,j);
            }else{
                // else we take the duplicated vertex
                FgDupl(i,j) = mapDupl(newFg(i,j)) + offset;
            }
            if(newFg(i,j)==320)cout<<FgDupl(i,j)<<" the orig row"<< mapDupl(newFg(i,j)) <<endl;
        }
    }
    MatrixXd fullVg(VgDupl.rows() + newVg.rows(), 3);
    fullVg<<newVg, VgDupl;
    // change the normal to align!
    VectorXi temp = FgDupl.col(1);
    FgDupl.col(1) = FgDupl.col(2);
    FgDupl.col(2) = temp;

    igl::writeOBJ("rightGarment.obj", fullVg, FgDupl);

    MatrixXi fullFg (FgDupl.rows()*2, 3);
    fullFg<<newFg, FgDupl;

    igl::writeOBJ("fullGarment.obj", fullVg, fullFg);

    MatrixXi FgDupl_pattern =  newFg_pattern;
    VectorXi tempF = FgDupl_pattern.col(1);
    FgDupl_pattern.col(1) = FgDupl_pattern.col(2);
    FgDupl_pattern.col(2) = tempF;
    MatrixXd VgDupl_pattern (newVg_pattern.rows(), 3);
    VgDupl_pattern = (rot*newVg_pattern.transpose()).transpose();

    T_sym = Vg_pattern.row(symVert1 ) - VgDupl_pattern.row(symVert2);
    VgDupl_pattern.rowwise() += T_sym.transpose();
    offset = newVg_pattern.rows();
    MatrixXi offsetM(FgDupl_pattern.rows(), FgDupl_pattern.cols()); offsetM.setConstant(offset);
    FgDupl_pattern+= offsetM;

    MatrixXi fullFg_pattern (FgDupl_pattern.rows() * 2, 3);
    fullFg_pattern << newFg_pattern, FgDupl_pattern;
    MatrixXd fullVg_pattern(newVg_pattern.rows() * 2, 3);
    fullVg_pattern << newVg_pattern, VgDupl_pattern;
    igl::writeOBJ("fullPattern.obj", fullVg_pattern, fullFg_pattern);

    Vg.resize(fullVg.rows(), 3);  Vg = fullVg;
    Fg.resize(fullFg.rows(), 3);  Fg = fullFg;
    Vg_pattern.resize(fullVg_pattern.rows(), 3);  Vg_pattern = fullVg_pattern;
    Fg_pattern.resize(fullFg_pattern.rows(), 3);  Fg_pattern = fullFg_pattern;


//


}