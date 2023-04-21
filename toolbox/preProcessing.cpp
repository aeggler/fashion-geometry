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
#include <igl/facet_components.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/boundary_loop.h>

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
                             map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& insertedIdxToPatternVert, VectorXi& isLeftVertPattern,MatrixXd& rightVert, string garment){
  /* filter the vertices on the left side and create maps form full pattern to half (if in wrong side it is not in the map),
   * and vice versa. */
   int n = Vg.rows();
   int m  = Fg.rows();
    VectorXi isLeftVert(n);
    VectorXi isRightVert(n);

    isLeftVert.setConstant(-1);
    isRightVert.setConstant(-1);
    double th = 0.2;
    if(garment == "skirt"){
        th = 0.8;
    }
    int leftCount = 0; int rightCount = 0 ;
    for(int i=0; i<n; i++){
        if(Vg(i, 0) <= th){
            isLeftVert(i) = 1;
            leftCount++;
        }
        if(Vg(i, 0) >= -th){
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
//        }else if (Vg(v0, 0)<0 ||Vg(v1, 0)<0 ||Vg(v2, 0)<0  ){
//            faceCount++;
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

    int patternHalfVert = isLeftVertPattern.sum();
    cout<<patternHalfVert<<" half vertices"<<endl;
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


igl::writeOBJ("halfPattern.obj", Vg_pattern_half,  Fg_pattern_half);

}
void insertPlane(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, string garment, string garmentEXT){
    /* we insert a plane along the x=0 axis, add vertices and introduce new faces
       this ensures we can split the garment safely in the pre processing
       we do the same for the pattern and duplicate x=0 vertices -> to make sure the patch is acutally disconnected and we can take only the half patch for symetry
    */
     map<std::pair<double, double>, std::pair<int, int>> yToFaceAndIdx;
    int vgsize = Vg_pattern.rows();
    int vgGarOrig = Vg.rows();
    std::vector< std::vector<int> > vfAdj,vfAdjG, vvAdj;
    createVertexFaceAdjacencyList(Fg, vfAdjG);

    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
    igl::adjacency_list(Fg_pattern, vvAdj);
    int addCount= 0 ;
    double precision = 1000000.;
    for(int i=0; i<Vg.rows(); i++){
        if(Vg(i ,0) == 0){
            addCount++;
        }
    }
    for(int i=0; i<Fg.rows(); i++){

        bool hasLeft = false;
        bool hasRight = false;
        Vector3i LR;
        for (int j=0; j<3; j++){
            double x = Vg(Fg(i, j), 0);
            if(x==0) {

                double el1 = std::floor( Vg(Fg(i, j), 1))*precision/ precision;
                double el2= std::floor( Vg(Fg(i, j), 2))*precision/ precision;
                pair<double, double> pair1 = make_pair(el1, el2);
                yToFaceAndIdx[ pair1] = std::make_pair(i,j );
//                addCount++;
//
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
            //todo
            VectorXd newPos1 = Vg.row(v1) + t2 * edge2.transpose();
            newPos1(1) = std::floor(newPos1(1) * precision)/ precision;
            newPos1(2) = std::floor(newPos1(2) * precision)/ precision;
            pair<double, double> pair1 = std::make_pair(newPos1(1), newPos1(2));
            newPos1(0)= 0;
            VectorXd newPos2 = Vg.row(v1) + t3 * edge3.transpose();
            newPos2(0) = 0;
            newPos2(1) = std::floor(newPos2(1) * precision)/ precision;
            newPos2(2) = std::floor(newPos2(2) * precision)/ precision;
            pair<double, double> pair2 = std::make_pair(newPos2(1), newPos2(2));


            int idx1 , idx2; int fac1, idfac1, fac2, idfac2;
            bool new1= false; bool new2 = false;
            bool extra1= false; bool extra2 = false;

            if(yToFaceAndIdx.find(pair1) == yToFaceAndIdx.end()){

                yToFaceAndIdx[pair1] = std::make_pair(i, (otherSide+1) % 3 );

                int vgrow =  Vg.rows();
                MatrixXd Vgnew( vgrow + 1, 3);
                Vgnew.block(0,0,vgrow, 3) = Vg;
                Vgnew.row(vgrow ) = newPos1;
                Vg.resize(vgrow+1, 3);
                Vg = Vgnew;
                idx1 = vgrow; new1= true;
            }else{
                fac1 = yToFaceAndIdx[pair1].first;
                idfac1 = yToFaceAndIdx[pair1].second;

                idx1 =Fg (fac1, idfac1 );
                if(i== 1128 && garment== "skirt_no2"){ extra1 = true;}
                if(i==102 && garment == "skirt") {
                    extra1 = true;
                }
                if(i==5596 && garmentEXT == "skirt_3") {
                    extra1 = true;
                }
                if(i == 6096 && garment == "hoodie"){
                    extra1 = true; //good
                }
                if(i==6427 && garment == "man_tshirt"){
                    extra1 = true;
                }

            }
            if(yToFaceAndIdx.find(pair2) == yToFaceAndIdx.end()){

                yToFaceAndIdx[pair2] = std::make_pair(i, (otherSide+2) % 3 );

                int vgrow =  Vg.rows();
                MatrixXd Vgnew( vgrow + 1, 3);
                Vgnew.block(0,0,vgrow, 3) = Vg;
                Vgnew.row(vgrow ) = newPos2;
                Vg.resize(vgrow+1, 3);
                Vg = Vgnew;
                idx2 = vgrow;
                new2= true;
            }else{
                fac2 = yToFaceAndIdx[pair2].first;
                idfac2 = yToFaceAndIdx[pair2].second;
                idx2 = Fg ( fac2, idfac2 );
                if(i== 1256 && garment== "skirt_no2"){ extra2 = true;}
                if(i==849 && garment == "skirt") {
                    extra2 = true;
                }
                if(i==5534 && garmentEXT == "skirt_3") {
                    extra2 = true;
                }
                if(i == 6433 && garment == "hoodie"){
                    extra2 = true;// good!
                }
                if(i == 6456 && garment == "hoodie"){
                    extra2 = true;// good
                }
                if(i == 2089 && garment == "hoodie"){
                    extra2 = true;//good
                }
                if(i==3872 && garment == "man_tshirt"){
                    extra2 = true;
                }

            }
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

            Fgnew(fgrow + 1, 0) = idx1;
            Fgnew(fgrow + 1, 1) = Fg_pattern(i, (otherSide + 2) % 3 );
            Fgnew(fgrow + 1, 2) = idx2;

            Fg_pattern.resize(Fgnew.rows(), 3);
            Fg_pattern= Fgnew;

        }
    }

    igl::writeOBJ("dress_3d.obj", Vg, Fg);
    Vg_pattern.col(2).setConstant(200);
    igl::writeOBJ("dress_2d.obj", Vg_pattern, Fg_pattern);
    int addedVert = Vg_pattern.rows()- vgsize;
    cout<<addCount<<" zero vertices and added vert "<<addedVert<<endl;
    int added = addedVert+ addCount;

    int vgnewsize = Vg_pattern.rows();
    MatrixXd Vgp (Vg_pattern.rows()+ added, 3);
    Vgp.block(0,0,vgnewsize, 3) = Vg_pattern;
    Vgp.block(vgnewsize, 0, added, 3) = Vg_pattern.block(vgsize, 0, addedVert, 3);
    int duplCount = 0;
    VectorXi addedDeja(Vg.rows()); addedDeja.setConstant(0);
    VectorXi duplVert(Vg_pattern.rows());
    for(int ll=0 ; ll<Fg.rows(); ll++){
        for(int lll=0; lll<3; lll++){
            int l = Fg(ll,lll);

            if(Vg(l, 0)==0 && addedDeja(l)==0 && l<vgGarOrig){
                addedDeja(l)++;
                cout<<"adding "<<l<<endl;
                Vgp.row(vgnewsize + addedVert+ duplCount)= Vg_pattern.row(Fg_pattern(ll,lll));
                duplVert(Fg_pattern(ll,lll))= vgnewsize + addedVert+ duplCount;
                duplCount++;
            }

        }
    }
    cout<<"continue"<<endl;

    // we duplicate the new vertices to split them
    duplCount=0;
    for (int i=0; i<Fg_pattern.rows(); i++){
        for(int j=0; j<3; j++){
            if(Fg_pattern(i,j)>= vgsize){

                int other = (j+1)%3;
                int co=0;

                while(Vg(Fg(i, other), 0) == 0 && co<6){
                    co++;
                    other++;
                    other %= 3;
                }// find one that is not 0
                if(co>=3) {cout<<i<<" face, vertex issues  "<<Vg(Fg(i, other), 0)<<" "<<Vg(Fg(i, (other+1)%3 ), 0)<<" "<<
                    Vg(Fg(i, (other+2)%3 ), 0)<<endl;
                    cout<<Fg(i, other)<<" "<<Fg(i, (other+1)%3 ) <<" "<<Fg(i, (other+2)%3 )<<endl;
                    cout<<Vg.row(Fg(i, (other+1)%3 ))<<endl;
                }
                bool isLeft = false;
                if(Vg(Fg(i, other), 0)<0){
                    isLeft= true;
                }

                if(isLeft){
                    Fg_pattern(i,j) +=added;
                }

            }
            else if (Vg(Fg(i,j),0)==0){
                if(Vg(Fg(i,(j+1)%3 ),0)<0|| Vg(Fg(i, (j+2)%3 ),0)<0){
                    Fg_pattern(i,j) =duplVert(Fg_pattern(i,j));
                }
            }
        }
    }
    igl::writeOBJ("dress_2d_dupl.obj", Vgp, Fg_pattern);
    Vg_pattern.resize(Vgp.rows(), Vgp.cols()); Vg_pattern = Vgp;

}
void laplFilter(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern){
    /* Vanilla laplacian smoothing to get better vertex positions.
     * */
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
bool isCorner(int id, MatrixXi& Fg,set<int>& cornersOfGar,  MatrixXi& Fg_pattern, MatrixXd& Vg, string garment){
    if(cornersOfGar.find(id) != cornersOfGar.end()){
        return true;
    }return false;
}
void edgeCollapse(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, string garment, string garmentEXT) {
    set<int> cornersOfGar;
    set<int> freecorners;// corners in the outer of the gament in3D.not sure about them
    double theresh = 7.;

    if (garment == "skirt_no2") {
        freecorners.insert(704);
        freecorners.insert(707);
        freecorners.insert(670);

        cornersOfGar.insert(704);
        cornersOfGar.insert(698);
        cornersOfGar.insert(707);
        cornersOfGar.insert(702);
        cornersOfGar.insert(0);
        cornersOfGar.insert(23);
        cornersOfGar.insert(679);

        cornersOfGar.insert(670);
        cornersOfGar.insert(708);
        cornersOfGar.insert(711);

    }
    else if(garment=="skirt" && garmentEXT !="skirt_3"){
        theresh= 9;
        cornersOfGar.insert(404);
        cornersOfGar.insert(400);
        cornersOfGar.insert(405);
        cornersOfGar.insert(406);

        cornersOfGar.insert(410);
        cornersOfGar.insert(409);
        cornersOfGar.insert(456);
        cornersOfGar.insert(448);

        freecorners.insert(404);
        freecorners.insert(406);
        freecorners.insert(456);
        freecorners.insert(409);
    }
    else if (garmentEXT == "skirt_3"){
        cornersOfGar.insert(1617);
        cornersOfGar.insert(733);
        cornersOfGar.insert(756);
        cornersOfGar.insert(1610);

        cornersOfGar.insert(43);
        cornersOfGar.insert(1528);
        cornersOfGar.insert(1507);
        cornersOfGar.insert(8);

        cornersOfGar.insert(1713);
        cornersOfGar.insert(1709);
        cornersOfGar.insert(1712);
        cornersOfGar.insert(1707);

//        freecorners.insert(1713);
//        freecorners.insert(1617);
//        freecorners.insert(1507);
//        freecorners.insert(1707);
    }
    else if (garment == "hoodie"){
        cornersOfGar.insert(2221);
        cornersOfGar.insert(2226);
        cornersOfGar.insert(1617);
        cornersOfGar.insert(1624);
        cornersOfGar.insert(1629);

        cornersOfGar.insert(2205);
        cornersOfGar.insert(1514);
        cornersOfGar.insert(1520);
        cornersOfGar.insert(1522);
        cornersOfGar.insert(2202);

        cornersOfGar.insert(2232);
        cornersOfGar.insert(1722);
        cornersOfGar.insert(3124);
        cornersOfGar.insert(2236);

        cornersOfGar.insert(11);
        cornersOfGar.insert(14);
        cornersOfGar.insert(2106);
        cornersOfGar.insert(2105);

        cornersOfGar.insert(1440);
        cornersOfGar.insert(1458);
        cornersOfGar.insert(1455);
        cornersOfGar.insert(1443);

        cornersOfGar.insert(723);
        cornersOfGar.insert(790);
        cornersOfGar.insert(768);
        cornersOfGar.insert(745);

        cornersOfGar.insert(2124);
        cornersOfGar.insert(155);
        cornersOfGar.insert(148);
        cornersOfGar.insert(136);
        cornersOfGar.insert(115);
        cornersOfGar.insert(2123);

        cornersOfGar.insert(63);
        cornersOfGar.insert(60);
        cornersOfGar.insert(2112);
        cornersOfGar.insert(2113);

    }
    else if (garment == "man_tshirt"){
        cornersOfGar.insert(1583);
        cornersOfGar.insert(1588);
        cornersOfGar.insert(1613);
        cornersOfGar.insert(1618);

        cornersOfGar.insert(1835);
        cornersOfGar.insert(64);
        cornersOfGar.insert(40);
        cornersOfGar.insert(28);
        cornersOfGar.insert(19);
        cornersOfGar.insert(1840);

        cornersOfGar.insert(2017);
        cornersOfGar.insert(2019);
        cornersOfGar.insert(1566);
        cornersOfGar.insert(1582);

        cornersOfGar.insert(1920);
        cornersOfGar.insert(798);
        cornersOfGar.insert(774);
        cornersOfGar.insert(761);
        cornersOfGar.insert(752);
        cornersOfGar.insert(1918);

        cornersOfGar.insert(1831);
        cornersOfGar.insert(1833);
        cornersOfGar.insert(6);
        cornersOfGar.insert(5);
        for(auto it: cornersOfGar){
            if(it!= 64 && it != 798 )
            freecorners.insert(it);
        }

    }
    int count = 0;
    int countItems = 0 ;

    vector<vector<int>> boundaryLoop;
    igl::boundary_loop(Fg_pattern, boundaryLoop);
    for(int iterations = 0; iterations < 5 ; iterations++){
    for (int i = 0; i < boundaryLoop.size(); i++) {
        for (int l = 0; l < boundaryLoop[i].size(); l++) {
//            freecorners.insert (-4);// should it always have at least one element?
            vector<vector<int>> vvAdj, vvAdjGar, vfAdj, vfAdjGar;
            createVertexFaceAdjacencyList(Fg_pattern, vfAdj);
            createVertexFaceAdjacencyList(Fg, vfAdjGar);
            igl::adjacency_list(Fg_pattern, vvAdj);
            igl::adjacency_list(Fg, vvAdjGar);

            int id0 = boundaryLoop[i][l];
            int id1 = boundaryLoop[i][(l + 1) % boundaryLoop[i].size()];
//            int face = adjacentFaceToEdge(id0, id1, -1, vfAdj);
            bool corner0 = isCorner(id0, Fg, cornersOfGar, Fg_pattern, Vg, garment);
            bool corner1 = isCorner(id1, Fg, cornersOfGar, Fg_pattern, Vg, garment);

            int v0Bound = isBoundaryVertex(Vg_pattern, id0, vvAdj, vfAdj);
            int v1Bound = isBoundaryVertex(Vg_pattern, id1, vvAdj, vfAdj);
            if (!v0Bound || !v1Bound)continue; // for now just look at two adjacent boundary vertices!
            double dist = (Vg_pattern.row(id0) - Vg_pattern.row(id1)).norm();

            if (dist < theresh) {

                int face = adjacentFaceToEdge(id0, id1, -1, vfAdj);

                int idx0 = 0;
                while (Fg_pattern(face, idx0) != id0) { idx0++; }
                int id0G = Fg(face, idx0);
                int idx1 = 0;
                while (Fg_pattern(face, idx1) != id1) { idx1++; }
                int id1G = Fg(face, idx1);

                RowVectorXd newpos, newposGar;
                if (corner0) {
                    if(freecorners.find(id0)==freecorners.end()){
                        continue;
                    }

                    newpos = Vg_pattern.row(id0);
                    newposGar = Vg.row(id0G);
                } else if (corner1) {
                    if(freecorners.find(id1)==freecorners.end()){
                        continue;
                    }
                    newpos = Vg_pattern.row(id1);
                    newposGar = Vg.row(id1G);

                } else {

                    newpos = (Vg_pattern.row(id0) + Vg_pattern.row(id1)) / 2;
                    newposGar = (Vg.row(id0G) + Vg.row(id1G)) / 2;

                }
                Vg_pattern.row(id0) = newpos;
                Vg.row(id0G) = newposGar;

                for (int j = 0; j < vfAdj[id1].size(); j++) {
                    int ii = 0;
                    while (Fg_pattern(vfAdj[id1][j], ii) != id1) {
                        ii++;
                    }
                    Fg_pattern(vfAdj[id1][j], ii) = id0;
                }

                for (int j = 0; j < vfAdjGar[id1G].size(); j++) {
                    int ii = 0;
                    while (Fg(vfAdjGar[id1G][j], ii) != id1G) {
                        ii++;
                    }
                    Fg(vfAdjGar[id1G][j], ii) = id0G;
                }

                MatrixXd VgNew(Vg_pattern.rows() - 1, 3);
                VgNew.block(0, 0, id1, 3) = Vg_pattern.block(0, 0, id1, 3);
                VgNew.block(id1, 0, VgNew.rows() - id1, 3) = Vg_pattern.block(id1 + 1, 0, VgNew.rows() - id1, 3);
                MatrixXd VgGarNew(Vg.rows() - 1, 3);
                VgGarNew.block(0, 0, id1G, 3) = Vg.block(0, 0, id1G, 3);
                VgGarNew.block(id1G, 0, VgGarNew.rows() - id1G, 3) = Vg.block(id1G + 1, 0, VgGarNew.rows() - id1G, 3);

                for (int j = 0; j < Fg_pattern.rows(); j++) {
                    for (int ii = 0; ii < 3; ii++) {
                        if (Fg_pattern(j, ii) >= id1) {
                            Fg_pattern(j, ii)--;
                        }
                    }
                }

                for (int j = 0; j < Fg.rows(); j++) {
                    for (int ii = 0; ii < 3; ii++) {
                        if (Fg(j, ii) >= id1G) {
                            Fg(j, ii)--;
                        }
                    }
                }
                set<int> cornersNew, freecornersNew;
                freecornersNew.clear();
                for (int it: cornersOfGar) {
                    if (it > id1) {
                        cornersNew.insert(it - 1);
                        if(freecorners.find(it)!= freecorners.end()){
                            freecornersNew.insert(it-1);
                        }

                    } else if (it == id1) {
                        cornersNew.insert(id0);
                        if(freecorners.find(it)!= freecorners.end()){
                            freecornersNew.insert(id0);
                        }
                    } else {
                        cornersNew.insert(it);
                        if(freecorners.find(it)!= freecorners.end()){
                            freecornersNew.insert(it);
                        }
                    }
                }
                cornersOfGar.clear();

                freecorners.clear();
                cornersOfGar = cornersNew;
                freecorners = freecornersNew;
                int ps = Fg_pattern.rows() - 1;
                MatrixXi FgNew_pat(ps, 3);
                FgNew_pat.block(0, 0, face, 3) = Fg_pattern.block(0, 0, face, 3);
                FgNew_pat.block(face, 0, FgNew_pat.rows() - face, 3) = Fg_pattern.block(face + 1, 0,
                                                                                        FgNew_pat.rows() - face, 3);

                MatrixXi FgGarNew(Fg.rows() - 1, 3);
                FgGarNew.block(0, 0, face, 3) = Fg.block(0, 0, face, 3);
                FgGarNew.block(face, 0, FgGarNew.rows() - face, 3) = Fg.block(face + 1, 0, Fg.rows() - face, 3);

                Vg_pattern.resize(VgNew.rows(), 3);
                Vg_pattern = VgNew;
                Fg_pattern.resize(FgNew_pat.rows(), 3);
                Fg_pattern = FgNew_pat;
//                igl::writeOBJ("leftPatternBeforeSmooth.obj", Vg_pattern, Fg_pattern);

                Vg.resize(VgGarNew.rows(), 3);
                Vg = VgGarNew;
                Fg.resize(FgGarNew.rows(), 3);
                Fg = FgGarNew;
//                igl::writeOBJ("leftGarmentBeforeSmooth.obj", Vg, Fg);
                vector<vector<int>> boundaryLoopNew;
                igl::boundary_loop(Fg_pattern, boundaryLoopNew);
                boundaryLoop.clear();
                boundaryLoop = boundaryLoopNew;
                count ++;

                countItems++;
            }
        }

    }
}
    cout<<countItems<<" items collapsed "<<endl;
    igl::writeOBJ("collapsed.obj", Vg_pattern, Fg_pattern);
    igl::writeOBJ("collapsed3d.obj", Vg, Fg);


}
void splitAndSmooth(MatrixXd& Vg,MatrixXi& Fg,MatrixXd& Vg_pattern,MatrixXi& Fg_pattern,
                    MatrixXd& VgPatternRet,MatrixXi& FgPatternRet,
                    MatrixXd& VgRet, MatrixXi& FgRet, string garment, string garmentEXT ){

    if(garment=="skirt"){
        VectorXi vertComp;
        igl::vertex_components( Fg_pattern, vertComp);
        for(int i=0; i< Vg_pattern.rows(); i++){
            if(vertComp(i)==0){
                Vg_pattern(i, 0) -=300;
            }
        }
    }
    if(garmentEXT =="skirt_3"){
        VectorXi vertComp;
        igl::vertex_components( Fg_pattern, vertComp);
        for(int i=0; i< Vg_pattern.rows(); i++){
            if(vertComp(i)==0){
                Vg_pattern(i, 0) +=300;
            }
        }
    }
    if(garment == "tshirt"){
        for(int i=0; i< Vg.rows(); i++){
            if(i<2)cout<<Vg(i,0)<<" for i "<<i<<endl;
            if(abs(Vg(i, 0)) <= 1.1){
                cout<<i<<endl;
                Vg(i, 0) = 0;
            }
        }
    }

    insertPlane(Vg, Fg, Vg_pattern, Fg_pattern, garment, garmentEXT);

    VectorXd leftFaces = VectorXd::Zero(Fg.rows());
    VectorXd leftVert = VectorXd::Zero(Vg.rows());
    VectorXd leftVert_pattern = VectorXd::Zero(Vg_pattern.rows());
    cout<<"finished insertion"<<endl;
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

    if(garment =="skirt"&& garmentEXT != "skirt_3"){
        VectorXi componentIdPerVert;
        igl::vertex_components(newFg_pattern, componentIdPerVert);
        for(int i=0; i<newVg_pattern.rows(); i++){
            if(componentIdPerVert(i) == 1){
                newVg_pattern(i,0) -= 50;
                newVg_pattern(i,1) -= 50;
            }
        }
    }
    igl::writeOBJ("leftPatternBeforecoll.obj", newVg_pattern, newFg_pattern);
    cout<<"starting edge collapse"<<endl;
    edgeCollapse (newVg, newFg, newVg_pattern, newFg_pattern, garment, garmentEXT);

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

void preProcessGarment(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, bool insPlane, int symVert1, int symVert2 ,VectorXd& T_sym, string garment, string garmentEXT){
    MatrixXd newVg, origVg;
    MatrixXd newVg_pattern;
    MatrixXi newFg;
    MatrixXi newFg_pattern;
    VectorXd onSeam;
    if(!insPlane) {
        // if there is no plane ot be inserted, e.g. if the seams are on the plane already and do not need to be split
        // just take all vertices on the left side in 3D, and the same vertices via *facial corresponance* from the pattern
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
                   if (abs(Vg(Fg(i, j), 0)) <= 0.4) {
                       onSeam(mapVert(Fg(i, j))) = 1;
//                       cout<<mapVert(Fg(i, j))<<" ";
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
       if(garment == "leggins"){
           VectorXi comp;
           igl::vertex_components(newFg_pattern, comp);
           for(int i=0; i<newVg_pattern.rows(); i++){
               if(comp(i)== 1 || comp(i)==3){
                   newVg_pattern(i, 0) += 100;
               }
           }
       }

       igl::writeOBJ("leftPattern.obj", newVg_pattern, newFg_pattern);
   }else{
        // if we have to insert a plane we did this in the previous step and smoothed the boundary of the cut. Now we apply lapalcian
        // smoothing to ensure the vvertices have a good positon
       igl::readOBJ("leftPattern_SmoothBound.obj", newVg_pattern, newFg_pattern);
       igl::readOBJ("leftGarmentBeforeSmooth.obj", newVg, newFg);
       cout<<newFg_pattern.rows()<<" pattern rows and garment rows"<<newFg.rows()<<endl;
        origVg = newVg;
       laplFilter(newVg, newFg , newVg_pattern, newFg_pattern);

        /* after changing the pattern by smoothing we have to ensure the garment is adapted similarly to
         * get good target stretch. We apply deformation transfer
         * */
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
        newVg.resize(initGuess.rows(), 3);
        newVg= initGuess;

   }

    /* finished the split, now duplicate to make it one again
   first duplicate all vertices and faces, then remove the ones that are on the seam.
     we iterate first to count how many vertices we need, then we create the new matrix

     If we inserted the plane we have to iterate over the faces instead of the vertices, because we have to find out in the garment which
     face are on the symmetry plane!
    */
    VectorXi mapDupl(newVg.rows());
    int count = 0;
    VectorXd onSeam2 = VectorXd::Zero(newVg.rows());
    if(!insPlane){
        onSeam2= onSeam;
        for(int i =0; i<newVg.rows(); i++){
            if( onSeam(i) == 0){
                mapDupl(i) = count;
                count++;
//                onSeam2(i)= onSeam(i);
            }
        }
    }else{
        VectorXi checked= VectorXi::Zero(newVg.rows());

        for(int i =0; i< newFg.rows(); i++){
            for(int j =0; j<3; j++){

                if(checked(newFg(i,j))!=0)continue;
                 if(abs( origVg(newFg(i,j),0)) >= 0.0005){
                    mapDupl(newFg(i,j)) = count;
                    count++;
                    onSeam2(newFg(i,j))= 0;
//                     if(newFg(i,j)==618)cout<<onSeam2(newFg(i,j))<<" on seam? and map "<<  mapDupl(newFg(i,j))<<endl;
                 }else{// if( newVg_pattern(newFg_pattern(i,j),0) == 0){
                    onSeam2(newFg(i,j))= 1;

//                     if(newFg(i,j)==618){
////                         cout<<onSeam2(newFg(i,j))<<" on seam? from face "<< i<<endl;
////                         cout<<origVg.row(618)<<endl;
//                     }

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
    /* We duplicate and rotate the garment , to get a full one. */
    VgDupl = (rot * VgDupl.transpose()).transpose();
    MatrixXi FgDupl( newFg.rows(), 3);
    int offset = newVg.rows();

    for(int i=0; i<newFg.rows(); i++){
        for(int j=0; j<3; j++){
            // if the vertex is on the seam, we take the same one!
            if( onSeam2(newFg(i,j )) != 0){

                FgDupl(i,j) = newFg(i,j);
//                if(i==1239){
//                    cout<<FgDupl(i,j)<<" from on seam "<<endl;
//                }
            }else{
                // else we take the duplicated vertex
                FgDupl(i,j) = mapDupl(newFg(i,j)) + offset;
//                if(FgDupl(i,j)<0){
//                    cout<<FgDupl(i,j)<<" map Dupl"<<endl;
//                }
            }
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
    /* the translation between the initial patten (wihtout symmetry) and the duplicated (without offset)
    * Pick the reference vertices wisely to ensure no overlaps in the pattern computation */
    if(garment =="skirt" && garmentEXT != "skirt_3"){
        T_sym.resize(3);
        T_sym(0)=-100;
        T_sym(1)=0;
        T_sym(2)=0;
    }else{
        T_sym = Vg_pattern.row(symVert1 ) - VgDupl_pattern.row(symVert2);
        if(garment =="top"){
            T_sym(0)+=50;
        }else if (garmentEXT == "skirt_3"){
            T_sym(0) = 300;
            T_sym(1) = 0;
            T_sym(2) = 0;
        }else if (garment =="hoodie"){
            T_sym (0) = 200;
        }

    }
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

}