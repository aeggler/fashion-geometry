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
                             map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& insertedIdxToPatternVert ){
   int n = Vg.rows();
   int m  = Fg.rows();
    VectorXi isLeftVert(Vg.rows());
    isLeftVert.setConstant(-1);
    int leftCount = 0;
    for(int i=0; i<n; i++){
        if(Vg(i, 0) <= 0){
            isLeftVert(i) = 1;
            leftCount++;
        }
    }
//    cout<<"Left count "<<leftCount<<"/ "<<n<<endl;

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
    VectorXi isLeftVertPattern(Vg_pattern.rows());isLeftVertPattern.setConstant(-1);
    VectorXi isLeftFace(m);
    isLeftFace.setConstant(0);
    int faceCount = 0;
    for(int i=0; i<m; i++) {
        int v0 = Fg(i, 0);
        int v1 = Fg(i, 1);
        int v2 = Fg(i, 2);
        isLeftVertPattern(Fg_pattern(i, 0)) = (isLeftVert(v0) == 1) ? 1 : 0;
        isLeftVertPattern(Fg_pattern(i, 1)) = (isLeftVert(v1) == 1) ? 1 : 0;
        isLeftVertPattern(Fg_pattern(i, 2)) = (isLeftVert(v2) == 1) ? 1 : 0;

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
        }
    }
    idx=0;
    int patternHalfVert = isLeftVertPattern.sum();
    Vg_pattern_half.resize(patternHalfVert, 3);
    for(int i=0; i< Vg_pattern.rows(); i++){
        if(isLeftVertPattern(i)){
            Vg_pattern_half.row(idx) = Vg_pattern.row(i);
            halfPatternVertToFullPatternVert[idx] = i;
            fullPatternVertToHalfPatternVert[i] = idx;
            idx++;

        }
    }



//    cout<<faceCount<<"/"<< m <<" left face count and middle face count "<<endl;
    vector<VectorXd> addedVert;
    int newIdx = leftCount;
    Fg_pattern_half.resize(faceCount, 3);
//    cout<<Fg_pattern_half.rows()<<" half pattern faces "<<endl ;


    idx = 0;
    for(int i = 0; i<m; i++){
        if(isLeftFace(i) == 1){
            Fg_pattern_half(idx, 0)= fullPatternVertToHalfPatternVert[Fg_pattern(i, 0)];
            Fg_pattern_half(idx, 1)= fullPatternVertToHalfPatternVert[Fg_pattern(i, 1)];
            Fg_pattern_half(idx, 2)= fullPatternVertToHalfPatternVert[Fg_pattern(i, 2)];
            idx++;

        }
//        else if(isLeftFace(i)==0){
//            int v0 = Fg(i, 0);
//            int v1 = Fg(i, 1);
//            int v2 = Fg(i, 2);
//            cout<<i<<" intermediate face "<<isLeftVert(v0)<<" "<<isLeftVert(v1)<<" "<<isLeftVert(v2)<<endl;
//
//            Fg_pattern_half(idx, 0) = fullPatternVertToHalfPatternVert[v0];
//            Fg_pattern_half(idx, 1) = fullPatternVertToHalfPatternVert[v1];
//            Fg_pattern_half(idx, 2) = fullPatternVertToHalfPatternVert[v2];
//
//            // at least one face is on the right side , it has to be replaced
//            if(isLeftVert(v0) == -1){
//                VectorXd vertPos = Vg_pattern.row(v0);
//                vertPos(0)=0;
//                // replace it and change the face
//                addedVert.push_back(vertPos);
//                Fg_pattern_half(idx, 0) = newIdx;
//                insertedIdxToPatternVert[newIdx]= v0;
//                newIdx++;
//            }
//            if(isLeftVert(v1) == -1){
//                VectorXd vertPos = Vg_pattern.row(v1);
//                vertPos(0)=0;
//                // replace it and change the face
//                addedVert.push_back(vertPos);
//                Fg_pattern_half(idx, 1) = newIdx;
//                insertedIdxToPatternVert[newIdx]= v1;
//                newIdx++;
//            }
//            if(isLeftVert(v2) == -1){
//                VectorXd vertPos = Vg_pattern.row(v2);
//                vertPos(0)=0;
//                // replace it and change the face
//                addedVert.push_back(vertPos);
//                Fg_pattern_half(idx, 2) = newIdx;
//                insertedIdxToPatternVert[newIdx]= v2;
//                newIdx++;
//            }
//
//            idx++;
//
//        }
    }

//    MatrixXd Vg_temp = Vg_pattern_half;
//    Vg_pattern_half.resize(Vg_pattern_half.rows() + addedVert.size(), 3);
//    Vg_pattern_half.block(0, 0, Vg_temp.rows(), 3) = Vg_temp;
//    idx = Vg_pattern_half.rows();
//    for(auto vert : addedVert){
//        Vg_pattern_half.row(idx) = vert;
//        idx ++;
//    }

//    cout<<Vg_pattern_half.rows()<<" half pattern vertices "<<endl;
}



//
//void oneShotLengthOpt(MatrixXi& Fg, int verts ){
//    cout<<endl<<endl<<endl;
//
//    int rowCount = 0;
//    std::vector<T> tripletList;
//    //UPDATE 2.12.22 we add a weighting , see https://online.stat.psu.edu/stat501/lesson/13/13.1
//    // we force the symmetry by having weight 1, the jacobian to have weight 0.5 as intial guess
//
//    for(int i =0; i<1; i++){//Fg.rows()
//        for(int j=0; j<3; j++){
//            int v0 = Fg(i, j);
//            v0=j;
//            tripletList.push_back(T(rowCount, 3 * v0, 1));
//            tripletList.push_back(T(rowCount+1, 3 * v0 + 1, 1));
//            tripletList.push_back(T(rowCount+2, 3 * v0 + 2, 1));
//
//            rowCount+=3;
//            if(i == 0) cout<<3 * v0<<" vert of face "<<v0<<" "<<rowCount-3<<endl;
//        }
//    }
//    A.resize( 1*9, 3* 4 );//Fg.rows()verts
//    A.setFromTriplets(tripletList.begin(), tripletList.end());
//    cout<<MatrixXd(A)<<endl;
////    cout<<MatrixXd(A)(1, 424) <<endl;
////    cout<<MatrixXd(A)(2, 425) <<endl;
////
//    SpMat LHS = (A.transpose() * A);
//    MatrixXd LHSd= MatrixXd(LHS);
//    for(int i=0; i<LHSd.rows(); i++){
//        for(int j=0; j<LHSd.cols(); j++){
//            if(LHSd(i,j)!=0){
//                cout<<i<<" i "<<j<<") = "<<LHSd(i,j)<<endl;
//            }
//        }
//    }
//    cholSolver.analyzePattern(LHS);
//    cholSolver.factorize(LHS);
//
//
//}
//
//void oneShotLengthSolve(MatrixXd& p_adaption, MatrixXi& Fg_pattern_curr, MatrixXd& baryCoordsUPattern, MatrixXd& baryCoordsVPattern, MatrixXd& mapFromVg, MatrixXi& mapFromFg){
//    int rowCounter = 0;
//
//    if(numVertPattern!= p_adaption.rows()){
//        numVertPattern = p_adaption.rows();
//        oneShotLengthOpt( Fg_pattern_curr, p_adaption.rows());
//    }
//
//    VectorXd b(1 * 9);
//    for(int i=0; i<1; i++){
//        Eigen::MatrixXd targetPositions(2, 3);
//        targetPositions.col(0)=  p_adaption.row(Fg_pattern_curr(i, 0)).leftCols(2).transpose() ;
//        targetPositions.col(1)=  p_adaption.row(Fg_pattern_curr(i, 1)).leftCols(2).transpose() ;
//        targetPositions.col(2)=  p_adaption.row(Fg_pattern_curr(i, 2)).leftCols(2).transpose() ;
//
//        VectorXd thisFaceU = baryCoordsUPattern(i,0) * targetPositions.col(0) +  baryCoordsUPattern(i,1) * targetPositions.col(1) + baryCoordsUPattern(i,2) * targetPositions.col(2) ;
//        VectorXd thisFaceV = baryCoordsVPattern(i,0) * targetPositions.col(0) +  baryCoordsVPattern(i,1) * targetPositions.col(1) + baryCoordsVPattern(i,2) * targetPositions.col(2) ;
//        VectorXd bary = targetPositions.rowwise().mean();
//
//        Eigen::MatrixXd Jnorm(2, 2);
//        Jnorm.col(0)= (thisFaceU - bary).normalized();
//        Jnorm.col(1)= (thisFaceV - bary).normalized();
//
//        double angle = acos((Jnorm.col(0)).dot(Jnorm.col(1)));
//        double deg = angle*180/M_PI;
//        double delta = abs(90-deg)/2;
//        delta = delta/180 * M_PI;
//        double DiagStiffness = abs(90-deg)/90;
//        Eigen::Matrix2d newRot= Eigen::MatrixXd::Identity(2, 2);
//        newRot(0, 0) = cos(DiagStiffness * delta);
//        newRot(1, 1) = newRot(0, 0);
//        newRot(0, 1) = - sin(DiagStiffness * delta);
//        newRot(1, 0) = sin( DiagStiffness * delta);
//
//        if(deg<=90){
//            Jnorm.col(0) = newRot.transpose() * Jnorm.col(0);
//            Jnorm.col(1) = newRot * Jnorm.col(1);
//        }else{
//            Jnorm.col(1) = newRot.transpose() * Jnorm.col(1);
//            Jnorm.col(0) = newRot * Jnorm.col(0);
//        }
//
//        Eigen::MatrixXd patternCoords(2, 3);
//        patternCoords.col(0) = mapFromVg.row(mapFromFg(i, 0)).leftCols(2).transpose();
//        patternCoords.col(1) = mapFromVg.row(mapFromFg(i, 1)).leftCols(2).transpose();
//        patternCoords.col(2) = mapFromVg.row(mapFromFg(i, 2)).leftCols(2).transpose();
//
//        if(rowCounter != 9*i) cout<<rowCounter<<" i "<<i<<endl;
//        Eigen::MatrixXd jacobiStretchedPattern = Jnorm * patternCoords;
//        Eigen::VectorXd pb = jacobiStretchedPattern.rowwise().mean();
//        Eigen::VectorXd qb = targetPositions.rowwise().mean();
//
//        VectorXd T_est = qb - pb;
//        Eigen::MatrixXd refTargetPos = jacobiStretchedPattern.colwise() + T_est;
//
//        b.block(rowCounter, 0, 2, 1) = refTargetPos.col(0);
//        rowCounter += 2;
//        b(rowCounter) =  mapFromVg(mapFromFg(i, 0), 2);
//        rowCounter++;
//
//        b.block(rowCounter, 0, 2, 1) = refTargetPos.col(1);
//        rowCounter += 2;
//        b(rowCounter) =  mapFromVg(mapFromFg(i, 0), 2);
//        rowCounter++;
//
//        b.block(rowCounter, 0, 2, 1) = refTargetPos.col(2);
//        rowCounter += 2;
//        b(rowCounter) =  mapFromVg(mapFromFg(i, 0), 2);
//        rowCounter++;
//    }
//
//
////    b= VectorXd::Ones(1 * 9);
//    MatrixXd RHS = A.transpose() * b;
//    for(int i=0; i<RHS.rows(); i++){
//        for(int j=0; j<RHS.cols(); j++){
//            if(RHS(i,j)!=0){
//                cout<<i<<" i "<<" = "<<RHS(i,j)<<endl;
//            }
//        }
//    }
//    cout<<RHS.rows()<<" RHS rows and cols "<<RHS.cols()<<endl;
//
//    auto p_asVec = cholSolver.solve(RHS);
//    if(A.transpose()*A*p_asVec != A.transpose()*b){
//        cout<<"solution wrong"<<endl ;
//    }
//    cout<<"  after solver "<<p_asVec.rows() <<endl ;
//    cout<<p_asVec(423)<<endl;
//    for(int i=0; i<p_asVec.rows(); i++){
//        if(p_asVec(i)!=0){
//            cout<<i<<" "<<p_asVec(i)<<endl ;
//        }
//    }
//
////    for(int i=0; i<p_adaption.rows(); i++){
////        p_adaption.row(i) = p_asVec.block(3*i, 0, 3, 1).transpose();
////        if(i> 161 && i<164) cout<<p_adaption.row(i)<<" i after "<<i <<" "<<p_asVec(3*i)<<endl ;
////        if(i==141) cout<<p_adaption.row(i)<<" i after "<< i<<endl ;
////
////    }
//}