//
// Created by Anna Maria Eggler on 30.03.23.
//

#include "exportFunctions.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>  // include the fstream header file
#include "../Clipper/clipper.hpp"
#include <map>
#include <igl/hsv_to_rgb.h>
#include <cmath>
#include <igl/barycentric_coordinates.h>
#include <igl/doublearea.h>
#include "postProcessing.h"
#include <igl/vertex_components.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>

using namespace Eigen;
using namespace std;
using namespace ClipperLib;
void writeMTL(MatrixXd& Ka, MatrixXd& Ks, MatrixXd& Kd, MatrixXd& Vg, MatrixXi& Fg, string garment, string avName, double interp, string dir){
    cout<<" write file "<<endl;
    string fileName = "testing_"+avName + "_"+garment+"_stress_in_"+dir+"_interpolate_"+to_string(interp);
    string objName= fileName+".obj";
    string mtlName = fileName +".mtl";

    ofstream outfile(objName);
    outfile << "mtllib " << mtlName << ".mtl" << endl;
    for (int i=0; i<Vg.rows(); i++) {
        outfile << "v " << Vg(i,0) << " " << Vg(i,1) << " " << Vg(i,2) << endl;
    }
    outfile << endl;
    for (int i = 0; i < Fg.rows(); i++) {
        outfile << "usemtl face" << i << endl;
        outfile << "f " << Fg(i,0)+1  << "/" << i+1 << " " << Fg(i,1)+1 << "/" << i+1 << " " <<Fg(i,2)+1 << "/" << i+1 << endl;
    }
    outfile.close();

    ofstream mtlFile(mtlName);
    for (int i = 0; i < Ka.rows(); i++) {
        mtlFile << "newmtl face" << i << endl;
        mtlFile << "Ka " << Ka(i,0)  << " " << Ka(i,1) << " " << Ka(i,2) << endl;
        mtlFile << "Kd " << Kd(i,0)  << " " << Kd(i,1) << " " << Kd(i,2) << endl;
        mtlFile << "Ks " << Ks(i,0)  << " " << Ks(i,1) << " " << Ks(i,2) << endl;
        mtlFile <<" Tr "<< 1.000000<<endl;
        mtlFile <<" illum "<< 2<<endl;
        mtlFile <<" Ns "<< 0.000000<<endl;

    }
    mtlFile.close();

}
bool isPointInTriangle(Vector3d& curr,Vector3d& a,Vector3d& b,Vector3d& c ) {
    VectorXd bary;
    igl::barycentric_coordinates(curr.transpose(), a.transpose(), b.transpose(), c.transpose(), bary);
    if (bary(0) + bary(1) + bary(2) != 1) {
        return false;
    }
    if (bary(0) > 1 || bary(0) < 0) return false;
    if (bary(1) > 1 || bary(1) < 0) return false;
    if (bary(2) > 1 || bary(2) < 0) return false;
    return true;
}
double triangleArea(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ac = c - a;
    double area = 0.5 * (ab.cross(ac)).norm();
    return area;
}
bool is_ear( vector<VectorXd>& p, int a, int b, int c){
    Vector3d p_a = p[a];
    Vector3d p_b = p[b];
    Vector3d p_c = p[c];

    Vector3d ba = p_b-p_a;
    Vector3d bc = p_c - p_b;
    auto crossp = ba.cross(bc);
    if(crossp(2) <= 0) {
        return false;
    }

    for(int i=0; i<p.size(); i++){
        if(i==a|| i ==b || i==c){
            continue;
        }
        Vector3d curr = p[i];
        if(isPointInTriangle(curr, p_a, p_b, p_c)){
            return false;
        }
    }
    auto banew = p_a-p_b;
    auto angle = std::acos(banew.normalized().dot(bc.normalized()))* 180 / M_PI;
    double thereshAngle = 10;// smaller angle and merge for better pattern!
    if(angle <thereshAngle|| angle >360-thereshAngle ){

        return true;
    }
    return false;
}
void clipEar( vector<vector<VectorXd>>& returnVec){
    int count = 0;
    for(int i=0; i<returnVec.size(); i++){
        MatrixXd cliV; MatrixXi cliF;
        startRetriangulation(returnVec[i], cliV, cliF);
        VectorXd dblA;
        igl::doublearea(cliV, cliF, dblA);
        double sum=0;
        for(int i=0; i<dblA.rows(); i++){
            sum+= dblA(i);
        }
        cout<<"Sum of "<<i<<" is "<<sum<<endl;
        if(sum< 100){
            cout<<"sorry too small"<<endl;
            continue;
        }

        count++;

        int size = returnVec[i].size();
        for(int j=0; j<returnVec[i].size(); j++) {
            cout << returnVec[i][(j) ].transpose() << endl;
        }
        int count = 0;
        int currVert = 0;
        bool has_earF = true;

        while (has_earF){
            has_earF = false;
            for(int j=0; j<returnVec[i].size(); j++){

                if(is_ear(returnVec[i], (j-1+size) % size,(j) % size,(j+1) % size)){
                    has_earF = true;
                }
            }
            if(!has_earF) break;
            while(!is_ear(returnVec[i], (currVert-1+size)%size, (currVert)%size, (currVert+1)%size )) {
                currVert++;
                currVert %= returnVec[i].size();
            }

            returnVec[i].erase(returnVec[i].begin()+currVert);

        }
    }

}
void clipDifference(vector<vector<int>>& boundaryL_adaptedFromPattern,vector<vector<int>>& boundaryL_toPattern, MatrixXd &
    currPattern, MatrixXd& Vg_to,  vector<vector<VectorXd>>& returnVec){

    Paths subj_dbl(boundaryL_adaptedFromPattern.size()), clip_dbl(boundaryL_toPattern.size()), clip_bef(boundaryL_toPattern.size());
    double offset_distance = 10.4;
    // Create a Clipper object

    int mult = 10000;
    int toCut;
    bool isInverse = false;
    if(boundaryL_adaptedFromPattern.size() == boundaryL_toPattern.size()){
        // heuristic for inverse
        toCut = boundaryL_adaptedFromPattern.size()/2;
    }else{
        toCut= boundaryL_adaptedFromPattern.size();
        isInverse = true;
    }
    for(int i=0; i<toCut ; i++){
        Path p;
        for(int j=0; j<boundaryL_adaptedFromPattern[i].size(); j++){
            int var = boundaryL_adaptedFromPattern[i][j];
            double x = currPattern(var, 0); int xi = x*mult;
            double y = currPattern(var, 1); int yi = y*mult;
            p<< IntPoint(xi, yi);

        }
        // Add the input polygon to the Clipper object
        ClipperOffset co;
        co.AddPath(p, jtSquare, etClosedPolygon);
        // Perform the offset operation
        Paths offset_clip;
        co.Execute(offset_clip, offset_distance);
        Path pp;
        for (int j = 0; j <offset_clip[0].size(); j++) {
            pp << IntPoint( offset_clip[0][j].X , offset_clip[0][j].Y );

        }

        clip_dbl.push_back(pp);
    }

    for(int i=0; i< boundaryL_toPattern.size()/2;i++){
        Path pp;
        for(int j=0; j<boundaryL_toPattern[i].size(); j++){

            int var = boundaryL_toPattern[i][j];
            double x = Vg_to(var, 0); int xi = x*mult;
            double y = Vg_to(var, 1); int yi = y*mult;
            double minDist = 0.;
            for(int ii=0; ii<toCut ; ii++){
                for(int jj=0; jj<boundaryL_adaptedFromPattern[ii].size(); jj++){
                    int vvar = boundaryL_adaptedFromPattern[ii][jj];
                    double dist = (currPattern.row(vvar)  - Vg_to.row(var)).norm();
                    if(dist<minDist) {
                        minDist = dist;
                        x =  currPattern(vvar, 0); xi = x*mult;
                         y = currPattern(vvar, 1); yi = y*mult;
                    }
                }
            }

            pp<<IntPoint(xi, yi );
        }
        subj_dbl.push_back(pp);
    }



    // create a Clipper object and clip the polygons
    Clipper c;
    Paths solution;
    c.AddPaths(subj_dbl, ptSubject, true);
    c.AddPaths(clip_dbl, ptClip, true);
    c.Execute(ctDifference, solution, pftNonZero, pftNonZero);

    // print the resulting solution
    double th = 0.001;
    for (int i = 0; i < solution.size(); i++) {
//        cout << "Polygon " << i << ": ";
        vector<VectorXd> poly;
        set<int> xmap, ymap;
        for (int j = 0; j < solution[i].size(); j++) {
//            cout << "(" << solution[i][j].X/ mult << "," << solution[i][j].Y/ mult << ") ";
            VectorXd point(3);
            point(0)= solution[i][j].X; point(0)/= double(mult);
            point(1) = solution[i][j].Y/ double(mult);
            point(2) = 200;
            bool ignore = false;
            for(auto it: poly){
                if((it-point).norm()<th){// should go!!
                    cout<<"too close"<<endl;
                    point(0)+= 0.01;
                    ignore = true;
                    continue;
                }
            }if(!ignore){
                poly.push_back(point);
            }

        }
        if(poly.size()<3){
            cout<<"skip this poly, degenerate"<<endl;
            continue;
        }
        returnVec.push_back(poly);
        cout << endl;
    }


    clipEar(returnVec);
}

// Convert HSV color to RGB color
void hsv_to_rgb(float h, float s, float v, float& r, float& g, float& b)
{
    if (s == 0.0) {
        r = g = b = v;
    }
    else {
        h *= 6.0;
        int i = static_cast<int>(h);
        float f = h - static_cast<float>(i);
        float p = v * (1.0 - s);
        float q = v * (1.0 - s * f);
        float t = v * (1.0 - s * (1.0 - f));

        switch (i) {
            case 0:
                r = v; g = t; b = p; break;
            case 1:
                r = q; g = v; b = p; break;
            case 2:
                r = p; g = v; b = t; break;
            case 3:
                r = p; g = q; b = v; break;
            case 4:
                r = t; g = p; b = v; break;
            default:
                r = v; g = p; b = q; break;
        }
    }
}

// Generate a set of distinct colors using the HSV color space
void generate_colors(int num,MatrixXd& cols )
{
    std::vector<std::vector<float>> colors;
    colors.reserve(num);

    float golden_ratio_conjugate = (1 + std::sqrt(5)) / 2;
    float hue = std::fmod(std::rand() * golden_ratio_conjugate, 1.0f);

    for (int i = 0; i < num; ++i) {
        float saturation = 0.5f + (float)i / (float)num / 2.0f;
        float value = 0.95f - (float)i / (float)num / 20.0f;

        float r, g, b;
        hsv_to_rgb(hue, saturation, value, r, g, b);

        colors.emplace_back(std::vector<float>{ r, g, b });

        hue += golden_ratio_conjugate;
        hue = std::fmod(hue, 1.0f);
    }

    std::random_shuffle(colors.begin(), colors.end());
    for(int i=0; i<num; i++){
        cout<<endl;
        for(int j=0; j<colors[i].size(); j++){
            cols(i,j)= colors[i][j];
            cout<<colors[i][j];
        }
    }

}


void computeCols(int num, MatrixXd& cols){
    cols.resize(num, 3);
    generate_colors(num, cols);


}

void duplicateInitPattern(MatrixXd& Vg ,MatrixXi& Fg){
    MatrixXd saveV = Vg;
    MatrixXi saveF = Fg;
    int vertSize = Vg.rows();
    int faceSize = Fg.rows();
// only for patch 2!!!!!!!!!!!!!!!!

    Eigen::VectorXi componentIdPerVert;
    igl::vertex_components(Fg, componentIdPerVert);
    MatrixXd R_symetry = MatrixXd::Identity(3, 3);
    R_symetry(0, 0) = -1;

    for(int i = 0; i<faceSize/2; i++){
        //bigger is the better oe
        for(int j = 0; j<3; j++){
            int id = Fg(i, j);
            if(componentIdPerVert(id)==2){
                Vg.row(id) = (R_symetry * Vg.row(id+vertSize/2).transpose()).transpose();
                Vg(id, 0) -=50;

            }
        }
//        int temp = Fg(i, 1);
//        Fg(i, 1) = Fg(i, 2);
//        Fg(i, 2) = temp;
    }
    igl::writeOBJ("patternSym.obj" , Vg, Fg);
//    Vg = saveV;
//    Fg = saveF;

}

void addedSquare(MatrixXi Fg, MatrixXd Vg) {
    string name1 ;
    if( Vg(300, 2) ==200){
        name1 = "addedSquare_2D.obj";
    }else{
        name1 = "addedSquare_3D.obj";
    }
    igl::writeOBJ(name1 , Vg, Fg);
    return;
    MatrixXd whichVg = Vg;
    MatrixXi whichFg = Fg;
    cout<<whichFg.rows()<<" prev rows"<<endl;
    MatrixXd addedHelperVg (whichVg.rows() + 4 ,3);
    MatrixXi addedHelperFg (whichFg.rows() + 2 ,3);
    addedHelperVg.block(0,0, whichVg.rows(), whichVg.cols()) = whichVg;
    addedHelperFg.block(0,0, whichFg.rows(), whichFg.cols()) = whichFg;
    cout<<addedHelperFg.rows()<<" new rows"<<endl;

    int offHelp = whichVg.rows();
    int offHelpF = whichFg.rows();
    addedHelperFg(offHelpF, 0)= offHelp;
    addedHelperFg(offHelpF, 1)= offHelp+1;
    addedHelperFg(offHelpF, 2)= offHelp+2;

    addedHelperFg(offHelpF + 1, 0)= offHelp+1;
    addedHelperFg(offHelpF + 1, 1)= offHelp+3;
    addedHelperFg(offHelpF + 1, 2)= offHelp+2;

    addedHelperVg(offHelp, 0)= 500;
    addedHelperVg(offHelp, 1)= 1200;
    addedHelperVg(offHelp, 2)= 200;
    offHelp++;
    addedHelperVg(offHelp, 0)= 510;
    addedHelperVg(offHelp, 1)= 1200;
    addedHelperVg(offHelp, 2)= 200;
    offHelp++;

    addedHelperVg(offHelp, 0)= 500;
    addedHelperVg(offHelp, 1)= 1210;
    addedHelperVg(offHelp, 2)= 200;
    offHelp++;

    addedHelperVg(offHelp, 0)= 510;
    addedHelperVg(offHelp, 1)= 1210;
    addedHelperVg(offHelp, 2)= 200;

    string name ;
    if( addedHelperVg(300, 2) ==200){
        name = "addedSquare_2D.obj";
    }else{
        name = "addedSquare_3D.obj";
    }
    igl::writeOBJ(name , addedHelperVg, addedHelperFg);
}
void movePatches(){
    MatrixXi Fg;
    MatrixXd Vg;
    igl::readOBJ("addedSquare_2D.obj", Vg, Fg);
    VectorXi comp;
    igl::vertex_components(Fg, comp);

    for(int i=0; i<comp.size(); i++){
        if(comp(i) >=4&& comp(i) <= 9){
            Vg(i, 0) += 800;
        }else if (comp(i)>13 ){
            Vg(i, 0 ) -= 800;
        }
    }
    igl::writeOBJ("patchMoved2D.obj", Vg, Fg) ;

}