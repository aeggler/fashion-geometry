//
// Created by Anna Maria Eggler on 30.03.23.
//

#include "exportFunctions.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>  // include the fstream header file
#include "../Clipper/clipper.hpp"
#include <map>


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
void clipDifference(vector<vector<int>>& boundaryL_adaptedFromPattern,vector<vector<int>>& boundaryL_toPattern, MatrixXd &
currPattern, MatrixXd& Vg_to,  vector<vector<VectorXd>>& returnVec){

    Paths subj_dbl(boundaryL_adaptedFromPattern.size()), clip_dbl(boundaryL_toPattern.size());
    double offset_distance = -2.20;
    // Create a Clipper object
    ClipperOffset co;
    for(int i=0; i< boundaryL_adaptedFromPattern.size();i++){
        Path p;
        for(int j=0; j<boundaryL_adaptedFromPattern[i].size(); j++){
            int var = boundaryL_adaptedFromPattern[i][j];
            double x = currPattern(var, 0); int xi = x*100;
            double y = currPattern(var, 1); int yi = y*100;
            p<< IntPoint(xi, yi);
        }
        // Add the input polygon to the Clipper object
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
            double x = Vg_to(var, 0); int xi = x*100;
            double y = Vg_to(var, 1); int yi = y*100;
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
double th = 0.1;
    for (int i = 0; i < solution.size(); i++) {
        cout << "Polygon " << i << ": ";
        vector<VectorXd> poly;
        set<int> xmap, ymap;
        for (int j = 0; j < solution[i].size(); j++) {
            cout << "(" << solution[i][j].X/ 100 << "," << solution[i][j].Y/ 100 << ") ";
            VectorXd point(3);
            point(0)= solution[i][j].X/ 100;
            point(1) = solution[i][j].Y/ 100;
            point(2) = 200;
            bool ignore = false;
            for(auto it: poly){
                if((it-point).norm()<th){// should go!!
                    cout<<"too close"<<endl;
                    ignore = true;
                    continue;
                }
            }if(!ignore){
                poly.push_back(point);
            }

        }
        if(poly.size()<3)continue;
        returnVec.push_back(poly);
        cout << endl;
    }
}