#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>
#include <igl/per_edge_normals.h>
#include <iostream>
#include <Eigen/Dense>
//#include "toolbox/TimeIntegration.h"
#include "toolbox/PositionBasedDynamics.h"
#include "toolbox/adjacency.h"
#include "toolbox/constraint_utils.h"
#include <igl/AABB.h>
#include "toolbox/Timer.h"
#include <igl/signed_distance.h>

using namespace std;
using namespace Eigen;
typedef double Real;
using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;

// The matrices of mesh and garment, original and modified
Eigen::MatrixXd Vg, Vm; // mesh for the garment and mannequin
Eigen::MatrixXi Fg, Fm, Fg_pattern;
Eigen::MatrixXd Vg_orig, Vm_orig; // original mesh for the garment and mannequin, restore for translation
Eigen::MatrixXd Vg_pattern, Vg_pattern_orig; // the pattern for the restshape, we might change this
Eigen::MatrixXi Fg_orig, Fm_orig;
Eigen::MatrixXi Eg; // garment edges
Eigen::Vector3d ambient, ambient_grey, diffuse, diffuse_grey, specular;
Eigen::Vector3f garment_translation (0., 0., 0.);// optimal for the given garment mesh
float garment_scale = 1.;
Eigen::Vector3f mannequin_translation (0., 0., 0.);
float mannequin_scale = 1.;
bool simulate= false;


//for the simulation

Real grav= 9.81;
Eigen::MatrixXi e4list;
int e4size, numVert, numFace;
Eigen::MatrixXd vel;
Eigen::Matrix<Matrix4r, Dynamic, 1> Q;
igl::AABB<Eigen::MatrixXd, 3> col_tree;
Eigen::MatrixXd edgeLengths;
MatrixXd C, N;
MatrixXi collisionVert;
Eigen::MatrixXd FN_m, VN_m, EN_m;	// vertices of the collision mesh
Eigen::MatrixXi E_m;				// triangles = faces of the garment mesh / faces of the collision mesh
Eigen::VectorXi EMAP_m;
Eigen::VectorXd w; // the particle weights
Eigen::MatrixXd p; // the proposed new positions
MatrixXd u1, u2; // precomputation for stretch
PositionBasedDynamics PBD;
Real timestep= 0.0005;
double stretchStiffness= 0.001;
double collisionStiffness = 1.;
Real bendingStiffness = 0.001;
double coll_EPS= 0.00; // like in Clo, 3 mm ? but for some reason this does not work well with the constraint function
int num_const_iterations = 5;
double blowFact= 0.003;// like in Clo, 3 mm ?
MatrixXd Vm_incr ;

bool pattern_loaded=false;
bool gar_loaded = false;
bool man_loaded = false;
MatrixXd faceAvg;
MatrixXd faceAvgWithU ;
MatrixXd faceAvgWithV ;
MatrixXd colU ;
MatrixXd colMixed;
MatrixXd colV ;
int whichStressVisualize= 0;


void setNewGarmentMesh(igl::opengl::glfw::Viewer& viewer);
void setNewMannequinMesh(igl::opengl::glfw::Viewer& viewer);
void showGarment(igl::opengl::glfw::Viewer& viewer);
void showMannequin(igl::opengl::glfw::Viewer& viewer);
void translateMesh(igl::opengl::glfw::Viewer& viewer, int which);

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers);
void reset(igl::opengl::glfw::Viewer& viewer);

void preComputeConstraintsForRestshape();
void dotimeStep(igl::opengl::glfw::Viewer& viewer);
void setCollisionMesh();
void setupCollisionConstraints();

void solveBendingConstraint();
void solveStretchConstraint();
void solveCollisionConstraint();
void preComputeStretch();
void computeStress(igl::opengl::glfw::Viewer& viewer);

bool pre_draw(igl::opengl::glfw::Viewer& viewer){
    viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;

    if(simulate){

        dotimeStep(viewer);
        computeStress(viewer);
        showGarment(viewer);// not sure if I actually need this, at least it breaks nothing
    }

    return false;
}
static bool noStress = true;
static bool StressU = false;
static bool StressV = false;
static bool StressMixed = false;
int main(int argc, char *argv[])
{
    // Init the viewer
    igl::opengl::glfw::Viewer viewer;

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    viewer.append_mesh();   // mesh for the garment
    viewer.append_mesh();   // mesh for the mannequin

    // set background
    viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 1);
    ambient_grey = Vector3d(0.4, 0.4, 0.4);
    ambient = Vector3d(0.26, 0.26, 0.26);
    diffuse_grey = Vector3d(0.5, 0.5, 0.5);
    diffuse = Vector3d(0.4, 0.57, 0.66);    // blue
    specular = Vector3d(0.01, 0.01, 0.01);

    // Load a mesh in OBJ format
    //string garment_file_name = igl::file_dialog_open();
    //for ease of use, for now let it be fixed

    string garment_file_name ="/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/dress2_3d_lowres/dress2_3d_lowres.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_merged_vertices_fixed.obj";
    igl::readOBJ(garment_file_name, Vg, Fg);
    igl::readOBJ(garment_file_name, Vg_orig, Fg_orig);

    string garment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/dress2_2d_lowres/dress2_2d_lowres.obj"; //"/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj";
    if(garment_pattern_file_name!=" "){
        pattern_loaded= true;
        igl::readOBJ(garment_pattern_file_name, Vg_pattern, Fg_pattern);
        Vg_pattern_orig= Vg_pattern;
    }

    preComputeConstraintsForRestshape();
    preComputeStretch();
    computeStress(viewer);
    setNewGarmentMesh(viewer);

    string avatar_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/avatar/avatar.obj";//avatar/avatar.obj";
    //string avatar_file_name = igl::file_dialog_open();
    igl::readOBJ(avatar_file_name, Vm, Fm);
    Vm_orig = Vm; Fm_orig = Fm;
    Vm_incr = (1+blowFact)*Vm;
    setNewMannequinMesh(viewer);
    cout<<" setting collision mesh "<<endl;
    setCollisionMesh();
    cout<<" collision mesh finished "<<endl;

    viewer.core().animation_max_fps = 100.;
    viewer.core().is_animating = false;

    //additional menu items
    menu.callback_draw_viewer_menu = [&]() {
        if (ImGui::CollapsingHeader("Garment", ImGuiTreeNodeFlags_OpenOnArrow)) {

            ImGui::InputFloat("Translation X", &garment_translation[0], 0, 0, "%0.4f");
            ImGui::InputFloat("Translation Y", &garment_translation[1], 0, 0, "%0.4f");
            ImGui::InputFloat("Translation Z", &garment_translation[2], 0, 0, "%0.4f");
            ImGui::InputFloat("Scaling factor X", &garment_scale, 0, 0, "%0.4f");
            if(ImGui::Button("Adjust garment", ImVec2(-1, 0))){
                translateMesh(viewer, 1 );
            }
        }
        if (ImGui::CollapsingHeader("Mannequin", ImGuiTreeNodeFlags_OpenOnArrow)) {

            ImGui::InputFloat("Translation X", &(mannequin_translation[0]),  0, 0, "%0.4f");
            ImGui::InputFloat("Translation Y", &(mannequin_translation[1]),  0, 0, "%0.4f");
            ImGui::InputFloat("Translation Z", &(mannequin_translation[2]),  0, 0, "%0.4f");
            ImGui::InputFloat("Scaling factor", &(mannequin_scale),  0, 0, "%0.4f");
            if(ImGui::Button("Adjust mannequin", ImVec2(-1, 0))){
                translateMesh(viewer, 2 );
            }
        }
        if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {

            if(ImGui::Button("Start Timestep", ImVec2(-1, 0))){
                dotimeStep(viewer);
                simulate= !simulate;
            }
            ImGui::InputDouble("Step size", &(timestep),  0, 0, "%0.4f");
            ImGui::InputDouble("Stretch Stiffness", &(stretchStiffness),  0, 0, "%0.4f");
            ImGui::InputDouble("Collision Stiffness", &(collisionStiffness),  0, 0, "%0.4f");
            ImGui::InputDouble("Bending Stiffness", &(bendingStiffness),  0, 0, "%0.4f");
            // figure out how that really works!!does not really do much
            ImGui::InputDouble("Collision thereshold", &(coll_EPS),  0, 0, "%0.6f");
            ImGui::InputInt("Number of constraint Iterations thereshold", &(num_const_iterations),  0, 0);
            ImGui::InputDouble("Mannequin blowup ", &(blowFact),  0, 0, "%0.4f");

            if(ImGui::Button("set new blowup", ImVec2(-1, 0))){
                translateMesh(viewer, 2 );
            }

            if(ImGui::Checkbox("Visualize no Stress", &noStress)){
                StressU = false;
                StressV=false;
                StressMixed = false;
                whichStressVisualize = 0;
                showGarment(viewer);
            }
            if(ImGui::Checkbox("Visualize U Stress", &StressU)){
                noStress = false;
                StressV = false;
                StressMixed = false;

                whichStressVisualize = 1;
                showGarment(viewer);
            }
            if(ImGui::Checkbox("Visualize V Stress ", &StressV)){
                noStress = false;
                StressU = false;
                StressMixed = false;
                whichStressVisualize = 2;
                showGarment(viewer);
           }
            if(ImGui::Checkbox("Visualize mixed Stress ", &StressMixed)){
                StressV= false;
                noStress = false;
                StressU = false;
                whichStressVisualize = 3;
                showGarment(viewer);
            }
            static bool remMan;
            if(ImGui::Checkbox("Remove Mannequin mesh ", &remMan)){
                viewer.selected_data_index = 1;
                viewer.data().show_faces = !viewer.data().show_faces;
            }

        }
        menu.draw_viewer_menu();
    };



    /*start new test of procrustes method */
    /*-----------------------------------*/

    cout<<"proc step started"<<endl;
    double rigidsum =0 ;
    for (int j =0; j<5; j++){
        int id0 = Fg_orig(j, 0);
        int id1 = Fg_orig(j, 1);
        int id2 = Fg_orig(j, 2);

        Matrix3d fromMat;// should be p after
        fromMat.row(0)= Vg.row(id0);
        fromMat.row(1)= Vg.row(id1);
        fromMat.row(2)= Vg.row(id2);

        Matrix3d toMat;
        toMat.row(0) = Vg_pattern.row(Fg_pattern(j, 0));
        toMat.row(1) = Vg_pattern.row(Fg_pattern(j, 1));
        toMat.row(2) = Vg_pattern.row(Fg_pattern(j, 2));


        Eigen::MatrixXd R_est;
        Eigen::VectorXd T_est;

        procrustes(fromMat, toMat, R_est, T_est);

        Eigen::MatrixXd appxToT = R_est * fromMat.transpose();
        appxToT = appxToT.colwise()+ T_est;
        MatrixXd appxTo = appxToT.transpose();
        cout<<" to mat "<<toMat<<endl;
        cout<<endl;
        cout<< " appx to mat "<<appxTo<<endl;

        Vector3d e1 = Vg_pattern.row(Fg_pattern(j, 1)) - Vg_pattern.row(Fg_pattern(j, 0));
        Vector3d e2 = Vg_pattern.row(Fg_pattern(j, 2)) - Vg_pattern.row(Fg_pattern(j, 0));

        Vector3d e1p = appxTo.row(1)- appxTo.row(0);
        Vector3d e2p = appxTo.row(2)- appxTo.row(0);

        double diff1 = e1(0) - e1p(0);
        double diff2 = e2(1) - e2p(1);
        rigidsum += (diff1 * diff1) + (diff2 * diff2)  ;
        cout<<rigidsum<<" rigidsum"<<endl ;
        cout<<endl;

    }
    cout<<"proc step finished"<<endl;





    // Add content to the default menu window
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &callback_key_down;

    viewer.launch();
}

void setNewGarmentMesh(igl::opengl::glfw::Viewer& viewer) {
    if (Vg.rows() == 0 || Fg.rows() == 0) {
        fprintf(stderr, "IOError: Could not load garment...\n");
        return;
    }
    igl::edges(Fg, Eg);
    showGarment(viewer);
    gar_loaded = true;
}
void showGarment(igl::opengl::glfw::Viewer& viewer) {
    viewer.selected_data_index = 0;
    viewer.data().clear();
    viewer.data().set_mesh(Vg, Fg);
    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_texture = false;
    viewer.data().set_face_based(false);
    //remove wireframe
    viewer.data().show_lines = false;
   // if 0 -> no face colour
    if(whichStressVisualize == 1){
        viewer.data().set_colors(colU);
    }else if (whichStressVisualize == 2){
        viewer.data().set_colors(colV);
    }else if (whichStressVisualize == 3 ){
        viewer.data().set_colors(colMixed);
    }

}
void setNewMannequinMesh(igl::opengl::glfw::Viewer& viewer) {
    if (Vm.rows() == 0 || Fm.rows() == 0) {
        fprintf(stderr, "IOError: Could not load garment...\n");
        return;
    }
    man_loaded=true;
    showMannequin(viewer);
}
void showMannequin(igl::opengl::glfw::Viewer& viewer) {
    viewer.selected_data_index = 1;
    viewer.data().clear();
    viewer.data().set_mesh(Vm, Fm);
    viewer.data().show_lines = false;
    viewer.data().uniform_colors(ambient_grey, diffuse_grey, specular);
    viewer.data().show_texture = false;
    viewer.data().set_face_based(false);
}
void translateMesh(igl::opengl::glfw::Viewer& viewer, int whichMesh ){
    if (whichMesh == 1){
        // we translate the garment, but based on the original mesh vg_orig in order to simplify gui
        for(int i=0; i<3; i++){
            auto temp = Vg_orig.col(i).array()+garment_translation[i];
            Vg.col(i) = temp;
        }
        Vg *= garment_scale;
        Vg_pattern *= garment_scale;
        // recompute the garment edge lengths
        for (int j=0; j<e4size; j++) {
            int id0 = e4list(j, 0);
            int id1 = e4list(j, 1);
            int id2 = e4list(j, 2);
            int id3 = e4list(j, 3);
            Vector3r pos0 = Vg.row(id0);
            Vector3r pos1 = Vg.row(id1);
            Vector3r pos2 = Vg.row(id2);
            Vector3r pos3 = Vg.row(id3);

            PBD.init_IsometricBendingConstraint(pos0, pos1, pos2, pos3, Q(j));
        }
        igl::edge_lengths(Vg, Fg, edgeLengths);
        preComputeConstraintsForRestshape();
        preComputeStretch();
        computeStress(viewer);
        showGarment(viewer);
    }else if (whichMesh == 2){
        // we translate the mannequin, hence the vertices Vm
        for(int i=0; i<3; i++){
            auto temp = Vm_orig.col(i).array()+mannequin_translation[i];
            Vm.col(i)=temp;
        }
        Vm *= mannequin_scale;
        showMannequin(viewer);
        Vm_incr = (1+blowFact)*Vm;
        setNewMannequinMesh(viewer);
        cout<<" setting collision mesh "<<endl;
        setCollisionMesh();
        cout<<" collision mesh finished "<<endl;
    }
}
bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers) {
    bool keyRecognition= false;

    if(key == 'S'){
        simulate= !simulate;
        viewer.core().is_animating= !viewer.core().is_animating;
        keyRecognition = true;
    }
    if (key == 'R')
    {
        reset(viewer);
        //reset simulation
        keyRecognition = true;
    }
    if (key == 'U')
    {
        whichStressVisualize = 1;
        StressV = false;
        StressU = true;
        noStress = false;
        showGarment(viewer);
        keyRecognition = true;
    }
    if (key == 'V')
    {
        whichStressVisualize = 2;
        StressV = true;
        StressU = false;
        noStress = false;
        showGarment(viewer);
        keyRecognition = true;
    }
    if (key == 'N')
    { // NONE stretch visualization
        whichStressVisualize = 0;
        StressV = false;
        StressU = false;
        noStress = true;
        showGarment(viewer);
        keyRecognition = true;
    }
    if(key == 'M'){
        whichStressVisualize = 3;
        StressV = false;
        StressU = false;
        noStress = false;
        StressMixed = true;
        showGarment(viewer);
        keyRecognition = true;
    }
    return keyRecognition;
}
//simulation part
void reset(igl::opengl::glfw::Viewer& viewer){
    cout<<" reset "<<endl;
    cout<<"---------"<<endl;
    Vg= Vg_orig;
    Vg_pattern = Vg_pattern_orig;
    vel = Eigen::MatrixXd::Zero(numVert, 3);
    Vm= Vm_orig;
    viewer.core().is_animating = false;
    simulate=false;
    StressV = false;
    StressU = false;
    StressMixed = false;
    noStress = true;

    whichStressVisualize = 0;
    translateMesh(viewer, 1);
    preComputeConstraintsForRestshape();
    preComputeStretch();
    computeStress(viewer);
    setNewGarmentMesh(viewer);
    //showGarment(viewer);
    showMannequin(viewer);
    setCollisionMesh();
}
void preComputeConstraintsForRestshape(){
    numVert= Vg.rows();
    numFace = Fg.rows();

    w = Eigen::VectorXd::Ones(numVert);
    vel = Eigen::MatrixXd::Zero(numVert, 3);
    vel.col(1) = w * (-1) * grav;

    // edge lengths from rest shape
    igl::edge_lengths(Vg_pattern, Fg_pattern, edgeLengths);

    createFacePairEdgeListWith4VerticeIDs(Fg, e4list);// from original since it has merged vertices adn thus more e4
    e4size= e4list.rows();
    // the pattern has fewer adjacent faces since the stitching does not count here but it does in the 3D case

    Q.resize(e4size, 1);

    for (int j=0; j<e4size; j++) {
        int id0 = e4list(j, 0);
        int id1 = e4list(j, 1);
        int id2 = e4list(j, 2);
        int id3 = e4list(j, 3);
        Vector3r pos0 = Vg.row(id0);
        Vector3r pos1 = Vg.row(id1);
        Vector3r pos2 = Vg.row(id2);
        Vector3r pos3 = Vg.row(id3);

        PBD.init_IsometricBendingConstraint(pos0, pos1, pos2, pos3, Q(j));
    }
}
void setCollisionMesh(){
    // the mesh the garment collides with -> Vm Fm the mannequin mesh
    //to have earlier detection, blow up the mannequin a bit and perform collision detection on this !

    col_tree.init(Vm_incr, Fm);
    igl::per_face_normals(Vm_incr, Fm, FN_m);
    igl::per_vertex_normals(Vm_incr, Fm, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_m, VN_m);
    igl::per_edge_normals(Vm_incr, Fm, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN_m, EN_m, E_m, EMAP_m);

}

void setupCollisionConstraints(){
    VectorXd S;
    VectorXi I;
    //MatrixXd C, N;
    collisionVert = Eigen::MatrixXi::Zero(numVert, 1);

    igl::signed_distance_pseudonormal(p, Vm_incr, Fm, col_tree, FN_m, VN_m, EN_m, EMAP_m, S, I, C, N);
    for(int i=0; i<numVert; i++){
        if(S(i)<coll_EPS){
            //if(S(i)<0){cout<<i<<" negative "<<S(i)<<endl; }
            collisionVert(i)=1;
        }
    }
}
void solveBendingConstraint(){
    // for each pair of adjacent triangles, precomputed from restshape (2D pattern TODO)
    for (int j = 0; j < e4size; j++) {
        Vector3r deltap0, deltap1, deltap2, deltap3;
        int id0 = e4list(j, 0);
        int id1 = e4list(j, 1);
        int id2 = e4list(j, 2);
        int id3 = e4list(j, 3);
        // now we got the correction delta p, add it to p
        PBD.solve_IsometricBendingConstraint(p.row(id0), w[id0], p.row(id1), w[id1], p.row(id2), w[id2],
                                             p.row(id3), w[id3], Q(j),
                                             bendingStiffness, deltap0, deltap1, deltap2, deltap3);
        p.row(id0) += deltap0;
        p.row(id1) += deltap1;
        p.row(id2) += deltap2;
        p.row(id3) += deltap3;
    }
}
void solveStretchConstraint(){
    /*each edges distance should remain, since we iterate over every face we iterate over every edge twice- but that should not be a problem */
    for (int j =0; j<numFace; j++){
        Vector3r deltap0, deltap1, deltap2;

        int id0 = Fg_orig(j, 0);
        int id1 = Fg_orig(j, 1);
        int id2 = Fg_orig(j, 2);

        // first edge
        PBD.solve_DistanceConstraint(p.row(id1), w(id1), p.row(id2), w(id2), edgeLengths(j, 0), stretchStiffness, deltap1, deltap2);
        p.row(id2) += deltap2;
        p.row(id1) += deltap1;

        //second edge
        PBD.solve_DistanceConstraint(p.row(id2), w(id2), p.row(id0), w(id0), edgeLengths(j, 1), stretchStiffness, deltap2, deltap0);
        p.row(id2)+= deltap2;
        p.row(id0)+= deltap0;

        // third edge
        PBD.solve_DistanceConstraint(p.row(id0), w(id0), p.row(id1), w(id1), edgeLengths(j, 2), stretchStiffness, deltap0, deltap1);
        p.row(id1)+= deltap1;
        p.row(id0)+= deltap0;
    }
}
void solveCollisionConstraint(){
    for (int j=0; j<numVert; j++){
        if(collisionVert(j)){
            Vector3r deltap0;
            // maybe I should compute the intersection instead of using the closest point C?
            PBD.solve_CollisionConstraint(p.row(j), w(j), C.row(j), N.row(j), deltap0, coll_EPS, vel.row(j));

            p.row(j) += collisionStiffness * deltap0;

        }
    }
}
void preComputeStretch(){
    faceAvg.resize(numFace, 3);
    faceAvgWithU.resize(numFace, 3);
    faceAvgWithV.resize(numFace, 3);
    colU = Eigen::MatrixXd ::Zero(numFace, 3);
    colV = Eigen::MatrixXd ::Zero (numFace, 3);
    colMixed = Eigen::MatrixXd ::Zero (numFace, 3);

    u1.resize(numFace, 2);
    u2.resize(numFace, 2);
    for(int j=0; j<numFace; j++){
        int id0 = Fg_pattern(j, 0);
        int id1 = Fg_pattern(j, 1);
        int id2 = Fg_pattern(j, 2);

        Vector2d u0, u1h, u2h;
        u0(0) = Vg_pattern(id0, 0);
        u0(1) = Vg_pattern(id0, 1);

        u1h( 0) = Vg_pattern(id1, 0);
        u1h(1) = Vg_pattern(id1, 1);

        u2h( 0) = Vg_pattern(id2, 0);
        u2h(1) = Vg_pattern(id2, 1);

        u1h -= u0;
        u2h -= u0;

        double det = u1h( 0)*u2h(1)- (u2h(0)*u1h(1));
        u1.row(j)= u1h/det;
        u2.row(j)= u2h/det;
    }
}
void computeStress(igl::opengl::glfw::Viewer& viewer){

    VectorXd normU (numFace);
    VectorXd normV (numFace);
    MatrixXd perFaceU (numFace, 3);
    MatrixXd perFaceV (numFace, 3);
    for(int j=0; j<numFace; j++){
        int id0 = Fg(j, 0);
        int id1 = Fg(j, 1);
        int id2 = Fg(j, 2);
        faceAvg.row(j) = (Vg.row(id0)+Vg.row(id1)+Vg.row(id2))/3;
        Vector3d p0 = Vg.row(id0);
        Vector3d p1 = Vg.row(id1);
        Vector3d p2 = Vg.row(id2);

        p1 -= p0;
        p2 -= p0;

        perFaceU.row(j) = p1*u2(j,1) - p2*u1(j,1);
        perFaceV.row(j) = p2 * u1(j, 0) - p1* u2(j, 0);
        double differenceIncrementFactor = 3.0;

        // deviation from 1 as the measure,
        /* large u stretch: norm > 1, thus y>0 , thus very red,little green => red
            compression = small u stretch: norm < 1, thus y<0 , thus little red, much green => green
            no stretch : y=0, thus very red , very green => yellow */
        normU(j) = perFaceU.row(j).norm();
        double y = (normU(j)-1) * differenceIncrementFactor; // to increase differences
        colU.row(j) = Vector3d(1.0 + y, 1. - y, 0.0);

        normV(j) = perFaceV.row(j).norm();
        y = (normV(j)-1) * differenceIncrementFactor;
        colV.row(j) = Vector3d ( 1. + y, 1.- y, 0.0);

        // this is an experiment
        y = (abs(normV(j)-1)+ abs(normU(j)-1))*3;
        colMixed.row(j) = Vector3d ( 1. + y, 1.- y, 0.0);

        faceAvgWithU.row(j) = faceAvg.row(j);
        faceAvgWithV.row(j) = faceAvg.row(j) ;

    }
    for(int i=0; i<numFace; i++){
        faceAvgWithU.row(i)+= (1 / normU.maxCoeff()) * normU(i) * 10 * perFaceU.row(i);
        faceAvgWithV.row(i)+= (1 / normV.maxCoeff()) * normV(i) * 10 * perFaceV.row(i);
    }
}

void dotimeStep(igl::opengl::glfw::Viewer& viewer){
    Timer t("Time step ");
    std::cout<<endl;
    std::cout<<"-------------- Time Step ------------"<<endl;
    // the stress is already computed, we can use it here
    Eigen::MatrixXd x_new = Vg;
    p = Vg;

    // line (5) of the algorithm https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
    // we only use it to add gravity to the system
    vel.col(1) += timestep * w*(-1)*grav*2;

    // (7)
    for (int i = 0; i<numVert; i++){
        p.row(i) = x_new.row(i).array()+ timestep*vel.row(i).array();
    }
    //t.printTime("do setup ");

    // detect collisions and solve for them in the loop
   // t.printTime(" setup collision constraints finished in ");
    setupCollisionConstraints();

    //(9)-(11), the loop should be repeated several times per timestep (according to Jan Bender)
    for(int i=0; i < num_const_iterations; i++){
        solveBendingConstraint();
        // t.printTime("computed bending");
        solveStretchConstraint();
        // t.printTime("computed stretch");


        /* we precomputed the normal and collision point of each vertex, now add the constraint (instead of checking collision in each iteration
         this is imprecise but much faster and acc to paper works fine in practice*/
//       setupCollisionConstraints();
        solveCollisionConstraint();

        // t.printTime("computed collision");
    }

    // (13) velocity and position update
    for(int i=0; i<numVert; i++){
        for(int j=0; j<3; j++){
            vel(i,j) = (p(i,j) - x_new(i,j)) / timestep;
        }
        x_new.row(i) = p.row(i);
    }
    //(14)
    Vg= x_new;
    // (16) Velocity update
    /*The velocity of each vertex for which a collision constraint has been generated is dampened perpendicular to the collision normal
     * and reflected in the direction of the collision normal.*/
    double collision_damping = 0.5;
    for(int i=0; i<numVert; i++){
        if(collisionVert(i)){
            Vector3d vel_i = vel.row(i);
            Vector3d normal = N.row(i).normalized();
            Vector3d n_vel = normal.dot(vel_i) * normal;
            Vector3d t_vel = vel_i-n_vel;
            vel.row(i)= (t_vel-collision_damping * n_vel);
        }
    }

    showGarment(viewer);
}
