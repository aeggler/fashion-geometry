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
#include "toolbox/PositionBasedDynamics.h"
#include "toolbox/adjacency.h"
#include "toolbox/constraint_utils.h"
#include <igl/AABB.h>
#include "toolbox/Timer.h"
#include "toolbox/body_interpolation.h"
#include "toolbox/garment_adaption.h"
#include <igl/signed_distance.h>

using namespace std;
using namespace Eigen;
typedef double Real;
using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;

// The matrices of mesh and garment, original and modified
Eigen::MatrixXd Vg, Vm, testMorph_V1, testMorph_V0; // mesh for the garment and mannequin
Eigen::MatrixXi Fg, Fm, Fg_pattern, testMorph_F1, testMorph_F0;
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
Eigen::MatrixXd procrustesPatternIn3D;
Real timestep= 0.0005;
double stretchStiffnessU= 0.001;
double stretchStiffnessV= 0.001;
double collisionStiffness = 1.;
Real bendingStiffness = 0.001;
double rigidStiffness = 0.001;
double coll_EPS= 0.000; // like in Clo, 3 mm ? but for some reason this does not work well with the constraint function
int num_const_iterations = 5;
double blowFact= 0.005;// like in Clo, 3 mm ?
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
VectorXd normU, normV;
MatrixXd perFaceU, perFaceV;
int whichStressVisualize= 0;
garment_adaption* gar_adapt;
int counter;
BodyInterpolator* body_interpolator;
bool bodyInterpolation= false;


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
void computeRigidMeasure();
void initProcrustesPatternTo3D();
void solveRigidEnergy();
void solveStretchUV();

double itCount =0;
bool pre_draw(igl::opengl::glfw::Viewer& viewer){
    viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;
        if(simulate){
            for(int i=0; i< 1; i++){
                computeStress(viewer);
                dotimeStep(viewer);
                counter++;
                showGarment(viewer);// not sure if I actually need this, at least it breaks nothing
            }
           // cout<<itCount<<endl;
            if(itCount>1000) itCount = 0;
            double p = itCount/1000.;
//            body_interpolator->interpolateMesh(p, Vm);
//            Vm_incr = (1+blowFact)*Vm;
//            setCollisionMesh();
            showMannequin(viewer);
            showGarment(viewer);// not sure if I actually need this, at least it breaks nothing

            itCount += 3.;

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
    counter = 0;
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

    string garment_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/dress2_3d_lowres/dress2_3d_lowres.obj"; //"/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_merged_vertices_fixed.obj"; //
    //string garment_file_name = "/Users/annaeggler/Desktop/aShapeDressBetter";
    igl::readOBJ(garment_file_name, Vg, Fg);
    igl::readOBJ(garment_file_name, Vg_orig, Fg_orig);

    string garment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/shrinkedGarment_2D.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj"; //
//    string garment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/dress2_2d_lowres/dress2_2d_lowres.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj"; //
    if(garment_pattern_file_name!=" "){
        pattern_loaded= true;
        igl::readOBJ(garment_pattern_file_name, Vg_pattern, Fg_pattern);
        Vg_pattern_orig= Vg_pattern;
    }

    preComputeConstraintsForRestshape();
    preComputeStretch();
    computeStress(viewer);
    setNewGarmentMesh(viewer);
// remember to adapt the collision constraint solving dep on avatar, sometimes normalization is needed, sometimes not for whatever magic
    string avatar_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/avatar/avatar.obj";
    // string avatar_file_name = "/Users/annaeggler/Desktop/custom_fit_garments/data/example_avatar/avatar_005.ply"; //// /avatar/avatar.obj";//
    //string avatar_file_name = igl::file_dialog_open();
    igl::readOBJ(avatar_file_name, Vm, Fm);
//    igl::readPLY(avatar_file_name, Vm, Fm);
    Vm_orig = Vm; Fm_orig = Fm;
    Vm_incr = (1+blowFact)*Vm;

   // string morphBody0 ="/Users/annaeggler/Desktop/custom_fit_garments/data/example_avatar/avatar_007.ply";
//    string morphBody1 = "/Users/annaeggler/Desktop/custom_fit_garments/data/example_avatar/avatar_001.ply";
//
//    igl::readPLY(morphBody1, testMorph_V1, testMorph_F1);
////    igl::readPLY(morphBody0, testMorph_V0, testMorph_F0);
//    cout<<Fm.row(5)<<endl;
//    cout<<testMorph_F1.row(5)<<endl;
//    cout<<endl;
//    body_interpolator = new BodyInterpolator(Vm, testMorph_V1, Fm);

//    Fm = testMorph_F0;
//    Vm = testMorph_V0;

    setNewMannequinMesh(viewer);
    cout<<" setting collision mesh "<<endl;
    setCollisionMesh();
    cout<<" collision mesh finished "<<endl;




    /*-----------------------*/
    /* Garment adaption update */
    string orig_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/dress2_2d_lowres/dress2_2d_lowres.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj"; //

//    string orig_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/shrinkedGarment.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj"; //
    Eigen::MatrixXd Vg_orig_pattern;
    igl::readOBJ(orig_file_name, Vg_orig_pattern, Fg_pattern);

    string orig_file_name_sugg = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj"; //

//    string orig_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/shrinkedGarment.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj"; //
    Eigen::MatrixXd computedPatten;
    igl::readOBJ(orig_file_name_sugg, computedPatten, Fg_pattern);
    cout<<computedPatten.col(2).maxCoeff()<<" max and min "<<computedPatten.col(2).minCoeff()<<" and diff "<<computedPatten.col(2).maxCoeff()-computedPatten.col(2).minCoeff()<<endl;

    double col0 = 0; double col1 = 0; double col2 = 0;
    for(int i=0; i< computedPatten.rows(); i++){

        col0 += computedPatten(i,0) - Vg_orig_pattern(i,0);
        col1 += computedPatten(i,1) - Vg_orig_pattern(i,1);
        col2 += computedPatten(i,2) - Vg_orig_pattern(i,2);

    }
    cout<< col0<<" zero, one : "<<col1<<" and two: "<<col2<<endl;

//    string skrinked_file_name_pattern = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/shrinkedGarment_2D.obj"; // "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_2D_2_fixed.obj"; //
//    Eigen::MatrixXd Vg_shrinked_pattern;
//    igl::readOBJ(skrinked_file_name_pattern, Vg_shrinked_pattern, Fg_pattern);

    gar_adapt = new garment_adaption(Vg, Fg,  Vg_orig_pattern, Fg_pattern);
    gar_adapt->computeJacobian();



    /*test */
    Eigen::MatrixXd Vg_test_pattern(5, 3);
    Vg_test_pattern <<
            0,0,00,
            1,0,00,
            0,1,00,
            -1,0,00,
            0,-1,00;
    Eigen::MatrixXi Fg_test(4, 3);
    Fg_test <<
            0,1,2,
            0,2,3,
            0,3,4,
            0,4,1;

    Eigen::MatrixXd Vg_test= Vg_test_pattern;
    Vg_test(2, 2)= 1;
    Vg_test(1, 2)= -1;
    cout<<" testmesh 3D "<<endl<<Vg_test<<endl<<endl;

    Eigen::MatrixXd Vg_test_shrinked = 0.5*Vg_test;
    cout<<" testmesh 3D shrinked  "<<endl<<Vg_test_shrinked<<endl<<endl;

    garment_adaption gar_adapt_test = *new garment_adaption(Vg_test, Fg_test,  Vg_test_pattern, Fg_test);
    gar_adapt_test.computeJacobian();

    gar_adapt_test.performJacobianUpdateAndMerge(Vg_test_shrinked);
    cout<<endl<<"shrinked: "<<endl;
    cout<<Vg_test_shrinked<<endl;
//    cout<<Vg_shrinked.row(50)<<endl;
//    igl::writeOBJ("suggested_shrinkedGarment_2D.obj", Vg_shrinked, Fg_pattern);

    /* end garment adaption update */
    /*-----------------------------*/

    viewer.core().animation_max_fps = 200.;
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
            ImGui::InputDouble("U Stretch Stiffness ", &(stretchStiffnessU),  0, 0, "%0.4f");
            ImGui::InputDouble("V Stretch Stiffness", &(stretchStiffnessV),  0, 0, "%0.4f");
            ImGui::InputDouble("Collision Stiffness", &(collisionStiffness),  0, 0, "%0.4f");
            ImGui::InputDouble("Bending Stiffness", &(bendingStiffness),  0, 0, "%0.4f");
            ImGui::InputDouble("Rigidity Stiffness", &(rigidStiffness),  0, 0, "%0.6f");

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
    if(bodyInterpolation){
        viewer.data().set_mesh(testMorph_V0, testMorph_F0);
    }else{
        viewer.data().set_mesh(Vm, Fm);
    }

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
//        igl::writeOBJ("shrinkedGarment_2D.obj", Vg_pattern, Fg_pattern);
//        std::cout<<" 2D Garment file written"<<endl;

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
//    if(key== 'W'){
//        keyRecognition=true;
//        igl::writeOBJ("shrinkedGarment.obj", Vg, Fg);
//          std::cout<<" Garment file written"<<endl;
//    }
    if(key == 'P'){       // Pattern
        keyRecognition = true;
        simulate = false;
        // we start computing the pattern for the current shape
        Eigen::MatrixXd computed_Vg_pattern= Vg;
        cout<<"start computing the pattern"<<endl;
        gar_adapt->performJacobianUpdateAndMerge(computed_Vg_pattern);
        igl::writeOBJ("patternComputed.obj", computed_Vg_pattern, Fg_pattern);
        cout<<"pattern written to *patternComputed*"<<endl;

    }
    if(key == 'B'){     // Bending
       // reset(viewer);
//        simulate= false;
//        Fm = testMorph_F0;
//        Vm = testMorph_V0;
//        cout<<Vm.rows()<<" and faces "<<Fm.rows()<<endl;
//        viewer.selected_data_index = 0;
//        viewer.data().clear();
//        bodyInterpolation = !bodyInterpolation;
//        setNewMannequinMesh(viewer);
//        cout<<"should be set"<<endl;
//        body_interpolator = new BodyInterpolator(testMorph_V0, testMorph_V1, testMorph_F0);
//       //
        keyRecognition = true;

    }
    return keyRecognition;
}
//simulation part
void reset(igl::opengl::glfw::Viewer& viewer){
    cout<<" reset "<<endl;
    cout<<"---------"<<endl;
    counter = 0;
    Vg= Vg_orig;
   // Vg_pattern = Vg_pattern_orig;
    vel = Eigen::MatrixXd::Zero(numVert, 3);
    Vm= Vm_orig; Vm_incr = (1+blowFact)*Vm;
    Fm= Fm_orig;
    viewer.core().is_animating = false;
    simulate=false;
    StressV = false;
    StressU = false;
    StressMixed = false;
    noStress = true;

    whichStressVisualize = 0;
    //translateMesh(viewer, 1);
    preComputeConstraintsForRestshape();
    preComputeStretch();
    computeStress(viewer);
    setNewGarmentMesh(viewer);
    //showGarment(viewer);
    showMannequin(viewer);
    setCollisionMesh();
    itCount=0;
}
void preComputeConstraintsForRestshape(){
    numVert= Vg.rows();
    numFace = Fg.rows();

    w = Eigen::VectorXd::Ones(numVert);
    vel = Eigen::MatrixXd::Zero(numVert, 3);
    //vel.col(1) = w * (-1) * grav;

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

// OLD VERSION, NOW UV DIFFERENTIATION EXISTS
void solveStretchConstraint(){
    /*each edges distance should remain, since we iterate over every face we iterate over every edge twice- but that should not be a problem */
    for (int j =0; j<numFace; j++){
        Vector3r deltap0, deltap1, deltap2;

        int id0 = Fg_orig(j, 0);
        int id1 = Fg_orig(j, 1);
        int id2 = Fg_orig(j, 2);

        // first edge
        PBD.solve_DistanceConstraint(p.row(id1), w(id1), p.row(id2), w(id2), edgeLengths(j, 0), stretchStiffnessU, deltap1, deltap2);
        p.row(id2) += deltap2;
        p.row(id1) += deltap1;

        //second edge
        PBD.solve_DistanceConstraint(p.row(id2), w(id2), p.row(id0), w(id0), edgeLengths(j, 1), stretchStiffnessU, deltap2, deltap0);
        p.row(id2)+= deltap2;
        p.row(id0)+= deltap0;

        // third edge
        PBD.solve_DistanceConstraint(p.row(id0), w(id0), p.row(id1), w(id1), edgeLengths(j, 2), stretchStiffnessU, deltap0, deltap1);
        p.row(id1)+= deltap1;
        p.row(id0)+= deltap0;
    }
}
Eigen::MatrixXd tarU, tarV;
void init_stretchUV(){
    tarU.resize(3*numFace, 3);
    tarV.resize(3*numFace, 3);
    for(int j = 0; j<numFace; j++){
        Eigen::MatrixXd patternCoords(2, 3);
        patternCoords(0,0) = Vg_pattern( Fg_pattern(j, 0), 0);
        patternCoords(1,0) = Vg_pattern( Fg_pattern(j, 0), 1);
        patternCoords( 0,1) = Vg_pattern( Fg_pattern(j, 1), 0);
        patternCoords(1,1) = Vg_pattern( Fg_pattern(j, 1), 1);
        patternCoords( 0,2) = Vg_pattern( Fg_pattern(j, 2), 0);
        patternCoords(1,2) = Vg_pattern( Fg_pattern(j, 2), 1);

        Eigen::MatrixXd targetPositions(3, 3);
        targetPositions.col(0)= p.row(Fg_orig(j, 0));
        targetPositions.col(1)= p.row(Fg_orig(j, 1));
        targetPositions.col(2)= p.row(Fg_orig(j, 2));

        Vector3r taru0 , taru1, taru2;

        PBD.init_UVStretch(normU(j), perFaceU.row(j), perFaceV.row(j), patternCoords, targetPositions,
                           stretchStiffnessU, taru0, taru1, taru2, 1);
        tarU.row(3*j)= taru0.transpose();
        tarU.row(3*j+1) = taru1.transpose();
        tarU.row(3*j+2 )= taru2.transpose();

        Vector3r tarv0, tarv1, tarv2 ;
        PBD.init_UVStretch( normV(j),perFaceU.row(j),perFaceV.row(j), patternCoords, targetPositions,
                        stretchStiffnessU,tarv0, tarv1, tarv2, 2);
        tarV.row(3*j)= tarv0.transpose();
        tarV.row(3*j+1) = tarv1.transpose();
        tarV.row(3*j+2 )= tarv2.transpose();

    }
}
void solveStretchUV(){
    double stretchEPS= 0.001;
    cout<<normU.sum()<<" sum of the norm u, and v norm  "<<normV.sum()<<endl;

    for (int j =0; j<numFace; j++){

        Vector3r deltap0, deltap1, deltap2;
        Eigen::MatrixXd targetPositions(3, 3);
        targetPositions.col(0)= p.row(Fg_orig(j, 0));
        targetPositions.col(1)= p.row(Fg_orig(j, 1));
        targetPositions.col(2)= p.row(Fg_orig(j, 2));

        Vector3r tarV0= tarV.row(3*j +0);
        Vector3r tarV1= tarV.row(3*j +1);
        Vector3r tarV2= tarV.row(3*j +2);

        Vector3r tarU0= tarU.row(3*j +0);
        Vector3r tarU1= tarU.row(3*j +1);
        Vector3r tarU2= tarU.row(3*j +2);


        if(abs(normU(j)-1) > stretchEPS ){
            PBD.solve_UVStretch(normU(j), targetPositions, 1,
            tarU0, tarU1, tarU2, stretchStiffnessU,deltap0,  deltap1, deltap2);

            p.row(Fg_orig(j, 2)) += deltap2;
            p.row(Fg_orig(j, 1)) += deltap1;
            p.row(Fg_orig(j, 0)) += deltap0;
        }
        if(abs(normV(j)-1) > stretchEPS ){

            PBD.solve_UVStretch(normV(j), targetPositions,2,
                                tarV0, tarV1, tarV2,
                               stretchStiffnessV, deltap0, deltap1, deltap2);

            p.row(Fg_orig(j, 2)) += deltap2;
            p.row(Fg_orig(j, 1)) += deltap1;
            p.row(Fg_orig(j, 0)) += deltap0;
        }
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

     normU.resize (numFace);
     normV.resize (numFace);
     perFaceU.resize (numFace, 3);
     perFaceV.resize (numFace, 3);
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
Eigen::VectorXd rigidEnergy;
void computeRigidMeasure(){
    rigidEnergy.resize(numFace);

    for (int j =0; j<numFace; j++){
        int id0 = Fg_orig(j, 0);
        int id1 = Fg_orig(j, 1);
        int id2 = Fg_orig(j, 2);

        Matrix3d fromMat;// should be p after
        fromMat.row(0)= p.row(id0);
        fromMat.row(1)= p.row(id1);
        fromMat.row(2)= p.row(id2);

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

        Vector3d e1 = Vg_pattern.row(Fg_pattern(j, 1)) - Vg_pattern.row(Fg_pattern(j, 0));
        Vector3d e2 = Vg_pattern.row(Fg_pattern(j, 2)) - Vg_pattern.row(Fg_pattern(j, 0));

        Vector3d e1p = appxTo.row(1)- appxTo.row(0);
        Vector3d e2p = appxTo.row(2)- appxTo.row(0);

        double diff1 = e1(0) - e1p(0);
        double diff2 = e2(1) - e2p(1);
        rigidEnergy(j) = (diff1 * diff1) + (diff2 * diff2);

    }
}
void solveRigidEnergy(){
    double rigideps= 0.000003;
    for(int j=0; j<numFace; j++){
        Vector3r delta0, delta1, delta2;

        PBD.solve_RigidEnergy(rigidEnergy(j), rigideps, rigidStiffness,
                             procrustesPatternIn3D.row(3*j+ 0), procrustesPatternIn3D.row(3*j+ 1), procrustesPatternIn3D.row(3*j+ 2),
                             p.row(Fg(j, 0)), p.row(Fg(j, 1)), p.row(Fg(j, 2)), delta0, delta1, delta2
                             );

        p.row(Fg(j, 0)) += delta0;
        p.row(Fg(j, 1)) += delta1;
        p.row(Fg(j, 2)) += delta2;
    }
}

void dotimeStep(igl::opengl::glfw::Viewer& viewer){
    Timer t("Time step ");
    std::cout<<endl;
    std::cout<<"-------------- Time Step ------------"<<counter<<endl;
    // the stress is already computed, we can use it here
    Eigen::MatrixXd x_new = Vg;
    p = Vg;

    // line (5) of the algorithm https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
    // we only use it to add gravity to the system
    vel.col(1) += timestep * w*(-1)*grav*0.01;

    // (7)
    for (int i = 0; i<numVert; i++){
        p.row(i) = x_new.row(i).array()+ timestep*vel.row(i).array();
    }
    // detect collisions and solve for them in the loop
    setupCollisionConstraints();
//    t.printTime(" collision constraint ");
    computeRigidMeasure();
//    t.printTime(" rigid measure");
//
    initProcrustesPatternTo3D(Vg_pattern, Fg_pattern, Fg_orig, p, procrustesPatternIn3D); // might be an imprecise but fast option to remove this from the loop
//    t.printTime(" pattern dimensions ");
    init_stretchUV();

    //(9)-(11), the loop should be repeated several times per timestep (according to Jan Bender)
    num_const_iterations= 1;
    for(int i=0; i < num_const_iterations; i++){
        solveBendingConstraint();
        t.printTime("computed bending");
        //solveStretchConstraint();
        solveStretchUV();
//        t.printTime("computed stretch");
        solveRigidEnergy();
//        t.printTime(" rigid energy ");

        /* we precomputed the normal and collision point of each vertex, now add the constraint (instead of checking collision in each iteration
         this is imprecise but much faster and acc to paper works fine in practice*/
        solveCollisionConstraint();
        t.printTime(" collision ");cout<<endl;
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
            vel.row(i) = (t_vel - collision_damping * n_vel);
        }
    }

    showGarment(viewer);
}
