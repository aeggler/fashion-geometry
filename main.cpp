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
Eigen::MatrixXi Fg, Fm;
Eigen::MatrixXd Vg_orig, Vm_orig; // original mesh for the garment and mannequin, restore for translation
Eigen::MatrixXi Fg_orig, Fm_orig;
Eigen::MatrixXi Eg; // garment edges
Eigen::Vector3d ambient, ambient_grey, diffuse, diffuse_grey, specular;
double EPS = 0.002f;							// in m (= 2 mm)
Eigen::Vector3f garment_translation (0., 2.2, 0.);// optimal for the given garment mesh
float garment_scale = 0.85;
Eigen::Vector3f mannequin_translation (0., 0., 0.);
float mannequin_scale = 1.;
bool simulate= false;

bool gar_loaded = false;// is a garment loaded
bool man_loaded = false; // is a mannequin loaded

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
Eigen::MatrixXi E_m;					// triangles = faces of the garment mesh / faces of the collision mesh
Eigen::VectorXi EMAP_m;
Eigen::VectorXd w; // the particle weights
Eigen::MatrixXd p; // the proposed new positions
PositionBasedDynamics PBD;
Real timestep= 0.0005;
double stretchStiffness= 0.8;
double collisionStiffness = 1.;
Real bendingStiffness = 0.1;
double coll_EPS= 0.0005;
int num_const_iterations = 5;
double blowFact= 0.001;
MatrixXd Vm_incr ;


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

bool pre_draw(igl::opengl::glfw::Viewer& viewer){
    viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;

    if(simulate){
        dotimeStep(viewer);
        showGarment(viewer);

    }
//    viewer.data().set_vertices(Vg);

    return false;
}
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
    string garment_file_name ="/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/garment/tshirt_1710.obj";
    igl::readOBJ(garment_file_name, Vg, Fg);
    igl::readOBJ(garment_file_name, Vg_orig, Fg_orig);

    preComputeConstraintsForRestshape();
    setNewGarmentMesh(viewer);

    string avatar_file_name ="/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/avatar/avatar_1710.obj";
    //string avatar_file_name = igl::file_dialog_open();
    igl::readOBJ(avatar_file_name, Vm, Fm);
    Vm_orig = Vm; Fm_orig = Fm;
    Vm_incr = (1+blowFact)*Vm;
    setNewMannequinMesh(viewer);
    cout<<" setting collision mesh "<<endl;
    setCollisionMesh();
    cout<<" collision mesh finished "<<endl;


    viewer.core().animation_max_fps = 200.;
    viewer.core().is_animating = false;

    //additional menu items
    menu.callback_draw_viewer_menu = [&]() {
        menu.draw_viewer_menu();
        if (ImGui::CollapsingHeader("Garment", ImGuiTreeNodeFlags_DefaultOpen)) {

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
            ImGui::InputDouble("Collision thereshold", &(coll_EPS),  0, 0, "%0.4f");
            ImGui::InputInt("Number of constraint Iterations thereshold", &(num_const_iterations),  0, 0);
            ImGui::InputDouble("Mannequin blowup ", &(blowFact),  0, 0, "%0.4f");

            if(ImGui::Button("set new blowup", ImVec2(-1, 0))){
                translateMesh(viewer, 2 );
            }


        }

    };

    // Add content to the default menu window
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &callback_key_down;

    viewer.launch();
}

void setNewGarmentMesh(igl::opengl::glfw::Viewer& viewer) {
    //Vg_initial = Vg;
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
    //viewer.core().align_camera_center(viewer.data().V, viewer.data().F);
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

        keyRecognition= true;
    }

    if (key == 'R')
    {
        reset(viewer);
        //reset simulation
        keyRecognition= true;
    }
    return keyRecognition;
}

//simulation part
void reset(igl::opengl::glfw::Viewer& viewer){
    cout<<" reset "<<endl;
    cout<<"---------"<<endl;
    cout<<Vg_orig.row(0)<<endl;
    Vg= Vg_orig;
    vel = Eigen::MatrixXd::Zero(numVert, 3);
    viewer.core().is_animating = false;
    simulate=false;
    translateMesh(viewer, 1);
    showGarment(viewer);
}

void preComputeConstraintsForRestshape(){
    numVert= Vg.rows();
    numFace = Fg.rows();

    w = Eigen::VectorXd::Ones(numVert);
    vel = Eigen::MatrixXd::Zero(numVert, 3);
    vel.col(1) =    w * (-1) * grav;
    igl::edge_lengths(Vg, Fg, edgeLengths);
    createFacePairEdgeListWith4VerticeIDs(Fg, e4list);
    e4size= e4list.rows();
    Q.resize(e4size, 1);
    cout<<" init triangle pairs " <<e4size<<endl;
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
    cout<<" setup finished"<<endl;
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
    double eps = coll_EPS; //0.005;
    for(int i=0; i<numVert; i++){
        if(S(i)<eps){
            collisionVert(i)=1;
        }
    }
}
void dotimeStep(igl::opengl::glfw::Viewer& viewer){
    Timer t("Time step ");
    std::cout<<endl;
    std::cout<<"-------------- Time Step ------------"<<endl;

    Eigen::MatrixXd x_new = Vg;
    p = Vg;

    // for fixing the vertices on top, set 32-40 to one
    Eigen::MatrixXd fixedVert= Eigen::MatrixXd::Zero(numVert, 1);

    // line (5) of the algorithm https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
    // we only use it to add gravity to the system
    vel.col(1) += timestep * w*(-1)*grav;

    // (7)
    for (int i = 0; i<numVert; i++){
        p.row(i) = x_new.row(i).array()+ timestep*vel.row(i).array();
    }
    // here we need to detect collisions and solve for them in the loop
    //t.printTime("do setup ");

    setupCollisionConstraints();
   // t.printTime(" setup collision constraints finished in ");

    //(9)-(11), the loop should be repeated several times per timestep (according to Jan Bender)
    for(int i=0; i < num_const_iterations; i++){
        // the bending, for each pair of adjacent triangles
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
               // if (!fixedVert(id0))
                    p.row(id0) += deltap0;
              //  if (!fixedVert(id1))
                    p.row(id1) += deltap1;
               // if (!fixedVert(id2))
                    p.row(id2) += deltap2;
               // if (!fixedVert(3))
                    p.row(id3) += deltap3;

            }
       // t.printTime("computed bending");

        //each edges distance should remain, now I have each edge covered twice (from both faces, but that should not be a problem)
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
       // t.printTime("computed stretch");

        // we precomputed the normal and collision point of each vertex, now add the constraint
        // this is imprecise but much faster and acc to paper works fine in practice

        for (int j=0; j<numVert; j++){
            if(collisionVert(j)){
                Vector3r deltap0;
                 PBD.solve_CollisionConstraint(p.row(j), w(j), C.row(j), N.row(j), deltap0, coll_EPS);

                 p.row(j) += collisionStiffness * deltap0;

            }
        }
       // t.printTime("computed collision");
    }
    // (13)
    for(int i=0; i<numVert; i++){
        for(int j=0; j<3; j++){
            vel(i,j)= (p(i,j)-x_new(i,j))/timestep;
        }
        x_new.row(i)=p.row(i);
    }
    //(14)
    Vg= x_new;
    // 16 Velocity update
    /*Friction and restitution can be handled by manipulating the velocities of colliding vertices in
     * step (16) of the algorithm. The velocity of each vertex for which a collision
     * constraint has been generated is dampened perpendicular to the collision normal
     * and reflected in the direction of the collision normal.*/
    double collision_damping = 0.5;
    for(int i=0; i<numVert; i++){
        if(collisionVert(i)){
            // friction update
            Vector3d vel_i = vel.row(i);
            Vector3d normal = N.row(i).normalized();
            Vector3d n_vel = normal.dot(vel_i)* normal;
            Vector3d t_vel = vel_i-n_vel;
            vel.row(i)= (t_vel-collision_damping *n_vel);
        }
    }


    showGarment(viewer);
}
