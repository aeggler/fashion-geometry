#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>
#include <igl/per_edge_normals.h>
#include <igl/adjacency_list.h>
#include <igl/facet_components.h>
#include <igl/vertex_components.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <iostream>
#include <Eigen/Dense>
#include "toolbox/PositionBasedDynamics.h"
#include "toolbox/adjacency.h"
#include "toolbox/constraint_utils.h"
#include <igl/AABB.h>
#include "toolbox/Timer.h"
#include "toolbox/body_interpolation.h"
#include "toolbox/garment_adaption.h"
#include "toolbox/MathFunctions.h"
#include <igl/signed_distance.h>
#include <map>
#include <string>
#include "toolbox/seam.h"
#include <fstream>

using namespace std;
using namespace Eigen;
typedef double Real;
using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;

// The matrices of mesh and garment, original and modified
Eigen::MatrixXd Vg, Vm, testMorph_V1; // mesh for the garment and mannequin
Eigen::MatrixXi Fg, Fm, Fg_pattern, testMorph_F1;
Eigen::MatrixXd Vg_orig, Vm_orig; // original mesh for the garment and mannequin, restore for translation
Eigen::MatrixXd Vg_pattern, Vg_pattern_orig; // the pattern for the restshape, we might change this
Eigen::MatrixXi Fg_orig, Fm_orig;
Eigen::MatrixXi Eg; // garment edges
Eigen::Vector3d ambient, ambient_grey, diffuse, diffuse_grey, specular;
Eigen::Vector3f garment_translation (0., 0., 0.);// optimal for the given garment mesh
float garment_scale = 1.;
Eigen::Vector3f mannequin_translation (0., 0., 0.);
float mannequin_scale = 1.;



//for the simulation
double timestep= 0.02;
bool simulate= false;
double grav= 9.81;
double stretchStiffnessU= 0.0010;
double stretchStiffnessV= 0.0010;
double stretchStiffnessD = 0.0080;
double collisionStiffness = 1.;
double boundaryStiffness = 0.9;
double bendingStiffness = 0.003;// smaller for better folds , bigger for smoother results
double edgeLengthStiffness = 0.00;
double coll_EPS= 4.500; // like in Clo, 3 mm ? but for some reason this does not work well with the constraint function
int num_const_iterations = 5;
double gravityfact =.0;
int localGlobalIterations= 2000;
int convergeIterations = 450;
int timestepCounter;

enum MouseMode { SELECTPATCH, SELECTBOUNDARY, NONE };
MouseMode mouse_mode = NONE;
// pre computations
Eigen::MatrixXi e4list;
int e4size, numVert, numFace;
Eigen::MatrixXd vel;
Eigen::Matrix<Matrix4r, Dynamic, 1> Q;
igl::AABB<Eigen::MatrixXd, 3> col_tree;
Eigen::MatrixXd edgeLengths;
MatrixXd C, N;
MatrixXi collisionVert;
vector<int> pureCollVert;

Eigen::MatrixXd FN_m, VN_m, EN_m;	// vertices of the collision mesh
Eigen::MatrixXi E_m;				// triangles = faces of the garment mesh / faces of the collision mesh
Eigen::VectorXi EMAP_m;
Eigen::VectorXd w; // the particle weights
Eigen::MatrixXd p; // the proposed new positions
MatrixXd u1, u2; // precomputation for stretch
PositionBasedDynamics PBD, PBD_adaption;

// colouring
MatrixXd colU, colJacDiff ;
MatrixXd colMixed;
MatrixXd colV ;
VectorXd normU, normV;
MatrixXd perFaceU, perFaceV;
int whichStressVisualize= 0;

garment_adaption* gar_adapt;
vector<seam*> seamsList;
vector<int> constrainedVertexIds;
VectorXi closestFaceId;

std::vector<std::pair<Eigen::Vector3d, int>> constrainedVertexBarycentricCoords;
std::vector<double> constrainedVertexDistance;
std::vector<std::pair<Eigen::Vector3d, int>> allVertexBarycentricCoords;
Eigen::MatrixXd tarU, tarV;// tarD1;
Eigen::MatrixXd baryCoords1, baryCoords2;
Eigen::MatrixXd baryCoordsd1, baryCoordsd2;
std::vector<std::vector<int> > boundaryL;


static bool noStress = true;
static bool StressU = false;
static bool StressV = false;
static bool StressDiffJac = false;
static bool StressJac = false;
int whichPatchMove=0;
std::vector<std::pair<double,double>> perFaceTargetNorm;
bool jacFlag=false;// not used anymore

// pattern adaption
bool adaptionFlag = false;
MatrixXd fromPattern, currPattern;
MatrixXd toPattern;
Eigen::MatrixXd p_adaption; // the proposed new positions
MatrixXd baryCoordsUPattern, baryCoordsVPattern;
vector<vector<pair<int, int>>> cornerPerBoundary;


MatrixXd patternPreInterpol,patternPreInterpol_temp ;
MatrixXd garmentPreInterpol,garmentPreInterpol_temp ;
MatrixXd mannequinPreInterpol, mannequinPreInterpol_temp;
Eigen::VectorXi componentIdPerFace,componentIdPerFaceNew, componentIdPerVert;

//test
//Eigen::SparseMatrix<double> L;
MatrixXd perFaceD2, perFaceD1;
VectorXd edgeVertices;

void computeBaryCoordsGarOnNewMannequin(igl::opengl::glfw::Viewer& viewer);
void setNewGarmentMesh(igl::opengl::glfw::Viewer& viewer);
void setNewMannequinMesh(igl::opengl::glfw::Viewer& viewer);
void showGarment(igl::opengl::glfw::Viewer& viewer);
void showMannequin(igl::opengl::glfw::Viewer& viewer);
void translateMesh(igl::opengl::glfw::Viewer& viewer, int which);

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers);
void reset(igl::opengl::glfw::Viewer& viewer);
void preComputeConstraintsForRestshape();
void dotimeStep(igl::opengl::glfw::Viewer& viewer);
void doAdaptionStep(igl::opengl::glfw::Viewer& viewer);
void setCollisionMesh();
void setupCollisionConstraints();
void computeBoundaryVertices();
void solveBendingConstraint();
void solveStretchConstraint();
void solveCollisionConstraint();
void preComputeStretch();
void computeStress(igl::opengl::glfw::Viewer& viewer);
void solveStretchUV();

// nice clicky interface
bool computePointOnMesh(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& Vuse, Eigen::MatrixXi& Fuse, Eigen::Vector3d& position, int& fid);
int computeClosestVertexOnMesh(Vector3d& b, int& fid, MatrixXi& F);
bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);

bool pre_draw(igl::opengl::glfw::Viewer& viewer){
    viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;
    if(simulate){
            Timer t("Draw loop");

            computeStress(viewer);
            t.printTime(" stress computed ");cout<<endl;
            dotimeStep(viewer);
            t.printTime(" timestep finished  ");cout<<endl;
            showGarment(viewer);// not sure if I actually need this, at least it breaks nothing
            t.printTime(" showing   ");cout<<endl;
            if(timestepCounter % convergeIterations == 10){
                gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, Vg_pattern,  seamsList, boundaryL);
                cout<<"after adaption"<<endl;
                preComputeConstraintsForRestshape();
                preComputeStretch();
                computeStress(viewer);
                t.printTime(" pattern computed ");cout<<endl;
            }
            timestepCounter++;
    }
     if (adaptionFlag){
//         cout<<" adaption flag on, arrived in loop "<<endl;

        doAdaptionStep(viewer);
    }

    return false;
}
void saveData(string fileName, VectorXi  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}
void saveData(string fileName, MatrixXi  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}
void saveData(string fileName, MatrixXd  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}
MatrixXi openData(string fileToOpen)
{

    // the inspiration for creating this function was drawn from here (I did NOT copy and paste the code)
    // https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix

    // the input is the file: "fileToOpen.csv":
    // a,b,c
    // d,e,f
    // This function converts input file data into the Eigen matrix format



    // the matrix entries are stored in this variable row-wise. For example if we have the matrix:
    // M=[a b c
    //    d e f]
    // the entries are stored as matrixEntries=[a,b,c,d,e,f], that is the variable "matrixEntries" is a row vector
    // later on, this vector is mapped into the Eigen matrix format
    vector<int> matrixEntries;

    // in this object we store the data from the matrix
    ifstream matrixDataFile(fileToOpen);

    // this variable is used to store the row of the matrix that contains commas
    string matrixRowString;

    // this variable is used to store the matrix entry;
    string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;


    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

        while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }

    // here we convet the vector variable into the matrix and return the resulting object,
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    return Map<Matrix<int, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}
MatrixXd openDataD(string fileToOpen)
{

    // the inspiration for creating this function was drawn from here (I did NOT copy and paste the code)
    // https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix

    // the input is the file: "fileToOpen.csv":
    // a,b,c
    // d,e,f
    // This function converts input file data into the Eigen matrix format



    // the matrix entries are stored in this variable row-wise. For example if we have the matrix:
    // M=[a b c
    //    d e f]
    // the entries are stored as matrixEntries=[a,b,c,d,e,f], that is the variable "matrixEntries" is a row vector
    // later on, this vector is mapped into the Eigen matrix format
    vector<double> matrixEntries;

    // in this object we store the data from the matrix
    ifstream matrixDataFile(fileToOpen);

    // this variable is used to store the row of the matrix that contains commas
    string matrixRowString;

    // this variable is used to store the matrix entry;
    string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;


    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

        while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }

    // here we convet the vector variable into the matrix and return the resulting object,
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}

void saveDataM(string fileName, vector<vector<int>>&  matrix)
{
    fstream file;
    file.open(fileName,ios_base::out);

    file<<matrix.size()<<endl;

    for(int i=0;i<matrix.size();i++)
    {
        file<<matrix[i].size()<<endl;
        for(int j=0; j< matrix[i].size(); j++){
            file<<matrix[i][j]<<endl;
        }
    }

    file.close();

}
void readDataM(string fileName, vector<vector<int>>& matrix ){
    fstream file;
    file.open(fileName,ios_base::in);
    int rows; file>> rows;
    matrix.resize(rows);
    for(int i=0; i<rows; i++){
        int cols; file>>cols;
        matrix[i].resize(cols);
        for(int j=0; j<cols; j++){
            file>>matrix[i][j];
        }
    }
}

int main(int argc, char *argv[])
{
    // Init the viewer
    igl::opengl::glfw::Viewer viewer;
    timestepCounter = 0;
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

    //string garment_file_name = igl::file_dialog_open();
//    string garment_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins/leggins_3d/leggins_3d_merged.obj"; //smaller collision thereshold to make sure it is not "eaten" after intiial step , 3.5 instead of 4.5
    string garment_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed3D_converged.obj";// smaller collision thereshold to make sure it is not "eaten" after intiial step , 3.5 instead of 4.5 is ok
//    string garment_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/dress_3d_lowres/dress_3d_lowres_merged_inlay.obj";// for the dress

    igl::readOBJ(garment_file_name, Vg, Fg);
    igl::readOBJ(garment_file_name, Vg_orig, Fg_orig);
    cout<<"loaded garment "<<endl;
    Timer t("Setup");
    garmentPreInterpol = Vg;
    Vg_orig = Vg; Fg_orig= Fg;

    string garment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins/leggins_2d/leggins_2d.obj"; //
//    string garment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/dress_2/dress_2d_lowres/dress_2d_lowres.obj"; //dress

    igl::readOBJ(garment_pattern_file_name, Vg_pattern, Fg_pattern);
    Vg_pattern_orig= Vg_pattern;
    patternPreInterpol= Vg_pattern;
    edgeVertices = VectorXd::Zero(Vg_pattern.rows());
    t.printTime(" init");
    preComputeConstraintsForRestshape();
    t.printTime(" preComputeConstraintsForRestshape");

    preComputeStretch();
    t.printTime( " preComputeStretch");
    jacFlag=false;

    setNewGarmentMesh(viewer);

// TODO remember to adapt the collision constraint solving dep on avatar, sometimes normalization is needed, sometimes not for whatever magic
    string avatar_file_name =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins/avatar/avatar_one_component.ply";
//    string avatar_file_name =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins_petite/avatar/avatar_one_component.ply";

    //string avatar_file_name = igl::file_dialog_open();
    igl::readPLY(avatar_file_name, Vm, Fm);
    Vm_orig = Vm; Fm_orig = Fm;
    mannequinPreInterpol = Vm;

    string morphBody1 =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins_petite/avatar/avatar_one_component.ply";
    igl::readPLY(morphBody1, testMorph_V1, testMorph_F1);

    if(Fm != testMorph_F1){
        cout<<"the faces are not the same!"<<endl;
    }

    setNewMannequinMesh(viewer);
    t.printTime( " set collison mesh ");

    std::map<int,int> vertexMapPattToGar;
    std::map<std::pair<int, int>,int> vertexMapGarAndIdToPatch;
    vertexMapPatternToGarment(Fg, Fg_pattern,vertexMapPattToGar);
    t.printTime( " vertexMapPatternToGarment ");

    igl::boundary_loop(Fg_pattern, boundaryL);
//    saveDataM("boundaryL_dress2_lowres.txt", boundaryL);
//    vector<vector<int>> newBoundaryL;
//    readDataM("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/boundaryL_dress2_lowres.txt", boundaryL);
//

    igl::facet_components(Fg_pattern, componentIdPerFace);
//    saveData("componentIdPerFace_dress2_lowres.csv", componentIdPerFace);
//    componentIdPerFace=  openData("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/componentIdPerFace_dress2_lowres.csv");

    igl::vertex_components(Fg_pattern, componentIdPerVert);
//    saveData("componentIdPerVert_dress2_lowres.csv", componentIdPerVert);
//    componentIdPerVert = openData("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/componentIdPerVert_dress2_lowres.csv");

    vertexMapGarmentAndPatchIdToPattern(Fg, Fg_pattern, componentIdPerVert, vertexMapGarAndIdToPatch);


    // use adjacentFacesToEdge of the 3D
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg, vfAdj);
    edgeVertices = VectorXd::Zero(Vg_pattern.rows());

    t.printTime( " before seams list  ");
    computeAllSeams( boundaryL,  vertexMapPattToGar, vertexMapGarAndIdToPatch, vfAdj, componentIdPerFace,
                     componentIdPerVert,edgeVertices, cornerPerBoundary,seamsList);
    t.printTime( " after seams list  ");

    gar_adapt = new garment_adaption(Vg, Fg,  Vg_pattern, Fg_pattern, seamsList, boundaryL); //none have been altered at this stage
    t.printTime( " garment init  ");
    gar_adapt->computeJacobian();
    t.printTime( " jacobian ");
    perFaceTargetNorm = gar_adapt->perFaceTargetNorm;
    Vg_orig = Vg;
    jacFlag = true;// not needed anymore...  was when we computed stress without reference jacobian

    // read constrained vertex ids and compute them as barycentric coordinates of the nearest face

    // save time
    setCollisionMesh();
//    col_tree.init(Vm, Fm);
//    t.printTime( " tree ");
//    FN_m = openDataD("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/FN_m_dress2.csv");
//    saveData("FN_m_dress2.csv", FN_m);
//    VN_m = openDataD("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/VN_m_dress2.csv");
//    saveData("VN_m_dress2.csv", VN_m);
//    EN_m =  openDataD("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/EN_m_dress2.csv");
//    saveData("EN_m_dress2.csv", EN_m);
//    E_m = openData("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/E_m_dress2.csv");
//    saveData("E_m_dress2.csv", E_m);
//    EMAP_m = openData("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/EMAP_m_dress2.csv");
//    saveData("EMAP_m_dress2.csv", EMAP_m);
// end save time

//    computeBoundaryVertices();
    t.printTime( " boundary ");
    computeBaryCoordsGarOnNewMannequin(viewer);// contains boundary vertices now
//    Vg = Vg_orig;
    t.printTime( " bary ");
    Vm = testMorph_V1;
    Vm_orig = testMorph_V1;
    showGarment(viewer);
    showMannequin(viewer);


    preComputeConstraintsForRestshape();
    preComputeStretch();
    computeStress(viewer);

//    setCollisionMesh();
    // save time
    setCollisionMesh();
//    col_tree.init(Vm, Fm);
//    t.printTime( " tree again  ");
//    FN_m = openDataD("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/FN_m_dress2_second.csv");
//    saveData("FN_m_dress2_second.csv", FN_m);
//    VN_m = openDataD("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/VN_m_dress2_second.csv");
//    saveData("VN_m_dress2_second.csv", VN_m);
//    EN_m =  openDataD("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/EN_m_dress2_second.csv");
//    saveData("EN_m_dress2_second.csv", EN_m);
//    E_m = openData("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/E_m_dress2_second.csv");
//    saveData("E_m_dress2_second.csv", E_m);
//    EMAP_m = openData("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/EMAP_m_dress2_second.csv");
//    saveData("EMAP_m_dress2_second.csv", EMAP_m);
// end save time


    viewer.core().animation_max_fps = 200.;
    viewer.core().is_animating = false;
    int whichSeam = 0;
    //additional menu items

    int whichBound=0;
    float movePatternX=0; float movePatternY=0;
    bool showPattern= false;
    menu.callback_draw_viewer_menu = [&]() {
        if (ImGui::CollapsingHeader("Garment", ImGuiTreeNodeFlags_OpenOnArrow)) {

            ImGui::InputFloat("Translation X", &garment_translation[0], 0, 0, "%0.4f");
            ImGui::InputFloat("Translation Y", &garment_translation[1], 0, 0, "%0.4f");
            ImGui::InputFloat("Translation Z", &garment_translation[2], 0, 0, "%0.4f");
            ImGui::InputFloat("Scaling factor X", &garment_scale, 0, 0, "%0.4f");
            if(ImGui::Button("Adjust garment", ImVec2(-1, 0))){
                translateMesh(viewer, 1 );
            }

            ImGui::InputInt("Vis Seam No", &whichSeam, 0, 0);
            if(ImGui::Button("Visualize Seam", ImVec2(-1, 0))){
                MatrixXd testCol= MatrixXd::Zero(Vg_pattern.rows(), 3);
                if(jacFlag){
                    std::vector<std::vector<int> > boundaryL;
                    igl::boundary_loop(Fg_pattern, boundaryL);
                    // testCol.col(0)= edgeVertices;
                    for(int j=whichSeam; j<whichSeam+1; j++){
                        seam* firstSeam = seamsList[j];
                        auto stP1 = firstSeam-> getStartAndPatch1();
                        auto stP2 = firstSeam-> getStartAndPatch2ForCorres();

                        int len = firstSeam -> seamLength();
                        int boundLen1 = boundaryL[stP1.second].size();
                        int boundLen2 = boundaryL[stP2.second].size();

                        /*
                int firstSide = boundaryL[stP1.second][(stP1.first + seamVert) % boundLen1];
                // what is the start of one side is the end for the other. (different traversal direction). ensure the counter does not get negative
                int secAccess = (stP2.first - seamVert) % boundLen2;
                if(secAccess < 0){
                    secAccess = boundLen2 + (stP2.first - seamVert);
                }
                if(seamsList[j]->inverted) secAccess = (stP2.first + seamVert) % boundLen2;

                int secondSide = boundaryL[stP2.second][secAccess];
                         * */
                        for(int i=0; i<=len; i++){
                            testCol(boundaryL[stP1.second][(stP1.first+i)% boundLen1],0) = 1.;
                            if(i==0)testCol(boundaryL[stP1.second][(stP1.first+i)% boundLen1],1) = 1.;
                            if(i==len)testCol(boundaryL[stP1.second][(stP1.first+i)% boundLen1],2) = 1.;

                            int setAccess = (stP2.first-i)% boundLen2;
                            if(setAccess < 0) {
                                setAccess +=boundLen2;
//                                testCol(boundaryL[stP2.second][(stP2.first-i)% boundLen2], 0) = 1.;
                            }
                            if(seamsList[j]->inverted) setAccess = (stP2.first + i) % boundLen2;
                            testCol(boundaryL[stP2.second][setAccess], 0) = 1.;

                        }
                    }

                }
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(Vg_pattern, Fg_pattern);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                //remove wireframe
                viewer.data().show_lines = false;
                // if 0 -> no face colour
                viewer.data().set_colors(testCol);
            }
            if(ImGui::Button("Visualize Corner", ImVec2(-1, 0))){
                MatrixXd testCol= MatrixXd::Zero(Vg_pattern.rows(), 3);
                for(int i=0; i<Vg_pattern.rows(); i++){
                    if(edgeVertices(i)){
                        testCol(i,0)=1;
                    }
                }
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(Vg_pattern, Fg_pattern);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                //remove wireframe
                viewer.data().show_lines = false;
                // if 0 -> no face colour
                viewer.data().set_colors(testCol);
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
        if (ImGui::CollapsingHeader("Pattern Computation", ImGuiTreeNodeFlags_DefaultOpen)) {
            if(ImGui::Checkbox("Show Pattern", &showPattern)){

                viewer.selected_data_index = 0;
                viewer.data().clear();
                mouse_mode = SELECTPATCH;
                if(showPattern){
                    viewer.data().set_mesh(Vg_pattern, Fg_pattern);
                }else{
                    viewer.data().set_mesh(Vg, Fg);
                }
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                //remove wireframe
                viewer.data().show_lines = true;

                // TODO we want the chance to click on a patch and make it constrained

            }
            ImGui::InputInt("Number of local global iterations", &(localGlobalIterations),  0, 0);
            // same as key down with key =='P'
            if (ImGui::CollapsingHeader("Alter Pattern ", ImGuiTreeNodeFlags_OpenOnArrow)){
                ImGui::InputInt("Which patch to move ", &(whichPatchMove),  0, 0);
                ImGui::InputFloat("Translation X", & movePatternX , 0, 0, "%0.4f");
                ImGui::InputFloat("Translation Y", & movePatternY , 0, 0, "%0.4f");

                if(ImGui::Button("Move ", ImVec2(-1, 0))){
                    MatrixXd testCol= MatrixXd::Zero(Vg_pattern.rows(), 3);
                    for(int i=0; i<Vg_pattern.rows(); i++){
                        if( componentIdPerVert(i) == whichPatchMove){
                            Vg_pattern(i, 0) += movePatternX;
                            Vg_pattern(i, 1) += movePatternY;
                            testCol(i,0)=1;
                        }
                    }
                    viewer.selected_data_index = 0;
                    viewer.data().clear();
                    viewer.data().set_mesh(Vg_pattern, Fg_pattern);
                    viewer.data().uniform_colors(ambient, diffuse, specular);
                    viewer.data().show_texture = false;
                    viewer.data().set_face_based(false);
                    //remove wireframe
                    viewer.data().show_lines = false;
                    // if 0 -> no face colour
                    viewer.data().set_colors(testCol);

                    igl::writeOBJ("patternComputed_translated.obj", Vg_pattern, Fg_pattern);
                    cout<<"pattern written to *patternComputed_translated*"<<endl;
                }
                ImGui::InputInt("Which Bound vis  ", &(whichBound),  0, 0);
                if(ImGui::Button("Vis Bound ", ImVec2(-1, 0))){
                    MatrixXd testCol= MatrixXd::Zero(Vg_pattern.rows(), 3);
                    for(int i=0; i<boundaryL[whichBound].size(); i++){
                            testCol(boundaryL[whichBound][i],0)=1;
                    }
                    viewer.selected_data_index = 0;
                    viewer.data().clear();
                    viewer.data().set_mesh(Vg_pattern, Fg_pattern);
                    viewer.data().uniform_colors(ambient, diffuse, specular);
                    viewer.data().show_texture = false;
                    viewer.data().set_face_based(false);
                    //remove wireframe
                    viewer.data().show_lines = false;
                    // if 0 -> no face colour
                    viewer.data().set_colors(testCol);
                }

            }
            if(ImGui::Button("Compute pattern", ImVec2(-1, 0))){
                simulate = false;
                // we start computing the pattern for the current shape
                Eigen::MatrixXd computed_Vg_pattern;//= Vg;
                cout<<"start computing the pattern with "<<localGlobalIterations<<" local global iterations"<<endl;
                gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, computed_Vg_pattern, seamsList, boundaryL);
                igl::writeOBJ("patternComputed.obj",Vg_pattern, Fg_pattern);
                cout<<"pattern written to *patternComputed*"<<endl;
            }
            if(ImGui::Button("Compute and Visualize stress of new pattern", ImVec2(-1, 0))){
                simulate = false;
                // we start computing the pattern for the current shape
                Eigen::MatrixXd computed_Vg_pattern;//= Vg;
                cout<<"start computing the pattern with "<<localGlobalIterations<<" local global iterations"<<endl;
                gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, computed_Vg_pattern, seamsList, boundaryL);
                igl::writeOBJ("patternComputed.obj", computed_Vg_pattern, Fg_pattern);
                cout<<"pattern written to *patternComputed*"<<endl;

                Vg_pattern_orig = computed_Vg_pattern;
                Vg_pattern = computed_Vg_pattern;

                preComputeConstraintsForRestshape();
                preComputeStretch();
                computeStress(viewer);

            }
        }
        if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {

            ImGui::InputInt("Converge Iterations ", &(convergeIterations), 0, 0);

            ImGui::InputDouble("Step size", &(timestep),  0, 0, "%0.4f");
            ImGui::InputDouble("U Stretch Stiffness ", &(stretchStiffnessU),  0, 0, "%0.4f");
            ImGui::InputDouble("V Stretch Stiffness", &(stretchStiffnessV),  0, 0, "%0.4f");
            ImGui::InputDouble("Diag Stretch Stiffness= Shear", &(stretchStiffnessD),  0, 0, "%0.6f");
            ImGui::InputDouble("Boundary Stiffness", &(boundaryStiffness), 0, 0, "%0.4f");
            ImGui::InputDouble("Collision Stiffness", &(collisionStiffness),  0, 0, "%0.4f");
            ImGui::InputDouble("Bending Stiffness", &(bendingStiffness),  0, 0, "%0.4f");
            ImGui::InputDouble("Edge length Stiffness", &(edgeLengthStiffness), 0, 0, "%0.6f");
            ImGui::InputDouble("Gravity factor", &(gravityfact),  0, 0, "%0.6f");

            ImGui::InputDouble("Collision thereshold", &(coll_EPS),  0, 0, "%0.2f");
            ImGui::InputInt("Number of constraint Iterations thereshold", &(num_const_iterations),  0, 0);

            if(ImGui::Checkbox("Visualize no Stress", &noStress)){
                StressU = false;
                StressV=false;
                StressDiffJac = false;
                whichStressVisualize = 0;
                showGarment(viewer);
            }
            if(ImGui::Checkbox("Visualize U Stress", &StressU)){
                noStress = false;
                StressV = false;
                StressDiffJac = false;
                StressJac = false;
                whichStressVisualize = 1;
                showGarment(viewer);
            }
            if(ImGui::Checkbox("Visualize V Stress ", &StressV)){
                noStress = false;
                StressU = false;
                StressDiffJac = false;
                StressJac = false;
                whichStressVisualize = 2;
                showGarment(viewer);
           }
            if(ImGui::Checkbox("Visualize diffFrom Jacobian ", &StressDiffJac)){
                StressV= false;
                noStress = false;
                StressU = false;
                StressJac = false;
                whichStressVisualize = 3;
                showGarment(viewer);
            }
            if(ImGui::Checkbox("Visualize Jacobian V", &StressJac)){
                StressV= false;
                noStress = false;
                StressU = false;
                StressDiffJac = false;

                Eigen::MatrixXd colJac(Fg.rows(), 3);
              for(int i=0; i<Fg.rows(); i++){
                  auto jac = gar_adapt->jacobians[i];
                  double y = (jac.col(1).norm()-1) * 3; // to increase differences
                  colJac.row(i) = Vector3d(1.0 + y, 1. - y, 0.0);
              }
                viewer.selected_data_index = 0;
                viewer.data().set_colors(colJac);
            }
            static bool remMan;
            if(ImGui::Checkbox("Original View ", &remMan)){
                if(remMan){
                    patternPreInterpol_temp= Vg_pattern;
                    garmentPreInterpol_temp= Vg;
                    mannequinPreInterpol_temp = Vm;
                    Vg_pattern= patternPreInterpol;
                    Vg = garmentPreInterpol;
                    Vm = mannequinPreInterpol;
                }else{
                    Vg_pattern= patternPreInterpol_temp;
                    Vg = garmentPreInterpol_temp;
                    Vm = mannequinPreInterpol_temp;
                }
                preComputeConstraintsForRestshape();
                preComputeStretch();
                computeStress(viewer);
                showGarment(viewer);
                showMannequin(viewer);
            }
        }
        if (ImGui::CollapsingHeader("Pattern adaption", ImGuiTreeNodeFlags_DefaultOpen)){
            if(ImGui::Button("Compute adaptation", ImVec2(-1, 0))){
                simulate = false;
                stretchStiffnessU = 0.8;
                stretchStiffnessD *= 2;
                boundaryStiffness = 0.99;
                cout<<" request to start the pattern adaptation"<<endl;
                cout<<"begin by mapping the original to the computed pattern "<<endl;
                if(cornerPerBoundary.size()==0){
                    cout<<" there are no corners to map"<<endl;

                }
                if(Vg_pattern.rows()!= Vg_pattern_orig.rows()){
                    cout<<" the initial number of vertices does not match ! "<<endl;
                }


                // copy the matrices to not mess with them
                string fromPatternFile = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed.obj";
                // quick
                igl::readOBJ(fromPatternFile, fromPattern, Fg_pattern);

//                 fromPattern = Vg_pattern;NOOOO
                 currPattern = fromPattern;
                 toPattern= Vg_pattern_orig;

                // use with PBD_adaption
                baryCoordsUPattern.resize(Fg_pattern.rows(), 3);
                baryCoordsVPattern.resize(Fg_pattern.rows(), 3);
                for(int i=0; i<Fg_pattern.rows(); i++){
                    int id0 = Fg_pattern(i, 0);
                    int id1 = Fg_pattern(i, 1);
                    int id2 = Fg_pattern(i, 2);

                    Vector2d Gu, Gv, G;
                    Vector2d p0, p1, p2;
                    p0 = fromPattern.block(id0, 0, 1, 2).transpose();
                    p1 = fromPattern.block(id1, 0, 1, 2).transpose();
                    p2 = fromPattern.block(id2, 0, 1, 2).transpose();

                    G = (1./3.) * p0 + (1./3.) * p1 + (1./3.) * p2;

                    Gu = G; Gu(0) += 1;
                    Gv = G; Gv(1) += 1;
                    Vector3d uInBary, vInBary;
                    MathFunctions mathFun;
                    mathFun.Barycentric(Gu, p0, p1, p2, uInBary);
                    mathFun.Barycentric(Gv, p0, p1, p2, vInBary);

                    baryCoordsUPattern.row(i) = uInBary;
                    baryCoordsVPattern.row(i) = vInBary;

                }
                viewer.core().is_animating = true;
                adaptionFlag = true;
            }
            if(ImGui::Button("Compute first Tear", ImVec2(-1, 0))){
                simulate = false;

//later these matrices shouold be fixed from the pattern adaption , ,this is precomputation for simplification

                // copy the matrices to not mess with them, , // quick
                string fromPatternFile = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed.obj";
                igl::readOBJ(fromPatternFile, fromPattern, Fg_pattern);

                string mappedFile = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/mappedPattern.obj";
                igl::readOBJ(mappedFile, currPattern, Fg_pattern);

                toPattern= Vg_pattern_orig;

                // use with PBD_adaption
                baryCoordsUPattern.resize(Fg_pattern.rows(), 3);
                baryCoordsVPattern.resize(Fg_pattern.rows(), 3);
                for(int i=0; i<Fg_pattern.rows(); i++){
                    int id0 = Fg_pattern(i, 0);
                    int id1 = Fg_pattern(i, 1);
                    int id2 = Fg_pattern(i, 2);

                    Vector2d Gu, Gv, G;
                    Vector2d p0, p1, p2;
                    p0 = fromPattern.block(id0, 0, 1, 2).transpose();
                    p1 = fromPattern.block(id1, 0, 1, 2).transpose();
                    p2 = fromPattern.block(id2, 0, 1, 2).transpose();

                    G = (1./3.) * p0 + (1./3.) * p1 + (1./3.) * p2;

                    Gu = G; Gu(0) += 1;
                    Gv = G; Gv(1) += 1;
                    Vector3d uInBary, vInBary;
                    MathFunctions mathFun;
                    mathFun.Barycentric(Gu, p0, p1, p2, uInBary);
                    mathFun.Barycentric(Gv, p0, p1, p2, vInBary);

                    baryCoordsUPattern.row(i) = uInBary;
                    baryCoordsVPattern.row(i) = vInBary;

                }
                viewer.core().is_animating = true;
                adaptionFlag = true;
            }

        }
        menu.draw_viewer_menu();
    };

    // Add content to the default menu window
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &callback_key_down;

    viewer.selected_data_index = 0;
    viewer.callback_mouse_down = &callback_mouse_down;

    t.printTime( " fin ");
    viewer.launch();
}

bool computePointOnMesh(igl::opengl::glfw::Viewer& viewer, MatrixXd& V, MatrixXi& F, Vector3d& b, int& fid) {
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    return igl::unproject_onto_mesh(Vector2f(x, y), viewer.core().view, viewer.core().proj, viewer.core().viewport, V, F, fid, b);
}
int computeClosestVertexOnMesh(Vector3d& b, int& fid, MatrixXi& F) {
    // get the closest vertex in that face
    int v_id;
    if (b(0) > b(1) && b(0) > b(2))
        v_id = F(fid, 0);
    else if (b(1) > b(0) && b(1) > b(2))
        v_id = F(fid, 1);
    else
        v_id = F(fid, 2);
    return v_id;
}
bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier){
    if(mouse_mode == SELECTBOUNDARY || mouse_mode == SELECTPATCH ){
        if (button == (int)igl::opengl::glfw::Viewer::MouseButton::Right)
            return false;

        int fid;
        Eigen::Vector3d b;
        MatrixXd Vrs = Vg_pattern;

        if (computePointOnMesh(viewer, Vrs, Fg_pattern, b, fid)) {
            int v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern);
//        cout<<v_id<<"computed closest pattern id "<<endl ;
            viewer.data().set_points(Vrs.row(v_id), RowVector3d(1.0, 0.0, 0.0));
            whichPatchMove = componentIdPerVert(v_id);
            // TODO now we could constrain the whole patch or the boundary
            return true;
        }
        return false;
    }
    return false;
}
void computeBaryCoordsGarOnNewMannequin(igl::opengl::glfw::Viewer& viewer){
    VectorXd S;
    VectorXd distVec(Vg.rows());

    constrainedVertexIds.clear();
    vector<vector<int> > vvAdj, vfAdj;
    igl::adjacency_list(Fg,vvAdj);
    createVertexFaceAdjacencyList(Fg, vfAdj);
    int boundarycount = 0;

    igl::signed_distance_pseudonormal(Vg, Vm, Fm, col_tree, FN_m, VN_m, EN_m, EMAP_m, S, closestFaceId, C, N);
    for(int i=0; i<Vg.rows(); i++){
        int closestFace = closestFaceId(i);
        Vector3d a = Vm.row(Fm(closestFace, 0));
        Vector3d b = Vm.row(Fm(closestFace, 1));
        Vector3d c = Vm.row(Fm(closestFace, 2));

        Vector3d bary = (a+b+c)/3;
        Vector3d vvec = Vg.row(i).transpose() - bary;
        Vector3d normalVec = (b-a).cross(c-a);
        normalVec = normalVec.normalized();
        distVec(i) = vvec.dot(normalVec);

        Vector3d currVert = Vg.row(i).transpose()- distVec(i)*normalVec;

        Vector3d currInBary;
        MathFunctions mathFun;
        mathFun.Barycentric3D(currVert, a, b, c, currInBary);
        if(isBoundaryVertex(Vg, i, vvAdj, vfAdj)){
            constrainedVertexIds.emplace_back(i); // (i)= 1;
            boundarycount++;
            constrainedVertexBarycentricCoords.emplace_back(std::make_pair(currInBary, closestFace));
            constrainedVertexDistance.push_back(distVec(i));
        }
        allVertexBarycentricCoords.emplace_back(std::make_pair(currInBary, closestFace));

        a = testMorph_V1.row(Fm(closestFace, 0));
        b = testMorph_V1.row(Fm(closestFace, 1));
        c = testMorph_V1.row(Fm(closestFace, 2));
        Vector3d newPos = currInBary(0) * a + currInBary(1) * b + currInBary(2) * c;
        normalVec = (b-a).cross(c-a);
        normalVec = normalVec.normalized();
        Vg.row(i) = newPos + distVec(i) * normalVec;
    }

}
void setNewGarmentMesh(igl::opengl::glfw::Viewer& viewer) {
    if (Vg.rows() == 0 || Fg.rows() == 0) {
        fprintf(stderr, "IOError: Could not load garment...\n");
        return;
    }
    igl::edges(Fg, Eg);
    showGarment(viewer);
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
        viewer.data().set_colors(colJacDiff);
    }

}
void setNewMannequinMesh(igl::opengl::glfw::Viewer& viewer) {
    if (Vm.rows() == 0 || Fm.rows() == 0) {
        fprintf(stderr, "IOError: Could not load garment...\n");
        return;
    }
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
        Vg_pattern = Vg_pattern_orig * garment_scale;
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
        setNewMannequinMesh(viewer);
        cout<<" setting collision mesh "<<endl;
        setCollisionMesh();
        cout<<" collision mesh finished "<<endl;
    }
}
bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers) {
    bool keyRecognition= false;

    if(key == 'A'){
        adaptionFlag = !adaptionFlag ;
        keyRecognition = true;
        viewer.core().is_animating = !viewer.core().is_animating;
    }
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
        StressDiffJac = true;
        showGarment(viewer);
        keyRecognition = true;
    }
    if(key== 'W'){
        keyRecognition=true;
        igl::writeOBJ("writtenPattern"+ to_string(garment_scale)+".obj", Vg_pattern, Fg_pattern);
          std::cout<<" Garment file written"<<endl;
    }
    if(key == 'P'){       // Pattern
        keyRecognition = true;
        simulate = false;
        // we start computing the pattern for the current shape
        Eigen::MatrixXd computed_Vg_pattern= Vg;
        cout<<"start computing the pattern with "<<localGlobalIterations<<" local global iterations"<<endl;
        gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, computed_Vg_pattern, seamsList, boundaryL);

        igl::writeOBJ("patternComputed.obj", computed_Vg_pattern, Fg_pattern);
        igl::writeOBJ("patternComputed3D.obj", Vg, Fg);

        cout<<"pattern written to *patternComputed*"<<endl;

    }
    if(key == 'B'){     // Bending
        cout<<" writing mapped pattern "<<endl;
        igl::writeOBJ("mappedPattern.obj", currPattern, Fg_pattern);
       // reset(viewer);
//        simulate= false;
//        Fm = testMorph_F0;
//        Vm = testMorph_V0;
//        cout<<Vm.rows()<<" and faces "<<Fm.rows()<<endl;
//        viewer.selected_data_index = 0;
//        viewer.data().clear();
 //       setNewMannequinMesh(viewer);
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
    timestepCounter = 0;
    Vg= Vg_orig;
    Vg_pattern = Vg_pattern_orig;
    vel = Eigen::MatrixXd::Zero(numVert, 3);
    Vm= Vm_orig;
    Fm= Fm_orig;
    viewer.core().is_animating = false;
    simulate=false;
    StressV = false;
    StressU = false;
    StressDiffJac = false;
    StressJac = false;
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
    // the pattern has fewer adjacent faces since the stitching does not count here, but it does in the 3D case

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

    col_tree.init(Vm, Fm);
    igl::per_face_normals(Vm, Fm, FN_m);
    igl::per_vertex_normals(Vm, Fm, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_m, VN_m);
    igl::per_edge_normals(Vm, Fm, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN_m, EN_m, E_m, EMAP_m);

}
void setupCollisionConstraints(){
    VectorXd S;
    //MatrixXd C, N;
    collisionVert = Eigen::MatrixXi::Zero(numVert, 1);
    pureCollVert.clear();

    igl::signed_distance_pseudonormal(p, Vm, Fm, col_tree, FN_m, VN_m, EN_m, EMAP_m, S, closestFaceId, C, N);
    int collCount=0;
    for(int i=0; i<numVert; i++){
        if(S(i)<coll_EPS){
            collCount++;
            collisionVert(i)=1;
            pureCollVert.push_back(i);

        }
    }
    if(pureCollVert.size()!= collCount){
        cout<<" size problem"<<endl;
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

//  UV DIFFERENTIATION EXISTS
void solveStretchConstraint(){
    /*each edges distance should remain, since we iterate over every face we iterate over every edge twice- but that should not be a problem */
    for (int j =0; j<numFace; j++){
        Vector3r deltap0, deltap1, deltap2;

        int id0 = Fg_orig(j, 0);
        int id1 = Fg_orig(j, 1);
        int id2 = Fg_orig(j, 2);

        // first edge
        PBD.solve_DistanceConstraint(p.row(id1), w(id1), p.row(id2), w(id2), edgeLengths(j, 0), edgeLengthStiffness, deltap1, deltap2);
        p.row(id2) += deltap2;
        p.row(id1) += deltap1;

        //second edge
        PBD.solve_DistanceConstraint(p.row(id2), w(id2), p.row(id0), w(id0), edgeLengths(j, 1), edgeLengthStiffness, deltap2, deltap0);
        p.row(id2)+= deltap2;
        p.row(id0)+= deltap0;

        // third edge
        PBD.solve_DistanceConstraint(p.row(id0), w(id0), p.row(id1), w(id1), edgeLengths(j, 2), edgeLengthStiffness, deltap0, deltap1);
        p.row(id1)+= deltap1;
        p.row(id0)+= deltap0;
    }
}
void init_stretchUV(){
    tarU.resize(3*numFace, 3);
    tarV.resize(3*numFace, 3);
    for(int j = 0; j<numFace; j++){
        Eigen::MatrixXd patternCoords(2, 3);
        patternCoords(0,0) = Vg_pattern( Fg_pattern(j, 0), 0);
        patternCoords(1,0) = Vg_pattern( Fg_pattern(j, 0), 1);
        patternCoords(0,1) = Vg_pattern( Fg_pattern(j, 1), 0);
        patternCoords(1,1) = Vg_pattern( Fg_pattern(j, 1), 1);
        patternCoords(0,2) = Vg_pattern( Fg_pattern(j, 2), 0);
        patternCoords(1,2) = Vg_pattern( Fg_pattern(j, 2), 1);

        Eigen::MatrixXd targetPositions(3, 3);
        targetPositions.col(0)= p.row(Fg_orig(j, 0));
        targetPositions.col(1)= p.row(Fg_orig(j, 1));
        targetPositions.col(2)= p.row(Fg_orig(j, 2));

        Vector3r tar0, tar1, tar2;

        PBD.init_UVStretch(perFaceU.row(j), perFaceV.row(j), patternCoords, targetPositions, tar0, tar1, tar2, 1,stretchStiffnessD );
        tarU.row(3*j) = tar0.transpose();
        tarU.row(3*j+1) = tar1.transpose();
        tarU.row(3*j+2) = tar2.transpose();

        PBD.init_UVStretch( perFaceU.row(j),perFaceV.row(j), patternCoords, targetPositions,tar0, tar1, tar2, 2, stretchStiffnessD);
        tarV.row(3*j) = tar0.transpose();
        tarV.row(3*j+1) = tar1.transpose();
        tarV.row(3*j+2) = tar2.transpose();

    }
}
void solveStretchUV(){
    for (int j =0; j<numFace; j++){

        Vector3r deltap0, deltap1, deltap2;
        Eigen::MatrixXd targetPositions(3, 3);
        targetPositions.col(0)= p.row(Fg_orig(j, 0));
        targetPositions.col(1)= p.row(Fg_orig(j, 1));
        targetPositions.col(2)= p.row(Fg_orig(j, 2));

        Vector3r tar0= tarU.row(3*j +0);
        Vector3r tar1= tarU.row(3*j +1);
        Vector3r tar2= tarU.row(3*j +2);
        targetPositions.col(0)= p.row(Fg_orig(j, 0));
        targetPositions.col(1)= p.row(Fg_orig(j, 1));
        targetPositions.col(2)= p.row(Fg_orig(j, 2));
        PBD.solve_Stretch(targetPositions,
                          tar0, tar1, tar2, stretchStiffnessU, deltap0, deltap1, deltap2);

            // only update if a vertex is not constrained!
        p.row(Fg_orig(j, 2)) += deltap2;
        p.row(Fg_orig(j, 1)) += deltap1;
        p.row(Fg_orig(j, 0)) += deltap0;

        tar0= tarV.row(3 * j + 0);
        tar1= tarV.row(3*j +1);
        tar2= tarV.row(3*j +2);
        targetPositions.col(0)= p.row(Fg_orig(j, 0));
        targetPositions.col(1)= p.row(Fg_orig(j, 1));
        targetPositions.col(2)= p.row(Fg_orig(j, 2));
        PBD.solve_Stretch(targetPositions, tar0, tar1, tar2,
                          stretchStiffnessV, deltap0, deltap1, deltap2);

        p.row(Fg_orig(j, 2)) += deltap2;
        p.row(Fg_orig(j, 1)) += deltap1;
        p.row(Fg_orig(j, 0)) += deltap0;
    }
}
void solveCollisionConstraint(){
    for(int i=0; i<pureCollVert.size(); i++){
        int j = pureCollVert[i];
        Vector3r deltap0;
        // maybe I should compute the intersection instead of using the closest point C?
        PBD.solve_CollisionConstraint(p.row(j),  C.row(j), N.row(j), deltap0, coll_EPS, vel.row(j));
        p.row(j) += collisionStiffness * deltap0;
    }
}

void preComputeStretch(){
    colU = Eigen::MatrixXd ::Zero(numFace, 3);
    colJacDiff = Eigen::MatrixXd ::Zero(numFace, 3);

    colV = Eigen::MatrixXd ::Zero (numFace, 3);
    colMixed = Eigen::MatrixXd ::Zero (numFace, 3);
    baryCoords1.resize(numFace, 3); // for Gu and Gv
    baryCoords2.resize(numFace, 3); // for Gu and Gv

    baryCoordsd1.resize(numFace, 3); //diags
    baryCoordsd2.resize(numFace, 3);

    u1.resize(numFace, 2);
    u2.resize(numFace, 2);
    for(int j=0; j<numFace; j++){
        int id0 = Fg_pattern(j, 0);
        int id1 = Fg_pattern(j, 1);
        int id2 = Fg_pattern(j, 2);

        Vector2d  gU, gV,centralG;

        // new part: we compute the barycentric coordinates for the vectors u, v. They are the reference for the stress
        centralG(0) = (Vg_pattern(id0, 0) + Vg_pattern(id1, 0) + Vg_pattern(id2, 0)) / 3.;
        centralG(1) = (Vg_pattern(id0, 1) + Vg_pattern(id1, 1) + Vg_pattern(id2, 1)) / 3.;
        gU = centralG;
        gV = centralG;

        gU(0) += 1;
        gV (1) += 1;

//        double det = gU( 0) * gV(1) - (gV(0)*gU(1));// 90 deg, should be 0 - no they are not vectors but pints
        u1(j,0)= 1; u1(j, 1)= 0;//gU/det;
        u2(j, 0)= 0; u2(j, 1) = 1; //gV/det;

        MathFunctions mathFun;
        Vector2d p0, p1, p2;
        p0(0) = Vg_pattern(id0, 0);
        p0(1) = Vg_pattern(id0, 1);

        p1(0) = Vg_pattern(id1, 0);
        p1(1) = Vg_pattern(id1, 1);

        p2(0) = Vg_pattern(id2, 0);
        p2(1) = Vg_pattern(id2, 1);
        Vector3d u1InBary, u2InBary;
        Vector3d d1InBary, d2InBary;

        mathFun.Barycentric(gU, p0, p1, p2, u1InBary);
        mathFun.Barycentric(gV, p0, p1, p2, u2InBary);
//
        baryCoords1.row(j) = u1InBary;
        baryCoords2.row(j) = u2InBary;
//
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

        double differenceIncrementFactor = 3.0;

        // deviation from 1 as the measure,
        /* large u stretch: norm > 1, thus y>0 , thus very red,little green => red
            compression = small u stretch: norm < 1, thus y<0 , thus little red, much green => green
            no stretch : y=0, thus very red , very green => yellow */
        Vector3d Gu, G, Gv;
        Vector3d Gd1, Gd2;

        Gu= baryCoords1(j, 0)*Vg.row(id0) + baryCoords1(j, 1)*Vg.row(id1) + baryCoords1(j, 2)*Vg.row(id2);
        Gv= baryCoords2(j, 0)*Vg.row(id0) + baryCoords2(j, 1)*Vg.row(id1) + baryCoords2(j, 2)*Vg.row(id2);
        G = (1./3)*Vg.row(id0) +(1./3)*Vg.row(id1) + (1./3)*Vg.row(id2);


        perFaceU.row(j) = (Gu-G);//*u2(j,1) - (Gv-G)*u1(j,1);
        perFaceV.row(j) = (Gv-G);// * u1(j, 0) - (Gu-G)* u2(j, 0);

        normU(j)= (Gu-G).norm();
        double y = (normU(j)-1) * differenceIncrementFactor; // to increase differences
        colU.row(j) = Vector3d(1.0 + y, 1. - y, 0.0);

        normV(j) = (Gv-G).norm();
        y = (normV(j)-1) * differenceIncrementFactor;
        colV.row(j) = Vector3d ( 1. + y, 1.- y, 0.0);

        if(jacFlag){
            double diffU = (normU(j)-perFaceTargetNorm[j].first)/ perFaceTargetNorm[j].first;
            double diffV = (normV(j)-perFaceTargetNorm[j].second)/ perFaceTargetNorm[j].second;
           y = diffU + diffV ;
           y*= 3; // to better see the difference
            colJacDiff.row(j)=  Vector3d ( 1. + y, 1.- y, 0.0);
        }else{
            colJacDiff.row(j)=  Vector3d ( 1. , 1., 0.0);
        }

        // this is an experiment
        y = (abs(normV(j)-1)+ abs(normU(j)-1))*3;
        colMixed.row(j) = Vector3d ( 1. + y, 1.- y, 0.0);
    }
}
void solveConstrainedVertices(){
    for(int i=0; i<constrainedVertexIds.size(); i++){

            int closestFace = get<1> (constrainedVertexBarycentricCoords[i]);
            Vector3d baryCoeff = get<0>(constrainedVertexBarycentricCoords[i]);
            Vector3d newSuggestedPos = baryCoeff(0)* Vm.row(Fm(closestFace, 0));
            newSuggestedPos += baryCoeff(1)*Vm.row(Fm(closestFace, 1));
            newSuggestedPos += baryCoeff(2)*Vm.row(Fm(closestFace, 2));

            // account for offset from body!otherwise it alternates between this and collision force .  works smoothly so far,
            Vector3d e1 = Vm.row(Fm(closestFace, 1)) - Vm.row(Fm(closestFace, 0));
            Vector3d e2 = Vm.row(Fm(closestFace, 2)) - Vm.row(Fm(closestFace, 0));
            Vector3d normal = e1.cross(e2);
            normal = normal.normalized();
            newSuggestedPos += constrainedVertexDistance[i]* normal;

            Vector3d dir = newSuggestedPos - p.row(constrainedVertexIds[i]).transpose();

            p.row(constrainedVertexIds[i])+= boundaryStiffness * dir;

    }
}
void solveStretchAdaption(MatrixXd& perFaceU_adapt,MatrixXd& perFaceV_adapt){
    // force that pulls back to the original position in fromPattern
    for(int i=0; i< Fg_pattern.rows(); i++){
        Eigen::MatrixXd patternCoords(2, 3);
        patternCoords.col(0) = fromPattern.row(Fg_pattern(i, 0)).leftCols(2).transpose();
        patternCoords.col(1) = fromPattern.row(Fg_pattern(i, 1)).leftCols(2).transpose();
        patternCoords.col(2) = fromPattern.row(Fg_pattern(i, 2)).leftCols(2).transpose();

        // where they would go to if no stretch in u
        Vector2r tarUV0, tarUV1 , tarUV2;

        Eigen::MatrixXd targetPositions(2, 3);
        targetPositions.col(0)=  p_adaption.row(Fg_pattern(i, 0)).leftCols(2).transpose() ;
        targetPositions.col(1)=  p_adaption.row(Fg_pattern(i, 1)).leftCols(2).transpose() ;
        targetPositions.col(2)=  p_adaption.row(Fg_pattern(i, 2)).leftCols(2).transpose() ;

        int uOrv = 1;
//        TODO STIFFNESS PARAMETER
        PBD_adaption.init_UVStretchPattern( perFaceU_adapt.row(i),  perFaceV_adapt.row(i), patternCoords,targetPositions,
                                                tarUV0, tarUV1,tarUV2, uOrv,  stretchStiffnessD);

        Vector2d dir0 = tarUV0 - p_adaption.row(Fg_pattern(i, 0)).leftCols(2).transpose() ;
        Vector2d dir1 = tarUV1 - p_adaption.row(Fg_pattern(i, 1)).leftCols(2).transpose() ;
        Vector2d dir2 = tarUV2 - p_adaption.row(Fg_pattern(i, 2)).leftCols(2).transpose() ;
//  TODO STIFFNESS PARAMETER
        p_adaption.row(Fg_pattern(i,0)).leftCols(2) += ( stretchStiffnessU * dir0);
        p_adaption.row(Fg_pattern(i,1)).leftCols(2) += ( stretchStiffnessU * dir1);
        p_adaption.row(Fg_pattern(i,2)).leftCols(2) += ( stretchStiffnessU * dir2);
//
    }
}
MatrixXd colPatternU, colPatternV;
void computePatternStress(MatrixXd& perFaceU_adapt,MatrixXd& perFaceV_adapt ){
    // we compute the stress between the current pattern and the one we started from
    colPatternU.resize(Fg_pattern.rows(), 3);
    colPatternV.resize(Fg_pattern.rows(), 3);

    for(int j=0; j<Fg_pattern.rows(); j++) {
        int id0 = Fg_pattern(j, 0);
        Vector2d p0 = p_adaption.block(id0, 0, 1, 2).transpose();
        int id1 = Fg_pattern(j, 1);
        Vector2d p1 = p_adaption.row(id1).leftCols(2);
        int id2 = Fg_pattern(j, 2);
        Vector2d p2 = p_adaption.block(id2, 0, 1, 2).transpose();

        Vector2d Gu, Gv, G;
        G = (1. / 3.) * p0 + (1. / 3.) * p1 + (1. / 3.) * p2;
        Gu = baryCoordsUPattern(j, 0) * p0 + baryCoordsUPattern(j, 1) * p1 + baryCoordsUPattern(j, 2) * p2;
        Gv = baryCoordsVPattern(j, 0) * p0 + baryCoordsVPattern(j, 1) * p1 + baryCoordsVPattern(j, 2) * p2;

        perFaceU_adapt.row(j) = (Gu - G).transpose();
        perFaceV_adapt.row(j) = (Gv - G).transpose();

        double normU= perFaceU_adapt.row(j).norm();
        double normV = perFaceV_adapt.row(j).norm();

        double differenceIncrementFactor = 5.;
        double y = (normU-1) * differenceIncrementFactor; // to increase differences
        double yV = (normV-1) * differenceIncrementFactor; // to increase differences

        colPatternU.row(j) = Vector3d((1.0 + y), (1. - y), 0.0);
        colPatternV.row(j) = Vector3d((1.0 + yV), (1. - yV), 0.0);

    }

}
void solveCornerMappedVertices(){

    for(int i=0; i< cornerPerBoundary.size(); i++){
        for(int j=0; j< cornerPerBoundary[i].size(); j++){
            int vertIdx = get<0>(cornerPerBoundary[i][j]);
            Vector2d newSuggestedPos = toPattern.row(vertIdx).leftCols(2);
            Vector2d dir = newSuggestedPos - p_adaption.row(vertIdx).leftCols(2).transpose();
// TODO PARAMETER
            p_adaption.row(vertIdx).leftCols(2) += boundaryStiffness * dir;
        }
    }

}
int adaptioncount = 0 ;
void doAdaptionStep(igl::opengl::glfw::Viewer& viewer){
    Timer t(" Adaption time step ");

    adaptioncount++;

    std::cout<<"-------------- Time Step ------------"<<adaptioncount<<endl;
    // we have no velocity or collision
    // but we do have stretch, constrainedStretch and bending
    // for the stretch, the current pattern is the reference, then its corners are mapped to another position
// the stretch as a simple solve stretch of the rest position and the stretched position?
// constrained corners as constrained condition: new suggested position= mapped position, then the direction is added to p
// bending: not needed if we don;t touch the z coordinate

    p_adaption.resize(currPattern.rows(), 3);
    p_adaption = currPattern;
    p_adaption.col(2)= Eigen::VectorXd::Ones(currPattern.rows());
    p_adaption.col(2)*= 200;

    MatrixXd perFaceU_adapt(Fg_pattern.rows(), 2);
    MatrixXd perFaceV_adapt(Fg_pattern.rows(), 2);

    t.printTime(" init ");
    for(int i=0; i<5; i++){
//        t.printTime(" iteration");
        solveCornerMappedVertices();
        // now we treat the stretch
        computePatternStress(perFaceU_adapt, perFaceV_adapt);
        solveStretchAdaption(perFaceU_adapt, perFaceV_adapt);

    }
    currPattern = p_adaption;

    viewer.selected_data_index = 1;
    viewer.data().clear();
    viewer.selected_data_index = 0;
    viewer.data().clear();
    viewer.data().set_mesh(currPattern, Fg_pattern);
    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_texture = false;
    viewer.data().set_face_based(false);
    //remove wireframe
    viewer.data().show_lines = true;

//    MatrixXd patternCol ( Fg_pattern.rows(), 3);
//    double differenceIncrementFactor = 5;
//    double maxnorm  = -1;
//    double maxnormV = -1;
//    for(int i=0; i<Fg_pattern.rows(); i++){
//
//        double normU= perFaceU_adapt.row(i).norm();
//        double normV = perFaceV_adapt.row(i).norm();
//        maxnorm = max(maxnorm ,normU);
//        maxnormV = max(maxnormV, normV);
////        double y = (normU-1) * differenceIncrementFactor; // to increase differences
////        double yV = (normV-1) * differenceIncrementFactor; // to increase differences
//
////        patternCol.row(i) = Vector3d(min(1.0 + y, 1.), max(1. - y, 0.), 0.0);
////        patternCol.row(i) = Vector3d(min(1.0 + yV, 1.), max(1. - yV, 0.), 0.0);
//
//    }
    cout<<(colPatternV.col(0).maxCoeff()-1)/5 +1 <<" max norm v direction "<<endl;
    cout<<(colPatternU.col(0).maxCoeff()-1)/5 +1<<" max norm in u direction"<<endl;
    viewer.data().set_colors(colPatternV);

}
void dotimeStep(igl::opengl::glfw::Viewer& viewer){
    Timer t("Time step ");
    std::cout<<endl;
    std::cout<<"-------------- Time Step ------------"<<timestepCounter<<endl;
    // the stress is already computed, we can use it here
    Eigen::MatrixXd x_new = Vg;
    p = Vg;
    // line (5) of the algorithm https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
    // we only use it to add gravity to the system
    vel.col(1) += timestep * w*(-1)*grav*gravityfact;

    // (7)
    for (int i = 0; i<numVert; i++){
        p.row(i) = x_new.row(i).array()+ timestep*vel.row(i).array();
    }

    // detect collisions and solve for them in the loop
    setupCollisionConstraints();
    t.printTime(" setup collision constraints ");

    init_stretchUV();
    t.printTime(" setup uv stretch ");
    //(9)-(11), the loop should be repeated several times per timestep (according to Jan Bender)
    for(int i = 0; i < num_const_iterations; i++){
            solveBendingConstraint();
            solveStretchConstraint();
            solveStretchUV();
            solveConstrainedVertices();
            /* we precomputed the normal and collision point of each vertex, now add the constraint (instead of checking collision in each iteration
             this is imprecise but much faster and acc to paper works fine in practice*/
            solveCollisionConstraint();
//        t.printTime(" collision ");
    }

        // (13) velocity and position update
    double collision_damping = 0.5;
    for(int i=0; i<numVert; i++){
         for(int j=0; j<3; j++){
            vel(i,j) = (p(i,j) - x_new(i,j)) / timestep;
         }
         x_new.row(i) = p.row(i);

        if(collisionVert(i)){
            Vector3d vel_i = vel.row(i);
            Vector3d normal = N.row(i).normalized();
            Vector3d n_vel = normal.dot(vel_i) * normal;
            Vector3d t_vel = vel_i-n_vel;
            vel.row(i) = (t_vel - collision_damping * n_vel);
        }
    }

    //(14)
    Vg= x_new;
    // (16) Velocity update
    /*The velocity of each vertex for which a collision constraint has been generated is dampened perpendicular to the collision normal
     * and reflected in the direction of the collision normal.*/
    cout<<normU.sum()<<" sum of the norm u, and v norm  "<<normV.sum()<<endl;

    showGarment(viewer);
}
