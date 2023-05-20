#include <iostream>
#include <Eigen/Dense>
#include <map>
#include <set>
#include <string>
#include <fstream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writePLY.h>
#include <igl/edges.h>
#include <igl/per_edge_normals.h>
#include <igl/adjacency_list.h>
#include <igl/facet_components.h>
#include <igl/vertex_components.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/jet.h>
#include <igl/AABB.h>
#include <igl/signed_distance.h>
#include <igl/exact_geodesic.h>
#include <igl/exact_geodesic.h>
#include <igl/parula.h>
#include <igl/isolines_map.h>

#include "toolbox/PositionBasedDynamics.h"
#include "toolbox/adjacency.h"
#include "toolbox/constraint_utils.h"
#include "toolbox/Timer.h"
#include "toolbox/body_interpolation.h"
#include "toolbox/garment_adaption.h"
#include "toolbox/MathFunctions.h"
#include "toolbox/patternAdaption.h"
#include "toolbox/postProcessing.h"
#include "toolbox/preProcessing.h"
#include "toolbox/seam.h"
#include "toolbox/exportFunctions.h"
using namespace std;
using namespace Eigen;
typedef double Real;
using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;

// The matrices of mesh and garment, original and modified
Eigen::MatrixXd Vg, Vm, testMorph_V1; // mesh for the garment and mannequin
MatrixXd testMorph_V1left, testMorph_V1right;
Eigen::MatrixXi Fg, Fm, Fg_pattern,Fg_pattern_orig, testMorph_F1;
MatrixXi testMorph_F1left, testMorph_F1right; // separate for case of wider legs
map<int, int> leftHalfToFullFaceMap, rightHalfToFullFaceMap; // when wider legs we need two different collision detection
Eigen::MatrixXd Vg_orig, Vm_orig; // original mesh for the garment and mannequin, restore for translation
Eigen::MatrixXd Vg_pattern, Vg_pattern_orig; // the pattern for the rest shape, we might change this
Eigen::MatrixXi Fg_orig, Fm_orig;
Eigen::MatrixXi Eg; // garment edges
Eigen::Vector3d ambient, ambient_grey, diffuse, diffuse_grey, specular;
Eigen::Vector3f garment_translation (0., 0., 0.);// optimal for the given garment mesh
float garment_scale = 1.;
float mannequin_scale = 1.;
Eigen::Vector3f mannequin_translation (0., 0., 0.);

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
bool LShapeAllowed = true;
vector<int> changeFitVert;

enum MouseMode { SELECTPATCH, SELECTBOUNDARY, NONE, SELECTVERT, SELECTAREA, SELECTAREAPATCHES, FINALSTITCH, CHANGEFIT, SELECTBOUNDSEAM };
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
MatrixXd u1, u2; // pre computation for stretch
PositionBasedDynamics PBD, PBD_adaption;

// colouring
MatrixXd colU, colJacDiff ;
MatrixXd colMixed;
MatrixXd colV ;
VectorXd normU, normV;
MatrixXd perFaceU, perFaceV;
int whichStressVisualize= 0;

garment_adaption* gar_adapt;
vector<seam*> seamsList ;
vector<minusOneSeam*> minusOneSeamsList;
map<int, int> mapCornerToCorner;// a map between the corners of the original pattern and the one after re triangulation

vector<int> constrainedVertexIds;
VectorXi closestFaceId;

std::vector<std::pair<Eigen::Vector3d, int>> constrainedVertexBarycentricCoords;
std::vector<double> constrainedVertexDistance;
std::vector<std::pair<Eigen::Vector3d, int>> allVertexBarycentricCoords;
Eigen::MatrixXd tarU, tarV;// tarD1;
Eigen::MatrixXd baryCoords1, baryCoords2;
Eigen::MatrixXd baryCoordsd1, baryCoordsd2;
std::vector<std::vector<int> > boundaryL, boundaryL_toPattern;
double setTheresholdlMid = 1.05;// 5% stress level allowed
double setTheresholdBound = 1.05;


static bool noStress = true;
static bool StressU = false;
static bool StressV = false;
static bool StressDiffJac = false;
static bool StressJac = false;
static bool prioInner = false;
static bool prioOuter = false;
int whichPatchMove=0;
bool prevTearFinished;// indicating if a cut is finished or not to make sure we sort only after a completed cut
bool preferManySmallCuts= false;
bool startTri = false;
std::vector<std::pair<double,double>> perFaceTargetNorm;

// pattern adaption
bool adaptionFlag = false;
MatrixXd fromPattern, currPattern;
MatrixXd toPattern;
Eigen::MatrixXd p_adaption; // the proposed new positions
MatrixXd baryCoordsUPattern, baryCoordsVPattern;
vector<vector<pair<int, int>>> cornerPerBoundary;// for each patch, for each corner it contains the vertex id and the loop id of the vertex
MatrixXd removePatchVg;
MatrixXi removePatchFg;

MatrixXd patternPreInterpol,patternPreInterpol_temp ;
MatrixXd garmentPreInterpol,garmentPreInterpol_temp ;
MatrixXd mannequinPreInterpol, mannequinPreInterpol_temp;
Eigen::VectorXi componentIdPerFace, componentIdPerVert;
vector<VectorXd> polylineSelected;
vector<int> polylineIndex;
vector<int> polyLineMeshIndicator;
//test
//Eigen::SparseMatrix<double> L;
MatrixXd Vg_retri;// re triangulated area
MatrixXi Fg_retri;// re triangulated face
vector<vector<VectorXd>> connectedVert; vector<int> newFaces;// indices of new faces (for coloring), and boundary vertices that are now duplicated
vector<int> isAscVert; vector<bool> isAscMerge;
VectorXd cornerVertices;
vector<cutVertEntry*> cutPositions;
map<int, pair<int, int>>  releasedVert; // all positions that need not be mapped to the boundary anymore from at least one side. keep track which side is released
set<int> toPattern_boundaryVerticesSet; // the boundary vertices of the toPattern, for visualization purposes
vector<int> startAndEnd; // start and end to do the smoothing
double taylor_lazyness = 1;
bool inverseMap, symetry, geoDistU, geoDistV;
int pos; VectorXd T_sym_pattern;
MatrixXd patternEdgeLengths_orig;
MatrixXi mapFromFg, mapToFg, Fg_pattern_curr; //from  = the rest shape of the garment
MatrixXd mapFromVg, mapToVg; //to = the target shape
int changedPos;
map<int, int> halfPatternFaceToFullPatternFace, fullPatternFaceToHalfPatternFace, halfPatternVertToFullPatternVert, fullPatternVertToHalfPatternVert, insertedIdxToPatternVert;
string avName, garment; // avatar name and garment name to write files properly
double geoDistMax = 30.; // radius for fit change
double geoDistChange = 0.97;// values <1 make it bigger. Conter inutitive
VectorXd geoDistDist; // geodesic distances
vector<vector<int>> boundaryLFrom;
string garmentExt;
set<int> nonSymSeam;
void preComputeAdaption();
void computeBaryCoordsGarOnNewMannequin(igl::opengl::glfw::Viewer& viewer,bool otherUsed, MatrixXd& otherM);
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
bool jacobianChanged= false;
// nice clicky interface
bool computePointOnMesh(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& Fuse, Eigen::Vector3d& position, int& fid);
int computeClosestVertexOnMesh(Vector3d& b, int& fid, MatrixXi& F);
bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
void updateChangedBaryCoordinates(int changedPosition, vector<vector<int>>& vfFromPatt);
bool midFractureForbidden= true ;
bool forceCut= false; bool forceClosed = false;
bool showOnly = false;BodyInterpolator* body_interpolator;
pair<int, int> constrainedSeamsSingle;
set<pair<int, int>> constrainedSeamsSet;
void visualizeSeam(pair<int, int> which, igl::opengl::glfw::Viewer& viewer);
int patchCounter = 0;int adaptioncount =0;
void smoothGarment();
void smoothGarmentOutline();
void smoothOutline(MatrixXd & Vpattern ,MatrixXi Fpattern);
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
            if(timestepCounter % convergeIterations == (convergeIterations-1)){
                gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, Vg_pattern,  seamsList, boundaryL, nonSymSeam);
                cout<<"after adaption"<<endl;
                preComputeConstraintsForRestshape();
                preComputeStretch();
                computeStress(viewer);
                t.printTime(" pattern computed ");cout<<endl;
            }
            timestepCounter++;
    }
     if (adaptionFlag){

        doAdaptionStep(viewer);
//         pos = 2354;
//         viewer.selected_data_index = 2;
//         viewer.data().clear();
//         MatrixXd st(1, 3);
//         st.row(0) = currPattern.row(pos);
//         MatrixXd tt(1, 3);
//         tt(0, 0) = -66.1871 ;
//         tt(0, 1) = 1104.35 ;
//         tt(0, 2) = 200;
//         const RowVector3d red(0.8,0.2,0.2);
//         viewer.data().add_edges(st, tt, red);
//        if(pos>0)viewer.data().set_points(currPattern.row(pos), RowVector3d(1.0, 0.0, 0.0));
     }

    return false;
}
bool patternExists;int zipIIdToClose=0;
VectorXd normUPattern, normVPattern;
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
    viewer.append_mesh(); // for other visualizations

    // set background
    viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 1);
    ambient_grey = Vector3d(0.4, 0.4, 0.4);
    ambient = Vector3d(0.26, 0.26, 0.26);
    diffuse_grey = Vector3d(0.5, 0.5, 0.5);
    diffuse = Vector3d(0.4, 0.57, 0.66);    // blue
    specular = Vector3d(0.01, 0.01, 0.01);

//    cout<<"choose garment 3D"<<endl;
//    string garment_file_name = igl::file_dialog_open();
//    cout<<garment_file_name<<" chosen file, thanks. "<<endl;

    string prefix = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/";
//    garment = "leggins";
    cout<<"inverse?  "<<endl;
    string inv ;
    std::getline(std::cin, inv);
    if(inv=="yes"){
        inverseMap  =true;
    }else{
        inverseMap = false;
    }

    cout<<"new pattern ?  "<<endl;
    string pattEx ;
    std::getline(std::cin, pattEx);
    if(pattEx=="no"){
        patternExists = true;
    }else{
        patternExists = false;
    }
//    string fromGarment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/girl24_patternSmoothed.obj"; //patternComputed_"+patternFromName+"_"+garmentExt+".obj";
//    string fromGarment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed_CLO_to_MH_woman_3_3_";
    string fromGarment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed_CLO_avatar_to_bodyScan_Raphael_rem_";

//    fromGarment_pattern_file_name = fromGarment_pattern_file_name+garmentExt+".obj";
//    string patternFromName ="5x5Morphed/CLO_to_MH_woman_3_3";// onnly needed if inter person flag is on
//    string patternFromName =  "girl24";//;
    string patternFromName = "CLO_avatar_to_bodyScan_Raphael_rem";
    bool interPersonFlag = false;
//    bool interPersonFlag= false;
    cout<<"garment?  "<<endl;

    string garInput ;
    std::getline(std::cin, garInput);
    if(garInput=="top"){
        garment  ="top";
        garmentExt = garment+ "_1";
    }
    else if (garInput=="skirt"){
        garment = "skirt";
        garmentExt = garment+ "_2";
    }
    else if (garInput=="hoodie"){
        garment = "hoodie";
        garmentExt = garment;
    }
    else if (garInput=="skirt1"){
        garment = "skirt";
        garmentExt = garment+"_1";
    }
    else if (garInput=="man_tshirt"){
        garment = "man_tshirt";
        garmentExt = garment;
    }
    else if (garInput=="man_tshirt2"){
        garment = "man_tshirt2";
        garmentExt = garment;
    }
    else if (garInput=="leggins"){
        garment = "leggins";
        garmentExt = garment;
    }
    else if (garInput=="man_pants"){
        garment = "man_pants";
        garmentExt = garment;
    }
    else if(garInput=="dress"){
        garment = garInput;
        garmentExt = garment+"_4";
    }

//    garment = "tshirt";
//    garment = "top";
//    garment = "top_fromAnna";
//    garment = "leggins";
//    garment = "man_pants";
//    string garment_file_name = prefix+ "leggins/leggins_3d/leggins_3d_merged.obj"; //smaller collision thereshold to make sure it is not "eaten" after intial step , 3.5 instead of 4.5
//    garment = "dress";
//    garment = "top";
//    garment = "hoodie"; // attention also needs flat starter
//    garment = "man_tshirt2";
//    garmentExt = garment;
    fixRafaPattern();
//     garmentExt = garment+ "_1";// skirt 3 needs another avatar!!// 2 is pencil skirt, 1 is mini skirt
//     cout<<fromGarment_pattern_file_name<<" fromGarment_pattern_file_name"<<endl;
    fromGarment_pattern_file_name = fromGarment_pattern_file_name + garmentExt+".obj";

    string garment_file_name = prefix + "moreGarments/"+ garmentExt+"/"+garment+"_3d.obj";

    igl::readOBJ(garment_file_name, Vg, Fg);
    Timer t("Setup");

//    string garment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed_"+patternFromName+"_"+garment+".obj";
//    string garment_pattern_file_name = prefix +"leggins/leggins_2d/leggins_2d.obj"; //
    string garment_pattern_file_name = prefix +"moreGarments/"+garmentExt+"/"+garment+"_2d.obj";

    //orig down here
//    string avatar_file_name =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins/avatar/avatar_Girl24.ply";
    string avatar_file_name;
    if(garment == "hoodie"){
        avatar_file_name =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins/avatar/avatar_flat.ply";

    }else if(garment == "man_tshirt"|| garment == "man_tshirt2"|| garment == "man_pants"){
        avatar_file_name =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins/avatar/male_avatar_rem_20.ply";

    }else{
        avatar_file_name =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/leggins/avatar/avatar_one_component.ply";

    }

    igl::readOBJ(garment_pattern_file_name, Vg_pattern, Fg_pattern);
//    garment = "skirt_no2";
    symetry = true;
    if(symetry){
        bool insertPlane = true;
        int symVert1 ,symVert2;
//
        if (garment == "leggins"){
            insertPlane = false;
            symVert1 = 1467;
            symVert2 = 23;
        }else if(garment == "top"){
            insertPlane = false;
            symVert1 = 353;// from orig
            symVert2 = 0;//from dupl
        }else if(garment == "top"){
            symVert1 = 2097;
            symVert2 = 1034;
        }else if (garmentExt == "skirt_2"){
            symVert1= 340;
            symVert2= 23;
        }else if (garmentExt == "skirt_1"){
            symVert1= 66;
            symVert2= 66;
        }
        else if(garmentExt == "man_pants"){
            insertPlane = false;
            symVert1= 0; symVert2= 0;

        }
        else{
            symVert1= 0; symVert2= 0;
        }
        if(insertPlane){

            MatrixXd VgPatternRet; MatrixXi FgPatternRet;
            MatrixXd VgRet; MatrixXi FgRet;
            splitAndSmooth(Vg, Fg, Vg_pattern, Fg_pattern, VgPatternRet, FgPatternRet, VgRet, FgRet, garment, garmentExt);
            startAndEnd.clear();
            if (garmentExt == "skirt_2"){
                startAndEnd.push_back(698);
                startAndEnd.push_back(695);
                startAndEnd.push_back(693);
            } else if (garmentExt == "skirt_3"){
                startAndEnd.push_back(1566);
                startAndEnd.push_back(1587);
                startAndEnd.push_back(1561);

            }
            else if (garmentExt == "skirt_1" ){
                startAndEnd.push_back(405);
                startAndEnd.push_back(406);
                startAndEnd.push_back(404);
            }
            else if (garment == "hoodie"){
                startAndEnd.push_back(2173);
                startAndEnd.push_back(2171);
                startAndEnd.push_back(2177);

            }
            else if (garment == "man_tshirt"){
                startAndEnd.push_back(1834);
                startAndEnd.push_back(1835);
                startAndEnd.push_back(1832);
            }
            else if ( garment == "man_tshirt2"){
                startAndEnd.push_back(1840);
                startAndEnd.push_back(1841);
                startAndEnd.push_back(1047);
            }else if ( garment == "dress"){
                startAndEnd.push_back(992);
                startAndEnd.push_back(993);
                startAndEnd.push_back(991);
            }
            smoothBetweenVertices(VgPatternRet, FgPatternRet, startAndEnd);
            igl::writeOBJ("leftPattern_SmoothBoundAfter1.obj", VgPatternRet, FgPatternRet);

            startAndEnd.clear();
            if (garmentExt == "skirt_2") {
                startAndEnd.push_back(700);
                startAndEnd.push_back(699);
                startAndEnd.push_back(696);
            }else if (garmentExt == "skirt_3"){
                startAndEnd.push_back(1616);
                startAndEnd.push_back(1622);
                startAndEnd.push_back(1620);
//                igl::writeOBJ("leftPattern_SmoothBoundAfter2.obj", VgPatternRet, FgPatternRet);
            }
            else if  (garmentExt == "skirt_1"){
                startAndEnd.push_back(425);
                startAndEnd.push_back(424);
                startAndEnd.push_back(427);
            }
            else if (garment == "hoodie"){
                startAndEnd.push_back(2164);
                startAndEnd.push_back(2163);
                startAndEnd.push_back(2161);

            }
            else if (garment == "man_tshirt"){
                startAndEnd.push_back(1876);
                startAndEnd.push_back(1878);
                startAndEnd.push_back(1877);
            }
            else if ( garment == "man_tshirt2"){
                startAndEnd.push_back(1786);
                startAndEnd.push_back(1791);
                startAndEnd.push_back(1787);
            }
            else if ( garment == "dress"){
                startAndEnd.push_back(0);
                startAndEnd.push_back(958);
                startAndEnd.push_back(918);
            }
            smoothBetweenVertices(VgPatternRet, FgPatternRet, startAndEnd);

            startAndEnd.clear();
            if (garmentExt == "skirt_2") {
                startAndEnd.push_back(677);
                startAndEnd.push_back(669);
                startAndEnd.push_back(651);
            }
            else if (garmentExt == "skirt_3"){
                startAndEnd.push_back(1619);
                startAndEnd.push_back(1618);
                startAndEnd.push_back(1615);
//                igl::writeOBJ("leftPattern_SmoothBoundAfter3.obj", VgPatternRet, FgPatternRet);

            }
            else if  (garmentExt == "skirt_1"){
                startAndEnd.push_back(400);
                startAndEnd.push_back(399);
                startAndEnd.push_back(398);
            }
            else if (garment == "hoodie"){
                startAndEnd.push_back(2180);
                startAndEnd.push_back(2179);
                startAndEnd.push_back(2183);

            }
            if(garment !="man_tshirt" && garment !="man_tshirt2"&& garment!="dress" ) smoothBetweenVertices(VgPatternRet, FgPatternRet, startAndEnd);
            startAndEnd.clear();
            if (garmentExt == "skirt_2") {
                startAndEnd.push_back(704);
                startAndEnd.push_back(703);
                startAndEnd.push_back(701);
            }
            else if (garmentExt == "skirt_3"){
                startAndEnd.push_back(1505);
                startAndEnd.push_back(1504);
                startAndEnd.push_back(1518);
//                igl::writeOBJ("leftPattern_SmoothBoundAfter4.obj", VgPatternRet, FgPatternRet);


            }
            else if  (garmentExt == "skirt_1"){
                startAndEnd.push_back(402);
                startAndEnd.push_back(397);
                startAndEnd.push_back(401);
            }
            else if (garment == "hoodie"){
                startAndEnd.push_back(2105);
                startAndEnd.push_back(2101);
                startAndEnd.push_back(2104);
            }
            if(garment !="man_tshirt" && garment !="man_tshirt2") smoothBetweenVertices(VgPatternRet, FgPatternRet, startAndEnd);
            if(garment == "hoodie"){
                startAndEnd.clear();
                startAndEnd.push_back(2118);
                startAndEnd.push_back(2114);
                startAndEnd.push_back(2117);
                smoothBetweenVertices(VgPatternRet, FgPatternRet, startAndEnd);

                startAndEnd.clear();
                startAndEnd.push_back(2111);
                startAndEnd.push_back(2107);
                startAndEnd.push_back(2110);
                smoothBetweenVertices(VgPatternRet, FgPatternRet, startAndEnd);

            }

            igl::writeOBJ("leftPattern_SmoothBound.obj", VgPatternRet, FgPatternRet);

        }
        preProcessGarment(Vg, Fg, Vg_pattern, Fg_pattern, insertPlane, symVert1, symVert2, T_sym_pattern, garment, garmentExt);
        Vg_orig = Vg; Fg_orig= Fg;
    }
    if(garment == "man_tshirt2"){
        nonSymSeam.insert(0);
        nonSymSeam.insert(2);
        nonSymSeam.insert(5);
        nonSymSeam.insert(8);
    }else if (garment == "hoodie"){
        nonSymSeam.insert(2);
        nonSymSeam.insert(6);
        nonSymSeam.insert(9);
        nonSymSeam.insert(19);
        nonSymSeam.insert(23);
        nonSymSeam.insert(26);

    }
    else if (garment == "man_tshirt"){
        nonSymSeam.insert(2);
        nonSymSeam.insert(3);
        nonSymSeam.insert(5);
        nonSymSeam.insert(6);
        nonSymSeam.insert(12);
        nonSymSeam.insert(14);
        nonSymSeam.insert(16);
        nonSymSeam.insert(17);
    }
    garmentPreInterpol = Vg;
    Vg_orig = Vg; Fg_orig= Fg;

    Vg_pattern_orig = Vg_pattern;
    Fg_pattern_orig = Fg_pattern;
    patternPreInterpol = Vg_pattern;

    cornerVertices = VectorXd::Zero(Vg_pattern.rows());
    t.printTime(" init");
    preComputeConstraintsForRestshape();
    preComputeStretch();

    setNewGarmentMesh(viewer);

// TODO remember to adapt the collision constraint solving dep on avatar, sometimes normalization is needed, sometimes not for whatever magic
//    string avatar_file_name =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/CLO_avatars_oneComponent/CLO_avatar_to_bodyScan_Anna_rem 2.ply";

    //string avatar_file_name = igl::file_dialog_open();
    igl::readPLY(avatar_file_name, Vm, Fm);
    Vm_orig = Vm; Fm_orig = Fm;
    mannequinPreInterpol = Vm;
    string folderName = "";
//     avName = "avatar_missy_straight_05_OC";// good for skirt
//     avName = "avatar_petite_curvy_01_OC";
//    avName = "avatar_maternity_05_OC";
//    avName = "avatar_missy_straight_09_OC";//new
//    avName = "avatar_plus_straight_05_OC";

// 5x5 av

//
    cout<<"enter which woman_ ? "<<endl;
    avName = "CLO_to_MH_woman_";
//    avName = "CLO_avatar_to_bodyScan_Raphael_rem";
    avName = "avatar_maternity_05_OC";
    string whichName ;
    std::getline(std::cin, whichName);
    if(whichName =="missy str"){
        avName = "avatar_missy_straight_09_OC";
    }else if (whichName =="petite"){
        avName = "avatar_petite_curvy_01_OC";
    }else if (whichName =="petite str"){
        avName = "avatar_petite_straight_01_OC";
    }
    else if (whichName =="missy 5"){
        avName = "avatar_missy_straight_05_OC";
    }
    else if (whichName =="normal"){
        avName = "avatar_oneComponent";
    }else if (whichName =="mat"){
        avName = "avatar_maternity_05_OC";
    }else if (whichName=="man"){
            avName = "CLO_avatar_to_male_large_avatar_rem_20";
    }
    else if (whichName=="man orig"){
        avName = "male_avatar_rem_20";

    }
    else if (whichName=="plus"){
        avName = "avatar_plus_straight_05_OC";

    }else if (whichName == "Paola"){
        avName= "CLO_avatar_to_bodyScan_Paola_rem";

    }else if (whichName == "Luka"){
        avName= "CLO_avatar_to_bodyScan_Luka_rem";

    }else if (whichName == "Rapha"){
        avName= "CLO_avatar_to_bodyScan_Raphael_rem";

    }
    else if (whichName == "Anna"){
        avName= "CLO_avatar_to_bodyScan_Anna_rem 2";

    }
    else{
//        folderName = "TeaserAvatarsRes/";
        folderName = "5x5Morphed/";
        avName = "CLO_to_MH_woman_" +whichName;
    }
//    avName = whichName;
//    avName += whichName;

    cout<<"entered "<<avName<<endl;
    string morphBody1 =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/CLO_avatars_oneComponent/"+folderName+ avName +".ply";//
    string morphBody1left =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/CLO_avatars_oneComponent/"+ folderName + avName +"_left.ply";
    string morphBody1right =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/data/CLO_avatars_oneComponent/"+folderName+ avName +"_right.ply";

    igl::readPLY(morphBody1, testMorph_V1, testMorph_F1);
    body_interpolator = new BodyInterpolator(Vm_orig, testMorph_V1, testMorph_F1);
    igl::readPLY(morphBody1left, testMorph_V1left, testMorph_F1left);
    igl::readPLY(morphBody1right, testMorph_V1right, testMorph_F1right);

    createHalfAvatarMap( testMorph_V1, testMorph_F1, testMorph_V1left, testMorph_F1left,
                         testMorph_V1right, testMorph_F1right, leftHalfToFullFaceMap, rightHalfToFullFaceMap);

    if(Fm != testMorph_F1){
        cout<<"the faces are not the same!"<<endl;
    }

    setNewMannequinMesh(viewer);
    std::map<int,int> vertexMapPattToGar;
    std::map<std::pair<int, int>,int> vertexMapGarAndIdToPatch;
    vertexMapPatternToGarment(Fg, Fg_pattern,vertexMapPattToGar);

    igl::boundary_loop(Fg_pattern, boundaryL);

    igl::facet_components(Fg_pattern, componentIdPerFace);
    vertex_componentsBasedOnFacet(Fg_pattern, componentIdPerFace, componentIdPerVert, Vg_pattern.rows());
    vertexMapGarmentAndPatchIdToPattern(Fg, Fg_pattern, componentIdPerVert, vertexMapGarAndIdToPatch);

    // use adjacentFacesToEdge of the 3D
    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg, vfAdj);
    cornerVertices = VectorXd::Zero(Vg_pattern.rows());// 1 for each corner, 0 for all other vertices

    map<int, vector<pair<int, int>>> seamIdPerCorner;    // contains corner id and a list of which seams start here (max 2),
    // each vector element  is a pair where first is if it's a -1 seam, and second is which index in corresponding List. If negative it's a backside ,i.e. part 2 of the seam

    computeAllSeams(boundaryL, vertexMapPattToGar, vertexMapGarAndIdToPatch, vfAdj, componentIdPerFace,
                    componentIdPerVert, cornerVertices, cornerPerBoundary, seamsList, minusOneSeamsList, seamIdPerCorner, garment, garmentExt);

    set<int> cornerSet;// a set containing all corner vertices

    for(auto cpbi : cornerPerBoundary){
        for(auto cpbij : cpbi){
            cornerSet.insert(cpbij.first);
        }
    }

    gar_adapt = new garment_adaption(Vg, Fg,  Vg_pattern, Fg_pattern, seamsList, boundaryL); //none have been altered at this stage
    gar_adapt->computeJacobian();
    perFaceTargetNorm = gar_adapt->perFaceTargetNorm;
    Vg_orig = Vg;

    setCollisionMesh();
    MatrixXd dummy;
    computeBaryCoordsGarOnNewMannequin(viewer, false,dummy );// contains boundary vertices now, needed for simulation, colsestFaceId
//    Vg = Vg_orig;
    Vm = testMorph_V1;
    Vm_orig = testMorph_V1;
    showGarment(viewer);
    showMannequin(viewer);

    preComputeConstraintsForRestshape();
    preComputeStretch();
    computeStress(viewer);
    setCollisionMesh();

    MatrixXd perfPattVg, perfPattVg_orig, addedFabricPatternVg;
    MatrixXi perfPattFg, perfPattFg_orig, addedFabricPatternFg;

    string startFile = "finished_retri_writtenPattern_"+avName+"_"+garmentExt+".obj";
    int vertIdOf0InDuplicated = 0;
    if(garmentExt =="top_1" ){
        vertIdOf0InDuplicated = 353;
    }else if (garment =="leggins" ){
        vertIdOf0InDuplicated = 1444; //-> check!
    }
    vector<seam*> seamsListDupl = seamsList;
    if(patternExists){
        cout<<" pattern exists"<<endl;
        if(interPersonFlag){

//            string fromGarment_pattern_file_name = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/girl24_patternSmoothed.obj"; //patternComputed_"+patternFromName+"_"+garmentExt+".obj";
            igl::readOBJ(fromGarment_pattern_file_name, Vg_pattern, Fg_pattern);
            Vg_pattern_orig.resize(Vg_pattern.rows(), 3);
            Fg_pattern_orig.resize(Fg_pattern.rows() ,3);
            Vg_pattern_orig = Vg_pattern;
            Fg_pattern_orig = Fg_pattern;
        }

        string prefPattern = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/";
        string perfPatternFile = prefPattern+ "patternComputed_"+avName+"_"+garmentExt+".obj";

        igl::readOBJ(perfPatternFile, perfPattVg_orig, perfPattFg_orig);
        if(garment == "top"){
            duplicateInitPattern( perfPattVg_orig, perfPattFg_orig);
        }
        perfPattVg_orig.col(2).setConstant(200);

        // copy the matrices to not mess with them
        if(inverseMap){
//            string helperToLocate = "/Users/annaeggler/Desktop/"+startFile;
            string helperToLocate = prefPattern + startFile;

            igl::readOBJ(helperToLocate, addedFabricPatternVg, addedFabricPatternFg);
//            mapFromVg = perfPattVg_orig;
//            mapFromFg = perfPattFg_orig;
            string addedFile = prefPattern+ "adaption2D_"+avName+"_"+garment+".obj";
            MatrixXd addedFabricPatternVgi; MatrixXi addedFabricPatternFgi;
            igl::readOBJ(addedFile, addedFabricPatternVgi, addedFabricPatternFgi);
            mapFromVg = addedFabricPatternVgi;
            mapFromFg = addedFabricPatternFgi;
//        Fg_pattern_curr = mapFromFg;
            mapToVg =  Vg_pattern_orig ;// curr = the current shape of the garment, something in between
            mapToFg = Fg_pattern_orig ;// the stress is computed between the rest shape and the current, ie mapFromVg and currPattern

            perfPattVg = perfPattVg_orig;
            perfPattFg = perfPattFg_orig;

        }else{
            mapToVg.resize(perfPattVg_orig.rows(), 3);
            mapToVg = perfPattVg_orig;
            cout<< mapToVg.rows()<<endl;
            mapToFg.resize(perfPattFg_orig.rows(), 3);
            mapToFg =  perfPattFg_orig;
            mapFromVg = Vg_pattern;
            mapFromFg = Fg_pattern;
            VectorXd pattComp;
            igl::vertex_components(Fg_pattern, pattComp);
            cout<<" no inverse mapping"<<endl;
            preComputeStretch();
            cout<<" no inverse mapping 2"<<endl;


        }
    }

    viewer.core().animation_max_fps = 200.;
    viewer.core().is_animating = false;
    int whichSeam = 0;
    changedPos = -1;
    //additional menu items
    int whichBound=0;
    float movePatternX=0; float movePatternY=0;
    bool showPattern= false;
    set<int> handledVerticesSet;
    pos = -1;

    MatrixXi Fg_pattern_half;
    MatrixXd Vg_pattern_half, rightVert;
    VectorXi isLeftVertPattern;
    int numFacesOneSide ;
    if(patternExists) {
        if (symetry && !inverseMap) {
            createHalfSewingPattern(Vg_orig, Fg_orig, Vg_pattern, Fg_pattern, Vg_pattern_half, Fg_pattern_half,
                                    halfPatternFaceToFullPatternFace, fullPatternFaceToHalfPatternFace,
                                    halfPatternVertToFullPatternVert,
                                    fullPatternVertToHalfPatternVert, insertedIdxToPatternVert, isLeftVertPattern,rightVert, garment);
            numFacesOneSide = Fg_pattern_half.rows();

        } else if (symetry && inverseMap) {
            // map from is already split in two sides,
            // do the same for map to
            createHalfSewingPattern(Vg_orig, Fg_orig, mapToVg, mapToFg, Vg_pattern_half, Fg_pattern_half,
                                    halfPatternFaceToFullPatternFace, fullPatternFaceToHalfPatternFace,
                                    halfPatternVertToFullPatternVert,
                                    fullPatternVertToHalfPatternVert, insertedIdxToPatternVert, isLeftVertPattern,rightVert, garment);
            numFacesOneSide = Fg_pattern_half.rows();

        } else {
            for (int i = 0; i < Vg_pattern.rows(); i++) {
                halfPatternVertToFullPatternVert[i] = i;
                fullPatternVertToHalfPatternVert[i] = i;
            }
            for (int i = 0; i < Fg_pattern.rows(); i++) {
                halfPatternFaceToFullPatternFace[i] = i;
                fullPatternFaceToHalfPatternFace[i] = i;
            }
        }
    }
    int patchcount=0;
    bool showPatchBoundary = false ;
    ofstream out("setupFile"+avName+"_"+garmentExt+ +".txt");
    out<<inverseMap<<" inverseMap"<<endl;
    out<<patternExists<<" pattern exists"<<endl;
    out<<garment<<"  garment"<<endl;

    bool trialFlag= false;
//    garmentExt = "top_fromAnna";
    menu.callback_draw_viewer_menu = [&]() {
        if (ImGui::CollapsingHeader("Garment", ImGuiTreeNodeFlags_OpenOnArrow)) {

            if(ImGui::Button(" Trial stress after ", ImVec2(-1, 0))){
                trialFlag = true;
                MatrixXd newPatt, newGar;
                MatrixXi newPattF, newGarF;
                string fileName = "adaption3D_"+avName + "_"+garment+".obj";
                igl::readOBJ(fileName, newGar, newGarF);
                fileName = "adaption2D_"+avName + "_"+garment+".obj";
                igl::readOBJ(fileName, newPatt, newPattF);

                Vg.resize(newGar.rows() ,3);
                Vg= newGar;
                Fg.resize(newGarF.rows() ,3);
                Fg= newGarF;

                Vg_pattern.resize(newPatt.rows() ,3);
                Vg_pattern= newPatt;
                Fg_pattern.resize(newPattF.rows() ,3);
                Fg_pattern= newPattF;
                numFace = Fg.rows();
                preComputeStretch();
                computeStress(viewer);

//                smoothGarmentOutline();
//                smoothGarmentOutline();
//                smoothGarmentOutline();
//                smoothGarment();
//                smoothGarment();
                showGarment(viewer);


            }
            if(ImGui::Button("Export Color .mtl ", ImVec2(-1, 0))){
                string dir;
                if(patternExists){
                    dir="afterPattern_";
                }
                if(garment=="top"||garment =="hoodie"){
                    smoothGarment();
//                    smoothGarmentOutline();
                }else if (garment == "skirt"){
                    smoothGarment();
                    smoothGarmentOutline();
                    smoothGarment();
                    smoothGarmentOutline();
                    smoothGarmentOutline();
                }
                for(int i =0;i<3; i++){
                    if(i==0){
                        dir = "v";
                        igl::jet(normV, 0., 2., colU);
                    }else if (i==1){
                        dir = "u";
                        igl::jet(normU, 0., 2., colU);
                    }else{
                        dir = "stressMax";
                        VectorXd maxNorm(normU.rows());
                        for(int i=0; i<normU.rows(); i++){
                            maxNorm(i)= max(normU(i) ,normV(i));
                        }
                        igl::jet(maxNorm, 0., 2., colU);
                    }
                    // make it symmetric
                    colU.block(colU.rows()/2, 0,colU.rows()/2, colU.cols())= colU.block(0, 0,colU.rows()/2, colU.cols());
                    viewer.data().set_colors(colU);
                    MatrixXd Ka, Ks, Kd;
                    Ka = viewer.data().F_material_ambient;
                    Kd = viewer.data().F_material_diffuse;
                    Ks = viewer.data().F_material_specular;

                    double itVal = 0.;
                    if(patternExists==false){
                        itVal = 20.;
                    }
                    string modifiedPattern  = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/garmentPattern_"+ avName +"_"+ garmentExt+"_backIn3d" +".obj";
                    MatrixXd addedFabricPatternV, adaptedPatternIn3d;
                    MatrixXi addedFabricPatternFg, adaptedPatternIn3d_faces;
                    igl::readOBJ(modifiedPattern, adaptedPatternIn3d, adaptedPatternIn3d_faces);

                    if(trialFlag){
                        cout<<"using the trial flag"<<endl;
                        writeMTL(Ka,Ks, Kd, Vg, Fg, garment, avName, 111, dir );

                    }else{
                        cout<<"not using trial flag"<<endl;
                        writeMTL(Ka,Ks, Kd, adaptedPatternIn3d, adaptedPatternIn3d_faces, garment, avName, 0, dir );

                    }
                    dir = dir+"_pattern";
                    if(patternExists){
                        cout<<"writing with pattern"<<endl;
                        if(trialFlag){
                            cout<<"again with trial flag"<<endl;
                            writeMTL(Ka,Ks, Kd, Vg_pattern, Fg_pattern, garment, avName,222, dir );

                        }
                        else
                            writeMTL(Ka,Ks, Kd, perfPattVg_orig, perfPattFg_orig, garment, avName, 0, dir );

                    }else{
                        writeMTL(Ka,Ks, Kd, Vg_pattern_orig, Fg_pattern_orig, garment, avName, 0, dir );

                    }

                }

            }
            ImGui::InputInt("Vis Seam No", &whichSeam, 0, 0);
            if(ImGui::Button("Visualize Seam", ImVec2(-1, 0))){
                MatrixXd testCol= MatrixXd::Zero(Vg_pattern.rows(), 3);

                for(int j=whichSeam; j<whichSeam+1; j++){
                        seam* firstSeam = seamsList[j];
                        auto stP1 = firstSeam-> getStartAndPatch1();
                        auto stP2 = firstSeam-> getStartAndPatch2ForCorres();
//
                        int len = firstSeam -> seamLength();
                        int boundLen1 = boundaryL[stP1.second].size();
                        int boundLen2 = boundaryL[stP2.second].size();
                        MatrixXi edgesMat (2*(len), 2);
                        cout<<"start iterating"<<endl;
                        for(int i=0; i<=len; i++){
                            testCol(boundaryL[stP1.second][(stP1.first+i)% boundLen1],0) = 1.;
                            int next = (stP1.first+i)% boundLen1;
//                            cout<<boundaryL[stP1.second][next]<<" next"<<endl;
                            if (i!= 0) edgesMat(2*(i-1), 1) = boundaryL[stP1.second][next];
                            if(i!=len)edgesMat(2*i, 0) = boundaryL[stP1.second][next];
//                            if (i!= 0) edgesMat((i-1), 1) = boundaryL[patch][(firstSeam->getStartIdx() +i)];
//                            cout<<boundaryL[patch][(firstSeam->getStartIdx() +i)]<<endl;
//                            edgesMat(i, 0) = boundaryL[patch][(firstSeam->getStartIdx() +i) ];
//
                            if(i==0)testCol(boundaryL[stP1.second][(stP1.first+i)% boundLen1],1) = 1.;
                            if(i==len)testCol(boundaryL[stP1.second][(stP1.first+i)% boundLen1],2) = 1.;

                            int setAccess = (stP2.first-i)% boundLen2;
                            if(setAccess < 0) {
                                setAccess +=boundLen2;
                                testCol(boundaryL[stP2.second][(stP2.first-i)% boundLen2], 0) = 1.;
                            }
                            if(seamsList[j]->inverted) setAccess = (stP2.first + i) % boundLen2;
                            testCol(boundaryL[stP2.second][setAccess], 0) = 1.;
                            if (i!= 0) edgesMat(2*i-1, 1) = boundaryL[stP2.second][setAccess];
                            if (i!= len)edgesMat(2*i+1, 0) = boundaryL[stP2.second][setAccess];
//
                        }
                        viewer.selected_data_index = 1;
                        viewer.data().clear();
                        viewer.data().set_edges(Vg_pattern, edgesMat, Eigen::RowVector3d(1, 0, 0));
                    }

//                viewer.selected_data_index = 1;
//                viewer.data().clear();

                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(Vg_pattern, Fg_pattern);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                //remove wireframe
                viewer.data().show_lines = false;
                // if 0 -> no face colour
//                viewer.data().set_colors(testCol);

            }
            if(ImGui::Button("Visualize Corner", ImVec2(-1, 0))){
                MatrixXd testCol= MatrixXd::Zero(Vg_pattern.rows(), 3);

                int coutn=0;
                for(int i=0; i<Vg_pattern.rows(); i++){
                    if(cornerVertices(i)){
//                        corners.row(coutn) = Vg_pattern.row(i);
                        coutn++;
                        testCol(i,0)=1;
                    }
                }
                MatrixXd corners(coutn, 3);
                coutn =0;
                for(int i=0; i<Vg_pattern.rows(); i++){
                    if(cornerVertices(i)){
                        corners.row(coutn) = Vg_pattern.row(i);
                        coutn++;
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
//                viewer.data().set_colors(testCol);
                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.data().point_size = 8.f;
                viewer.data().set_points(corners, RowVector3d(1.0, 0.0, 0.0));

            }
        }
        if (ImGui::CollapsingHeader("Pattern Computation", ImGuiTreeNodeFlags_OpenOnArrow)) {
            if(ImGui::Button("Smooth Garment ", ImVec2(-1, 0))){
                smoothGarment();
                showGarment(viewer);
                out<<"garment smoothed"<<endl;
                computeStress(viewer);


            }
            if(ImGui::Button("Smooth Outline ", ImVec2(-1, 0))){
                smoothGarmentOutline();
                showGarment(viewer);
                out<<"outline garment smoothed"<<endl;
                computeStress(viewer);

            }
            if(ImGui::Checkbox("Show Pattern", &showPattern)){
                cout<<Vg_pattern.rows()<<" "<<Fg_pattern.rows()<<endl;
                viewer.selected_data_index = 1;
                viewer.data().clear();
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
                    MatrixXd testCol= MatrixXd::Ones(Vg_pattern.rows(), 3);
                    testCol.col(0).setConstant(0.01);
                    for(int i=0; i<Vg_pattern.rows(); i++){
                        if( componentIdPerVert(i) == whichPatchMove){
                            Vg_pattern(i, 0) += movePatternX;
//                            Vg_pattern(i, 1) += movePatternY;
//                            testCol(i,0)=1; testCol(i,1)=0; testCol(i,2)=0;

                        }
                    }
                    viewer.selected_data_index = 1;
                    viewer.data().clear();
                    viewer.data().set_mesh(Vg_pattern, Fg_pattern);
                    viewer.data().uniform_colors(ambient, diffuse, specular);
                    viewer.data().show_texture = false;
                    viewer.data().set_face_based(false);
                    viewer.data().show_lines = true;

                    viewer.selected_data_index = 0;
                    viewer.data().clear();

                    igl::writeOBJ("patternComputed_translated.obj", Vg_pattern, Fg_pattern);
                    cout<<"pattern written to *patternComputed_translated*"<<endl;
                }
                ImGui::InputInt("Which Bound vis  ", &(whichBound),  0, 0);
                if(ImGui::Button("Vis Bound ", ImVec2(-1, 0))){
                    MatrixXd testCol= MatrixXd::Zero(Vg_pattern.rows(), 3);
                    for(auto bli: boundaryL[whichBound]){
                            testCol(bli,0)=1;
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
                gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, computed_Vg_pattern, seamsList, boundaryL, nonSymSeam);
                if(!jacobianChanged ){
                    igl::writeOBJ("patternComputed_"+avName+"_"+garmentExt+".obj",Vg_pattern, Fg_pattern);
                    cout<<"pattern written to *patternComputed*"<<endl;
                }else{
                    string dir = (geoDistU)? "U" : "V";
                    igl::writeOBJ("patternComputed_changedFit_" + avName + "_" + garmentExt + "_" +to_string(geoDistMax)+ "_" + to_string(geoDistChange) + "_" + dir + ".obj",Vg_pattern, Fg_pattern);
                    cout<<"pattern written to * patternComputed changed fit *"<<endl;
                }
            }
            if(ImGui::Button("Compute and Visualize stress of new pattern", ImVec2(-1, 0))){
                simulate = false;
                // we start computing the pattern for the current shape
                Eigen::MatrixXd computed_Vg_pattern;//= Vg;
                cout<<"start computing the pattern with "<<localGlobalIterations<<" local global iterations"<<endl;
                gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, computed_Vg_pattern, seamsList, boundaryL, nonSymSeam);
                igl::writeOBJ("patternComputed_"+avName+"_"+garmentExt+".obj", computed_Vg_pattern, Fg_pattern);
                cout<<"pattern written to *patternComputed*"<<endl;

                Vg_pattern_orig = computed_Vg_pattern;
                Vg_pattern = computed_Vg_pattern;

                preComputeConstraintsForRestshape();
                preComputeStretch();
                computeStress(viewer);

            }
        }
        if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_OpenOnArrow)) {

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
                StressJac = false;
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
                StressV= false;
                noStress = false;
                StressU = false;
                StressDiffJac = false;
                StressJac = false;

                if(remMan){

                    patternPreInterpol_temp = Vg_pattern;
                    garmentPreInterpol_temp = Vg;
                    mannequinPreInterpol_temp = Vm;
                    Vg_pattern = patternPreInterpol;
                    Vg = garmentPreInterpol;
                    Vm = mannequinPreInterpol;
                }else{
                    Vg_pattern = patternPreInterpol_temp;
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
        if (ImGui::CollapsingHeader("Smooth Pattern", ImGuiTreeNodeFlags_OpenOnArrow)){
            bool startSmooth = false;
            bool choosePatchArea = false;
            bool choosePatches = false;
            if(ImGui::Checkbox("init ", &startSmooth)) {
//                string smoothPatternFile ="final_addedFace_patternNotMerged.obj";
//                string smoothPatternFile ="finished_tear_writtenPattern_CLO_avatar_to_bodyScan_Paola_rem_top.obj";
                string smoothPatternFile =  "patternComputed_"+avName+"_"+garmentExt+".obj";
                igl::readOBJ(smoothPatternFile, currPattern, Fg_pattern_curr);
                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.selected_data_index = 2;
                viewer.data().clear();
                viewer.selected_data_index = 0;
                viewer.data().clear();

                viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                viewer.data().show_lines = true;
            }

            if(ImGui::Checkbox("Start smooth", &startSmooth)) {
                cout<<"Please choose 3 points to smooth between.  "<<endl;
                startAndEnd.clear();
                simulate = false;
                adaptionFlag = false;
                mouse_mode= SELECTVERT;
            }
            if(ImGui::Button("Confirm smooth", ImVec2(-1, 0))) {
                if (startAndEnd.size() == 3) {
                    cout << "Great, you selected " << startAndEnd[0] << " and " << startAndEnd[2]
                         << " as endpoints. Let's get to work on smoothing. " << endl;
                    out<<"Garment smoothed" << startAndEnd[0]<<", "<<startAndEnd[1] << " and " << startAndEnd[2]<<endl;

                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                }
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(currPattern, Fg_pattern_curr);

                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                viewer.data().show_lines = true;
            }
            if(ImGui::Button("Predefined smooth", ImVec2(-1, 0))) {
                out<<"predefined smoothed"<<endl;

                if(garment =="leggins"){
                    simulate = false;
                    adaptionFlag = false;
                    startAndEnd.clear();
                    startAndEnd.push_back(3019);
                    startAndEnd.push_back(3020);
                    startAndEnd.push_back(3027);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();


                    startAndEnd.push_back(2968);
                    startAndEnd.push_back(2969);
                    startAndEnd.push_back(2976);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(2987);
                    startAndEnd.push_back(2986);
                    startAndEnd.push_back(2978);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1524);
                    startAndEnd.push_back(1525);
                    startAndEnd.push_back(1533);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(2363);
                    startAndEnd.push_back(2364);
                    startAndEnd.push_back(2372);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(2372);
                    startAndEnd.push_back(2373);
                    startAndEnd.push_back(2307);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1524);
                    startAndEnd.push_back(1638);
                    startAndEnd.push_back(1592);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(2312);
                    startAndEnd.push_back(2311);
                    startAndEnd.push_back(2307);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1592);
                    startAndEnd.push_back(1591);
                    startAndEnd.push_back(1586);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();
// right side
                    startAndEnd.push_back(1452);
                    startAndEnd.push_back(1451);
                    startAndEnd.push_back(1444);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1503);
                    startAndEnd.push_back(1502);
                    startAndEnd.push_back(1495);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1484);
                    startAndEnd.push_back(1485);
                    startAndEnd.push_back(1493);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1454);
                    startAndEnd.push_back(1455);
                    startAndEnd.push_back(1463);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();
//main leg
                    startAndEnd.push_back(9);
                    startAndEnd.push_back(8);
                    startAndEnd.push_back(0);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(848);
                    startAndEnd.push_back(847);
                    startAndEnd.push_back(839);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(848);
                    startAndEnd.push_back(849);
                    startAndEnd.push_back(783);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(0);
                    startAndEnd.push_back(114);
                    startAndEnd.push_back(68);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(783);
                    startAndEnd.push_back(784);
                    startAndEnd.push_back(788);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(62);
                    startAndEnd.push_back(63);
                    startAndEnd.push_back(68);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();
//only once
                    startAndEnd.push_back(788);
                    startAndEnd.push_back(789);
                    startAndEnd.push_back(827);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(23);
                    startAndEnd.push_back(24);
                    startAndEnd.push_back(62);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1586);
                    startAndEnd.push_back(1585);
                    startAndEnd.push_back(1547);

                    startAndEnd.push_back(2351);
                    startAndEnd.push_back(2313);
                    startAndEnd.push_back(2312);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(9);
                    startAndEnd.push_back(10);
                    startAndEnd.push_back(23);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1533);
                    startAndEnd.push_back(1546);
                    startAndEnd.push_back(1547);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();
                    viewer.selected_data_index = 0;
                    viewer.data().clear();
                    viewer.data().set_mesh(currPattern, Fg_pattern_curr);

                    viewer.data().uniform_colors(ambient, diffuse, specular);
                    viewer.data().show_texture = false;
                    viewer.data().set_face_based(false);
                    viewer.data().show_lines = true;

                }else if (garmentExt =="skirt_3" ){
                    simulate = false;
                    adaptionFlag = false;
                    startAndEnd.clear();
                    startAndEnd.push_back(3128);
                    startAndEnd.push_back(3127);
                    startAndEnd.push_back(3141);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1566);
                    startAndEnd.push_back(1564);
                    startAndEnd.push_back(1561);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(3189);
                    startAndEnd.push_back(3187);
                    startAndEnd.push_back(3184);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();
                    startAndEnd.push_back(1505);
                    startAndEnd.push_back(1504);
                    startAndEnd.push_back(1518);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                }
                else if (garmentExt == "skirt_02"){
                    simulate = false;
                    adaptionFlag = false;
                    startAndEnd.clear();
                    startAndEnd.push_back(1392);
                    startAndEnd.push_back(1391);
                    startAndEnd.push_back(1418);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(704);
                    startAndEnd.push_back(703);
                    startAndEnd.push_back(701);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(1445);
                    startAndEnd.push_back(1444);
                    startAndEnd.push_back(1442);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(677);
                    startAndEnd.push_back(669);
                    startAndEnd.push_back(651);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();
                }
                else if (garmentExt == "skirt_1"){
                    simulate = false;
                    adaptionFlag = false;
                    startAndEnd.clear();
                    startAndEnd.push_back(864);
                    startAndEnd.push_back(863);
                    startAndEnd.push_back(866);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(405);
                    startAndEnd.push_back(406);
                    startAndEnd.push_back(404);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(844);
                    startAndEnd.push_back(845);
                    startAndEnd.push_back(843);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();

                    startAndEnd.push_back(425);
                    startAndEnd.push_back(424);
                    startAndEnd.push_back(427);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                    startAndEnd.clear();
                }
                else {
                    simulate = false;
                    adaptionFlag = false;
                    int count =0;
                    int maxSize = cornerSet.size();
                    startAndEnd.clear();
                    vector<vector<int>>bd;
                    igl::boundary_loop(Fg_pattern_curr, bd);
                    insertToStartEnd(startAndEnd, cornerSet, currPattern, Fg_pattern_curr, bd);
//                        smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
//                        smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);


                    startAndEnd.clear();
                }
            }
            if(ImGui::Button("End smooth", ImVec2(-1, 0))) {
                mouse_mode = NONE;
                cout<<"End smoothing selection"<<endl;
                string smoothPatternFile =  "patternComputed_"+avName+"_"+garmentExt+".obj";
                igl::writeOBJ(smoothPatternFile, currPattern, Fg_pattern_curr);
            }


        }
        if (ImGui::CollapsingHeader("Pattern adaption", ImGuiTreeNodeFlags_OpenOnArrow)){
            if(ImGui::Button("Compute adaptation", ImVec2(-1, 0))){

                currPattern = mapFromVg;
                Fg_pattern_curr = mapFromFg;
//                cout<<endl<<currPattern.rows()<<" curr pattern rows, faces  "<<Fg_pattern_curr.rows()<<endl;
                preComputeAdaption();
                if(inverseMap){
                    initialGuessAdaption(currPattern, mapToVg, perfPattVg, Fg_pattern_curr, perfPattFg, symetry, cornerSet,
                                         mapCornerToCorner, halfPatternVertToFullPatternVert.size(), halfPatternVertToFullPatternVert, garment);
                }else if (!symetry){
                    currPattern= mapToVg;
                }else{
                    int halfVert = Vg_pattern_half.rows();
                    currPattern.resize(halfVert, 3);
                    int usedIdx = 0;
                    for(int i=0; i< mapFromVg.rows(); i++){
                        if(fullPatternVertToHalfPatternVert.find(i) != fullPatternVertToHalfPatternVert.end())
                            currPattern.row(fullPatternVertToHalfPatternVert[i]) = mapToVg.row(i);
                            usedIdx++;
                    }
                    Fg_pattern_curr = Fg_pattern_half;
                    cout<<" finished half map  "<<endl;

                }
//                cout<<endl<<currPattern.rows()<<" curr pattern rows, second  faces  "<<Fg_pattern_curr.rows()<<endl;
                viewer.selected_data_index = 2;
                viewer.data().clear();
                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.selected_data_index = 0;
                viewer.data().clear();
//                viewer.data().set_mesh(Vg_pattern_half, Fg_pattern_half);
//
                viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                //remove wireframe
                viewer.data().show_lines = true;

                std::vector<std::vector<int> > boundaryLnew;
                igl::boundary_loop(Fg_pattern_curr, boundaryLnew);
                boundaryL.clear();
                boundaryL = boundaryLnew;
                boundaryLFrom.clear();
                igl::boundary_loop(mapFromFg, boundaryLFrom);

                if(inverseMap) {
                    // corner per boundary contains only corners of the half pattern, but their original id -> always apply map when using these constructs
                    createMapCornersToNewCorner(currPattern, mapToVg, cornerPerBoundary, mapCornerToCorner, boundaryL, halfPatternVertToFullPatternVert, fullPatternVertToHalfPatternVert,
                    symetry, halfPatternFaceToFullPatternFace, fullPatternFaceToHalfPatternFace);
                    cornerVertices = VectorXd::Zero(currPattern.rows());
                    if(symetry){
                        updateCornerUtilsInverse(cornerSet, cornerPerBoundary, seamIdPerCorner, mapCornerToCorner,
                                          cornerVertices, fullPatternVertToHalfPatternVert);
                    }else {
                        updateCornerUtils(cornerSet, cornerPerBoundary, seamIdPerCorner, mapCornerToCorner,
                                          cornerVertices, fullPatternVertToHalfPatternVert);
                    }
                    updateSeamCorner(seamsList, minusOneSeamsList, mapCornerToCorner, boundaryL);
                }
                cout<<"ready for adaption"<<endl;
                igl::writeOBJ("finished_tear_Special_"+avName+"_"+garmentExt+".obj", currPattern, Fg_pattern_curr);
                cout<<"File written to special tear "<<endl;
//                viewer.core().is_animating = true;
//                adaptionFlag = true;
            }
            ImGui::InputDouble("Taylor Lazyness ", &(taylor_lazyness),  0, 0, "%0.2f");
            ImGui::InputDouble("Thereshold Mid  ", &(setTheresholdlMid),  0, 0, "%0.4f");
            ImGui::InputDouble("Thereshold Boundary  ", &(setTheresholdBound),  0, 0, "%0.4f");
            if (ImGui::CollapsingHeader("Constrain seam", ImGuiTreeNodeFlags_OpenOnArrow)){
                if(ImGui::Button("Select boundary vertex ", ImVec2(-1, 0))){
                    mouse_mode = SELECTBOUNDSEAM;
                    viewer.selected_data_index = 1;
                    viewer.data().clear();
                    viewer.selected_data_index = 0;
                    viewer.data().clear();
                    viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                    viewer.data().uniform_colors(ambient, diffuse, specular);
                    viewer.data().show_texture = false;
                    viewer.data().set_face_based(false);
                    //remove wireframe
                    viewer.data().show_lines = true;

                }
                if(ImGui::Button("Confirm Selection ", ImVec2(-1, 0))){
                    mouse_mode = NONE;
                    cout<<constrainedSeamsSingle.first<<" is the seam "<<endl;
                    constrainedSeamsSet.insert(constrainedSeamsSingle);
                    viewer.selected_data_index = 0;
                    viewer.data().clear();
                    viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                    viewer.data().uniform_colors(ambient, diffuse, specular);
                    viewer.data().show_texture = false;
                    viewer.data().set_face_based(false);
                    //remove wireframe
                    viewer.data().show_lines = true;

                }
                if(ImGui::Button("Clear Selection ", ImVec2(-1, 0))){
                    constrainedSeamsSet.clear();
                }


            }
            if(ImGui::Button("Compute first Tear", ImVec2(-1, 0))){
                out<<" init cut "<<endl;
                out<<"taylor_lazyness "<<taylor_lazyness<<endl;
                out<<"setTheresholdlMid "<<setTheresholdlMid<<endl;
                out<<"setTheresholdBound "<<to_string(setTheresholdBound)<<endl;
                out<<"constrained seams: ";

                if(!midFractureForbidden){
                    for(int i =0;i<3; i++){
                        string dir="frac_";
                        if(i==0){
                            dir += "v";
                            igl::jet(normVPattern, 0., 2., colU);
                        }else if (i==1){
                            dir += "u";
                            igl::jet(normUPattern, 0., 2., colU);
                        }else{
                            dir += "stressMax";
                            VectorXd maxNorm(normUPattern.rows());
                            for(int i=0; i<normUPattern.rows(); i++){
                                maxNorm(i)= max(normUPattern(i) ,normVPattern(i));
                            }
                            igl::jet(maxNorm, 0., 2., colU);
                        }

                        viewer.data().set_colors(colU);
                        MatrixXd Ka, Ks, Kd;
                        Ka = viewer.data().F_material_ambient;
                        Kd = viewer.data().F_material_diffuse;
                        Ks = viewer.data().F_material_specular;

                        dir = dir;
                        writeMTL(Ka,Ks, Kd, currPattern, Fg_pattern_curr, garment, avName, 100, dir );

                    }
                }
                for(auto it :  constrainedSeamsSet) out<<it.first<<","<<it.second<<" ";
                out<<endl;

                simulate = false;
                adaptionFlag = false;
                viewer.core().is_animating = false;

                bool fin = false;
                auto copyPattern = mapFromVg;
                set<int> tipVert;
                if(garment == "top"){
                    tipVert.insert(958);
                    tipVert.insert( 726);
                    tipVert.insert( 354);
                    tipVert.insert( 267);
                    tipVert.insert( 35);
                    tipVert.insert( 1045);

                }else if (garmentExt == "skirt_2"){
                    tipVert.insert(31);
                    tipVert.insert(369);
                }
                 pos = computeTear(inverseMap, mapFromVg, currPattern, Fg_pattern_curr, patternEdgeLengths_orig, seamsList ,
                            minusOneSeamsList, boundaryL, fin, cornerPerBoundary, // updated in adaption
                            seamIdPerCorner,
                            cornerVertices, cutPositions, releasedVert, toPattern_boundaryVerticesSet, cornerSet,
                            handledVerticesSet, prevTearFinished, LShapeAllowed,
                            prioInner, prioOuter, taylor_lazyness, mapFromFg, setTheresholdlMid,
                                 setTheresholdBound, fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace,
                                 symetry, tipVert, midFractureForbidden, constrainedSeamsSet);

                 changedPos = pos;
                 cout<<pos<<" Pos was changed"<<endl;

                if( copyPattern != mapFromVg){// rest shape changes if we insert middle cuts
                    vector<vector<int> > vvAdjPatt;
                    igl::adjacency_list( mapFromFg, vvAdjPatt);
                    vector<vector<int>> vfAdjPatternFrom;
                    createVertexFaceAdjacencyList(mapFromFg,vfAdjPatternFrom );
                    for(auto adjPosi: vvAdjPatt[halfPatternVertToFullPatternVert[pos]]){
                        if(copyPattern.row(adjPosi) != mapFromVg.row(adjPosi)){
                            changedPos = adjPosi;
                            cout<<copyPattern.row(adjPosi)<<" old one "<<changedPos <<endl<< mapFromVg.row(adjPosi)<<" new pos"<<endl;
                            updateChangedBaryCoordinates(changedPos, vfAdjPatternFrom);
//
                        }
                    }
                    MatrixXd lengthsOrig;
                    igl::edge_lengths(mapFromVg, mapFromFg, lengthsOrig);
                    patternEdgeLengths_orig = lengthsOrig;
                    cout<<"---------------------"<<endl;
                    cout<<"updated edge lengths"<<endl;
                    viewer.selected_data_index = 0;
                    viewer.data().clear();
                    viewer.data().uniform_colors(ambient, specular, diffuse);
                    viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                }

                if(pos!=-1){
                    viewer.selected_data_index = 2;
                    viewer.data().uniform_colors(ambient, diffuse, specular);

                    viewer.data().set_points(currPattern.row(pos), RowVector3d(1.0, 0.0, 0.0));

                }
                std::vector<std::vector<int> > boundaryLnew;
                igl::boundary_loop(Fg_pattern_curr, boundaryLnew);
                if(boundaryLnew.size() != boundaryL.size()){
                    // TODO CHECK
                    updatePatchId(cutPositions, boundaryLnew , seamsList, minusOneSeamsList,fullPatternVertToHalfPatternVert );
                }
                boundaryL.clear();
                boundaryL = boundaryLnew;

                viewer.core().is_animating = true;
                adaptionFlag = true;

            }

            if(ImGui::Button("Compute further Tear", ImVec2(-1, 0))){
                out<<" further "<<endl;
                simulate = false;
                adaptionFlag = false;
                viewer.core().is_animating = false;
                auto copyPattern = mapFromVg;


                pos = tearFurther(cutPositions, currPattern, Fg_pattern_curr, seamsList, minusOneSeamsList, releasedVert,
                            toPattern_boundaryVerticesSet, boundaryL, cornerSet, handledVerticesSet, prevTearFinished,
                            preferManySmallCuts, LShapeAllowed, patternEdgeLengths_orig, mapFromVg, mapFromFg, prioInner, prioOuter, setTheresholdlMid, setTheresholdBound,
                                  fullPatternVertToHalfPatternVert, halfPatternVertToFullPatternVert, halfPatternFaceToFullPatternFace,
                                  showOnly, forceCut
                            );

                if(pos!=-1){
                    viewer.selected_data_index = 2;
                    viewer.data().clear();
                    viewer.data().set_points(currPattern.row(pos), RowVector3d(1.0, 0.0, 0.0));
                }else{
                    // we are done
                    igl::writeOBJ("adaptation_"+to_string(inverseMap)+"_" + avName+ "_"+garmentExt+".obj" , currPattern, Fg_pattern_curr);

                }

                std::vector<std::vector<int> > boundaryLnew;
                igl::boundary_loop(Fg_pattern_curr, boundaryLnew);

                if(boundaryLnew.size() != boundaryL.size()){
                    cout<<seamsList[1]->getStart1()<<" and end "<<seamsList[1] ->getEndCornerIds().first<<endl;
                    updatePatchId(cutPositions, boundaryLnew,seamsList, minusOneSeamsList, fullPatternVertToHalfPatternVert );
                    cout<<seamsList[1]->getStart1()<<" and after "<<seamsList[1] ->getEndCornerIds().first<<endl;
                    igl::writeOBJ("debugPattern.obj", currPattern, Fg_pattern_curr);
                }
                boundaryL.clear();
                boundaryL= boundaryLnew;

                // also adapt the new edge lengths since we might have changed the rest position when cutting in middl
                bool changeFlag = false;
                if( copyPattern != mapFromVg){// rest shape changes if we insert middle cuts
                    vector<vector<int> > vvAdjPatt;
                    igl::adjacency_list( mapFromFg, vvAdjPatt);
                    vector<vector<int>> vfAdjPatternFrom;
                    createVertexFaceAdjacencyList(mapFromFg,vfAdjPatternFrom );
                    for(auto adjPosi: vvAdjPatt[halfPatternVertToFullPatternVert[pos]]){
                        if(copyPattern.row(adjPosi) != mapFromVg.row(adjPosi)){
                            changedPos = adjPosi;
//                            cout<<copyPattern.row(adjPosi)<<" old one "<<changedPos <<endl<< mapFromVg.row(adjPosi)<<" new pos"<<endl;
                            updateChangedBaryCoordinates(changedPos, vfAdjPatternFrom);
                            // todo uvv should also change

                        }
                    }
                    if(prevTearFinished){
                        // maybe we just cut through, update also all neigh neigh
                        for(auto adjPos: vvAdjPatt[halfPatternVertToFullPatternVert[pos]]){
                            for(auto adjPosi : vvAdjPatt[adjPos]){
                                if(copyPattern.row(adjPosi) != mapFromVg.row(adjPosi)){
                                    changedPos = adjPosi;
//                                    cout<<copyPattern.row(adjPosi)<<" additional change old one "<<changedPos <<endl<< mapFromVg.row(adjPosi)<<" new pos"<<endl;
                                    updateChangedBaryCoordinates(changedPos, vfAdjPatternFrom);

                                }
                            }
                        }
                    }
                    MatrixXd lengthsOrig;
                    igl::edge_lengths(mapFromVg, mapFromFg, lengthsOrig);
                    patternEdgeLengths_orig = lengthsOrig;
                    cout<<"---------------------"<<endl;

                    cout<<"updated edge lengths"<<endl;
                    changeFlag= true;
                    viewer.selected_data_index = 0;
                    viewer.data().clear();
                    viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                }
                viewer.core().is_animating = true;
                adaptionFlag = true;
                showOnly = false;

            }

            if(ImGui::Button("show top ", ImVec2(-1, 0))){
                MatrixXd topV;
                MatrixXi topF;
                igl::readOBJ("top_2d_single.obj", topV, topF);
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(topV, topF);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                //remove wireframe
                viewer.data().show_lines = true;
            }

            if(ImGui::Button("write stress file ", ImVec2(-1, 0))){
                if(!midFractureForbidden){
                    for(int i =0;i<3; i++){
                        string dir="frac_"+to_string(adaptioncount);
                        if(i==0){
                            dir += "v";
                            igl::jet(normVPattern, 0., 2., colU);
                        }else if (i==1){
                            dir += "u";
                            igl::jet(normUPattern, 0., 2., colU);
                        }else{
                            dir += "stressMax";
                            VectorXd maxNorm(normUPattern.rows());
                            for(int i=0; i<normUPattern.rows(); i++){
                                maxNorm(i)= max(normUPattern(i) ,normVPattern(i));
                            }
                            igl::jet(maxNorm, 0., 2., colU);
                        }

                        viewer.data().set_colors(colU);
                        MatrixXd Ka, Ks, Kd;
                        Ka = viewer.data().F_material_ambient;
                        Kd = viewer.data().F_material_diffuse;
                        Ks = viewer.data().F_material_specular;

                        dir = dir;
                        writeMTL(Ka,Ks, Kd, currPattern, Fg_pattern_curr, garment, avName, 100, dir );

                    }
                }
            }
            if(ImGui::Button("Show next ", ImVec2(-1, 0))){
                showOnly = true;
            }

            if(ImGui::Checkbox("Force next cut", &forceCut)){
                out<<" force next cut"<<endl;
            }

            if(ImGui::Button("Finished Tear", ImVec2(-1, 0))){
                igl::writeOBJ("finished_tear_writtenPattern_"+avName+"_"+garmentExt+".obj", currPattern, Fg_pattern_curr);
                cout<<"File written to finished_tear_writtenPattern_ "<<endl;
            }

            if (ImGui::CollapsingHeader("Zip Tears ", ImGuiTreeNodeFlags_OpenOnArrow)){

                if(ImGui::Button("Zip Tears:Check or Next", ImVec2(-1, 0))){
                    out<<"zip"<<endl;
                    zipTears( cutPositions, currPattern, Fg_pattern_curr, mapFromFg, mapFromVg, halfPatternFaceToFullPatternFace, inverseMap, forceClosed);
                    forceClosed = false;
                }
                if(ImGui::Checkbox("Force keep closed ", &forceClosed)){
                    out<<"force closed"<<endl;
                }
                ImGui::InputInt("Force keep closed id", &zipIIdToClose, 0, 0);
                if(ImGui::Button("Force close cut position ", ImVec2(-1, 0))){
                    out<<" forcing zip id "<<zipIIdToClose<<endl;
                    forceZipId(cutPositions, currPattern, Fg_pattern_curr, mapFromFg, mapFromVg,
                               halfPatternFaceToFullPatternFace, inverseMap, true, zipIIdToClose);
                }
            }

            if(ImGui::Checkbox("Allow L-shaped fabric insertion", &LShapeAllowed)){
                out<<"Allow L-shaped"<<endl;
            }
            if(ImGui::Checkbox("Prioritize Inner Cuts", &prioInner)){
                prioOuter = false;
            }
            if(ImGui::Checkbox("Prioritize Outer Cuts", &prioOuter)){
                prioInner = false;
            }
            if(ImGui::Checkbox("Forbid mid fracture ", &midFractureForbidden)){
                out<<"forbid mid frac "<<endl;
            }

            if(ImGui::Button("Remove priorities ", ImVec2(-1, 0))){
                prioInner = false;
                prioOuter = false;
            }

        }
        if (ImGui::CollapsingHeader("Modify adapted Pattern ", ImGuiTreeNodeFlags_OpenOnArrow)){
            bool startSmooth = false;
            bool choosePatchArea = false;
            bool choosePatches = false;
            if(ImGui::Checkbox("Start triangulating", &startTri)) {
                simulate = false;
                adaptionFlag = false;
//                string modifiedPattern = "/Users/annaeggler/Desktop/mappedPattern.obj"; //
//                string modifiedPattern = "/Users/annaeggler/Desktop/writtenPattern_currentFractured.obj"; //
                string modifiedPattern  = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/finished_tear_writtenPattern_"+avName+"_"+garmentExt+".obj";
//                string modifiedPattern  = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/finished_retri_writtenPattern_"+avName+"_"+garment+".obj";
//
                igl::readOBJ(modifiedPattern, currPattern, Fg_pattern_curr);

                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.selected_data_index = 2;
                viewer.data().clear();
                viewer.selected_data_index = 0;
                viewer.data().clear();

                viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                viewer.data().show_lines = true;

            }
            if(ImGui::Checkbox("Start smooth", &startSmooth)) {
                cout<<"Please choose 3 points to smooth between.  "<<endl;
                startAndEnd.clear();
                simulate = false;
                adaptionFlag = false;
                mouse_mode= SELECTVERT;
            }
            if(ImGui::Button("Confirm smooth", ImVec2(-1, 0))) {
                if (startAndEnd.size() == 3) {
                    cout << "Great, you selected " << startAndEnd[0] << " and " << startAndEnd[2]
                         << " as endpoints. Let's get to work on smoothing. " << endl;

                    smoothBetweenVertices(currPattern, Fg_pattern_curr, startAndEnd);
                }
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(currPattern, Fg_pattern_curr);

                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                viewer.data().show_lines = true;
            }
            if(ImGui::Button("End smooth", ImVec2(-1, 0))) {
                mouse_mode = NONE;
                cout<<"End smoothing selection"<<endl;
            }
            if(ImGui::Button("Clip area", ImVec2(-1, 0))){
                std::vector<std::vector<int> > boundaryL_adaptedFromPattern;
                igl::boundary_loop( Fg_pattern_curr, boundaryL_adaptedFromPattern );
                MatrixXi Fg_toPattern = mapToFg; // todo not always
                igl::boundary_loop(Fg_toPattern, boundaryL_toPattern);
                vector<vector<VectorXd>> returnVec;
                clipDifference(boundaryL_adaptedFromPattern,
                               boundaryL_toPattern, currPattern, mapToVg, returnVec);
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
                    if(sum>=100){
                        if(garment=="top" && avName == "CLO_avatar_to_bodyScan_Paola_rem" && patchCounter == 4){
                            cliV.row(6) = currPattern.row(374);
                            cout<<" in ACTION"<<endl;
                        }else  if(garment=="top" && avName == "CLO_avatar_to_bodyScan_Paola_rem" && patchCounter == 5) {
                            cliV.row(3) = currPattern.row(418);
                        }else if (i==5) cout<<"NOT TAKEN"<<endl;
                        igl::writeOBJ("clipper_"+to_string(patchCounter)+".obj", cliV, cliF);
                        patchCounter++;
                        connectedVert.emplace_back(returnVec[i]);
                    }
                }
                cout<<boundaryL_adaptedFromPattern.size()<<" the boundary sizes "<<boundaryL_toPattern.size()<<endl;

                igl::readOBJ("clipper_"+to_string(0)+".obj", Vg_retri, Fg_retri);

                mouse_mode = NONE;
                for(int patchi=0; patchi< patchCounter; patchi++) {
                    MatrixXd Vg_retrii;
                    MatrixXi Fg_retrii;
                    igl::readOBJ("clipper_" + to_string(patchi) + ".obj", Vg_retrii, Fg_retrii);
                    Fg_retri.resize(Fg_retrii.rows(), Fg_retrii.cols());
                    Fg_retri= Fg_retrii;

                    Vg_retri.resize(Vg_retrii.rows(), Vg_retrii.cols());
                    Vg_retri= Vg_retrii;

                    if (!inverseMap ) {
                        mergeTriagulatedAndPattern(connectedVert[patchi],Vg_retri, Fg_retri,
                                                   currPattern,
                                                   Fg_pattern_curr, newFaces, avName, garment, garmentExt);
                    } else {
                        int offset = removePatchVg.rows();
                        for (int i = 0; i < Fg_retri.rows(); i++) {
                            for (int j = 0; j < 3; j++) {
                                Fg_retri(i, j) += offset;
                            }
                        }
                        MatrixXd dupl = removePatchVg;
                        removePatchVg.resize(offset + Vg_retri.rows(), 3);
                        removePatchVg.block(0, 0, offset, 3) = dupl;
                        removePatchVg.block(offset, 0, Vg_retri.rows(), 3) = Vg_retri;

                        MatrixXi duplF = removePatchFg;
                        removePatchFg.resize(duplF.rows() + Fg_retri.rows(), 3);
                        removePatchFg.block(0, 0, duplF.rows(), 3) = duplF;
                        removePatchFg.block(duplF.rows(), 0, Fg_retri.rows(), 3) = Fg_retri;
                        igl::writeOBJ("removedAreas_" + avName + "_" + garmentExt + ".obj", removePatchVg, removePatchFg);

                    }
                }
                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                MatrixXd C = MatrixXd::Zero(Fg_pattern_curr.rows(), 3);
                C.col(1).setConstant(1);
                C.col(0).setConstant(1);
                for(auto i: newFaces){
                    C(i, 0) =0;
                }
                viewer.data().set_colors(C);
                patchcount++;

            }

            if(ImGui::Button("Visualize to be removed Area", ImVec2(-1, 0))) {
                // or use removedPatchVg & Fg
                string prefPattern = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/";
                string addedAreaFile = prefPattern+"removedAreas_"+avName+"_"+garmentExt+".obj";
//                 addedAreaFile =  "/Users/annaeggler/Desktop/mappedPatternWithSmoothedCuts.obj";
                igl::readOBJ(addedAreaFile, currPattern, Fg_pattern_curr);
                if(symetry){
                    MatrixXd duplRemV; MatrixXi duplRemF;
                    duplicatePattern(duplRemV, duplRemF,currPattern, Fg_pattern_curr, T_sym_pattern);
                    igl::writeOBJ("duplicate_PatternRemoved_final_of_"+avName+"_"+garmentExt+".obj" , duplRemV, duplRemF);
                    currPattern.resize(duplRemV.rows(), 3); currPattern = duplRemV;
                    Fg_pattern_curr.resize(duplRemF.rows(), 3); Fg_pattern_curr = duplRemF ;
                }



                viewer.selected_data_index = 0;

                viewer.data().clear();
                viewer.data().show_lines = false;
                viewer.data().set_mesh(currPattern, Fg_pattern_curr);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_texture = false;
                viewer.data().set_face_based(false);
                MatrixXd C = MatrixXd::Zero(Fg_pattern_curr.rows(), 3);

                C.col(0).setConstant(0.8);
                C.col(1).setConstant(0.1);
                viewer.data().set_colors(C);

                MatrixXd Ka, Ks, Kd;
                Ka = viewer.data().F_material_ambient;
                Kd = viewer.data().F_material_diffuse;
                Ks = viewer.data().F_material_specular;
                string dir = "u";
                // 10 for inverse red area
                writeMTL(Ka,Ks, Kd, currPattern, Fg_pattern_curr, garment, avName, 10, dir );

                // make sure the outline pattern is symmetric
                MatrixXi halfFaces (mapToFg.rows()/2, 3);
                halfFaces = mapToFg.block(0,0,mapToFg.rows()/2, 3);
                MatrixXd duplRemV; MatrixXi duplRemF;
                duplicatePattern(duplRemV, duplRemF,mapToVg, halfFaces, T_sym_pattern);
//                igl::writeOBJ("testSym.obj" ,duplRemV ,duplRemF);
                mapToVg.resize(duplRemV.rows(), 3); mapToVg = duplRemV;
                mapToFg.resize(duplRemF.rows(), 3); mapToFg = duplRemF;
                viewer.selected_data_index = 2;
                viewer.data().clear();
                C.resize(mapToFg.rows(), 3);
                C.col(2).setConstant(0);
                C.col(1).setConstant(1);
                C.col(0).setConstant(1);
                viewer.data().set_mesh(mapToVg, mapToFg);
                viewer.data().set_colors(C);
                MatrixXd Kas, Kss, Kds;
                Kas = viewer.data().F_material_ambient;
                Kds = viewer.data().F_material_diffuse;
                Kss = viewer.data().F_material_specular;
                // 11 for inverse base
                writeMTL(Kas , Kss, Kds, mapToVg, mapToFg, garment, avName, 11, dir );

                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.data().show_lines = false;


                MatrixXi Fg_toPattern = mapToFg; //mapFromFg; // todo not always
                MatrixXd Vg_toPattern = mapToVg; // mapFromVg; // todo not always

                MatrixXi boundaryOfToPattern;
                computeBoundaryEdges(Fg_toPattern, boundaryOfToPattern);

                viewer.data().set_edges(Vg_toPattern, boundaryOfToPattern, Eigen::RowVector3d(0, 0, 1));
                string objfileName = "outlineRemoval2D_"+avName + "_"+garment + ".obj";
                std::ofstream obj_file(objfileName);

                // Write the vertex positions to the OBJ file
                for (int i = 0; i < Vg_toPattern.rows(); i ++) {
                    obj_file << "v " << Vg_toPattern(i,0) << " " <<  Vg_toPattern(i,1) << " " <<  Vg_toPattern(i,2) << std::endl;
                }

                for (int j = 0; j < boundaryOfToPattern.rows(); j++) {
                    obj_file << "l "<<boundaryOfToPattern(j,0)+1<<" "<<boundaryOfToPattern(j,1)+1<<endl;
                }

                // Close the output file
                obj_file.close();
//
            }
            if(ImGui::Button("End Modification", ImVec2(-1, 0))){
                igl::writeOBJ("finished_retri_writtenPattern_"+avName+"_"+garmentExt +".obj", currPattern, Fg_pattern_curr);
                cout<<" End modification "<<endl;

            }
            if(ImGui::Button("Manual Modification", ImVec2(-1, 0))){
                MatrixXd modif; MatrixXi modifF;
                igl::readOBJ("finished_retri_writtenPattern_"+avName+"_"+garmentExt +".obj", modif, modifF);
                string id ;
                cout<<"Start manual modification"<<endl;
                std::getline(std::cin, id);
                while (id != "q"){
                    std::stringstream ss(id);

                    int face;
                    int idx;
                    int shouldBe;
                    ss >> face;
                    ss>> idx;
                    ss>>shouldBe;
                    int was = modifF(face, idx);
                    modifF(face, idx) = shouldBe;
                    cout<<was<<" is now "<<modifF(face, idx)<<" ok ?"<<endl;
                    string acc;
                    std::getline(std::cin, acc);
                    if(acc== "no"){
                        modifF(face, idx) = was;
                    }
                    std::getline(std::cin, id);
                }
                cout<<" End modification "<<endl;
                igl::writeOBJ("finished_retri_writtenPattern_"+avName+"_"+garmentExt +".obj", modif, modifF);

            }
            if(ImGui::Button("Manual face insertion", ImVec2(-1, 0))){
                MatrixXd modif; MatrixXi modifF;
                igl::readOBJ("finished_retri_writtenPattern_"+avName+"_"+garmentExt +".obj", modif, modifF);
                string id ;
                cout<<"Start manual face insertion"<<endl;
                cout<<"face to be split and the two indices, new vert gets positioin"<<endl;
                int idx1, idx0, newvert;
                std::getline(std::cin, id);

                std::stringstream ss(id);
                int face;
                ss >> face;
                ss>> idx0;
                ss>>idx1;
                ss>>newvert;
                int fsize = modifF.rows();
                int vsize = modif.rows();
                MatrixXd modifNew(vsize+1, 3);
                modifNew.block(0,0,vsize, 3)= modif;
                MatrixXi modifFNew(fsize+1, 3);
                modifFNew.block(0,0,fsize, 3)= modifF;
                modifFNew.row(fsize) = modifF.row(face);
                modifNew.row(vsize) = modif.row(newvert);
                modifFNew(fsize,idx0) = vsize;
                modifFNew(face,idx1) = vsize;

                cout<<" End modification "<<endl;
                igl::writeOBJ("finished_retri_writtenPattern_"+avName+"_"+garmentExt +".obj", modifNew, modifFNew);

            }

        }
        if (ImGui::CollapsingHeader("Manual stitch", ImGuiTreeNodeFlags_OpenOnArrow)) {
            if(ImGui::Button("full  ", ImVec2(-1, 0))){
                MatrixXd modif; MatrixXi modifF;
                vector<vector<int> >  vfAdjmodi;
                igl::readOBJ("stitched3d2.obj", modif, modifF);
                createVertexFaceAdjacencyList(modifF, vfAdjmodi);

                string id ;
                cout<<"Start manual modification"<<endl;
                std::getline(std::cin, id);
                while (id != "q"){
                    std::stringstream ss(id);

                    int face;
                    int idx;
                    int shouldBe;
                    ss >> face;
                    ss>> idx;
                    ss>>shouldBe;
                    int was = modifF(face, idx);
                    modifF(face, idx) = shouldBe;
                    cout<<was<<" is now "<<modifF(face, idx)<<" ok ?"<<endl;
                    string acc;
                    std::getline(std::cin, acc);
                    if(acc== "no"){
                        modifF(face, idx) = was;
                    }else{
                        for(int i = 0; i<vfAdjmodi[was].size(); i++){
                            for(int j=0; j<3; j++){
                                if(modifF(vfAdjmodi[was][i],j)== was){
                                    modifF(vfAdjmodi[was][i],j)= shouldBe;
                                }
                            }
                        }
                    }
                    igl::writeOBJ("stitched3d2.obj", modif, modifF);

                    std::getline(std::cin, id);
                }
                cout<<" End modification "<<endl;
                igl::writeOBJ("stitched3d2save.obj", modif, modifF);
            }
            if(ImGui::Button("single  ", ImVec2(-1, 0))){
                MatrixXd modif; MatrixXi modifF;
                vector<vector<int> >  vfAdjmodi;
                igl::readOBJ("stitched3d2.obj", modif, modifF);
                createVertexFaceAdjacencyList(modifF, vfAdjmodi);

                string id ;
                cout<<"Start manual modification"<<endl;
                std::getline(std::cin, id);
                while (id != "q"){
                    std::stringstream ss(id);

                    int face;
                    int idx;
                    int shouldBe;
                    ss >> face;
                    ss>> idx;
                    ss>>shouldBe;
                    int was = modifF(face, idx);
                    modifF(face, idx) = shouldBe;
                    cout<<was<<" is now "<<modifF(face, idx)<<" ok ?"<<endl;
                    string acc;
                    std::getline(std::cin, acc);
                    if(acc== "no"){
                        modifF(face, idx) = was;
                    }
                    igl::writeOBJ("stitched3d2.obj", modif, modifF);

                    std::getline(std::cin, id);
                }
                cout<<" End modification "<<endl;
                igl::writeOBJ("stitched3d2save.obj", modif, modifF);
            }
        }
        if (ImGui::CollapsingHeader("Inverse direction", ImGuiTreeNodeFlags_OpenOnArrow)) {
            if(ImGui::Button("Map back for pattern ", ImVec2(-1, 0))){
                string modifiedPattern  = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/finished_tear_Special_"+ avName +"_"+ garmentExt +".obj";
                MatrixXd addedFabricPatternVg; MatrixXi addedFabricPatternFg;
                igl::readOBJ(modifiedPattern, addedFabricPatternVg, addedFabricPatternFg);

                mouse_mode = NONE;
                simulate = false;
                adaptionFlag = false; // this is the pattern after the second mapping direction, it is in shape of mapFrom
                currPattern.resize(addedFabricPatternVg.rows(), addedFabricPatternVg.cols());
                currPattern = addedFabricPatternVg;
                Fg_pattern_curr.resize(addedFabricPatternFg.rows(), addedFabricPatternFg.cols());
                Fg_pattern_curr = addedFabricPatternFg;

                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
                backTo3Dmapping(currPattern, Fg_pattern_curr, perfPattVg_orig, perfPattFg_orig, Vg,
                                Fg, adaptedPatternIn3d, adaptedPatternIn3d_faces ,symetry, garment);

                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = true;
                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                igl::writeOBJ("garmentPattern_"+avName+"_"+garmentExt+"_backIn3d.obj", adaptedPatternIn3d, adaptedPatternIn3d_faces);
            }

            if(ImGui::Button("Map back ", ImVec2(-1, 0))){

//                smoothGarmentOutline();
//                smoothGarmentOutline();
//                smoothGarmentOutline();
//                smoothGarment();
//                smoothGarment();

                string modifiedPattern  = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/finished_retri_writtenPattern_"+ avName +"_"+ garmentExt +".obj";
                MatrixXd addedFabricPatternVg; MatrixXi addedFabricPatternFg;
                igl::readOBJ(modifiedPattern, addedFabricPatternVg, addedFabricPatternFg);

                mouse_mode = NONE;
                simulate = false;
                adaptionFlag = false; // this is the pattern after the second mapping direction, it is in shape of mapFrom
                currPattern.resize(addedFabricPatternVg.rows(), addedFabricPatternVg.cols());
                currPattern = addedFabricPatternVg;
                Fg_pattern_curr.resize(addedFabricPatternFg.rows(), addedFabricPatternFg.cols());
                Fg_pattern_curr = addedFabricPatternFg;

                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
                backTo3Dmapping(currPattern, Fg_pattern_curr, perfPattVg_orig, perfPattFg_orig, Vg,
                                Fg, adaptedPatternIn3d, adaptedPatternIn3d_faces ,symetry, garment);

                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = true;
                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                igl::writePLY("finalGarmentPattern_"+avName+"_"+garmentExt+"_backIn3d.ply", adaptedPatternIn3d, adaptedPatternIn3d_faces);
                igl::writeOBJ("finalGarmentPattern_"+avName+"_"+garmentExt+"_backIn3d.obj", adaptedPatternIn3d, adaptedPatternIn3d_faces);

                showMannequin(viewer);
                if(symetry){
                    duplicatePattern(currPattern, Fg_pattern_curr,addedFabricPatternVg, addedFabricPatternFg, T_sym_pattern);
                    igl::writePLY("duplicate_Pattern_final_of_"+avName+"_"+garmentExt+".ply" , currPattern, Fg_pattern_curr);
                    igl::writeOBJ("duplicate_Pattern_final_of_"+avName+"_"+garmentExt+".obj" , currPattern, Fg_pattern_curr);
                    MatrixXd Vg_notMerged ,Vg_notMerged_dupl; MatrixXi Fg_notMerged, Fg_notMerged_dupl;
                    igl::readOBJ("notMergedPatches.obj", Vg_notMerged, Fg_notMerged);

                    duplicatePattern(Vg_notMerged_dupl, Fg_notMerged_dupl,Vg_notMerged, Fg_notMerged, T_sym_pattern);

                    igl::writeOBJ("notMergedPatches_dupl.obj", Vg_notMerged_dupl, Fg_notMerged_dupl);
                    addedSquare(Fg_notMerged_dupl, Vg_notMerged_dupl);

                }
            }
            if(ImGui::Button("Stitch 3D", ImVec2(-1, 0))){
                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
                igl::readPLY("finalGarmentPattern_"+avName+"_"+garmentExt+"_backIn3d.ply", adaptedPatternIn3d, adaptedPatternIn3d_faces);

                string perfPatternFile = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/patternComputed_"+avName+"_"+garmentExt+".obj";
                igl::readOBJ(perfPatternFile, perfPattVg_orig, perfPattFg_orig);

                string mapFromFile = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/finished_retri_writtenPattern_"+avName+"_"+garmentExt+".obj";
                igl::readOBJ(mapFromFile, mapFromVg, mapFromFg);
                currPattern = mapFromVg;
                Fg_pattern_curr = mapFromFg;

                mapToVg =  Vg_pattern_orig ;// curr = the current shape of the garment, something in between
                mapToFg = Fg_pattern_orig ;// the stress is computed between the rest shape and the current, ie mapFromVg and currPattern

                initialGuessAdaption(currPattern, mapToVg, perfPattVg_orig, Fg_pattern_curr, perfPattFg_orig, symetry, cornerSet,
                                         mapCornerToCorner, halfPatternVertToFullPatternVert.size(), halfPatternVertToFullPatternVert, garment);

                igl::readOBJ("duplicate_Pattern_final_of_"+avName+"_"+garmentExt+".obj" , currPattern, Fg_pattern_curr);
                stitchAdapted3D(adaptedPatternIn3d, adaptedPatternIn3d_faces,Fg_pattern_orig, seamsListDupl, mapCornerToCorner, halfPatternVertToFullPatternVert, currPattern, Fg_pattern_curr);// compute adaptation first
                igl::writeOBJ("finalGarmentPattern_"+avName+"_"+garmentExt+"_Merged_In3d.obj", adaptedPatternIn3d, adaptedPatternIn3d_faces);
                igl::writeOBJ("finalGarmentPattern_"+avName+"_"+garmentExt+"_Merged_FaceCorrespWith3DMerged.obj", currPattern, Fg_pattern_curr );
                addedSquare(adaptedPatternIn3d_faces, adaptedPatternIn3d);
                movePatches();

                MatrixXd cornersMat (mapCornerToCorner.size(), 3); int count = 0;
                for(auto it: mapCornerToCorner){
                    cornersMat.row(count) = adaptedPatternIn3d.row(it.second);
                    count ++;
                }
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = true;
                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().set_points(cornersMat, RowVector3d(1.0, 0.0, 0.0));
            }
            if(ImGui::Button("Not Stitched 3D corresp Faces", ImVec2(-1, 0))){
                string finalFacesFile = "finalGarmentPattern_"+avName+"_"+garmentExt+"_Merged_FaceCorrespWith3DMerged.obj";
                MatrixXd addedFabricPatternVg; MatrixXi addedFabricPatternFg;
                igl::readOBJ(finalFacesFile, addedFabricPatternVg, addedFabricPatternFg);
                mouse_mode = NONE;
                simulate = false;
                adaptionFlag = false; // this is the pattern after the second mapping direction, it is in shape of mapFrom
                currPattern.resize(addedFabricPatternVg.rows(), addedFabricPatternVg.cols());
                currPattern = addedFabricPatternVg;
                Fg_pattern_curr.resize(addedFabricPatternFg.rows(), addedFabricPatternFg.cols());
                Fg_pattern_curr = addedFabricPatternFg;

                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
                backTo3Dmapping(currPattern, Fg_pattern_curr, perfPattVg_orig, perfPattFg_orig, Vg,
                                Fg, adaptedPatternIn3d, adaptedPatternIn3d_faces ,false, garment);

                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = true;
                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
                viewer.data().uniform_colors(ambient, diffuse, specular);
                igl::writeOBJ("finalFaces_3D_notMerged_"+avName+"_"+garmentExt+"_Merged_FaceCorrespWith3DMerged.obj",adaptedPatternIn3d, adaptedPatternIn3d_faces) ;

            }
            bool origIn3D = false;
            if(ImGui::Button("Show inserted in 3D ", ImVec2(-1, 0))){
                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
                cout<< newFaces.size()<<" with color??? "<<endl;
                int size; vector<vector<int>> perFaceNewFaces;
                string witCol ;
                std::getline(std::cin, witCol);
                if(witCol!="no"){
                    cout<<"reading new faces from file"<<endl;
//                if(garment == "leggins"|| "top") patchcount = 5;

                string filename = "newFacesAfterPatch_"+avName+"_"+garmentExt+"_"+ to_string(patchcount) +".txt";
                     filename  = "newFacesAfterPatch_"+avName+"_"+garmentExt+"_final" +".txt";
                ifstream in("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/" + filename);
                in>>size;
                for(int i=0; i<size; i++ ){
                    int faceSize; in>>faceSize;
                    vector<int> currF;
                    for(int j=0;j<faceSize; j++){
                        int newFace;
                        in>>newFace;
                        currF.push_back(newFace);
                    }
                    perFaceNewFaces.push_back(currF);
                }
                }
//

                igl::readPLY("finalGarmentPattern_"+avName+"_"+garmentExt+"_backIn3d.ply", adaptedPatternIn3d, adaptedPatternIn3d_faces);
                MatrixXd C = MatrixXd::Zero(adaptedPatternIn3d_faces.rows(), 3);
//                MatrixXd colrScatter;
//                computeCols(size, colrScatter);
                C.col(1).setConstant(1);
                C.col(0).setConstant(1);
                int offset = C.rows()/2;
                Eigen::RowVector3d dodger_blue(30./255., 144./255., 1);

                for(int i=0; i<size; i++){
                    for(int j=0; j<perFaceNewFaces[i].size(); j++){
                        double val = (((double)i)/((double)(size+1)));
                        C.row(perFaceNewFaces[i][j]) = dodger_blue; //colrScatter.row(i);

                        if(symetry){
                            C.row(perFaceNewFaces[i][j]+offset) =  C.row(perFaceNewFaces[i][j]);
                        }
                    }
                }

                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = false;
                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
                viewer.data().set_colors(C);
                // colors are set, export them!
                MatrixXd Ka, Ks, Kd;
                Ka = viewer.data().F_material_ambient;
                Kd = viewer.data().F_material_diffuse;
                Ks = viewer.data().F_material_specular;
                string dir = "u";

// 1 for 3D
                writeMTL(Ka,Ks, Kd, adaptedPatternIn3d, adaptedPatternIn3d_faces, garment, avName, 1, dir );
            }
            if(ImGui::Button("Show initial in 2D", ImVec2(-1, 0))){
                MatrixXd C = MatrixXd::Zero(Fg_pattern_orig.rows(), 3);
                C.col(1).setConstant(1);
                C.col(0).setConstant(1);

                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = false;
                viewer.data().set_mesh(Vg_pattern_orig, Fg_pattern_orig);
                viewer.data().set_colors(C);

                vector<vector<int>> boundaryL_toPattern;
                igl::boundary_loop(Fg_pattern_orig, boundaryL_toPattern);
                viewer.selected_data_index = 2;
                viewer.data().clear();
                int boundVert=0;
                for (auto bli: boundaryL_toPattern) {
                    boundVert+= bli.size();
                }


                MatrixXi boundaryOfToPattern(boundVert, 2);
                MatrixXd vertPoints(boundVert, 3);
                int curr = 0;
                for (auto bli: boundaryL_toPattern) {
                    for (int j = 0; j < bli.size(); j++) {
                        boundaryOfToPattern(curr, 0) = bli[j];
                        boundaryOfToPattern(curr, 1) = bli[(j + 1) % (bli.size())];
                        curr++;
                    }
                }
                viewer.data().set_edges(Vg_pattern_orig, boundaryOfToPattern, Eigen::RowVector3d(0, 0, 1));


            }
            if(ImGui::Button("Show inserted in 2D", ImVec2(-1, 0))){
                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
                cout<< newFaces.size()<<"with color???"<<endl;
                int size; vector<vector<int>> perFaceNewFaces;
                string witCol ;
                std::getline(std::cin, witCol);
                if(witCol!="no"){
                    cout<<"reading new faces from file"<<endl;
                    string filename = "newFacesAfterPatch_"+avName+"_"+garmentExt+"_"+ to_string(patchcount) +".txt";
                    filename  = "newFacesAfterPatch_"+avName+"_"+garmentExt+"_final" +".txt";
                    ifstream in("/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/" + filename);
                    in>>size;
                    for(int i=0; i<size; i++ ){
                        int faceSize; in>>faceSize;
                        vector<int> currF;
                        for(int j=0;j<faceSize; j++){
                            int newFace;
                            in>>newFace;
                            currF.push_back(newFace);
                        }
                        perFaceNewFaces.push_back(currF);
                    }
                }


                igl::readOBJ("duplicate_Pattern_final_of_"+avName+"_"+garmentExt+".obj", adaptedPatternIn3d, adaptedPatternIn3d_faces);
                MatrixXd C = MatrixXd::Zero(adaptedPatternIn3d_faces.rows(), 3);
                C.col(1).setConstant(1);
                C.col(0).setConstant(1);
//                MatrixXd colrScatter;
//                computeCols(size, colrScatter);
                if(witCol !="no"){
                    Eigen::RowVector3d dodger_blue(30./255., 144./255., 1);
                    int offset = C.rows()/2;
                    for(int i=0; i<size; i++){
                        for(int j=0; j<perFaceNewFaces[i].size(); j++){
                            double val = (((double)i)/((double)(size+1)));
                            C.row(perFaceNewFaces[i][j]) = dodger_blue;
//
                            if(symetry){
                                C.row(perFaceNewFaces[i][j]+offset) = C.row(perFaceNewFaces[i][j]);
                            }

                        }
                    }
                }


                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = false;
                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
                viewer.data().set_colors(C);
                // colors are set, export them!
                MatrixXd Ka, Ks, Kd;
                Ka = viewer.data().F_material_ambient;
                Kd = viewer.data().F_material_diffuse;
                Ks = viewer.data().F_material_specular;
                string dir = "u";

                // 3 for 2D
                writeMTL(Ka,Ks, Kd, adaptedPatternIn3d, adaptedPatternIn3d_faces, garment, avName, 3, dir );

                vector<vector<int>> boundaryL_toPattern;
                igl::boundary_loop(adaptedPatternIn3d_faces, boundaryL_toPattern);

                viewer.selected_data_index = 2;
                viewer.data().clear();
                int boundVert=0;
                for (auto bli: boundaryL_toPattern) {
                    boundVert+= bli.size();
                }

                MatrixXi boundaryOfToPattern(boundVert, 2);
                MatrixXd vertPoints(boundVert, 3);
                int curr = 0;


            }

            if(ImGui::Button("Visualize Patch boundary ", ImVec2(-1, 0))) {
                showPatchBoundary = !showPatchBoundary;
                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
//
                fixRafaPattern();

                igl::readPLY("finalGarmentPattern_"+avName+"_"+garmentExt+"_backIn3d.ply", adaptedPatternIn3d, adaptedPatternIn3d_faces);

                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = false;
//                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
//                viewer.data().set_colors(C);

                vector<vector<int>> boundaryL_toPattern;
                igl::boundary_loop(adaptedPatternIn3d_faces, boundaryL_toPattern);
                viewer.selected_data_index = 2;
                viewer.data().clear();
                int boundVert=0;
                for (auto bli: boundaryL_toPattern) {
                    boundVert+= bli.size();
                }

                if(showPatchBoundary){

                    MatrixXi boundaryOfToPattern(boundVert, 2);
                    MatrixXd vertPoints(boundVert, 3);
                    int curr = 0;
                    // Open the output file
                    string objfileName = "outline3D_"+avName + "_"+garment + ".obj";
                    std::ofstream obj_file(objfileName);

                    MatrixXd moved3d = adaptedPatternIn3d;
                    Eigen::VectorXi componentIdPerVert_pattMoved;
                    igl::vertex_components(adaptedPatternIn3d_faces, componentIdPerVert_pattMoved);


                    // Write the vertex positions to the OBJ file
                    for (int i = 0; i < adaptedPatternIn3d.rows(); i ++) {
                        obj_file << "v " << moved3d(i,0) << " " <<  moved3d(i,1) << " " <<  moved3d(i,2) << std::endl;
                    }

                    for (int i=0;i<boundaryL_toPattern.size() ;i++) {
                        vector<int> bli = boundaryL_toPattern[i];

                        for (int j = 0; j < bli.size(); j++) {
                            int next = (j+1)% bli.size();
                            next = bli[next];
                            if((adaptedPatternIn3d.row(next)-adaptedPatternIn3d.row(bli[j])).norm()>45){
                                cout<<bli[j]<<" skipping "<<next<<endl;
                                continue;
                            }
                            if(next==bli[j]){
                                //we connect them
                                next = boundaryL_toPattern[(i+1)% boundaryL_toPattern.size()][0];
                                if((adaptedPatternIn3d.row(next)-adaptedPatternIn3d.row(bli[j])).norm()>45){
                                    continue;
                                }
                            }

                            obj_file << "l "<<bli[j]+1<<" "<<next+1<<endl;
                            boundaryOfToPattern(curr, 0) = bli[j];
                            boundaryOfToPattern(curr, 1) = next;
                            curr++;
                        }
                    }

                    string id ;
                    cout<<"Start manual addition, insert to be connected vertices"<<endl;
                    std::getline(std::cin, id);
                    while (id != "q"){
                        std::stringstream ss(id);

                        int first;
                        int second;
                        ss >> first;
                        ss>> second;
                        cout<<first<<" and  "<<second<<" ok ?"<<endl;
                        string acc;
                        std::getline(std::cin, acc);
                        if(acc== "no"){
                            // do nothing
                        }else{
                            obj_file << "l "<<first+1<<" "<<second +1<<endl;

                        }
                        std::getline(std::cin, id);
                    }
                    cout<<" End modification "<<endl;

                    // Close the output file
                    obj_file.close();
                    viewer.data().set_edges(adaptedPatternIn3d, boundaryOfToPattern, Eigen::RowVector3d(0, 0, 1));

                }

            }
            if(ImGui::Checkbox("Perfect Pattern in 3D", &origIn3D)){
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = true;
                viewer.data().set_mesh(Vg, Fg);
                showMannequin(viewer);
            }
            if(ImGui::Button("Jacobian for target ", ImVec2(-1, 0))){
                string modifiedPattern  = "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/duplicate_Pattern_final_of_"+ avName +"_"+ garmentExt +".ply";
                MatrixXd addedFabricPatternVg; MatrixXi addedFabricPatternFg;
                igl::readPLY(modifiedPattern, addedFabricPatternVg, addedFabricPatternFg);
                MatrixXd adaptedPatternIn3d;
                MatrixXi adaptedPatternIn3d_faces;
                string adaptedPatternIn3dfile =  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/build/finalGarmentPattern_"+avName+"_"+garmentExt+"_backIn3d.ply";
                igl::readPLY(adaptedPatternIn3dfile, adaptedPatternIn3d, adaptedPatternIn3d_faces);
                VectorXd jacUAdapted, jacVAdapted, jacDiffAdapted;
                computeFinalJacobian(addedFabricPatternVg, addedFabricPatternFg, adaptedPatternIn3d, adaptedPatternIn3d_faces, jacUAdapted, jacVAdapted, jacDiffAdapted);

                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().set_mesh(adaptedPatternIn3d, adaptedPatternIn3d_faces);
//                viewer.data().uniform_colors(ambient, diffuse, specular);
                viewer.data().show_lines = false;
                MatrixXd colJac;
                igl::jet(jacUAdapted, 0., 2, colJac);

                viewer.data().set_colors(colJac);

            }

        }
        if (ImGui::CollapsingHeader("Change fit", ImGuiTreeNodeFlags_OpenOnArrow)){
            if(ImGui::Button("Clear vert ", ImVec2(-1, 0))){
                changeFitVert.clear();
                currPattern.resize(Vg.rows(), 3);
                currPattern = Vg;
                Fg_pattern_curr.resize(Fg.rows(), 3);
                Fg_pattern_curr = Fg;
                mouse_mode = CHANGEFIT;

            }
            ImGui::InputDouble("Max dist", &geoDistMax, 0, 0, "%0.2f");
            ImGui::InputDouble("Change ", &geoDistChange, 0, 0, "%0.4f");
            if(ImGui::Checkbox("Change in U", &geoDistU)){}
            if(ImGui::Checkbox("Change in V", &geoDistV)){}

            if(ImGui::Button("Change in Jacobian ", ImVec2(-1, 0))){
                mouse_mode = NONE;
                VectorXd affectedFaces;
                computeAffection(geoDistDist, geoDistMax, Fg_pattern_curr, affectedFaces);
                gar_adapt ->changeFitViaJacobian( geoDistU, geoDistV, geoDistChange, affectedFaces);
                jacobianChanged= true;
                perFaceTargetNorm = gar_adapt->perFaceTargetNorm;
                computeStress(viewer);
            }
            if(ImGui::Button("Compare Patterns", ImVec2(-1, 0))){
                string dir = (geoDistU)? "U" : "V";
                string fileName = "patternComputed_changedFit_" + avName + "_" + garmentExt + "_" +to_string(geoDistMax)+ "_" + to_string(geoDistChange) + "_" + dir + ".obj";
                MatrixXd changedFitGarV;
                MatrixXi changedFitGarF;
                igl::readOBJ(fileName, changedFitGarV, changedFitGarF);
                cout<<"read file previously  written to * patternComputed changed fit *"<<endl;

                viewer.selected_data_index = 1;
                viewer.data().clear();
                viewer.data().show_lines = true;
                MatrixXi boundaryOfToPattern;
                computeBoundaryEdges(changedFitGarF, boundaryOfToPattern);
                viewer.data().set_edges(changedFitGarV, boundaryOfToPattern, Eigen::RowVector3d(1, 0, 1));

                /* visualize the original pattern for comparison */
                igl::readOBJ("patternComputed_"+avName+"_"+garmentExt+".obj",currPattern, Fg_pattern_curr);
                viewer.selected_data_index = 0;
                viewer.data().clear();
                viewer.data().show_lines = true;
                MatrixXi boundaryOfToPatternOrig;
                computeBoundaryEdges(Fg_pattern_curr, boundaryOfToPatternOrig);
                viewer.data().set_edges(currPattern, boundaryOfToPatternOrig, Eigen::RowVector3d(0, 0, 1));

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
MatrixXi EdgesVisFromPattern;
void updateChangedBaryCoordinates(int changedPosition, vector<vector<int>>& vfFromPatt){
    vector<int> adjFace = vfFromPatt[changedPosition];
    for (auto i: adjFace){
        cout<<"updating face "<<i<<endl ;
        int id0 = mapFromFg(i, 0);
        int id1 = mapFromFg(i, 1);
        int id2 = mapFromFg(i, 2);

        Vector2d Gu, Gv, G;
        Vector2d p0, p1, p2;
        p0 = mapFromVg.block(id0, 0, 1, 2).transpose();
        p1 = mapFromVg.block(id1, 0, 1, 2).transpose();
        p2 = mapFromVg.block(id2, 0, 1, 2).transpose();

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

}
void preComputeAdaption(){
    simulate = false;
    stretchStiffnessU = 0.08;
    stretchStiffnessD *= 2;
    boundaryStiffness = 0.999;

    cout<<" ***** start pattern adaptation ****** "<<endl;
    cout<<"mapping the original to the computed pattern "<<endl;
    if(cornerPerBoundary.empty()){
        cout<<" there are no corners to map"<<endl;
    }

//    patternEdgeLengths.resize(Fg_pattern.rows() ,3);
    igl::edge_lengths(mapFromVg,mapFromFg, patternEdgeLengths_orig);

    igl::boundary_loop(mapToFg, boundaryL_toPattern);
    int boundVert = 0;
    for(auto bli : boundaryL_toPattern){
        boundVert += bli.size();
        for(auto blij : bli){
            toPattern_boundaryVerticesSet.insert(blij);
        }
    }
    EdgesVisFromPattern.resize(boundVert, 2);
    int curr = 0;
    for (auto bli: boundaryL_toPattern){
        for(int j=0; j<bli.size(); j++){
            EdgesVisFromPattern(curr, 0) = bli[j];
            EdgesVisFromPattern(curr, 1) = bli[(j + 1) % (bli.size())];
            curr++;
        }
    }

    // use with PBD_adaption
    baryCoordsUPattern.resize(mapFromFg.rows(), 3);
    baryCoordsVPattern.resize(mapFromFg.rows(), 3);
    for(int i=0; i<Fg_pattern_curr.rows(); i++){
        int id0 = mapFromFg(i, 0);
        int id1 = mapFromFg(i, 1);
        int id2 = mapFromFg(i, 2);

        Vector2d Gu, Gv, G;
        Vector2d p0, p1, p2;
        p0 = mapFromVg.block(id0, 0, 1, 2).transpose();
        p1 = mapFromVg.block(id1, 0, 1, 2).transpose();
        p2 = mapFromVg.block(id2, 0, 1, 2).transpose();

        G = (1./3.) * p0 + (1./3.) * p1 + (1./3.) * p2;
        Gu = G; Gu(0) += 1;
        Gv = G; Gv(1) += 1;

        Vector3d uInBary, vInBary;
        MathFunctions mathFun;
        mathFun.Barycentric(Gu, p0, p1, p2, uInBary);
        mathFun.Barycentric(Gv, p0, p1, p2, vInBary);

        baryCoordsUPattern.row(i) = uInBary.transpose();
        baryCoordsVPattern.row(i) = vInBary.transpose();

    }

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
void visualizeSeam(pair<int, int> which, igl::opengl::glfw::Viewer& viewer){
    int seamId = which.first;
    int seamType = which.second;
    MatrixXi edgesMat;
    vector<vector<int>> boundaryUsed;
    igl::boundary_loop( mapFromFg, boundaryUsed);
    if(seamType ==1 ) {
        int whichSeam = (seamId<0)? (seamId+1)*(-1): seamId;
        for(int j=whichSeam; j<whichSeam+1; j++){
            seam* firstSeam = seamsList[j];
            auto stP1 = firstSeam-> getStartAndPatch1();
            auto stP2 = firstSeam-> getStartAndPatch2ForCorres();
            int len = firstSeam -> seamLength();
            int boundLen1; bool es1 =  (stP1.second >=boundaryUsed.size() )? false : true;
            if(es1 ) {
                boundLen1 = boundaryUsed[stP1.second].size();
            }
            int boundLen2; bool es2 =  (stP2.second >=boundaryUsed.size() )? false : true;
            if(es2){
                boundLen2 = boundaryUsed[stP2.second].size();
            }
            if(!es1 && !es2)continue;
            int size ;
            if(!es1 || !es2) {size = len; }else size = 2*len;
            edgesMat.resize(size, 2);
            for(int i=0; i<=len; i++){
                if(es1) {
                    int next = (stP1.first + i) % boundLen1;
                    if (i != 0) edgesMat(2 * (i - 1), 1) = boundaryUsed[stP1.second][next];
                    if (i != len)edgesMat(2 * i, 0) = boundaryUsed[stP1.second][next];
                }
                if(es2){
                int setAccess = (stP2.first-i)% boundLen2;
                if(setAccess < 0) {
                    setAccess +=boundLen2;
                }
                if(seamsList[j]->inverted) setAccess = (stP2.first + i) % boundLen2;
                if (i!= 0) edgesMat(2*i-1, 1) = boundaryUsed[stP2.second][setAccess];
                if (i!= len)edgesMat(2*i+1, 0) = boundaryUsed[stP2.second][setAccess];
                }
            }
        }
    }else{
        int whichSeam = seamId;
        for(int j = whichSeam; j<whichSeam+1; j++) {
            minusOneSeam *currSeam = minusOneSeamsList[j];
            int patch = currSeam->getPatch();
            bool es1 =  (patch >=boundaryUsed.size() )? false : true;
            if(!es1) continue;
            int idx = currSeam->getStartIdx();
            int len = currSeam->getLength();
            edgesMat.resize(2 * (len), 2);
            int ps = boundaryUsed[patch].size();
            for (int i = 0; i <= len; i++) {

                int next = (idx + i) % ps;
                if (i != 0) edgesMat(2 * (i - 1), 1) = boundaryUsed[patch][next];
                if (i != len)edgesMat(2 * i, 0) = boundaryUsed[ps][next];
            }
        }
    }

    viewer.selected_data_index = 1;
    viewer.data().clear();
    viewer.data().set_edges(currPattern, edgesMat, Eigen::RowVector3d(1, 0, 0));

    viewer.selected_data_index = 0;
    viewer.data().clear();
    viewer.data().set_mesh(currPattern, Fg_pattern_curr);
    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_texture = false;
    viewer.data().set_face_based(false);
    //remove wireframe
    viewer.data().show_lines = false;
}
void findSeamOfVert(int v, igl::opengl::glfw::Viewer& viewer){
    bool found = false; int seamId, seamType;
    vector<vector<int>> boundaryUsed;
    igl::boundary_loop( mapFromFg, boundaryUsed);
    for(int j=0; j<seamsList.size(); j++ ) {
        seam* firstSeam = seamsList[j];
        auto stP1 = firstSeam-> getStartAndPatch1();
        auto stP2 = firstSeam-> getStartAndPatch2ForCorres();
        int len = firstSeam -> seamLength();

        int boundLen1 = boundaryUsed[stP1.second].size();
        int boundLen2 = boundaryUsed[stP2.second].size();

        for(int i=0; i<=len; i++){
            if(boundaryUsed[stP1.second][(stP1.first+i)% boundLen1] == v) {
                found = true;
                seamId = j; seamType = 1;
                constrainedSeamsSingle =make_pair(seamId, seamType);

            }
            int setAccess = (stP2.first-i)% boundLen2;
            if(setAccess < 0) {
                setAccess +=boundLen2;
            }
            if(boundaryUsed[stP2.second][setAccess] == v) {
                found = true;

                seamId = (j+1)*(-1);
                seamType = 1;
                constrainedSeamsSingle =make_pair(seamId, seamType);

            }
        }
    }
    for(int j=0; j<minusOneSeamsList.size(); j++){
        minusOneSeam *currSeam = minusOneSeamsList[j];
        int patch = currSeam->getPatch();
        int idx =  currSeam->getStartIdx();
        int len = currSeam->getLength();
        int ps = boundaryUsed[patch].size();
        for(int i=0; i<len; i++){
            if(boundaryUsed[patch][(idx+i)% ps] == v ){
                found = true;
                seamType = -1;
                seamId = j ;

                constrainedSeamsSingle = make_pair(seamId, seamType);

            }
        }
    }
    if(found){
        visualizeSeam(constrainedSeamsSingle, viewer);
    }
}
bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier){
    if(mouse_mode == SELECTBOUNDARY || mouse_mode == SELECTPATCH ){
        if (button == (int)igl::opengl::glfw::Viewer::MouseButton::Right)
            return false;

        int fid;
        Eigen::Vector3d b;
        MatrixXd Vrs = currPattern;

        if (computePointOnMesh(viewer, Vrs, Fg_pattern_curr, b, fid)) {
            int v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern_curr);
//        cout<<v_id<<"computed closest pattern id "<<endl ;
            viewer.data().set_points(Vrs.row(v_id), RowVector3d(1.0, 0.0, 0.0));
            whichPatchMove = componentIdPerVert(v_id);
            // TODO now we could constrain the whole patch or the boundary
            return true;
        }
        return false;
    }
    if(mouse_mode == SELECTVERT){
        int fid;
        Eigen::Vector3d b;
        MatrixXd Vrs = currPattern;

        if (computePointOnMesh(viewer, Vrs, Fg_pattern_curr, b, fid)) {
            int v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern_curr);
            viewer.data().set_points(Vrs.row(v_id), RowVector3d(1.0, 0.0, 0.0));
            cout<<"Selected vertex "<<v_id<<endl;
             startAndEnd.push_back(v_id);

            return true;
        }
    }
    if(mouse_mode == CHANGEFIT){
        int fid;
        Eigen::Vector3d b;
        // it is in the from mesh, thus snap to the closest vertex on the mesh
        if (computePointOnMesh(viewer, currPattern, Fg_pattern_curr, b, fid)) {
            int v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern_curr);

            changeFitVert.push_back(v_id);
            int half = Fg_pattern_curr.rows()/2;
            if(v_id == Fg_pattern_curr(fid, 0)){
                changeFitVert.push_back(Fg_pattern_curr((fid + half)% (2*half), 0) );
            }else if (v_id == Fg_pattern_curr(fid, 1)){
                changeFitVert.push_back(Fg_pattern_curr((fid + half)% (2*half), 2) );
            }else {
                changeFitVert.push_back(Fg_pattern_curr((fid + half)% (2*half), 1) );

            }
            MatrixXd setPointsMatrix (changeFitVert.size(), 3);
            int rowIdx = 0;
            Eigen::VectorXi VS,FS,VT,FT;
            VS.resize(changeFitVert.size());
            for( auto it: changeFitVert){
                VS(rowIdx) = it;
                setPointsMatrix.row(rowIdx) = currPattern.row(it);
                rowIdx++;
            }

            // All vertices are the targets
            VT.setLinSpaced(currPattern.rows(),0,currPattern.rows()-1);

            igl::exact_geodesic(currPattern,Fg_pattern_curr,VS,FS,VT,FT,geoDistDist);
            Eigen::MatrixXd CM;
            igl::parula(Eigen::VectorXd::LinSpaced(21,0,1).eval(),false,CM);
            igl::isolines_map(Eigen::MatrixXd(CM),CM);
            viewer.data().set_colormap(CM);
            viewer.data().set_data(geoDistDist);
        }
    }
    if(mouse_mode == SELECTAREA){
        int fid, v_id, whichMesh;
        Eigen::Vector3d b;
        MatrixXd Vrs = currPattern;
        // it is in the from mesh, thus snap to the closest vertex on the mesh
        if (computePointOnMesh(viewer, Vrs, Fg_pattern_curr, b, fid)) {
             v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern_curr);
            viewer.data().set_points(Vrs.row(v_id), RowVector3d(1.0, 1.0, 0.0));
            whichMesh=1;
            cout<<"found on inner mesh"<<endl;
            polylineSelected.emplace_back(Vrs.row(v_id));
            polylineIndex.push_back(v_id);

            // the vertex is not in our original fromMesh. Locate it in the toMesh and use the exackt mouse chosen position by the barycentric coordinates to set its position
        }else{
            Vrs = mapToVg; // TODO CHANGE
            MatrixXi Frs = mapToFg;
            if (computePointOnMesh(viewer, Vrs, Frs, b, fid)) {
                cout<<"on outer"<<endl;
                VectorXd chosen = b(0) * Vrs.row(Frs(fid ,0)) + b(1) * Vrs.row(Frs(fid, 1))+ b(2) * Vrs.row(Frs(fid, 2));

                viewer.data().set_points(chosen.transpose(), RowVector3d(.0, 1.0, 0.0));
                cout<<"Vertex from toPattern was chosen"<<endl;
                whichMesh = 2;
                polylineSelected.push_back(chosen);
                polylineIndex.push_back(fid);
            }
        }

        polyLineMeshIndicator.push_back(whichMesh);
        if(  polylineSelected.size() % 2 != 0){
            cout<<"please select endpoint from same mesh "<<endl;
        }else{
            cout<<"finished, please confirm "<<endl;
        }
    }
    if(mouse_mode == SELECTAREAPATCHES){
        int fid, v_id, whichMesh;
        Eigen::Vector3d b;
        MatrixXd Vrs = currPattern;
        // it is in the from mesh, thus snap to the closest vertex on the mesh
        if (computePointOnMesh(viewer, Vrs, Fg_pattern_curr, b, fid)) {
            v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern_curr);
//            viewer.data().set_points(Vrs.row(v_id), RowVector3d(1.0, 1.0, 0.0));
            whichMesh=1;
            cout<<"found on inner mesh"<<endl;
            polylineSelected.emplace_back(Vrs.row(v_id));
            polylineIndex.push_back(v_id);


            // the vertex is not in our original fromMesh. Locate it in the toMesh and use the exackt mouse chosen position by the barycentric coordinates to set its position
        }else{
            Vrs = mapToVg; // TODO CHANGE
            MatrixXi Frs = mapToFg;
            if (computePointOnMesh(viewer, Vrs, Frs, b, fid)) {
                v_id = computeClosestVertexOnMesh(b, fid, Frs);
                cout<<"on outer"<<endl;

//                viewer.data().set_points(Vrs.row(v_id), RowVector3d(.0, 1.0, 0.0));
                whichMesh = 2;
                polylineSelected.emplace_back(Vrs.row(v_id));
                polylineIndex.push_back(v_id);
            }
        }

        polyLineMeshIndicator.push_back(whichMesh);
        if(  polylineSelected.size() % 12 != 0){
            cout<<"please select "<<12 - (polylineSelected.size() % 12) <<" more points according to the structure "<<endl;
        }else{
            cout<<"finished, please confirm "<<endl;
        }
    }
    if(mouse_mode == FINALSTITCH){
        int fid;
        Eigen::Vector3d b;
        MatrixXd Vrs = currPattern;
        if (computePointOnMesh(viewer, Vrs, Fg_pattern_curr, b, fid)) {
            int v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern_curr);
            viewer.data().set_points(Vrs.row(v_id), RowVector3d(1.0, 0.0, 0.0));
            cout<<"Selected vertex "<<v_id<<endl;

            startAndEnd.push_back(v_id);
            return true;
        }

    }
    if(mouse_mode == SELECTBOUNDSEAM){
        int fid;
        Eigen::Vector3d b;
        MatrixXd Vrs = currPattern;

        if (computePointOnMesh(viewer, Vrs, Fg_pattern_curr, b, fid)) {
            int v_id = computeClosestVertexOnMesh(b, fid, Fg_pattern_curr);
            viewer.data().set_points(Vrs.row(v_id), RowVector3d(1.0, 0.0, 0.0));
            cout<<"Selected vertex "<<v_id<<endl;
            findSeamOfVert(v_id, viewer);

            return true;
        }
    }
    return false;
}

void computeBaryCoordsGarOnNewMannequin(igl::opengl::glfw::Viewer& viewer, bool useOtherModel, MatrixXd& otherM){
    VectorXd distVec(garmentPreInterpol.rows());
    MatrixXd targetM;
    if(useOtherModel){targetM = otherM;
    }else {targetM =  testMorph_V1; }
    constrainedVertexIds.clear();
    vector<vector<int> > vvAdj, vfAdj;
    igl::adjacency_list(Fg,vvAdj);
    createVertexFaceAdjacencyList(Fg, vfAdj);
    int boundarycount = 0;

    igl::signed_distance(garmentPreInterpol, mannequinPreInterpol, Fm, igl::SIGNED_DISTANCE_TYPE_UNSIGNED, distVec, closestFaceId, C, N);
    N.resize(Vg.rows(), 3);

    VectorXd vis (Vg.rows()); vis.setConstant(0);
    for(int ii=0; ii<Fg.rows(); ii++){
        for(int j=0; j<3; j++) {
            int i = Fg(ii, j);
            if(vis(i)!= 0 ) continue;
            vis(i) ++;
            int closestFace = closestFaceId(i);

            Vector3d a = mannequinPreInterpol.row(Fm(closestFace, 0));
            Vector3d b = mannequinPreInterpol.row(Fm(closestFace, 1));
            Vector3d c = mannequinPreInterpol.row(Fm(closestFace, 2));
//don't quite understand why it is not the normal of the new mannequin
            Vector3d normalVec = FN_m.row(closestFace);// N.row(i);

            Vector3d currVert = garmentPreInterpol.row(i) - (distVec(i) * normalVec).transpose();
            MatrixXd input(1, 3);
            input.row(0) = currVert;
            MatrixXd Bary;
            igl::barycentric_coordinates(input, mannequinPreInterpol.row(Fm(closestFace, 0)), mannequinPreInterpol.row(Fm(closestFace, 1)),
                                         mannequinPreInterpol.row(Fm(closestFace, 2)), Bary);

            Vector3d currInBary = Bary.row(0);
            if(!useOtherModel) {
                if (isBoundaryVertex(garmentPreInterpol, i, vvAdj, vfAdj)) {
                    constrainedVertexIds.emplace_back(i); // (i)= 1;
                    boundarycount++;
                    constrainedVertexBarycentricCoords.emplace_back(std::make_pair(currInBary, closestFace));
                    constrainedVertexDistance.push_back(distVec(i));
                }
                allVertexBarycentricCoords.emplace_back(std::make_pair(currInBary, closestFace));
            }
            a = targetM.row(Fm(closestFace, 0));
            b = targetM.row(Fm(closestFace, 1));
            c = targetM.row(Fm(closestFace, 2));
            N.row(i) = ((b - a).cross(c - a)).normalized();
            Vector3d newPos = currInBary(0) * a + currInBary(1) * b + currInBary(2) * c;

            double oldx = -1;//todo
//            if(abs(garmentPreInterpol(i, 0))<1){
//               oldx = garmentPreInterpol(i, 0);
//            }
            Vg.row(i) = newPos.transpose() + distVec(i) * N.row(i);

            if(oldx !=-1){
                Vg(i, 0) = oldx ;
            }
        }
    }

    MatrixXd Vgsmooth = Vg;
    vector<vector<int> > vvAdjGar, vfAdjGar;
    igl::adjacency_list(Fg, vvAdjGar);
    createVertexFaceAdjacencyList(Fg, vfAdjGar);

    for(int i = 0; i<garmentPreInterpol.rows(); i++){
        bool isBound = isBoundaryVertex(garmentPreInterpol, i, vvAdjGar, vfAdjGar );
        if(abs(Vg(i,0))<1 && !isBound){
            double zcoord = 0;
            for(int j=0; j<vvAdjGar[i].size(); j++){
                zcoord += Vgsmooth(vvAdjGar[i][j], 2);
            }
            Vg(i,2) = zcoord/vvAdjGar[i].size();
        }

    }

}
void smoothGarment() {
    MatrixXd Vg_dupl = Vg;
    vector<vector<int> > vvAdj, vfAdj;
    igl::adjacency_list(Fg, vvAdj);
    createVertexFaceAdjacencyList(Fg, vfAdj);
    double lamda = 0.1;
    double mu = -0.1;
    int iterations = 300;
    for(int it = 0; it<iterations; it++){
        //shrink
        Vg_dupl = Vg;
        for (int i = 0; i < garmentPreInterpol.rows(); i++) {
            bool isBound = isBoundaryVertex(garmentPreInterpol, i, vvAdj, vfAdj);
            if (!isBound) {// we can smooth
                VectorXd avg = VectorXd::Zero(3);
                int count = 0;
                for (int j = 0; j < vvAdj[i].size(); j++) {
                    avg.transpose() += ( Vg_dupl.row(vvAdj[i][j])- Vg_dupl.row(i));
                    count++;
                }
                avg /= count;
                Vg.row(i) += lamda * avg.transpose();

            }
        }
        Vg_dupl = Vg;
        // inflate
        for (int i = 0; i < garmentPreInterpol.rows(); i++) {
            bool isBound = isBoundaryVertex(garmentPreInterpol, i, vvAdj, vfAdj);
            if (!isBound) {// we can smooth
                VectorXd avg = VectorXd::Zero(3);
                int count = 0;
                for (int j = 0; j < vvAdj[i].size(); j++) {
                    avg.transpose() += ( Vg_dupl.row(vvAdj[i][j])- Vg_dupl.row(i));
                    count++;
                }
                avg /= count;
                Vg.row(i) += mu * avg.transpose();

            }
        }
    }
}
void smoothOutline(MatrixXd & Vpattern ,MatrixXi Fpattern){

    std::vector<std::vector<int> > boundaryLnew;
    igl::boundary_loop(Fpattern, boundaryLnew);
    double lamda = 0.1;
    double mu = -0.1;
    int iterations = 100;

    for(int it = 0; it<iterations; it++){
        //shrink
        MatrixXd Vg_dupl= Vpattern;
        for(int i =0; i<boundaryLnew.size(); i++) {
            Vector3d prev = Vg_dupl.row(boundaryLnew[i][0]);
            for (int j = 1; j <= boundaryLnew[i].size(); j++) {
                int id = j % boundaryLnew[i].size();
                int idn = (j + 1) % (boundaryLnew[i].size());
                Vector3d curr = Vg_dupl.row(boundaryLnew[i][id]);
                Vector3d next = Vg_dupl.row(boundaryLnew[i][idn]);
                auto diff = ((prev + next) / 2).transpose()- Vpattern.row(boundaryLnew[i][id]);
                Vpattern.row(boundaryLnew[i][id]) += lamda * diff;
                prev = curr;
            }
        }
        //inflate
        Vg_dupl= Vpattern;
        for(int i =0; i<boundaryLnew.size(); i++) {
            Vector3d prev = Vg_dupl.row(boundaryLnew[i][0]);
            for (int j = 1; j <= boundaryLnew[i].size(); j++) {
                int id = j % boundaryLnew[i].size();
                int idn = (j + 1) % (boundaryLnew[i].size());
                Vector3d curr = Vg_dupl.row(boundaryLnew[i][id]);
                Vector3d next = Vg_dupl.row(boundaryLnew[i][idn]);
                auto diff = ((prev + next) / 2).transpose()- Vpattern.row(boundaryLnew[i][id]);
                Vpattern.row(boundaryLnew[i][id]) += mu * diff;
                prev = curr;
            }
        }
    }
}
void smoothGarmentOutline(){
    std::vector<std::vector<int> > boundaryLnew;
    igl::boundary_loop(Fg, boundaryLnew);
    double lamda = 0.1;
    double mu = -0.1;
    int iterations = 100;
    for(int it = 0; it<iterations; it++){
        //shrink
        MatrixXd Vg_dupl= Vg;
        for(int i =0; i<boundaryLnew.size(); i++) {
            Vector3d prev = Vg_dupl.row(boundaryLnew[i][0]);
            for (int j = 1; j <= boundaryLnew[i].size(); j++) {
                int id = j % boundaryLnew[i].size();
                int idn = (j + 1) % (boundaryLnew[i].size());
                Vector3d curr = Vg_dupl.row(boundaryLnew[i][id]);
                Vector3d next = Vg_dupl.row(boundaryLnew[i][idn]);
                auto diff = ((prev + next) / 2).transpose()- Vg.row(boundaryLnew[i][id]);
                Vg.row(boundaryLnew[i][id]) += lamda * diff;
                prev = curr;
            }
        }
        //inflate
        Vg_dupl= Vg;
        for(int i =0; i<boundaryLnew.size(); i++) {
            Vector3d prev = Vg_dupl.row(boundaryLnew[i][0]);
            for (int j = 1; j <= boundaryLnew[i].size(); j++) {
                int id = j % boundaryLnew[i].size();
                int idn = (j + 1) % (boundaryLnew[i].size());
                Vector3d curr = Vg_dupl.row(boundaryLnew[i][id]);
                Vector3d next = Vg_dupl.row(boundaryLnew[i][idn]);
                auto diff = ((prev + next) / 2).transpose()- Vg.row(boundaryLnew[i][id]);
                Vg.row(boundaryLnew[i][id]) += mu * diff;
                prev = curr;
            }
        }
    }

    VectorXd distVecNew(Vg.rows());
    VectorXi closestFaceIdNew(Vg.rows());
    MatrixXd Nnew, Cnew;
    igl::signed_distance(Vg, Vm, Fm, igl::SIGNED_DISTANCE_TYPE_UNSIGNED, distVecNew, closestFaceIdNew, Cnew, Nnew);
    constrainedVertexIds.clear();
    constrainedVertexBarycentricCoords.clear();
    constrainedVertexDistance.clear();

    for(int i =0; i<boundaryLnew.size(); i++) {
        for (int j = 0; j < boundaryLnew[i].size(); j++) {
            int id = boundaryLnew[i][j];
            int closestFace = closestFaceIdNew(id);

            Vector3d a = Vm.row(Fm(closestFace, 0));
            Vector3d b = Vm.row(Fm(closestFace, 1));
            Vector3d c = Vm.row(Fm(closestFace, 2));

            Vector3d currVert = Vg.row(id) - (distVecNew(id) * N.row(id).transpose()).transpose();

            MatrixXd input(1, 3);
            input.row(0) = currVert;
            MatrixXd Bary;
            igl::barycentric_coordinates(input, Vm.row(Fm(closestFace, 0)), Vm.row(Fm(closestFace, 1)),
                                         Vm.row(Fm(closestFace, 2)), Bary);
            Vector3d currInBary = Bary.row(0);
            constrainedVertexIds.emplace_back(id); // (i)= 1;

            constrainedVertexBarycentricCoords.emplace_back(std::make_pair(currInBary, closestFace));
            constrainedVertexDistance.push_back(distVecNew(id));

        }
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
double interp=0;
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
//    MatrixXd C = MatrixXd::Zero(Vg.rows(), 3);
//    C.col(1).setConstant(1);
//    C.col(0).setConstant(1);
//            viewer.data().set_colors(C);
//
    if(whichStressVisualize == 1){
        igl::jet(normU, 0., 2., colU);
        viewer.data().set_colors(colU);
        MatrixXd Ka, Ks, Kd;
        Ka = viewer.data().F_material_ambient;
        Kd = viewer.data().F_material_diffuse;
        Ks = viewer.data().F_material_specular;
        string dir = "u";
//        smoothGarment();
//        smoothGarmentOutline();
        writeMTL(Ka,Ks, Kd, Vg, Fg, garment, avName, interp, dir );

        cout<<"finished writing "<<endl;
    }else if (whichStressVisualize == 2){
        igl::jet(normV, 0., 2., colV);
        viewer.data().set_colors(colV);
    }else if (whichStressVisualize == 3 ){
        VectorXd maxNorm(normU.rows());
        for (int i=0; i<normU.rows(); i++){
            maxNorm(i) = max(normU(i), normV(i));
        }
        MatrixXd diffCol;
        igl::jet(maxNorm, 0.0, 2., diffCol);
        viewer.data().set_colors(diffCol);
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
        StressDiffJac = false;
        showGarment(viewer);
        keyRecognition = true;
    }
    if (key == 'V')
    {
        whichStressVisualize = 2;
        StressV = true;
        StressU = false;
        noStress = false;
        StressDiffJac = false;
        showGarment(viewer);
        keyRecognition = true;
    }
    if (key == 'N')
    { // NONE stretch visualization
        whichStressVisualize = 0;
        StressV = false;
        StressU = false;
        noStress = true;
        StressDiffJac = false;
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
        VectorXd pattComp;
        igl::vertex_components(Fg_pattern_curr, pattComp);
//        for(int i=0; i<Vg_pattern.rows(); i++){
//            if(pattComp(i)==1|| pattComp(i)==5 ){
//                currPattern(i,0) -= 70;
//            }else if (pattComp(i)==3|| pattComp(i)==6){
//                currPattern(i,0) += 70;
//            }
//        }
        igl::writeOBJ("writtenPattern_"+avName+"_"+garmentExt+".obj", currPattern, Fg_pattern_curr);
          std::cout<<" Garment file written"<<endl;
    }
    if(key == 'P'){       // Pattern
        keyRecognition = true;
        simulate = false;
        smoothGarmentOutline();

        // we start computing the pattern for the current shape
        Eigen::MatrixXd computed_Vg_pattern= Vg;
        cout<<"start computing the pattern with "<<localGlobalIterations<<" local global iterations"<<endl;
        gar_adapt->performJacobianUpdateAndMerge(Vg, localGlobalIterations, baryCoords1, baryCoords2, computed_Vg_pattern, seamsList, boundaryL ,nonSymSeam);


        if(!jacobianChanged){
            cout<<"pattern written to *patternComputed*"<<endl;
            igl::writeOBJ("patternComputed_"+avName+"_"+garmentExt+".obj", computed_Vg_pattern, Fg_pattern);
            igl::writeOBJ("patternComputed3D_"+avName+"_"+garmentExt+".obj", Vg, Fg);
        }else{
            cout<<"pattern written to * patternComputed_changedFit *"<<endl;
            string dir = (geoDistU)? "U" : "V";
            string fileName = "patternComputed_changedFit_" + avName + "_" + garmentExt + "_" +to_string(geoDistMax)+ "_" + to_string(geoDistChange) + "_" + dir+ ".obj";

            igl::writeOBJ(fileName , computed_Vg_pattern, Fg_pattern);
            igl::writeOBJ("patternComputed3D_"+avName+"_"+garmentExt+".obj", Vg, Fg);// does not change with changed fit !
        }


    }
    if(key == 'B'){     // Bending
        cout<<" writing mapped pattern "<<endl;
 //       setNewMannequinMesh(viewer);

        MatrixXd V_updated;
        interp +=0.1;
        if(interp>1) interp -= 1;
        body_interpolator->interpolateMesh(interp, V_updated);
        viewer.selected_data_index = 0;
        viewer.data().clear();
        viewer.selected_data_index = 2;
        viewer.data().clear();
        viewer.selected_data_index = 1;
        viewer.data().clear();
        viewer.data().set_mesh(V_updated, Fm);
        viewer.data().uniform_colors(ambient_grey, diffuse_grey, specular);
        computeBaryCoordsGarOnNewMannequin(viewer, true, V_updated);
        computeStress(viewer);
        showGarment(viewer);
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
    initCollMeshCall(testMorph_V1left, testMorph_F1left,testMorph_V1right, testMorph_F1right);
    col_tree.init(Vm, Fm);
    igl::per_face_normals(Vm, Fm, FN_m);
    igl::per_vertex_normals(Vm, Fm, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN_m, VN_m);
    igl::per_edge_normals(Vm, Fm, igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN_m, EN_m, E_m, EMAP_m);

}
vector<VectorXd> CleftRight, NleftRight;
void setupCollisionConstraints(){
//    igl::writeOBJ("garment3DonNewAvatar.obj", p, Fg);
    setupCollisionConstraintsCall( collisionVert, pureCollVert, testMorph_V1left, testMorph_F1left, testMorph_V1right, testMorph_F1right,p, numVert, coll_EPS,
                               leftHalfToFullFaceMap, rightHalfToFullFaceMap,CleftRight, NleftRight, closestFaceId, Vm, Fm, Fg);
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
    // idea: have for each vertex a set of faces that it may intersect with
    // check collision for noth sides seperately
    // if closestFaceID not in set of allowed faces (computed in initial guess) (first trial with closestFaceId)

    for(int i=0; i<pureCollVert.size(); i++){
        int j = pureCollVert[i];
        Vector3r deltap0;
        PBD.solve_CollisionConstraint(p.row(j),  CleftRight[i], NleftRight[i], deltap0, coll_EPS, vel.row(j));

        // maybe I should compute the intersection instead of using the closest point C?
//        PBD.solve_CollisionConstraint(p.row(j),  C.row(j), N.row(j), deltap0, coll_EPS, vel.row(j));
        p.row(j) += collisionStiffness * deltap0;
    }
}

void preComputeStretch(){
    MatrixXd temp;
    if(mapToVg.rows()== Vg_pattern.rows()){
        cout<<"changed input pattern"<<endl;
        temp.resize(mapToVg.rows(), 3);
        temp = Vg_pattern;
        Vg_pattern = mapToVg;
    }else{
        cout<<mapToVg.rows()<<" the different number of vertiecs mapTo and Vg_pattern "<< Vg_pattern.rows()<<endl;
    }
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
    cout<< baryCoords1.row(0)<<" pattern row "<<endl;
    if(mapToVg.rows()== Vg_pattern.rows()){
        Vg_pattern = temp;
    }
}
void smoothStress(){
    VectorXd normDupl = normU;

    vector<vector<int>> ffAdj;
    createFaceFaceAdjacencyList(Fg, ffAdj);
    for(int i=0; i<normU.rows(); i ++){
        double summation = 0;
        int size = ffAdj[i].size();
        for(int j =0; j<size; j++){
            summation += normDupl(ffAdj[i][j]);
        }
        normU(i) = summation/size;
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
        normV(j) = (Gv-G).norm();

        double diffU = (normU(j)-perFaceTargetNorm[j].first)/ perFaceTargetNorm[j].first;
        double diffV = (normV(j)-perFaceTargetNorm[j].second)/ perFaceTargetNorm[j].second;
        double y = abs(diffU) + abs(diffV) ;
        colJacDiff.row(j)=  Vector3d (  y,  y, 0.0);

        // this is an experiment
        y = (abs(normV(j)-1)+ abs(normU(j)-1))*3;
        colMixed.row(j) = Vector3d ( 1. + y, 1.- y, 0.0);
    }
    smoothStress();

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
void solveStretchAdaption(){

    MatrixXd correctionTerm = MatrixXd::Zero(currPattern.rows(), 3);
    VectorXd itemCount = VectorXd::Zero(currPattern.rows());
    VectorXd dblA;
    igl::doublearea(currPattern, Fg_pattern_curr, dblA);
//    oneShotLengthSolve( p_adaption,  Fg_pattern_curr, baryCoordsUPattern, baryCoordsVPattern, mapFromVg, mapFromFg);
        // force that pulls back to the original position in fromPattern
    // it does not quite work after tthe 3rd cut. Jacobian seems to be fine but it messes up
    for(int j=0; j< Fg_pattern_curr.rows(); j++){
//        if(dblA(j)<4) {
//            for(int l=0; l<3; l++){
//                itemCount(Fg_pattern_curr(j,l))++;
//            }
//            continue;
//        }
        // check all three angles. if any of them is too small continue

        Eigen::MatrixXd patternCoords(2, 3);
        int i  = (inverseMap) ? j : halfPatternFaceToFullPatternFace[j];
        patternCoords.col(0) = mapFromVg.row(mapFromFg(i, 0)).leftCols(2).transpose();
        patternCoords.col(1) = mapFromVg.row(mapFromFg(i, 1)).leftCols(2).transpose();
        patternCoords.col(2) = mapFromVg.row(mapFromFg(i, 2)).leftCols(2).transpose();
        // where they would go to if no stretch in u
        vector<Vector2r> tar(3);
        vector<Vector2r> tarAngle(3);

        Eigen::MatrixXd targetPositions(2, 3);
        targetPositions.col(0)=  p_adaption.row(Fg_pattern_curr(j, 0)).leftCols(2).transpose() ;
        targetPositions.col(1)=  p_adaption.row(Fg_pattern_curr(j, 1)).leftCols(2).transpose() ;
        targetPositions.col(2)=  p_adaption.row(Fg_pattern_curr(j, 2)).leftCols(2).transpose() ;
        int uOrv = 1;

//        TODO STIFFNESS PARAMETER
        VectorXd thisFaceU = baryCoordsUPattern(i,0) * targetPositions.col(0) +  baryCoordsUPattern(i,1) * targetPositions.col(1) + baryCoordsUPattern(i,2) * targetPositions.col(2) ;
        VectorXd thisFaceV = baryCoordsVPattern(i,0) * targetPositions.col(0) +  baryCoordsVPattern(i,1) * targetPositions.col(1) + baryCoordsVPattern(i,2) * targetPositions.col(2) ;
        VectorXd bary = (targetPositions.col(0) + targetPositions.col(1) + targetPositions.col(2)) / 3;

        PBD_adaption.init_UVStretchPattern( thisFaceU- bary,  thisFaceV - bary, patternCoords,targetPositions,
                                            tar[0], tar[1],tar[2], uOrv,  stretchStiffnessD);
//        if(j==1235) uOrv =11;
        PBD_adaption.init_UVStretchPatternCorrectAngle( thisFaceU- bary,  thisFaceV - bary, patternCoords,targetPositions,
                                                        tarAngle[0], tarAngle[1],tarAngle[2], uOrv,  1);

        for(int l=0; l<3; l++){

            Vector2d dir0 = tarAngle[l] - p_adaption.row(Fg_pattern_curr(j, l)).leftCols(2).transpose() ;
            correctionTerm.row(Fg_pattern_curr(j,l)).leftCols(2) += ( stretchStiffnessD * dir0);
            itemCount(Fg_pattern_curr(j,l))++;
            if(dblA(j)<4) {
//                itemCount(Fg_pattern_curr(j,l))++;
//                continue;
            }
            Vector3d e1 =  p_adaption.row(Fg_pattern_curr(j,l) )- p_adaption.row(Fg_pattern_curr(j,(l + 1) % 3));
            Vector3d e2 =  p_adaption.row(Fg_pattern_curr(j,l)) - p_adaption.row(Fg_pattern_curr(j,(l + 2) % 3));
            auto dot = e1.dot(e2);
            dot /= (e1.norm() * e2.norm());
            //cos angle
            if(dot >= 0.95){
//                itemCount(Fg_pattern_curr(j,l))++;
//                continue;
            }
            // just do this if the pattern is not too small
             dir0 = tar[l] - p_adaption.row(Fg_pattern_curr(j, l)).leftCols(2).transpose() ;
            correctionTerm.row(Fg_pattern_curr(j,l)).leftCols(2) += ( stretchStiffnessU * dir0);
            itemCount(Fg_pattern_curr(j,l))++;

        }

    }
    for(int i=0; i<p_adaption.rows(); i++){
        if(itemCount(i) ==0) continue;
        p_adaption.row(i) += (correctionTerm.row(i))/itemCount(i);

    }

}
MatrixXd colPatternU, colPatternV;
void computePatternStress(MatrixXd& perFaceU_adapt,MatrixXd& perFaceV_adapt ){
    // we compute the stress between the current pattern and the one we started from
    colPatternU.resize(Fg_pattern.rows(), 3);
    colPatternV.resize(Fg_pattern.rows(), 3);
    for(int j=0; j<Fg_pattern.rows(); j++) {
        int id0 = Fg_pattern(j, 0);
        int id1 = Fg_pattern(j, 1);
        int id2 = Fg_pattern(j, 2);

        Vector2d p0 = p_adaption.row(id0).leftCols(2);
        Vector2d p1 = p_adaption.row(id1).leftCols(2);
        Vector2d p2 = p_adaption.row(id2).leftCols(2);

        Vector2d Gu, Gv, G;
        G = (1. / 3.) * p0 + (1. / 3.) * p1 + (1. / 3.) * p2;
        Gu = baryCoordsUPattern(j, 0) * p0 + baryCoordsUPattern(j, 1) * p1 + baryCoordsUPattern(j, 2) * p2;
        Gv = baryCoordsVPattern(j, 0) * p0 + baryCoordsVPattern(j, 1) * p1 + baryCoordsVPattern(j, 2) * p2;

        perFaceU_adapt.row(j) = (Gu - G).transpose();

        perFaceV_adapt.row(j) = (Gv - G).transpose();

        double uNorm = perFaceU_adapt.row(j).norm();
        double vNorm = perFaceV_adapt.row(j).norm();

        double differenceIncrementFactor = 5.;
        double y = (uNorm-1) * differenceIncrementFactor; // to increase differences
        double yV = (vNorm-1) * differenceIncrementFactor; // to increase differences

        colPatternU.row(j) = Vector3d((1.0 + y), (1. - y), 0.0);
        colPatternV.row(j) = Vector3d((1.0 + yV), (1. - yV), 0.0);

    }

}

void solveCornerMappedVertices(){

    for(auto cpb: cornerPerBoundary){
        for(auto cpbj :  cpb ){
            int vertIdx = get<0>(cpbj);
            if(vertIdx < 0){
                vertIdx*= -1;
            }else{
                if(inverseMap && (get<1>(cpbj)==-1)){
                    continue;
                }else if (!inverseMap && fullPatternVertToHalfPatternVert.find(vertIdx) == fullPatternVertToHalfPatternVert.end()){
                    continue;
                }
            }
            int ppos = (get<0>(cpbj) < 0)? get<0>(cpbj)*(-1) : fullPatternVertToHalfPatternVert[vertIdx];
            if(releasedVert.find(ppos)!= releasedVert.end()){
                continue;
            }
            Vector2d newSuggestedPos;
            if(!inverseMap || !symetry){
                newSuggestedPos = mapToVg.row(vertIdx).leftCols(2);
            }else{
                newSuggestedPos = mapToVg.row(get<1>(cpbj)).leftCols(2);

            }

            Vector2d dir = newSuggestedPos - p_adaption.row(ppos).leftCols(2).transpose();
            // TODO PARAMETER

            p_adaption.row(ppos).leftCols(2) += boundaryStiffness * dir;
        }
    }

}
map<int,int> IdMap;

void doAdaptionStep(igl::opengl::glfw::Viewer& viewer){
    Timer t(" Adaption time step ");
    if(adaptioncount ==0){
        for(int i=0; i<currPattern.rows(); i++){
            IdMap[i]=i;
        }
        setUpMap(boundaryL_toPattern , fullPatternVertToHalfPatternVert);

    }


    adaptioncount++;
//    if(adaptioncount>1)return;

//    std::cout<<"-------------- Time Step ------------"<<adaptioncount<<endl;

    // we have no velocity or collision, but we do have stretch, constrainedStretch and bending
    // for the stretch, the current pattern is the reference, then its corners are mapped to another position
// the stretch as a simple solve stretch of the rest position and the stretched position?
// constrained corners as constrained condition: new suggested position= mapped position, then the direction is added to p
// bending: not needed if we don;t touch the z coordinate

    p_adaption.resize(currPattern.rows(), 3);
    p_adaption = currPattern;
    p_adaption.col(2)= Eigen::VectorXd::Ones(currPattern.rows());
    p_adaption.col(2)*= 200;

    changedPos = -1;
//    t.printTime(" init ");
    for(int i=0; i<12; i++){
        solveStretchAdaption();
//        t.printTime(" stretch ");


        // before cutting the boundaries should be the same
        map<int, int> mapUsed = fullPatternVertToHalfPatternVert;
        map<int, int> extFHV = fullPatternVertToHalfPatternVert;
//        t.printTime(" maps ");

        if(symetry && inverseMap&& false) {
            mapUsed = IdMap;
            for(auto item : mapCornerToCorner){
                if(item.second<0){
//                    cout<<"adding negative"<<endl;
                    extFHV[item.second]= (-1)* item.second;
                }
            }
        }
//        t.printTime(" maps made ");

        projectBackOnBoundary( mapToVg, p_adaption, seamsList, minusOneSeamsList, boundaryL_toPattern,
                                boundaryLFrom, releasedVert ,inverseMap,  mapUsed, extFHV);
//        t.printTime(" proj ");

//        ensurePairwiseDist(p_adaption, toPattern, Fg_pattern);
        solveCornerMappedVertices();
//        t.printTime(" corner ");

        // this causes the weired ange issue
//        ensureAngle(p_adaption, mapFromVg, Fg_pattern_curr, mapFromFg);
    }


    currPattern = p_adaption;
    viewer.selected_data_index = 1;
    viewer.data().clear();
    viewer.data().line_width= 30.f;
    MatrixXd visToPattern = mapToVg;
    visToPattern.col(2) *= 1.001;
    viewer.data().set_edges(visToPattern, EdgesVisFromPattern, Eigen::RowVector3d(0, 0, 1));

    viewer.selected_data_index = 0;
    viewer.data().clear();
    viewer.data().set_mesh(currPattern, Fg_pattern_curr);

    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_texture = false;
    viewer.data().set_face_based(false);
    viewer.data().show_lines = true;

    vector<vector<int> > vfAdj;
    createVertexFaceAdjacencyList(Fg_pattern_curr, vfAdj);

    MatrixXd startPerEdge(Fg_pattern_curr.rows(), 3);
    MatrixXd uPerEdge = MatrixXd::Zero(Fg_pattern_curr.rows(), 3);
    MatrixXd vPerEdge = MatrixXd::Zero(Fg_pattern_curr.rows(), 3);
    MatrixXd colvPerEdge(Fg_pattern_curr.rows(), 3);
    MatrixXd coluPerEdge(Fg_pattern_curr.rows(), 3);
    VectorXd area2; igl::doublearea(mapFromVg, mapFromFg, area2);
    normUPattern.resize(Fg_pattern_curr.rows());
    normVPattern.resize(Fg_pattern_curr.rows());
    for(int i=0; i<Fg_pattern_curr.rows(); i++){
        Vector3d v0new = currPattern.row(Fg_pattern_curr(i, 0)).transpose();
        Vector3d v1new = currPattern.row(Fg_pattern_curr(i, 1)).transpose();
        Vector3d v2new = currPattern.row(Fg_pattern_curr(i, 2)).transpose();
        startPerEdge.row(i) = ( (v0new + v1new + v2new)/3).transpose();
        int idx = (inverseMap)? i: halfPatternFaceToFullPatternFace[i];

        /*trial */
        int id0 = mapFromFg(i, 0);
        int id1 = mapFromFg(i, 1);
        int id2 = mapFromFg(i, 2);

        Vector2d Gu, Gv, G;
        Vector2d p0, p1, p2;
        p0 = mapFromVg.block(id0, 0, 1, 2).transpose();
        p1 = mapFromVg.block(id1, 0, 1, 2).transpose();
        p2 = mapFromVg.block(id2, 0, 1, 2).transpose();

        G = (1./3.) * p0 + (1./3.) * p1 + (1./3.) * p2;
        double maxX = max(p0(0), max(p1(0), p2(0) ) );
        double minX = min(p0(0), min(p1(0), p2(0) ) );
        double maxY = max(p0(1), max(p1(1), p2(1) ) );
        double minY = min(p0(1), min(p1(1), p2(1) ) );
        Gu = G; auto GGu = G; GGu (0)+=1; Gu(0) += area2(i); //(maxX - minX);
        Gv = G; auto GGv = G; GGv(1)+=1;  Gv(1) += area2(i); //(maxY - minY);
        double origLenU = area2(i);//(maxX-minX);
        double origLenV = area2(i); // (maxY-minY);
        Vector3d uInBary, vInBary,uInBaryG, vInBaryG;
        MathFunctions mathFun;
        mathFun.Barycentric(Gu, p0, p1, p2, uInBary);
        mathFun.Barycentric(Gv, p0, p1, p2, vInBary);
        mathFun.Barycentric(GGu, p0, p1, p2, uInBaryG);
        mathFun.Barycentric(GGv, p0, p1, p2, vInBaryG);

        Vector3d ubary = baryCoordsUPattern.row(idx );
        Vector3d vbary = baryCoordsVPattern.row(idx);
//
//        uPerEdge.row(i) = (ubary(0) * v0new + ubary(1) * v1new + ubary(2) * v2new).transpose() -  startPerEdge.row(i);
//        vPerEdge.row(i) = (vbary(0) * v0new + vbary(1) * v1new + vbary(2) * v2new).transpose() -  startPerEdge.row(i);
        uPerEdge.row(i) = (uInBary(0) * v0new + uInBary(1) * v1new + uInBary(2) * v2new).transpose() - startPerEdge.row(i);
        vPerEdge.row(i) = (vInBary(0) * v0new + vInBary(1) * v1new + vInBary(2) * v2new).transpose() - startPerEdge.row(i);
        Vector3d red; red(0)=G(0); red(1) = G(1); red(2) = 200;
        normUPattern(i) = ((uInBaryG(0) * v0new + uInBaryG(1) * v1new + uInBaryG(2) * v2new).transpose()- startPerEdge.row(i)).norm();
        normVPattern(i) = ((vInBaryG(0) * v0new + vInBaryG(1) * v1new + vInBaryG(2) * v2new).transpose()- startPerEdge.row(i)).norm();
//        double ulen = uPerEdge.row(i).norm();
        double ulen = uPerEdge.row(i).norm()/ origLenU;
        uPerEdge.row(i) = uPerEdge.row(i).normalized()*(ulen*ulen);
        double vlen = vPerEdge.row(i).norm() / origLenV;
        vPerEdge.row(i) = vPerEdge.row(i).normalized()*(vlen*vlen);

        coluPerEdge(i,0) = ulen;
        coluPerEdge(i, 1) = vlen;

    }

    VectorXd uHelp = coluPerEdge.col(0);
    VectorXd vHelp = coluPerEdge.col(1);
    igl::jet(uHelp, 0.5, 1.5, coluPerEdge);
    igl::jet(vHelp, 0.5, 1.5, colvPerEdge);

    viewer.data().add_edges(startPerEdge, startPerEdge + 3 * uPerEdge, coluPerEdge);
    viewer.data().add_edges(startPerEdge, startPerEdge - 3 * uPerEdge, coluPerEdge);

    viewer.data().add_edges(startPerEdge, startPerEdge + 3 * vPerEdge, colvPerEdge);
    viewer.data().add_edges(startPerEdge, startPerEdge - 3 * vPerEdge, colvPerEdge);

//    adaptionFlag = false;
}
void dotimeStep(igl::opengl::glfw::Viewer& viewer){
    Timer t("Time step ");
    std::cout<<endl;
    std::cout<<"-------------- Time Step ------------"<<timestepCounter<<endl;
    if( !patternExists &&timestepCounter >= 15){
        simulate = false;
        return;
    }
    if( patternExists &&timestepCounter >= 10){
        simulate = false;
        return;
    }
    // the stress is already computed, we can use it here
    Eigen::MatrixXd x_new = Vg;
    p = Vg;

    t.printTime(" set p ");
    // line (5) of the algorithm https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
    // we only use it to add gravity to the system
    vel.col(1) += timestep * w*(-1)*grav*gravityfact;
    t.printTime(" set vel ");

    // (7)
    for (int i = 0; i<numVert; i++){
        p.row(i) = x_new.row(i).array()+ timestep*vel.row(i).array();
    }

    t.printTime(" setup p again ");
    setupCollisionConstraints();
    t.printTime(" setup collision constraints ");
//
    init_stretchUV();
    t.printTime(" setup uv stretch ");
    //(9)-(11), the loop should be repeated several times per timestep (according to Jan Bender)
    for(int i = 0; i < num_const_iterations; i++){
        solveBendingConstraint();
//        t.printTime(" bend   ");
//        cout<<p.row(86)<<" 1"<<endl;

        solveStretchConstraint();
//        t.printTime("  stretch ");
//        cout<<p.row(86)<<" 2"<<endl;

        solveStretchUV();
//        t.printTime("  uv  ");
//        cout<<p.row(86)<<" 3"<<endl;

        solveConstrainedVertices();
//        t.printTime("  constr ");
//        cout<<p.row(86)<<" 4"<<endl;

        /* we precomputed the normal and collision point of each vertex, now add the constraint (instead of checking collision in each iteration
         this is imprecise but much faster and acc to paper works fine in practice*/
        solveCollisionConstraint();
//        cout<<p.row(86)<<" 5"<<endl;

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
