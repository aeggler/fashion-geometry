//
// Created by Anna Maria Eggler on 17.10.22.
//

#ifndef EXAMPLE_ADJACENCY_H
#define EXAMPLE_ADJACENCY_H

#include <vector>
#include <Eigen/Core>
#include <map>

// special methods for this project:
// for each non-boundary edge, we collect all 4 corresponding vertex ids in a specific order
//         x2
//         /\
//        /  \
//     e1/    \e3
//      /  t0  \
//     /        \
//    /    e0    \
//  x0------------x1
//    \          /
//     \   t1   /
//      \      /
//     e2\    /e4
//        \  /
//         \/
//         x3
//
// this one is for the bending energy
void createFacePairEdgeListWith4VerticeIDs(
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const std::vector< std::vector<int> >& vf_adj,
        Eigen::MatrixXi& E4,
        Eigen::MatrixXi& EF6,
        Eigen::MatrixXi& ef_adj
);

// this one is slightly different for the quadratic bending energy
void createFacePairEdgeListWith4VerticeIDs(
        const Eigen::MatrixXi& F,
        Eigen::MatrixXi& E4
);

// lists for each vertex the adjacent faces' indices
void createVertexFaceAdjacencyList(
        const Eigen::MatrixXi &F,
        std::vector< std::vector<int> > &adjecencyList
);

void createFaceEdgeAdjecencyList(
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXi & E,
        const std::vector< std::vector<int> > & vertexFaceAdjecencyList,
        std::vector< std::vector<int> > & faceEdgeAdjecencyList
);

void createVertexEdgeAdjecencyList(
        const Eigen::MatrixXi & E,
        std::vector< std::vector<int> > & vertexEdgeAdjecencyList
);

void createFaceFaceAdjacencyList(
        const Eigen::MatrixXi & F,
        std::vector< std::vector<int> > & faceFaceAdjecencyList
);

int adjacentFaceToEdge(
        const int v1,
        const int v2,
        const int old_face,
        const std::vector< std::vector<int> > & vertexFaceAdjecencyList
);

void adjacentFacesToEdge(
        const int v1,
        const int v2,
        const std::vector< std::vector<int> > & vertexFaceAdjecencyList,
        std::pair<int, int> & faces
);

int adjacentFaceToVertices(
        const int v1,
        const int v2,
        const int v3,
        const std::vector< std::vector<int> > & vertexFaceAdjecencyList
);

bool isBoundaryVertex(
        const Eigen::MatrixXd & V,
        int v,
        const std::vector< std::vector<int> > & vvAdj,
        const std::vector< std::vector<int> > &vfAdj
);

int edgeBetweenVertices(
        int v1,
        int v2,
        const std::vector< std::vector<int> > &veAdj
);
/* A helper function to create a map between a vertex on the pattern to the garment. Note that several vertices can map to the same 3D garment vertex, thus this
 * map is not bijective. The faces, however, have to correspond to each other, also in the order of the vertices for a triangle has to be the same*/
void vertexMapPatternToGarment(
        const Eigen::MatrixXi& Fg_test,
        const Eigen::MatrixXi& Fg_patternTest,
        std::map<int,int>& vertexMapPattToGar
        );
/* A helper function to create a map between a vertex on the garment to its corresponding vertex on a certain patch. Note that a vertex
 * can be on several patches, hence the additional patch id parameter is needed to get a injective map*/
void vertexMapGarmentAndPatchIdToPattern(
        const Eigen::MatrixXi& Fg,
        const Eigen::MatrixXi& Fg_pattern,
        Eigen::VectorXi& componentIdPerVert,
        std::map<std::pair<int, int>,int>& vertexMapGarAndIdToPatch);

#endif //EXAMPLE_ADJACENCY_H
