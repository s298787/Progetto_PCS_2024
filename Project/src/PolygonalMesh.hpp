# pragma once

#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace MeshTools {
struct PolygonalMesh {
    unsigned int NumberCell0d; // Numero di celle 0d
    vector<unsigned int> IdCell0d; // Array di id delle celle 0d
    vector<Vector3d> CoordCell0d; // Array con le coordinate delle celle 0d

    unsigned int NumberCell1d; // Numero di celle 1d
    vector<unsigned int> IdCell1d; // Array di id delle celle 1d
    vector<vector<unsigned int>> IdVerticesCell1d; // Array con gli id dei vertici delle celle 1d

    unsigned int NumberCell2d; // Numero di celle 2d
    vector<unsigned int> IdCell2d; // Array di id delle celle 2d
    vector<vector<unsigned int>> IdVerticesCell2d; // // Array con gli id dei vertici delle celle 2d
    vector<vector<unsigned int>> IdEdgesCell2d; // // Array con gli id dei lati delle celle 2d
};
}
