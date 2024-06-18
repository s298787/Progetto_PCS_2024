# pragma once

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace Polygons {
struct Fractures {
    unsigned int FracturesNumber; // Numero di fratture nel file
    vector<vector<Vector3d>> CoordVertices; // Array con le matrici dei vertici
    vector<unsigned int> FracturesId; // Array con gli id delle fratture   
    vector<Vector4d> Spheres; // Array che contiene il baricentro (posizioni 0,1,2) e la distanza dai vertici (posizione 3) dei poligoni
    vector<Vector3d> Normals; // Array  che contiene il vettore normale al poligono
};
}
