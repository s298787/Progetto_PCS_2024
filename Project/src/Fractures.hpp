# pragma once

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace Polygons {
struct Fractures {
    unsigned int FracturesNumber; //numero di fratture nel file
    vector<vector<Vector3d>> CoordVertices; //array con le matrici dei vertici
    vector<unsigned int> FracturesId; //array con gli id delle fratture
    map<unsigned int, list<unsigned int>> NumVertices; //mappa che associa il numero di vertici ai rispettivi id
    vector<Vector4d> Spheres; //array che contiene baricentri e le distanze massime dai vertici dei poligoni
    vector<Vector3d> Normals; //array  che contiene il vettore normale al poligono
};
}
