#include "Fractures.hpp"
#include "Traces.hpp"
#include <vector>
#include <Eigen/Eigen>
#include <iostream>

using namespace std;
using namespace Eigen;

namespace Polygons {
bool importdfn(const string& filename, Fractures& fractures); //funzione che apre e legge il file in input
}

namespace Analytics {
Vector3d trovaBaricentro(const vector<Vector3d>& vertices_data); //funzione che calcola le sfere
}
