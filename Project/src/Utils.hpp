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
Vector4d calcsphere(const vector<Vector3d>& vertex_data); //funzione che calcola baricentro e raggio
double distance(const Vector3d& point1, const Vector3d& point2); //calcola la distanza tra due punti 3d
Vector3d normal(const vector<Vector3d>& vertex_data); //funzione che calcola e normalizza il vettore normale al poligono
}

// namespace SortLibrary {
// template<typename T>
// void BubbleSort(vector<T>& data); //ordina il vettore
// }
