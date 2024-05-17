#include "Fractures.hpp"
#include "Traces.hpp"
#include <vector>
#include <Eigen/Eigen>
#include <iostream>

using namespace std;
using namespace Eigen;

namespace Polygons {
bool importdfn(const string& filename, Fractures& fractures); //funzione che apre e legge il file in input
list<vector<unsigned int>> checkspheres(const Fractures& fractures); //funzione che screma le fratture in base alla distanza reciproca
void tracesfinder(const Fractures& fractures, const list<vector<unsigned int>>& goodcouples, Traces& traces, const double& epsilon); //funzione che trova le tracce
}

namespace Analytics {
Vector4d calcsphere(const vector<Vector3d>& vertex_data); //funzione che calcola baricentro e raggio
double distance(const Vector3d& point1, const Vector3d& point2); //calcola la distanza tra due punti 3d
Vector3d normal(const vector<Vector3d>& vertex_data); //funzione che calcola e normalizza il vettore normale al poligono
Vector3d intersectrettaretta(const Vector3d& p1, const Vector3d& d1, const Vector3d& p2, const Vector3d& d2); //funzione che trova il punto di intersezione retta-segmento
int sign(double x); //funzione segno
bool positive(vector<int> checksigns);
}

// namespace SortLibrary {
// template<typename T>
// void BubbleSort(vector<T>& data); //ordina il vettore
// }
