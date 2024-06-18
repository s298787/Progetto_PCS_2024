#include "Fractures.hpp"
#include "Traces.hpp"
#include "PolygonalMesh.hpp"
#include <vector>
#include <Eigen/Eigen>
#include <iostream>

using namespace std;
using namespace Eigen;

// Parte 1
namespace Polygons {
bool importdfn(const string& filename, Fractures& fractures); //funzione che apre e legge il file in input
list<vector<unsigned int>> checkspheres(const Fractures& fractures); //funzione che screma le fratture in base alla distanza reciproca
void tracesfinder(const Fractures& fractures, const list<vector<unsigned int>>& goodcouples, Traces& traces, const double& epsilon); //funzione che trova le tracce
}

namespace Analytics {
Vector4d calcsphere(const vector<Vector3d>& vertex_data); //funzione che calcola baricentro e raggio
double distance(const Vector3d& point1, const Vector3d& point2); //calcola la distanza tra due punti 3d
Vector3d normal(const vector<Vector3d>& vertex_data); //funzione che calcola e normalizza il vettore normale al poligono
bool intersectrettaretta(const double &epsilon, const Vector3d& point, const Vector3d& dir,
                         const Vector3d& point1, const Vector3d& dir1, Vector3d& inter); //funzione che trova il punto di intersezione retta-segmento
bool intersectrettasemiretta(const double& epsilon, const Vector3d& point, const Vector3d& dir,
                             const Vector3d& point1, const Vector3d& dir1, Vector3d& control); //funzione che trova il punto di intersezione semiretta-segmento
double calcolangolo(const Vector3d& v1, const Vector3d& v2, const Vector3d& normal);
void antiorario(vector<Vector3d>& sottopol, const Vector3d& normal); //funzione che memorizza i vertici in senso antiorario
}

namespace OutputFileTools {
using namespace Polygons;
bool printtraces(const string& tracesfileout, const Traces& traces); //funzione che stampa su file i risultati delle tracce
bool printtips(const string& tipsfileout, Traces& traces, const Fractures& fractures); //funzione che stampa su file le tracce divise in passanti e non
}

// Parte 2
namespace MeshLibrary {
using namespace Polygons;
bool meshcalc(const double& epsilon, const Traces& traces, const Fractures& fractures, vector<PolygonalMesh>& mesh, const string& meshfolderout);
//funzione che calcola la mesh e stampa l'output su file
}

namespace Export{
using namespace Polygons;
using namespace MeshLibrary;
void exportMesh(const vector<PolygonalMesh>& mesh, const unsigned int& exportingFracture);
}
