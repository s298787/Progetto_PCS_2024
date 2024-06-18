#include "Utils.hpp"
#include "Fractures.hpp"
#include "Traces.hpp"
#include "PolygonalMesh.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace Polygons;
using namespace Analytics;
using namespace OutputFileTools;
using namespace MeshLibrary;
using namespace Export;

int main(int argc, char** argv)
{
    chrono::steady_clock::time_point t_begin = chrono::steady_clock::now();
    string filename = "";
    double epsilon = 10*numeric_limits<double>::epsilon();
    if(argc == 1)
        cerr << "Not file" << endl;
    else if (argc == 3)
    {
        istringstream str(argv[1]);
        str >> filename;

        istringstream strTol(argv[2]);
        double tol;
        strTol >> tol;
        epsilon = max(epsilon, tol);
    }

    Fractures dfn;
    Traces traces;
    vector<PolygonalMesh> mesh;

    //Parte 1
    // Lettura file input
    if(importdfn(filename, dfn))
    {
        cout << "File " << filename << " read" << endl;
    }
    if(dfn.CoordVertices.empty())
    {
        cout << "dfn è vuoto" << endl;
    }

    // Calcola tracce
    list<vector<unsigned int>> goodcouples = checkspheres(dfn);
    tracesfinder(dfn, goodcouples, traces, epsilon);

    // Scrittura file output
    string tracesfileout = "traces_" + to_string(dfn.FracturesNumber) + ".txt";
    if(printtraces(tracesfileout, traces)) {
        cout << "File " << tracesfileout << " written successfully" << endl;
    }
    string tipsfileout = "tips_" + to_string(dfn.FracturesNumber) + ".txt";
    if(printtips(tipsfileout, traces, dfn)) {
        cout << "File " << tipsfileout << " written successfully" << endl;
    }


    // Parte 2
    // Calcola mesh
    string meshfileout = "mesh_" + to_string(dfn.FracturesNumber) + ".txt";
    if (meshcalc(epsilon, traces, dfn, mesh, meshfileout)) {
        cout << "File " << meshfileout << " written successfully" << endl;
    }

    // exportMesh(mesh);

    chrono::steady_clock::time_point t_end = chrono::steady_clock::now();
    double tempoTrascorso = chrono::duration_cast<chrono::milliseconds>(t_end-t_begin).count();
    cout << "Tempo trascorso: " << tempoTrascorso << "ms. " << endl;
    return 0;
}
