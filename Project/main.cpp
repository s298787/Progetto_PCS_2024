#include "Utils.hpp"
#include "Fractures.hpp"
#include "Traces.hpp"
#include "PolygonalMesh.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;
using namespace Polygons;
using namespace Analytics;
using namespace OutputFileTools;
using namespace MeshLibrary;

int main()
{
    Fractures dfn;
    Traces traces;
    PolygonalMesh mesh;
    string filename = "/home/alberto/Documenti/Materiale scuola Alberto/Appunti/Programmazione e calcolo scientifico/Consegne/Progetto_PCS_2024/Project/DFN/FR3_data.txt";
    if(importdfn(filename, dfn))
    {
        cout << "File " << filename << " read" << endl;
    }

    // if(dfn.CoordVertices.empty())
    // {
    //     cout << "dfn è vuoto" << endl;
    // }
    // for (vector<Vector3d> matrix : dfn.CoordVertices)
    // {
    //     for (Vector3d coord : matrix)
    //     {
    //         for(unsigned int i = 0; i<3; ++i)
    //         {
    //             cout << coord(i) << ";";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    double epsilon = 10*numeric_limits<double>::epsilon();
    list<vector<unsigned int>> goodcouples = checkspheres(dfn);
    tracesfinder(dfn, goodcouples, traces, epsilon);
    // for (size_t i = 0; i < traces.TracesId.size(); ++i)
    // {
    //     cout << traces.TracesId[i] << ": " << traces.TracesFracturesId[i][0] << "&" << traces.TracesFracturesId[i][1] << endl;
    //     cout << traces.TracesExtremesCoord[i][0].transpose() << " ; " << traces.TracesExtremesCoord[i][1].transpose() << endl;
    // }
    // cout << "Passanti: " << endl;
    // for (size_t i = 0; i < traces.TipsTrue.size(); ++i) {
    //     cout << traces.TipsTrue[i][0] << " passa per " << traces.TipsTrue[i][1] << endl;
    // }
    // cout << "Non passanti: " << endl;
    // for (size_t i = 0; i < traces.TipsFalse.size(); ++i) {
    //     cout << traces.TipsFalse[i][0] << " non passa per " << traces.TipsFalse[i][1] << endl;
    // }
    string tracesfileout = "traces_" + to_string(dfn.FracturesNumber) + ".txt";
    if(printtraces(tracesfileout, traces)) {
        cout << "File " << tracesfileout << " written successfully" << endl;
    }
    string tipsfileout = "tips_" + to_string(dfn.FracturesNumber) + ".txt";
    if(printtips(tipsfileout, traces, dfn)) {
        cout << "File " << tipsfileout << " written successfully" << endl;
    }
    meshcalc(traces, dfn, mesh);
    return 0;
}
