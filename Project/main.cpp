#include "Utils.hpp"
#include "Fractures.hpp"
#include "Traces.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;
using namespace Polygons;
using namespace Analytics;

int main()
{
    Fractures dfn;
    Traces traces;
    string filename = "/home/alberto/Documenti/Materiale scuola Alberto/Appunti/Programmazione e calcolo scientifico/Consegne/Progetto_PCS_2024/Project/DFN/FR4_data.txt";
    if(importdfn(filename, dfn))
    {
        cout << "File read" << endl;
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

    // if(dfn.NumVertices.empty())
    // {
    //     cout << "La mappa NumVertices è vuota" << endl;
    // }

    // for (Vector4d sp : dfn.Spheres)
    // {
    //     for (unsigned int i = 0; i < 4; ++i)
    //     {
    //         cout << sp(i) << ";";
    //     }
    //     cout << endl;
    // }

    // for (Vector3d n : dfn.Normals)
    // {
    //     for (size_t i = 0; i < n.size(); ++i)
    //     {
    //         cout << n(i) << ";";
    //     }
    //     cout << endl;
    // }

    double epsilon = 10*numeric_limits<double>::epsilon();
    list<vector<unsigned int>> goodcouples = checkspheres(dfn);
    // for (vector<unsigned int> couple : goodcouples)
    // {
    //     cout << couple[0] << ";" << couple [1] << endl;
    // }
    // cout << goodcouples.size() << endl;
    tracesfinder(dfn, goodcouples, traces, epsilon);
    // for (size_t i = 0; i < traces.TracesId.size(); ++i)
    // {
    //     cout << traces.TracesId[i] << ": " << traces.TracesFracturesId[i][0] << "&" << traces.TracesFracturesId[i][1] << endl;
    //     cout << traces.TracesExtremesCoord[i][0].transpose() << " ; " << traces.TracesExtremesCoord[i][1].transpose() << endl;
    // }
    cout << "Passanti: " << endl;
    for (size_t i = 0; i < traces.TipsTrue.size(); ++i) {
        cout << traces.TipsTrue[i][0] << " passa per " << traces.TipsTrue[i][1] << endl;
    }
    cout << "Non passanti: " << endl;
    for (size_t i = 0; i < traces.TipsFalse.size(); ++i) {
        cout << traces.TipsFalse[i][0] << " non passa per " << traces.TipsFalse[i][1] << endl;
    }
    return 0;
}
