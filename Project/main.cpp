#include "Utils.hpp"
#include "Fractures.hpp"
#include "Traces.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace Polygons;
using namespace Analytics;

int main()
{
    Fractures dfn;
    string filename = "/home/alberto/Documenti/Materiale scuola Alberto/Appunti/Programmazione e calcolo scientifico/Consegne/Progetto_PCS_2024/Project/DFN/FR10_data.txt";
    if(importdfn(filename, dfn))
    {
        cout << "File read" << endl;
    }

    if(dfn.CoordVertices.empty())
    {
        cout << "dfn è vuoto" << endl;
    }
    for (vector<Vector3d> matrix : dfn.CoordVertices)
    {
        for (Vector3d coord : matrix)
        {
            for(unsigned int i = 0; i<3; ++i)
            {
                cout << coord(i) << ";";
            }
            cout << endl;
        }
        cout << endl;
    }

    if(dfn.NumVertices.empty())
    {
        cout << "La mappa NumVertices è vuota" << endl;
    }

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
    return 0;
}
