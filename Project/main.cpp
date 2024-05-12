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

    if(dfn.NumVertices.empty())
    {
        cout << "La mappa NumVertices è vuota" << endl;
    }

    for (vector<Vector3d> matrix : dfn.CoordVertices)
    {
        Vector3d baricentro = trovaBaricentro(matrix);
        for (unsigned int i = 0; i<3; ++i)
        {
            cout << baricentro(i) << ";";
        }
        cout << endl;
    }

    return 0;
}
