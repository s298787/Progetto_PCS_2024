#include "Utils.hpp"
#include "Fractures.hpp"
#include "Traces.hpp"
#include <Eigen/Eigen>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace Polygons;

int main()
{
    Fractures dfn;
    string filename = "/home/alberto/Documenti/Materiale scuola Alberto/Appunti/Programmazione e calcolo scientifico/Consegne/Progetto_PCS_2024/Project/DFN/FR10_data.txt";
    if(importdfn(filename, dfn))
    {
        cout << "File read" << endl;
    }
    return 0;
}
