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
    string filename = "DFN/";
    importdfn(filename, dfn);
    return 0;
}
