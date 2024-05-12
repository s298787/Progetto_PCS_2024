# pragma once

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace Polygons {
struct Traces {
    vector<Matrix3d> TracesExtremesCoord; //array con le coordinate degli estremi delle tracce
    vector<unsigned int> TracesId; //array con gli id delle tracce
    vector<vector<unsigned int>> TracesFracturesId; //array con gli id delle coppie di fratture che generano le tracce
    map<bool, list<unsigned int>> TipsFlag; //mappa che collega il bool che indica se la traccia è passante o meno alle rispettive tracce
};
}
