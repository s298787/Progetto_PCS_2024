# pragma once

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace Polygons {
struct Traces {
    unsigned int TracesNumber; //numero di tracce
    vector<vector<Vector3d>> TracesExtremesCoord; //array con le coordinate degli estremi delle tracce
    vector<unsigned int> TracesId; //array con gli id delle tracce
    vector<vector<unsigned int>> TracesFracturesId; //array con gli id delle coppie di fratture che generano le tracce
    vector<vector<unsigned int>> TipsTrue; //lista delle tracce passanti: in posizione 0 nel vector è salvato l'id della traccia mentre in posizione 1 l'id della frattura
    vector<vector<unsigned int>> TipsFalse; //lista delle tracce non passanti: in posizione 0 nel vector è salvato l'id della traccia mentre in posizione 1 l'id della frattura
};
}
