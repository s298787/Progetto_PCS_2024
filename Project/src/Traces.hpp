# pragma once

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace Polygons {
struct Traces {
    unsigned int TracesNumber; // Numero di tracce
    vector<array<Vector3d, 2>> TracesExtremesCoord; // Array con le coordinate degli estremi delle tracce
    vector<unsigned int> TracesId; // Array con gli id delle tracce
    vector<array<unsigned int, 2>> TracesFracturesId; // Array con gli id delle coppie di fratture che generano le tracce
    vector<double> TracesLengths; // Array con le lunghezze delle tracce
    vector<array<unsigned int, 2>> TipsTrue; // Lista delle tracce passanti: in posizione 0 nel vector è salvato l'id della traccia mentre in posizione 1 l'id della frattura
    vector<array<unsigned int, 2>> TipsFalse; // Lista delle tracce non passanti: in posizione 0 nel vector è salvato l'id della traccia mentre in posizione 1 l'id della frattura
};
}
