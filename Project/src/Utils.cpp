#include "Utils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>
#include <algorithm>

using namespace std;
using namespace Eigen;

namespace Polygons {
bool importdfn(const string& filename, Fractures& dfn)
{
    ifstream file(filename);
    if (file.fail())
    {
        cerr << "File not found" << endl;
        return false;
    }

    string header;
    string numberfractures;
    stringstream ss;
    getline(file, header);
    getline(file, numberfractures);
    ss << numberfractures;
    ss >> dfn.FracturesNumber;
    ss.clear();

    dfn.CoordVertices.reserve(dfn.FracturesNumber);
    dfn.Spheres.reserve(dfn.FracturesNumber);
    dfn.Normals.reserve(dfn.FracturesNumber);

    string line;
    unsigned int count = 1;
    while(getline(file, line))
    {
        string token;
        if (count == 1 || count == 3)
        {
            continue;
        }

        istringstream convert(line);
        if (count == 2)
        {
            unsigned int id;
            getline(line, token, '; ');
            ss << token;
            ss >> id;
            dfn.FracturesId.push_back(id);
        }

        count++;
    }

    return true;
}
}
