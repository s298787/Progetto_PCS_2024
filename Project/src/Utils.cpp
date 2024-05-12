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
bool importdfn(const string& filename, Fractures& fractures)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Errore nell'apertura del file" << endl;
        return false;
    }

    string header;
    string numberfractures;
    stringstream ss;
    getline(file, header);
    getline(file, numberfractures);
    ss << numberfractures;
    ss >> fractures.FracturesNumber;
    ss.clear();

    string line;
    // unsigned int currentFractureId = 0;
    while (getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        istringstream iss(line);
        string token;

        // Leggi il FractureId e il NumVertices
        getline(iss, token, ';');
        unsigned int fractureId = stoi(token);
        getline(iss, token);
        unsigned int numVertices = stoi(token);

        // Popola FracturesId
        fractures.FracturesId.push_back(fractureId);

        getline(file, line);

        // Popola CoordVertices
        vector<Vector3d> vertex_data;
        for (unsigned int i = 0; i < 3; ++i)
        {
            getline(file, line);
            istringstream vertex_iss(line);
            for (unsigned int k = 0; k < numVertices - 1; ++k)
            {
                getline(vertex_iss, token, ';');
                vertex_data[k](i) = stod(token);
            }
            getline(vertex_iss, token);
            vertex_data[numVertices - 1](i) = stod(token);
        }

        fractures.CoordVertices.push_back(vertex_data);

        // Popola NumVertices
        // fractures.NumVertices[currentFractureId] = list<unsigned int>(numVertices);

        // ++currentFractureId;
    }

    // Imposta FracturesNumber
    // fractures.FracturesNumber = fractures.FracturesId.size();

    return true;
}
}

// namespace Polygons {
// bool importdfn(const string& filename, Fractures& dfn)
// {
//     ifstream file(filename);
//     if (file.fail())
//     {
//         cerr << "File not found" << endl;
//         return false;
//     }

//     string header;
//     string numberfractures;
//     stringstream ss;
//     getline(file, header);
//     getline(file, numberfractures);
//     ss << numberfractures;
//     ss >> dfn.FracturesNumber;
//     ss.clear();

//     dfn.CoordVertices.reserve(dfn.FracturesNumber);
//     dfn.Spheres.reserve(dfn.FracturesNumber);
//     dfn.Normals.reserve(dfn.FracturesNumber);

//     string line;
//     unsigned int count = 1;
//     while(getline(file, line))
//     {
//         string token;
//         if (count == 1 || count == 3)
//         {
//             continue;
//         }

//         istringstream convert(line);
//         if (count == 2)
//         {
//             unsigned int id;
//             getline(convert, token, ';');
//             ss << token;
//             ss >> id;
//             dfn.FracturesId.push_back(id);
//             ss.clear();
//         }
//         else
//         {
//             double coord;
//             getline(convert, token, ';');
//             ss << token;
//             ss >> coord;
//         }

//         count++;
//     }

//     return true;
// }
// }
