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
    while (getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        istringstream iss(line);
        string token;

        // Legge il FractureId e il NumVertices
        unsigned int fractureId;
        getline(iss, token, ';');
        ss << token;
        ss >> fractureId;
        ss.clear();
        unsigned int numVertices;
        getline(iss, token);
        ss << token;
        ss >> numVertices;
        ss.clear();

        // Popola FracturesId
        fractures.FracturesId.push_back(fractureId);

        getline(file, line);

        // Popola CoordVertices
        vector<Vector3d> vertex_data;
        vertex_data.resize(numVertices);
        for (unsigned int i = 0; i < 3; ++i)
        {
            getline(file, line);
            istringstream vertex_iss(line);
            for (unsigned int k = 0; k < numVertices - 1; ++k)
            {
                getline(vertex_iss, token, ';');
                ss << token;
                ss >> vertex_data[k](i);
                ss.clear();
            }
            getline(vertex_iss, token);
            ss << token;
            ss >> vertex_data[numVertices - 1](i);
            ss.clear();
        }

        fractures.CoordVertices.push_back(vertex_data);

        // Aggiorniamo la mappa NumVertices
        if(fractures.NumVertices.find(numVertices) == fractures.NumVertices.end())
        {
            fractures.NumVertices.insert({numVertices, {fractureId}});
        }
        else
        {
            fractures.NumVertices[numVertices].push_back(fractureId);
        }
    }
    return true;
}
}

namespace Analytics {
Vector3d trovaBaricentro(const vector<Vector3d>& vertices_data)
{
    Vector3d baricentro{0.0,0.0,0.0};

    for (size_t k = 0; k < vertices_data.size(); ++k)
    {
        baricentro += vertices_data[k];

    }

    baricentro /= vertices_data.size();
    return baricentro;
}
}
