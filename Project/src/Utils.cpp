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
    if (!file.is_open())
    {
        cerr << "Errore nell'apertura del file" << endl;
        return false;
    }

    string header;
    string line;
    stringstream ss;
    getline(file, header);
    getline(file, line); //Leggiamo il numero di fratture nel file
    ss << line;
    ss >> fractures.FracturesNumber;
    ss.clear();

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

        // Legge le informazioni sui vertici:
        vector<Vector3d> vertex_data;
        vertex_data.resize(numVertices);
        for (unsigned int i = 0; i < 3; ++i) // Legge per 3 volte la righe con le coordinate
        {
            getline(file, line);
            istringstream vertex_iss(line);
            for (unsigned int k = 0; k < numVertices - 1; ++k)
            {
                getline(vertex_iss, token, ';'); // Separa le singole coordinate
                ss << token;
                ss >> vertex_data[k](i); // Iscatola le coordinate iterando sulle righe di vertex_data
                ss.clear();
            }
            getline(vertex_iss, token); // Legge l'ultima coordinata che ha separatore \n
            ss << token;
            ss >> vertex_data[numVertices - 1](i);
            ss.clear();
        }

        // Popola CordVertices
        fractures.CoordVertices.push_back(vertex_data);

        // Aggiorna la mappa NumVertices
        if(fractures.NumVertices.find(numVertices) == fractures.NumVertices.end())
        {
            fractures.NumVertices.insert({numVertices, {fractureId}});
        }
        else
        {
            fractures.NumVertices[numVertices].push_back(fractureId);
        }

        // Popola Spheres
        fractures.Spheres.push_back(Analytics::calcsphere(vertex_data));
        // Popola Normals
        // Calcola e normalizza il prodotto scalare tra i vettori che congiungono i primi tre vertici
        fractures.Normals.push_back(Analytics::normal(vertex_data));
    }
    return true;
}
list<vector<unsigned int>> checkspheres(const Fractures& fractures)
{
    list<vector<unsigned int>> goodcouples;
    for (unsigned int id1 = 0; id1 < fractures.FracturesNumber - 1; ++id1)
    {
        for (unsigned int id2 = id1 + 1; id2 < fractures.FracturesNumber; ++id2)
        {
            Vector3d point1(fractures.Spheres[id1](0), fractures.Spheres[id1](1), fractures.Spheres[id1](2));
            Vector3d point2(fractures.Spheres[id2](0), fractures.Spheres[id2](1), fractures.Spheres[id2](2));
            double r1 = fractures.Spheres[id1](3);
            double r2 = fractures.Spheres[id2](3);
            if(Analytics::distance(point1, point2) < r1+r2)
            {
                vector<unsigned int> ids = {id1, id2};
                goodcouples.push_back(ids);
            }
        }
    }
    return goodcouples;
}
}

namespace Analytics {
Vector4d calcsphere(const vector<Vector3d>& vertex_data)
{
    Vector4d sphere;
    Vector3d baricentro{0.0,0.0,0.0};

    // Calcola il baricentro come media delle posizioni dei vertici
    for (size_t k = 0; k < vertex_data.size(); ++k)
    {
        baricentro += vertex_data[k];
    }
    baricentro /= vertex_data.size();
    // Copia il baricentro nelle prime tre posizioni di sphere
    for (unsigned int i = 0; i < 3; ++i)
    {
        sphere(i) = baricentro(i);
    }

    // Salva al fondo di sphere il raggio (tutti i vertici hanno la stessa distanza dal baricentro)
    sphere(3) = distance(baricentro, vertex_data[0]);
    return sphere;
}
double distance(const Vector3d& point1, const Vector3d& point2)
{
    return (point1-point2).norm();
}
Vector3d normal(const vector<Vector3d>& vertex_data)
{
    Vector3d n = (vertex_data[1] - vertex_data[0]).cross(vertex_data[2]-vertex_data[0]);
    n.normalize();
    return n;
}
}

// namespace SortLibrary {
// template<typename T>
// void BubbleSort(vector<T>& data)
// {
//     size_t rem_size = data.size();
//     size_t last_seen = rem_size;
//     bool swapped = true;

//     while (swapped)
//     {
//         swapped = false;
//         for (size_t i = 1; i < rem_size; i++)
//         {
//             if (data[i-1] < data[i]) //ordina dal maggiore al minore
//             {
//                 swap(data[i-1], data[i]);
//                 swapped = true;
//                 last_seen = i;
//             }
//         }
//         //        rem_size = rem_size - 1;
//         rem_size = last_seen;
//     }
// }
// }
