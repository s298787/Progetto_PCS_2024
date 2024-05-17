#include "Utils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>
#include <algorithm>
#include <bits/stdc++.h>

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
void tracesfinder(const Fractures& fractures, const list<vector<unsigned int>>& goodcouples, Traces& traces, const double& epsilon)
{
    Vector3d limit(100,100,100);
    for(vector<unsigned int> couple : goodcouples)
    {
        unsigned int id1 = couple[0];
        unsigned int id2 = couple[1];
        Vector3d n1 = fractures.Normals[id1];
        Vector3d n2 = fractures.Normals[id2];
        Vector3d tangent = n1.cross(n2); // Vettore direzione della retta di intersezione
        if (tangent.norm() > epsilon) // Controlla che i piani non siano paralleli
        {
            // Calcola i termini noti del piano
            double d1 = n1.dot(fractures.CoordVertices[id1][0]);
            double d2 = n2.dot(fractures.CoordVertices[id2][0]);

            Matrix<double, 2, 3> A; // Matrice delle normali
            A.row(0) = n1.transpose();
            A.row(1) = n2.transpose();

            Vector2d b(d1,d2); // Vettore dei termini noti

            // Risolve il sistema di equazioni lineari per trovare un punto sulla retta di intersezione
            Vector3d point = A.fullPivLu().solve(b);

            vector<Vector3d> intersections;
            Vector3d inter;
            for (size_t i = 0; i < fractures.CoordVertices[id1].size() - 1; ++i)
            {
                // Cerca l'intersezione tra la retta tangent ed il segmento
                inter = Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id1][i], fractures.CoordVertices[id1][i+1] - fractures.CoordVertices[id1][i]);
                if (inter != limit)
                {
                    // Verifica se inter è interno al poligono id2
                    vector<int> checksigns;
                    for (size_t j = 0; j < fractures.CoordVertices[id2].size()-1;++j)
                    {

                        Vector3d u = fractures.CoordVertices[id2][j+1]-fractures.CoordVertices[id2][j];
                        Vector3d v = -(fractures.CoordVertices[id2][j]-inter);
                        double uv = u.dot(v)/(u.norm()*v.norm());
                        int alpha = Analytics::sign(sin(acos(uv)));
                        // cout << alpha << " ";
                        checksigns.push_back(alpha);
                    }

                    Vector3d l = -(fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1]-inter);
                    Vector3d m = fractures.CoordVertices[id2][0]-fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1];
                    double lm = l.dot(m)/(l.norm()*m.norm());
                    int alpha1 = Analytics::sign(sin(acos(lm)));
                    // cout << alpha1 << " " << endl;
                    checksigns.push_back(alpha1);
                    if(Analytics::positive(checksigns))
                    {
                        intersections.push_back(inter);
                    }
                }
            }
            inter = Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id1][fractures.CoordVertices[id1].size() - 1], fractures.CoordVertices[id1][0] - fractures.CoordVertices[id1][fractures.CoordVertices[id1].size() - 1]);
            if (inter != limit)
            {
                vector<int> checksigns;
                for (size_t j = 0; j < fractures.CoordVertices[id2].size()-1;++j)
                {

                    Vector3d u = fractures.CoordVertices[id2][j+1]-fractures.CoordVertices[id2][j];
                    Vector3d v = -(fractures.CoordVertices[id2][j]-inter);
                    double uv ;
                    uv = u.dot(v)/(u.norm()*v.norm());
                    int alpha = Analytics::sign(sin(acos(uv)));
                    // cout << alpha << " ";
                    checksigns.push_back(alpha);
                }
                Vector3d l = -(fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1]-inter);
                Vector3d m = fractures.CoordVertices[id2][0]-fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1];
                double lm = l.dot(m)/(l.norm()*m.norm());
                int alpha1 = Analytics::sign(sin(acos(lm)));
                // cout << alpha1 << " ";
                checksigns.push_back(alpha1);
                if(Analytics::positive(checksigns))
                {
                    intersections.push_back(inter);
                }
            }
            for (size_t i = 0; i < fractures.CoordVertices[id2].size() - 1; ++i)
            {
                inter = Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id2][i], fractures.CoordVertices[id2][i+1] - fractures.CoordVertices[id2][i]);
                if (inter != limit)
                {
                    vector<int> checksigns;
                    for (size_t j = 0; j < fractures.CoordVertices[id1].size()-1;++j)
                    {

                        Vector3d u = fractures.CoordVertices[id1][j+1]-fractures.CoordVertices[id1][j];
                        Vector3d v = -(fractures.CoordVertices[id1][j]-inter);
                        double uv ;
                        uv = u.dot(v)/(u.norm()*v.norm());
                        int alpha = Analytics::sign(sin(acos(uv)));
                        // cout <<alpha << " ";
                        checksigns.push_back(alpha);
                    }
                    Vector3d l = -(fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1]-inter);
                    Vector3d m = fractures.CoordVertices[id1][0]-fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1];
                    double lm = l.dot(m)/(l.norm()*m.norm());
                    int alpha1 = Analytics::sign(sin(acos(lm)));
                    // cout << alpha1 << " ";
                    checksigns.push_back(alpha1);

                    if(Analytics::positive(checksigns))
                    {
                        intersections.push_back(inter);
                    }
                }
            }
            inter = Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id2][fractures.CoordVertices[id2].size() - 1], fractures.CoordVertices[id2][0] - fractures.CoordVertices[id2][fractures.CoordVertices[id2].size() - 1]);
            if (inter != limit)
            {
                vector<int> checksigns;
                for (size_t j = 0; j < fractures.CoordVertices[id1].size()-1;++j)
                {

                    Vector3d u = fractures.CoordVertices[id1][j+1]-fractures.CoordVertices[id1][j];
                    Vector3d v = -(fractures.CoordVertices[id1][j]-inter);
                    double uv ;
                    uv = u.dot(v)/(u.norm()*v.norm());
                    int alpha = Analytics::sign(sin(acos(uv)));
                    // cout << alpha << " ";
                    checksigns.push_back(alpha);
                }
                Vector3d l = -(fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1]-inter);
                Vector3d m = fractures.CoordVertices[id1][0]-fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1];
                double lm = l.dot(m)/(l.norm()*m.norm());
                int alpha1 = Analytics::sign(sin(acos(lm)));
                // cout << alpha1 << " ";
                checksigns.push_back(alpha1);
                if(Analytics::positive(checksigns))
                {
                    intersections.push_back(inter);
                }
            }
            cout << intersections.size() << endl;
            cout << id1 << "&" << id2 << ": ";
            for (Vector3d v : intersections)
            {
                for (unsigned int i = 0; i < v.size(); ++i)
                {
                    cout << v(i) << ";";
                }
                cout << "   ";
            }
            cout << endl;
        }
    }
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

Vector3d intersectrettaretta(const Vector3d& point, const Vector3d& dir,
                             const Vector3d& point1, const Vector3d& dir1)
{
    Vector3d intersec; //= {1e10,1e10,1e10};

    Vector3d cross_product = dir.cross(dir1);
    double cross_product_norm = cross_product.norm();
    Vector3d diff = point1 - point;
    // cout << diff.transpose()<<endl;

    double t = (diff.cross(dir1)).dot(cross_product)/ pow(cross_product_norm,2);
    double s = (diff.cross(dir)).dot(cross_product) / pow(cross_product_norm,2);

    // double s = (-diff(1)/dir1(1) + diff(0)*dir(1)/(dir1(1)*dir(0)))/(1-dir(1)*dir1(0)/(dir1(1)*dir(0)));
    // double t = dir1(0)*s/dir(0) + diff(0)/dir(0);
    // cout << dir(0) << " ";
    Vector3d pointOnL1 = point + t * dir;
    Vector3d pointOnL2 = point1 + s * dir1;
    // cout << (pointOnL1 - pointOnL2).norm() << endl;
    if ((pointOnL1 - pointOnL2).norm() < 1e-10)
    {
        if (s >= 0.0 && s <= 1.0)
        {
            intersec = pointOnL1;
        }
    }
    else
    {
        intersec = {1e10,1e10,1e10};
    }
    return intersec;

}
int sign(double x)
{
    if (x > 0)
    {
        return 1;  // Positive
    }
    else if (x < 0)
    {
        return -1; // Negative
    }
    else
    {
        return 0;  // Zero
    }
}
bool positive(vector<int> checksigns)
{
    for (int s : checksigns)
    {
        if (s == -1)
        {
            return false;
        }
    }
    return true;
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
