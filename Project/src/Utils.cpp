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
using namespace Analytics;

namespace Polygons {
bool importdfn(const string& filename, Fractures& fractures)
{
    ifstream file(filename);
    if (!file.is_open()) {
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

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
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
        for (unsigned int i = 0; i < 3; ++i) {
            // Legge per 3 volte la righe con le coordinate
            getline(file, line);
            istringstream vertex_iss(line);
            for (unsigned int k = 0; k < numVertices - 1; ++k) {
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
        if(fractures.NumVertices.find(numVertices) == fractures.NumVertices.end()) {
            fractures.NumVertices.insert({numVertices, {fractureId}});
        }
        else {
            fractures.NumVertices[numVertices].push_back(fractureId);
        }

        // Popola Spheres
        fractures.Spheres.push_back(Analytics::calcsphere(vertex_data));
        // Calcola e normalizza il prodotto scalare tra i vettori che congiungono i primi tre vertici, popola Normals
        fractures.Normals.push_back(Analytics::normal(vertex_data));
    }
    return true;
}
list<vector<unsigned int>> checkspheres(const Fractures& fractures)
{
    list<vector<unsigned int>> goodcouples;
    for (unsigned int id1 = 0; id1 < fractures.FracturesNumber - 1; ++id1) {
        for (unsigned int id2 = id1 + 1; id2 < fractures.FracturesNumber; ++id2) {
            Vector3d point1(fractures.Spheres[id1](0), fractures.Spheres[id1](1), fractures.Spheres[id1](2));
            Vector3d point2(fractures.Spheres[id2](0), fractures.Spheres[id2](1), fractures.Spheres[id2](2));
            double r1 = fractures.Spheres[id1](3);
            double r2 = fractures.Spheres[id2](3);
            if(Analytics::distance(point1, point2) < r1+r2) {
                vector<unsigned int> ids = {id1, id2};
                goodcouples.push_back(ids);
            }
        }
    }
    return goodcouples;
}
void tracesfinder(const Fractures& fractures, const list<vector<unsigned int>>& goodcouples, Traces& traces, const double& epsilon)
{
    Vector3d inter;
    for(vector<unsigned int> couple : goodcouples) {
        unsigned int id1 = couple[0];
        unsigned int id2 = couple[1];
        Vector3d n1 = fractures.Normals[id1];
        Vector3d n2 = fractures.Normals[id2];
        Vector3d tangent = n1.cross(n2);
        if (tangent.norm() > epsilon) { // Controlla che i piani non siano paralleli
            double d1 = n1.dot(fractures.CoordVertices[id1][0]);
            double d2 = n2.dot(fractures.CoordVertices[id2][0]);

            Matrix<double, 2, 3> A;
            A.row(0) = n1.transpose();
            A.row(1) = n2.transpose();
            Vector2d b(d1,d2);

            // Risolvi il sistema di equazioni lineari per trovare il punto di intersezione
            Vector3d point = A.fullPivLu().solve(b);

            vector<Vector3d> intersections;
            for (size_t i = 0; i < fractures.CoordVertices[id1].size() - 1; ++i) {
                if (Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id1][i], fractures.CoordVertices[id1][i+1] - fractures.CoordVertices[id1][i], inter)) {
                    vector<Vector3d> verifica;
                    Vector3d control;
                    for (size_t j = 0; j<fractures.CoordVertices[id2].size()-1;++j) {
                        if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][j],fractures.CoordVertices[id2][j+1]-fractures.CoordVertices[id2][j], control)) {
                            verifica.push_back(control);
                        }
                    }
                    if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1],fractures.CoordVertices[id2][0]-fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1], control)) {
                        verifica.push_back(control);
                    }
                    if (verifica.size()%2 != 0) {
                        intersections.push_back(inter);
                    }
                    if (find(verifica.begin(),verifica.end(),inter)!=verifica.end() && verifica.size()==2) {
                        intersections.push_back(inter);
                    }
                }
            }

            if (Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id1][fractures.CoordVertices[id1].size() - 1], fractures.CoordVertices[id1][0] - fractures.CoordVertices[id1][fractures.CoordVertices[id1].size() - 1], inter)) {
                vector<Vector3d> verifica;
                Vector3d control;
                for (size_t j = 0; j<fractures.CoordVertices[id2].size()-1;++j) {
                    if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][j],fractures.CoordVertices[id2][j+1]-fractures.CoordVertices[id2][j], control)) {
                        verifica.push_back(control);
                    }
                }
                if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1],fractures.CoordVertices[id2][0]-fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1], control)) {
                    verifica.push_back(control);
                }
                if (verifica.size()%2 != 0) {
                    intersections.push_back(inter);
                }
                if (find(verifica.begin(),verifica.end(),inter)!=verifica.end() && verifica.size()==2) {
                    intersections.push_back(inter);
                }
            }

            for (size_t i = 0; i < fractures.CoordVertices[id2].size() - 1; ++i) {
                if (Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id2][i], fractures.CoordVertices[id2][i+1] - fractures.CoordVertices[id2][i], inter)) {
                    vector<Vector3d> verifica;
                    Vector3d control;
                    for (size_t j = 0; j<fractures.CoordVertices[id1].size()-1;++j) {
                        if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][j],fractures.CoordVertices[id1][j+1]-fractures.CoordVertices[id1][j], control)) {
                            verifica.push_back(control);
                        }
                    }
                    if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1],fractures.CoordVertices[id1][0]-fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1], control)) {
                        verifica.push_back(control);
                    }
                    if (verifica.size()%2 != 0) {
                        intersections.push_back(inter);
                    }
                    if (find(verifica.begin(),verifica.end(),inter)!=verifica.end() && verifica.size()==2) {
                        intersections.push_back(inter);
                    }
                }
            }

            if (Analytics::intersectrettaretta(point, tangent, fractures.CoordVertices[id2][fractures.CoordVertices[id2].size() - 1], fractures.CoordVertices[id2][0] - fractures.CoordVertices[id2][fractures.CoordVertices[id2].size() - 1], inter)) {
                vector<Vector3d> verifica;
                Vector3d control;
                for (size_t j = 0; j<fractures.CoordVertices[id1].size()-1;++j) {
                    if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][j],fractures.CoordVertices[id1][j+1]-fractures.CoordVertices[id1][j], control)) {
                        verifica.push_back(control);
                    }
                }
                if (Analytics::intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][fractures.CoordVertices[id2].size()-1],fractures.CoordVertices[id1][0]-fractures.CoordVertices[id1][fractures.CoordVertices[id2].size()-1], control)) {
                    verifica.push_back(control);
                }
                if (verifica.size()%2 != 0) {
                    intersections.push_back(inter);
                }
                if (find(verifica.begin(),verifica.end(),inter)!=verifica.end() && verifica.size()==2) {
                    intersections.push_back(inter);
                }
            }

            cout << id1 << "&" << id2 << ": ";
            for (Vector3d v : intersections) {
                for (unsigned int i = 0; i < v.size(); ++i) {
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
    vector<double> raggi;
    double raggio;

    // Calcola il baricentro come media delle posizioni dei vertici e per ogni vertice il raggio
    for (size_t k = 0; k < vertex_data.size(); ++k) {
        baricentro += vertex_data[k];
        raggio = distance(baricentro, vertex_data[k]);
        raggi.push_back(raggio);
    }
    baricentro /= vertex_data.size();
    // Copia il baricentro nelle prime tre posizioni di sphere
    for (unsigned int i = 0; i < 3; ++i) {
        sphere(i) = baricentro(i);
    }
    // Trova il raggio massimo
    double max_raggio = 0;
    for (size_t i = 0; i < raggi.size(); ++i) {
        if (raggi[i] > max_raggio) {
            max_raggio = raggi[i];
        }
    }
    // Salva al fondo di sphere il raggio
    sphere(3) = max_raggio;
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
bool intersectrettaretta(const Vector3d& point, const Vector3d& dir,
                         const Vector3d& point1, const Vector3d& dir1, Vector3d& inter)
{
    Vector3d cross_product = dir.cross(dir1);
    Vector3d diff = point1 - point;

    double t = (diff.cross(dir1)).dot(cross_product)/pow(cross_product.norm(),2);
    double s = (diff.cross(dir)).dot(cross_product)/pow(cross_product.norm(),2);

    Vector3d pl = point+t*dir;
    Vector3d pl1 = point1+s*dir1;

    if((pl-pl1).norm() < 1e-10) {
        if (s > -1e-10 && s < 1.0+1e-10) {
            inter = point + t * dir;
            return true;
        }
        else {
            inter = {1000,1000,1000};
            return false;
        }
    }
    else {
        inter = {1000,1000,1000};
        return false;
    }
}

bool intersectrettasemiretta(const Vector3d& point, const Vector3d& dir,
                             const Vector3d& point1, const Vector3d& dir1, Vector3d& control)
{
    Vector3d cross_product = dir.cross(dir1);
    Vector3d diff = point1 - point;

    double t = (diff.cross(dir1)).dot(cross_product)/pow(cross_product.norm(),2) ;
    double s = (diff.cross(dir)).dot(cross_product)/pow(cross_product.norm(),2);

    Vector3d pl = point+t*dir;
    Vector3d pl1 = point1+s*dir1;

    if((pl-pl1).norm() < 1e-10 && t >= 0) {
        if (s > -1e-10 && s < 1.0 + 1e-10) {
            control = point + t * dir;
            return true;
        }
        else {
            control = {1000,1000,1000};
            return false;
        }
    }
    else {
        control = {1000,1000,1000};
        return false;
    }
}
}
