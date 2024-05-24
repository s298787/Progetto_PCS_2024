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
        cerr << "Errore while opening " << filename << endl;
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
    file.close();
    return true;
}
list<vector<unsigned int>> checkspheres(const Fractures& fractures)
{
    list<vector<unsigned int>> goodcouples;
    double r1;
    double r2;
    vector<unsigned int> ids;
    for (unsigned int id1 = 0; id1 < fractures.FracturesNumber - 1; ++id1) {
        for (unsigned int id2 = id1 + 1; id2 < fractures.FracturesNumber; ++id2) {
            Vector3d point1(fractures.Spheres[id1](0), fractures.Spheres[id1](1), fractures.Spheres[id1](2));
            Vector3d point2(fractures.Spheres[id2](0), fractures.Spheres[id2](1), fractures.Spheres[id2](2));
            r1 = fractures.Spheres[id1](3);
            r2 = fractures.Spheres[id2](3);
            if(distance(point1, point2) < r1+r2) {
                ids = {id1, id2};
                goodcouples.push_back(ids);
            }
        }
    }
    return goodcouples;
}
void tracesfinder(const Fractures& fractures, const list<vector<unsigned int>>& goodcouples, Traces& traces, const double& epsilon)
{
    Vector3d inter;
    vector<Vector3d> intersections;
    unsigned int id1;
    unsigned int id2;
    Vector3d n1;
    Vector3d n2;
    Vector3d tangent;
    double d1;
    double d2;
    Matrix<double, 2, 3> A;
    Vector3d point;
    for(vector<unsigned int> couple : goodcouples) {
        id1 = couple[0];
        id2 = couple[1];
        n1 = fractures.Normals[id1];
        n2 = fractures.Normals[id2];
        tangent = n1.cross(n2);
        if (tangent.norm() > epsilon) { // Controlla che i piani non siano paralleli
            d1 = n1.dot(fractures.CoordVertices[id1][0]);
            d2 = n2.dot(fractures.CoordVertices[id2][0]);

            // Matrix<double, 2, 3> A;
            A.row(0) = n1.transpose();
            A.row(1) = n2.transpose();
            Vector2d b(d1,d2);

            // Risolvi il sistema di equazioni lineari per trovare il punto di intersezione
            point = A.fullPivLu().solve(b);

            // vector<Vector3d> intersections;
            for (size_t i = 0; i < fractures.CoordVertices[id1].size() - 1; ++i) {
                if (intersectrettaretta(point, tangent, fractures.CoordVertices[id1][i], fractures.CoordVertices[id1][i+1] - fractures.CoordVertices[id1][i], inter)) {
                    vector<Vector3d> verifica;
                    Vector3d control;
                    for (size_t j = 0; j<fractures.CoordVertices[id2].size()-1;++j) {
                        if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][j],fractures.CoordVertices[id2][j+1]-fractures.CoordVertices[id2][j], control)) {
                            verifica.push_back(control);
                        }
                    }
                    if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1],fractures.CoordVertices[id2][0]-fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1], control)) {
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

            if (intersectrettaretta(point, tangent, fractures.CoordVertices[id1][fractures.CoordVertices[id1].size() - 1], fractures.CoordVertices[id1][0] - fractures.CoordVertices[id1][fractures.CoordVertices[id1].size() - 1], inter)) {
                vector<Vector3d> verifica;
                Vector3d control;
                for (size_t j = 0; j<fractures.CoordVertices[id2].size()-1;++j) {
                    if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][j],fractures.CoordVertices[id2][j+1]-fractures.CoordVertices[id2][j], control)) {
                        verifica.push_back(control);
                    }
                }
                if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1],fractures.CoordVertices[id2][0]-fractures.CoordVertices[id2][fractures.CoordVertices[id2].size()-1], control)) {
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
                if (intersectrettaretta(point, tangent, fractures.CoordVertices[id2][i], fractures.CoordVertices[id2][i+1] - fractures.CoordVertices[id2][i], inter)) {
                    vector<Vector3d> verifica;
                    Vector3d control;
                    for (size_t j = 0; j<fractures.CoordVertices[id1].size()-1;++j) {
                        if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][j],fractures.CoordVertices[id1][j+1]-fractures.CoordVertices[id1][j], control)) {
                            verifica.push_back(control);
                        }
                    }
                    if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1],fractures.CoordVertices[id1][0]-fractures.CoordVertices[id1][fractures.CoordVertices[id1].size()-1], control)) {
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

            if (intersectrettaretta(point, tangent, fractures.CoordVertices[id2][fractures.CoordVertices[id2].size() - 1], fractures.CoordVertices[id2][0] - fractures.CoordVertices[id2][fractures.CoordVertices[id2].size() - 1], inter)) {
                vector<Vector3d> verifica;
                Vector3d control;
                for (size_t j = 0; j<fractures.CoordVertices[id1].size()-1;++j) {
                    if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][j],fractures.CoordVertices[id1][j+1]-fractures.CoordVertices[id1][j], control)) {
                        verifica.push_back(control);
                    }
                }
                if (intersectrettasemiretta(inter, tangent,fractures.CoordVertices[id1][fractures.CoordVertices[id2].size()-1],fractures.CoordVertices[id1][0]-fractures.CoordVertices[id1][fractures.CoordVertices[id2].size()-1], control)) {
                    verifica.push_back(control);
                }
                if (verifica.size()%2 != 0) {
                    intersections.push_back(inter);
                }
                if (find(verifica.begin(),verifica.end(),inter)!=verifica.end() && verifica.size()==2) {
                    intersections.push_back(inter);
                }
            }
            // Elimina gli elementi ripetuti da intersections
            for (size_t i = 0; i < intersections.size(); ++i) {
                for (size_t j = i + 1; j < intersections.size(); ++j) {
                    if ((intersections[i] - intersections[j]).norm() < 1e-10) {
                        intersections.erase(intersections.begin()+j);
                    }
                }
            }
            if (intersections.size() == 2) {
                // Popola TracesExtremesCoord
                traces.TracesExtremesCoord.push_back(intersections);
                // Popola TracesFracturesId
                traces.TracesFracturesId.push_back({id1, id2});
                // Popola TracesId
                unsigned int traceid = traces.TracesExtremesCoord.size() - 1;
                traces.TracesId.push_back(traceid);
                // Popola TracesLengths
                traces.TracesLengths.push_back(distance(intersections[0], intersections[1]));
                // Controlla se la traccia è passante per il poligono id1, per farlo vede se entrambi gli estremi giacciono su un lato di id1
                // Per verificare se un estremo giace su un lato controlla se la somma delle distanze tra l'estremo e i vertici è uguale alla lunghezza del lato
                unsigned int count = 0;
                // Controlla tutti i lati di id1 per il primo estremo della traccia
                for (size_t i = 0; i < fractures.CoordVertices[id1].size() - 1; ++i) {
                    if (distance(intersections[0], fractures.CoordVertices[id1][i]) + distance(intersections[0], fractures.CoordVertices[id1][i+1])
                        - distance(fractures.CoordVertices[id1][i+1], fractures.CoordVertices[id1][i]) <= 1e-6) {
                        count++;
                    }
                }
                if (distance(intersections[0], fractures.CoordVertices[id1][0]) + distance(intersections[0], fractures.CoordVertices[id1][fractures.CoordVertices.size() - 1])
                    - distance(fractures.CoordVertices[id1][0], fractures.CoordVertices[id1][fractures.CoordVertices.size()-1]) <= 1e-6) {
                    count++;
                }
                // Controlla tutti i lati di id1 per il primo estremo della traccia
                for (size_t i = 0; i < fractures.CoordVertices[id1].size() - 1; ++i) {
                    if (distance(intersections[1], fractures.CoordVertices[id1][i]) + distance(intersections[1], fractures.CoordVertices[id1][i+1])
                        - distance(fractures.CoordVertices[id1][i+1], fractures.CoordVertices[id1][i]) <= 1e-6) {
                        count++;
                    }
                }
                if (distance(intersections[1], fractures.CoordVertices[id1][0]) + distance(intersections[1], fractures.CoordVertices[id1][fractures.CoordVertices.size() - 1])
                    - distance(fractures.CoordVertices[id1][0], fractures.CoordVertices[id1][fractures.CoordVertices.size()-1]) <= 1e-6) {
                    count++;
                }
                // Popola TipsTrue o TipsFalse
                if (count == 2) {
                    traces.TipsTrue.push_back({traceid, id1});
                }
                else {
                    traces.TipsFalse.push_back({traceid, id1});
                }
                // Ripete per il poligono id2
                count = 0;
                for (size_t i = 0; i < fractures.CoordVertices[id2].size() - 1; ++i) {
                    if (distance(intersections[0], fractures.CoordVertices[id2][i]) + distance(intersections[0], fractures.CoordVertices[id2][i+1])
                        - distance(fractures.CoordVertices[id2][i+1], fractures.CoordVertices[id2][i]) <= 1e-6) {
                        count++;
                    }
                }
                if (distance(intersections[0], fractures.CoordVertices[id2][0]) + distance(intersections[0], fractures.CoordVertices[id2][fractures.CoordVertices.size() - 1])
                    - distance(fractures.CoordVertices[id2][0], fractures.CoordVertices[id2][fractures.CoordVertices.size()-1]) <= 1e-6) {
                    count++;
                }
                for (size_t i = 0; i < fractures.CoordVertices[id2].size() - 1; ++i) {
                    if (distance(intersections[1], fractures.CoordVertices[id2][i]) + distance(intersections[1], fractures.CoordVertices[id2][i+1])
                        - distance(fractures.CoordVertices[id2][i+1], fractures.CoordVertices[id2][i]) <= 1e-6) {
                        count++;
                    }
                }
                if (distance(intersections[1], fractures.CoordVertices[id2][0]) + distance(intersections[1], fractures.CoordVertices[id2][fractures.CoordVertices.size() - 1])
                    - distance(fractures.CoordVertices[id2][0], fractures.CoordVertices[id2][fractures.CoordVertices.size()-1]) <= 1e-6) {
                    count++;
                }
                if (count == 2) {
                    traces.TipsTrue.push_back({traceid, id2});
                }
                else {
                    traces.TipsFalse.push_back({traceid, id2});
                }
            }
        }
        intersections.clear(); // Svuota intersections
    }
    // Popola TracesNumber
    traces.TracesNumber = traces.TracesId.size();
}
}

namespace Analytics {
Vector4d calcsphere(const vector<Vector3d>& vertex_data)
{
    Vector4d sphere;
    Vector3d baricentro{0.0,0.0,0.0};
    vector<double> raggi;
    double raggio;

    // Calcola il baricentro come media delle posizioni dei vertici
    for (size_t k = 0; k < vertex_data.size(); ++k) {
        baricentro += vertex_data[k];
    }
    baricentro /= vertex_data.size();
    // Calcola la distanza tra il baricentro e ciascuno dei raggi
    for (size_t k = 0; k < vertex_data.size(); ++k) {
        raggio = distance(baricentro, vertex_data[k]);
        raggi.push_back(raggio);
    }
    // Copia il baricentro nelle prime tre posizioni di sphere
    for (unsigned int i = 0; i < 3; ++i) {
        sphere(i) = baricentro(i);
    }
    // Trova il raggio massimo
    raggio = 0;
    for (size_t i = 0; i < raggi.size(); ++i) {
        if (raggi[i] > raggio) {
            raggio = raggi[i];
        }
    }
    // Salva il raggio al fondo di sphere
    sphere(3) = raggio;
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

namespace OutputFileTools {
using namespace Polygons;
bool printtraces(const string& tracesfileout, const Traces& traces)
{
    ofstream fileout(tracesfileout);
    if (fileout.fail()) {
        cerr << "Error while creating/opening " << tracesfileout << endl;
        return false;
    }

    fileout << "# Number of Traces" << endl; // Stampa la riga di intestazione
    fileout << traces.TracesNumber << endl; // Stampa il numero di tracce
    for (size_t i = 0; i < traces.TracesNumber; ++i) {
        // Per ogni traccia stampa gli id e le coordinate degli estremi
        fileout << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
        fileout << traces.TracesId[i] << "; " << traces.TracesFracturesId[i][0] << "; " << traces.TracesFracturesId[i][1] << "; ";
        for (unsigned int j = 0; j < 2; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                fileout << traces.TracesExtremesCoord[i][j](k) << "; ";
            }
        }
        fileout << endl;
    }
    fileout.close();
    return true;
}
bool printtips(const string& tipsfileout, const Traces& traces, const Fractures& fractures)
{
    ofstream fileout(tipsfileout);
    if (fileout.fail()) {
        cerr << "Error while creating/opening " << tipsfileout << endl;
        return false;
    }
    unsigned int id_t;
    vector<unsigned int> confronto;
    vector<unsigned int> idtraces;
    // Scorre gli id delle fratture
    for (unsigned int id = 0; id < fractures.FracturesNumber; ++id) {        
        // Verifica quali tracce passano per la frattura id e le memorizza in idtraces
        for (size_t j = 0; j < traces.TracesNumber; ++j) {
            if (traces.TracesFracturesId[j][0] == id || traces.TracesFracturesId[j][1] == id) {
                idtraces.push_back(j);
            }
        }
        // Stampa su file le informazioni sulla frattura
        fileout << "# FractureId; NumTraces" << endl;
        fileout << id << "; " << idtraces.size() << endl;

        // Ordina idtraces in base alla lunghezza delle tracce
        size_t rem_size = idtraces.size();
        size_t last_seen = rem_size;
        bool swapped = true;
        while (swapped) {
            swapped = false;
            for (size_t i = 1; i < rem_size; i++) {
                if (traces.TracesLengths[idtraces[i-1]] < traces.TracesLengths[idtraces[i]]) {
                    swap(idtraces[i-1], idtraces[i]);
                    swapped = true;
                    last_seen = i;
                }
            }
            rem_size = last_seen;
        }

        // Stampa su file le informazioni sulle tracce passanti appartenenti a id
        for (size_t k = 0; k < idtraces.size(); ++k) {
            id_t = idtraces[k];
            confronto = {id_t, id};
            // Controlla se la coppia idtraccia - idfrattura è presente in TipsTrue e stampa
            if (find(traces.TipsTrue.begin(), traces.TipsTrue.end(), confronto) != traces.TipsTrue.end()) {
                fileout << "# TraceId; Tips; Length" << endl;
                fileout << idtraces[k] << "; " << "True" << "; " << traces.TracesLengths[id_t] << endl;
            }
        }
        // Stampa su file le informazioni sulle tracce non passanti appartenenti a id
        for (size_t k = 0; k < idtraces.size(); ++k) {
            id_t = idtraces[k];
            confronto = {id_t, id};
            // Controlla se la coppia idtraccia - idfrattura è presente in TipsFalse e stampa
            if (find(traces.TipsFalse.begin(), traces.TipsFalse.end(), confronto) != traces.TipsFalse.end()) {
                fileout << "# TraceId; Tips; Length" << endl;
                fileout << idtraces[k] << "; " << "False" << "; " << traces.TracesLengths[id_t] << endl;
            }
        }
        fileout << endl;
        idtraces.clear(); // Svuota idtraces
    }
    fileout.close();
    return true;
}
}
