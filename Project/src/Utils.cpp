#include "Utils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>
#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace Analytics;

// Parte 1
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
            if((point1-point2).squaredNorm() < (r1+r2)*(r1+r2)) {
                ids = {id1, id2};
                goodcouples.push_back(ids);
            }
        }
    }
    return goodcouples;
}
void tracesfinder(const Fractures& fractures, const list<vector<unsigned int>>& goodcouples, Traces& traces, const double& epsilon)
{
    Vector3d inter; // Vettore in cui viene salvata la singola intersezione
    vector<Vector3d> intersections; // Array delle intersezioni
    unsigned int id1; // Id della prima fattura
    unsigned int id2; // Id della seconda frattura
    Vector3d n1; // Normale a id1
    Vector3d n2;// Normale a id2
    Vector3d tangent; // Direttrice retta di intersezione
    double d1; // Termine noto relativo al piano id1
    double d2; // Termine noto relativo al piano id2
    Matrix<double, 2, 3> A; // Matrice del sistema lineare di piani
    Vector3d point;
    for(const vector<unsigned int>& couple : goodcouples) {
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
            point = A.colPivHouseholderQr().solve(b);

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
            // Elimina gli elementi ripetuti in intersections
            for (size_t i = 0; i < intersections.size(); ++i) {
                for (size_t j = i + 1; j < intersections.size(); ++j) {
                    if ((intersections[i] - intersections[j]).norm() < epsilon) {
                        intersections.erase(intersections.begin()+j);
                    }
                }
            }
            if (intersections.size() == 2 && distance(intersections[0], intersections[1]) > epsilon) {
                // Popola TracesExtremesCoord
                traces.TracesExtremesCoord.push_back({intersections[0], intersections[1]});
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
                bool tips = false;
                // Controlla tutti i lati di id1 per il primo estremo della traccia
                for (size_t i = 0; i < fractures.CoordVertices[id1].size() - 1; ++i) {
                    if (abs(distance(intersections[0], fractures.CoordVertices[id1][i]) + distance(intersections[0], fractures.CoordVertices[id1][i+1])
                            - distance(fractures.CoordVertices[id1][i+1], fractures.CoordVertices[id1][i])) < epsilon) {
                        count++;
                    }
                }                
                if (abs(distance(intersections[0], fractures.CoordVertices[id1][0]) + distance(intersections[0], fractures.CoordVertices[id1][fractures.CoordVertices.size() - 1])
                        - distance(fractures.CoordVertices[id1][0], fractures.CoordVertices[id1][fractures.CoordVertices.size()-1])) < epsilon) {
                    count++;
                }
                if (count >= 1) {
                    tips = true;
                }
                // Controlla tutti i lati di id1 per il secondo estremo della traccia
                count = 0;
                for (size_t i = 0; i < fractures.CoordVertices[id1].size() - 1; ++i) {
                    if (abs(distance(intersections[1], fractures.CoordVertices[id1][i]) + distance(intersections[1], fractures.CoordVertices[id1][i+1])
                            - distance(fractures.CoordVertices[id1][i+1], fractures.CoordVertices[id1][i])) < epsilon) {
                        count++;
                    }
                }
                if (abs(distance(intersections[1], fractures.CoordVertices[id1][0]) + distance(intersections[1], fractures.CoordVertices[id1][fractures.CoordVertices.size() - 1])
                        - distance(fractures.CoordVertices[id1][0], fractures.CoordVertices[id1][fractures.CoordVertices.size()-1])) < epsilon) {
                    count++;
                }
                if (count == 0) {
                    tips = false;
                }
                // Popola TipsTrue o TipsFalse
                if (tips) {
                    traces.TipsTrue.push_back({traceid, id1});
                }
                else {
                    traces.TipsFalse.push_back({traceid, id1});
                }
                // Ripete per il poligono id2
                count = 0;
                tips = false;
                for (size_t i = 0; i < fractures.CoordVertices[id2].size() - 1; ++i) {
                    if (distance(intersections[0], fractures.CoordVertices[id2][i]) + distance(intersections[0], fractures.CoordVertices[id2][i+1])
                        - distance(fractures.CoordVertices[id2][i+1], fractures.CoordVertices[id2][i]) < epsilon) {
                        count++;
                    }
                }
                if (distance(intersections[0], fractures.CoordVertices[id2][0]) + distance(intersections[0], fractures.CoordVertices[id2][fractures.CoordVertices.size() - 1])
                    - distance(fractures.CoordVertices[id2][0], fractures.CoordVertices[id2][fractures.CoordVertices.size()-1]) < epsilon) {
                    count++;
                }
                if (count >= 1) {
                    tips = true;
                }
                for (size_t i = 0; i < fractures.CoordVertices[id2].size() - 1; ++i) {
                    if (distance(intersections[1], fractures.CoordVertices[id2][i]) + distance(intersections[1], fractures.CoordVertices[id2][i+1])
                        - distance(fractures.CoordVertices[id2][i+1], fractures.CoordVertices[id2][i]) < epsilon) {
                        count++;
                    }
                }
                if (distance(intersections[1], fractures.CoordVertices[id2][0]) + distance(intersections[1], fractures.CoordVertices[id2][fractures.CoordVertices.size() - 1])
                    - distance(fractures.CoordVertices[id2][0], fractures.CoordVertices[id2][fractures.CoordVertices.size()-1]) < epsilon) {
                    count++;
                }
                if (count == 0) {
                    tips = false;
                }
                if (tips) {
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
    double crossnorm = cross_product.norm();
    crossnorm *= crossnorm;

    double t = (diff.cross(dir1)).dot(cross_product)/crossnorm;
    double s = (diff.cross(dir)).dot(cross_product)/crossnorm;

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
    double crossnorm = cross_product.norm();
    crossnorm *= crossnorm;

    double t = (diff.cross(dir1)).dot(cross_product)/crossnorm;
    double s = (diff.cross(dir)).dot(cross_product)/crossnorm;

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
double calcolangolo(const Vector3d& v1, const Vector3d& v2, const Vector3d& normal) {
    Vector3d cross_product = v1.cross(v2);
    double dot_product = v1.dot(v2);
    double angle = atan2(cross_product.norm(), dot_product);
    if (normal.dot(cross_product) < 0) {
        angle = 2 * M_PI - angle;
    }
    return angle;
}
void antiorario(vector<Vector3d>& sottopol, const Vector3d& normal)
{
    Vector3d baricentro = {0,0,0};
    for (const Vector3d& stt : sottopol) {
        baricentro += stt;
    }
    baricentro /= sottopol.size();
    sort(sottopol.begin(), sottopol.end(), [&baricentro, &normal](const Vector3d&  a,const Vector3d&  b)
         {
        return (calcolangolo(Vector3d(1,0,0),a-baricentro,normal) <
                calcolangolo(Vector3d(1,0,0),b-baricentro,normal));
         });
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
bool printtips(const string& tipsfileout, Traces& traces, const Fractures& fractures)
{
    // Ordina TracesId in base alla lunghezza delle tracce
    vector<double>& len = traces.TracesLengths;
    // sort(traces.TracesId.begin(), traces.TracesId.end(), [len](unsigned int  a,  unsigned int  b)
    //      {
    //     return len[a] < len[b];
    //      });

    sort(traces.TipsTrue.begin(), traces.TipsTrue.end(), [&len](array<unsigned int, 2>& a, array<unsigned int, 2>& b)
         {
        return len[a[0]] > len[b[0]];
         });

    sort(traces.TipsFalse.begin(), traces.TipsFalse.end(), [&len](array<unsigned int, 2>& a, array<unsigned int, 2>& b)
         {
        return len[a[0]] > len[b[0]];
         });

    ofstream fileout(tipsfileout);
    if (fileout.fail()) {
        cerr << "Error while creating/opening " << tipsfileout << endl;
        return false;
    }
    unsigned int id_t;
    array<unsigned int, 2> confronto;
    vector<unsigned int> idtraces;
    // Scorre gli id delle fratture
    for (unsigned int id = 0; id < fractures.FracturesNumber; ++id) {        
        // Verifica quali tracce passano per la frattura id e le memorizza in idtraces
        for (size_t j = 0; j < traces.TipsTrue.size(); ++j) {
            id_t = traces.TipsTrue[j][0];
            if (traces.TipsTrue[j][1] == id) {
            // if (traces.TracesFracturesId[id_t][0] == id || traces.TracesFracturesId[id_t][1] == id) {
                idtraces.push_back(id_t);
            }
        }
        for (size_t j = 0; j < traces.TipsFalse.size(); ++j) {
            id_t = traces.TipsFalse[j][0];
            if (traces.TipsFalse[j][1] == id) {
            // if (traces.TracesFracturesId[id_t][0] == id || traces.TracesFracturesId[id_t][1] == id) {
                idtraces.push_back(id_t);
            }
        }
        // Stampa su file le informazioni sulla frattura
        fileout << "# FractureId; NumTraces" << endl;
        fileout << id << "; " << idtraces.size() << endl;

        // Stampa su file le informazioni sulle tracce appartenenti a id
        for (size_t k = 0; k < idtraces.size(); ++k) {
            id_t = idtraces[k];
            confronto = {id_t, id};
            // Controlla se la coppia idtraccia - idfrattura è presente in TipsTrue e stampa
            if (find(traces.TipsTrue.begin(), traces.TipsTrue.end(), confronto) != traces.TipsTrue.end()) {
                fileout << "# TraceId; Tips; Length" << endl;
                fileout << id_t << "; " << "True" << "; " << traces.TracesLengths[id_t] << endl;
            }
            // Controlla se la coppia idtraccia - idfrattura è presente in TipsFalse e stampa
            if (find(traces.TipsFalse.begin(), traces.TipsFalse.end(), confronto) != traces.TipsFalse.end()) {
                fileout << "# TraceId; Tips; Length" << endl;
                fileout << id_t << "; " << "False" << "; " << traces.TracesLengths[id_t] << endl;
            }
        }
        // // Stampa su file le informazioni sulle tracce non passanti appartenenti a id
        // for (size_t k = 0; k < idtraces.size(); ++k) {
        //     id_t = idtraces[k];
        //     confronto = {id_t, id};
        //     // Controlla se la coppia idtraccia - idfrattura è presente in TipsFalse e stampa
            // if (find(traces.TipsFalse.begin(), traces.TipsFalse.end(), confronto) != traces.TipsFalse.end()) {
            //     fileout << "# TraceId; Tips; Length" << endl;
            //     fileout << id_t << "; " << "False" << "; " << traces.TracesLengths[id_t] << endl;
            // }
        // }
        fileout << endl;
        idtraces.clear(); // Svuota idtraces
    }
    fileout.close();
    return true;
}
}

// Parte 2
namespace MeshLibrary {
using namespace Polygons;
using namespace Analytics;
struct vector3d_hash {
    size_t operator()(const Vector3d& v) const {
        size_t seed = 0;
        for (int i = 0; i < v.size(); ++i) {
            seed ^= hash<double>()(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
struct pair_hash {
    size_t operator()(const pair<unsigned int, unsigned int>& p) const {
        auto hash1 = hash<unsigned int>{}(p.first);
        auto hash2 = hash<unsigned int>{}(p.second);
        return hash1 ^ (hash2 << 1);
    }
};
bool meshcalc(const double& epsilon, const Traces& traces, const Fractures& fractures, vector<PolygonalMesh>& mesh, const string& meshfileout)
{
    ofstream fileout(meshfileout);
    if (fileout.fail()) {
        cerr << "Error while creating/opening " << meshfileout << endl;
        return false;
    }
    for (unsigned int id = 0; id < fractures.FracturesNumber; ++id) {
        // Crea un array di tracce della singola frattura
        vector<unsigned int> idtraces;
        // Verifica quali tracce passano per la frattura id e le memorizza in idtraces
        for (size_t j = 0; j < traces.TipsTrue.size(); ++j) {
            unsigned int id_t = traces.TipsTrue[j][0];
            if (traces.TipsTrue[j][1] == id) {
                idtraces.push_back(id_t);
            }
        }
        for (size_t j = 0; j < traces.TipsFalse.size(); ++j) {
            unsigned int id_t = traces.TipsFalse[j][0];
            if (traces.TipsFalse[j][1] == id) {
                idtraces.push_back(id_t);
            }
        }
        // // Popola con le tracce passanti
        // vector<unsigned int> confronto;
        // for (unsigned int j = 0; j < traces.TracesNumber; ++j) {
        //     unsigned int id_tr = traces.TracesId[j];
        //     confronto = {id_tr, id};
        //     if ((traces.TracesFracturesId[id_tr][0] == id || traces.TracesFracturesId[id_tr][1] == id) &&
        //         find(traces.TipsTrue.begin(), traces.TipsTrue.end(), confronto) != traces.TipsTrue.end()) {
        //         idtraces.push_back(id_tr);
        //     }
        // }
        // // Popola con le tracce tracce non passanti
        // for (unsigned int j = 0; j < traces.TracesNumber; ++j) {
        //     unsigned int id_tr = traces.TracesId[j];
        //     confronto = {id_tr, id};
        //     if ((traces.TracesFracturesId[id_tr][0] == id || traces.TracesFracturesId[id_tr][1] == id) &&
        //         find(traces.TipsFalse.begin(), traces.TipsFalse.end(), confronto) != traces.TipsFalse.end()) {
        //         idtraces.push_back(id_tr);
        //     }
        // }

        vector<vector<Vector3d>> sottopoligoni; // Array in cui memorizzare i sottopoligoni
        vector<Vector3d> fracture = fractures.CoordVertices[id];
        sottopoligoni.reserve(50*idtraces.size());
        sottopoligoni.push_back(fracture); // Inizializza sottopoligoni uguale alla frattura
        vector<vector<Vector3d>> copia; // Array copia di sottopoligoni
        copia.reserve(2*sottopoligoni.size());
        // vector<Vector3d> verticesontracelist;
        // Cicla su tutte le tracce della frattura (idtraces)
        for (const unsigned int& id_t : idtraces) {
            // vector<vector<Vector3d>> copia; // Array copia di sottopoligoni
            // copia.reserve(2*sottopoligoni.size());
            array<Vector3d, 2> traceverts = traces.TracesExtremesCoord[id_t];
            // Cicla su tutti gli attuali sottopoligoni
            for (const vector<Vector3d>& currentpolygon : sottopoligoni) {
                Vector3d inter;

                // Salta al prossimo se currentpolygon e la traccia sono troppo lontani
                Vector4d currentpolygonsphere = calcsphere(currentpolygon);
                Vector4d tracevertsphere = calcsphere({traceverts[0], traceverts[1]});
                if (distance({currentpolygonsphere(0), currentpolygonsphere(1), currentpolygonsphere(2)},
                             {tracevertsphere(0), tracevertsphere(1), tracevertsphere(2)}) >
                    currentpolygonsphere(3) + tracevertsphere(3)) {
                    copia.push_back(currentpolygon);
                    continue;
                }

                // Vediamo se la traccia interseca due segmenti di currentpolygon
                bool doubleintersections = false;               
                unsigned int countdoubleintersections = 0;
                Vector3d vert = traceverts[0];
                Vector3d direction = traceverts[1] - traceverts[0];
                vector<Vector3d> newvertices;

                // Controlla se la retta della traccia interseca il lato e poi verifica
                // se l'intersezione è interna alla traccia
                for (size_t j = 0; j < currentpolygon.size() - 1; ++j) {
                    if (intersectrettaretta(vert, direction, currentpolygon[j],
                                            currentpolygon[j+1] - currentpolygon[j], inter)) {
                        newvertices.push_back(inter);
                        if (distance(inter, traceverts[0]) + distance(inter, traceverts[1]) -
                            distance(traceverts[0], traceverts[1]) < epsilon) {
                            countdoubleintersections++;
                        }
                    }
                }
                if (intersectrettaretta(vert, direction, currentpolygon[currentpolygon.size()-1],
                                        currentpolygon[0] - currentpolygon[currentpolygon.size()-1], inter)) {
                    newvertices.push_back(inter);
                    if (distance(inter, traceverts[0]) + distance(inter, traceverts[1]) -
                        distance(traceverts[0], traceverts[1]) < epsilon) {
                        countdoubleintersections++;
                    }
                }
                if (countdoubleintersections == 2) {
                    doubleintersections = true;
                }

                bool raycasting = false;
                if (!doubleintersections) {
                    // Usiamo il ray-casting sugli estremi della traccia
                    for (size_t k = 0; k < traceverts.size(); ++k) {
                        Vector3d vert = traceverts[k];
                        unsigned int countraycasting = 0;
                        unsigned int vertexinter = 0;
                        if (!raycasting) {
                            direction = traceverts[1 - k] - vert;
                            for (size_t j = 0; j < currentpolygon.size() - 1; ++j) {
                                if (intersectrettasemiretta(vert, direction, currentpolygon[j],
                                                            currentpolygon[j+1] - currentpolygon[j], inter)) {
                                    if (distance(inter, vert) < epsilon) {
                                        raycasting = true;
                                        break;
                                    }
                                    if (distance(inter,currentpolygon[j])) {
                                        vertexinter++;
                                    }
                                    countraycasting++;
                                }
                            }
                            if (raycasting) {
                                break;
                            }
                            if (intersectrettasemiretta(vert, direction, currentpolygon[currentpolygon.size()-1],
                                                        currentpolygon[0] - currentpolygon[currentpolygon.size()-1], inter)) {
                                if (distance(inter, vert) < epsilon) {
                                    raycasting = true;
                                    break;
                                }
                                if (distance(inter,currentpolygon[currentpolygon.size()])) {
                                    vertexinter++;
                                }
                                countraycasting++;
                            }
                            countraycasting -= vertexinter/2;
                            if (countraycasting == 1) {
                                raycasting = true;
                                break;
                            }
                        }
                    }
                }

                // cout << "Id frattura: " << id << "; Id traccia: " << id_t << "; raycasting: " << raycasting << "; doubleintersections: " << doubleintersections << endl;
                // Se nessuna delle condizioni è soddisfatta, salta il sottopoligono corrente
                if (!raycasting && !doubleintersections) {
                    copia.push_back(currentpolygon);                    
                    continue;
                }               

                // Elimina gli elementi ripetuti in newvertices
                for (size_t t = 0; t < newvertices.size()-1; ++t) {
                    for (size_t j = t + 1; j < newvertices.size(); ++j) {
                        if (distance(newvertices[t], newvertices[j]) < epsilon) {
                            newvertices.erase(newvertices.begin()+j);
                        }
                    }
                }

                if (newvertices.size() == 2) {
                    vector<Vector3d> sottopol1;
                    vector<Vector3d> sottopol2;
                    sottopol1.reserve(currentpolygon.size()+1);
                    sottopol2.reserve(currentpolygon.size()+1);
                    sottopol1.push_back(newvertices[0]);
                    sottopol1.push_back(newvertices[1]);
                    sottopol2.push_back(newvertices[0]);
                    sottopol2.push_back(newvertices[1]);
                    for (size_t l = 0; l< currentpolygon.size(); ++l) {
                        if (distance(newvertices[0], currentpolygon[l]) < epsilon ||
                            distance(newvertices[1], currentpolygon[l]) < epsilon) {
                            continue;
                        }
                        if(((newvertices[0]-currentpolygon[l]).cross(newvertices[1]-currentpolygon[l])).dot(fractures.Normals[id]) > 0) {
                            sottopol1.push_back(currentpolygon[l]);
                        }
                        else{
                            sottopol2.push_back(currentpolygon[l]);
                        }
                    }
                    antiorario(sottopol1, fractures.Normals[id]);
                    antiorario(sottopol2, fractures.Normals[id]);

                    // Memorizza i sottopoligoni in copia
                    copia.push_back(sottopol1);
                    copia.push_back(sottopol2);

                    // Aggiunge il vertice se si trova sulla traccia
                    // for (size_t tr = 0; tr < id_t; ++tr) {
                    //     traceverts = traces.TracesExtremesCoord[tr];
                    //     for (const Vector3d& vertex : newvertices) {
                    //         if (abs(distance(vertex,traceverts[0])+distance(vertex,traceverts[1])
                    //                 -distance(traceverts[0], traceverts[1])) <= 1e-6) {
                    //             for (size_t cp = 0; cp < sottopoligoni.size(); ++cp) {
                    //                 if (sottopoligoni[cp] != currentpolygon) {
                    //                     bool vertexadded = false;
                    //                     for (size_t l = 0; l < sottopoligoni[cp].size()-1; ++l) {
                    //                         if (abs(distance(vertex, sottopoligoni[cp][l])+distance(vertex, sottopoligoni[cp][l+1])
                    //                                 -distance(sottopoligoni[cp][l], sottopoligoni[cp][l+1])) >= 1e-6) {
                    //                             if (!vertexadded) {
                    //                                 sottopoligoni[cp].push_back(vertex);
                    //                                 vertexadded = true;
                    //                             }
                    //                         }
                    //                     }
                    //                     if (abs(distance(vertex, sottopoligoni[cp][sottopoligoni[cp].size()])+distance(vertex, sottopoligoni[cp][0])
                    //                             -distance(sottopoligoni[cp][sottopoligoni[cp].size()], sottopoligoni[cp][0])) >= 1e-6) {
                    //                         if (!vertexadded) {
                    //                             sottopoligoni[cp].push_back(vertex);
                    //                             vertexadded = true;
                    //                         }
                    //                     }
                    //                     if (vertexadded) {
                    //                         antiorario(sottopoligoni[cp], fractures.Normals[id]);
                    //                     }
                    //                 }
                    //             }
                    //         }
                    //     }
                    // }
                }
                else {
                    copia.push_back(currentpolygon);
                }
            }
            sottopoligoni = move(copia); // Memorizza copia in sottopoligoni
            copia.clear(); // Svuota copia
        }
        // Elimina gli elementi ripetuti in sottopoligoni
        // for (size_t t = 0; t < sottopoligoni.size(); ++t) {
        //     for (size_t n = 0; n < sottopoligoni[t].size()-1; ++n) {
        //         for (size_t j = n + 1; j < sottopoligoni[t].size(); ++j) {
        //             if (distance(sottopoligoni[t][n], sottopoligoni[t][j]) < 1e-6) {
        //                 sottopoligoni[t].erase(sottopoligoni[t].begin()+j);
        //             }
        //         }
        //     }
        // }

        // Popola mesh
        chrono::steady_clock::time_point t_begin = chrono::steady_clock::now();
        PolygonalMesh fracturemesh;
        unsigned int id0d = 0;
        for (unsigned int n = 0; n < sottopoligoni.size(); ++n) {
            for (unsigned int j = 0; j < sottopoligoni[n].size(); ++j)  {
                bool repetition = false;
                for (size_t k = 0; k < fracturemesh.CoordCell0d.size(); ++k) {
                    if (distance(sottopoligoni[n][j], fracturemesh.CoordCell0d[k]) < epsilon) {
                        repetition = true;
                    }
                }
                if (!repetition) {
                    fracturemesh.CoordCell0d.push_back(sottopoligoni[n][j]);
                    fracturemesh.IdCell0d.push_back(id0d);
                    id0d++;
                }
            }
        }
        fracturemesh.NumberCell0d = fracturemesh.IdCell0d.size();


        unsigned int id1d = 0;
        unordered_map<pair<unsigned int, unsigned int>, bool, pair_hash> checked_pairs;

        // Celle 1d
        for (const auto& poly : sottopoligoni) {
            for (size_t j = 0; j < poly.size(); ++j) {
                size_t next_j = (j + 1) % poly.size();
                for (unsigned int p = 0; p < fracturemesh.CoordCell0d.size(); ++p) {
                    if (distance(fracturemesh.CoordCell0d[p], poly[j]) < epsilon) {
                        for (unsigned int q = 0; q < fracturemesh.CoordCell0d.size(); ++ q) {
                            if (distance(fracturemesh.CoordCell0d[q], poly[next_j]) < epsilon) {
                                auto lato1 = make_pair(fracturemesh.IdCell0d[p], fracturemesh.IdCell0d[q]);
                                auto lato2 = make_pair(fracturemesh.IdCell0d[q], fracturemesh.IdCell0d[p]);
                                if (checked_pairs.find(lato1) == checked_pairs.end() && checked_pairs.find(lato2) == checked_pairs.end()) {
                                    fracturemesh.IdVerticesCell1d.push_back({fracturemesh.IdCell0d[p], fracturemesh.IdCell0d[q]});
                                    fracturemesh.IdCell1d.push_back(id1d);
                                    id1d++;
                                    checked_pairs[lato1] = true;
                                    checked_pairs[lato2] = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        fracturemesh.NumberCell1d = fracturemesh.IdCell1d.size();

        // Celle 2d
        fracturemesh.IdCell2d.reserve(sottopoligoni.size());
        fracturemesh.IdVerticesCell2d.reserve(sottopoligoni.size());

        unordered_map<Vector3d, unsigned int, vector3d_hash> coord_to_id;
        for (unsigned int i = 0; i < fracturemesh.CoordCell0d.size(); ++i) {
            coord_to_id[fracturemesh.CoordCell0d[i]] = fracturemesh.IdCell0d[i];
        }

        for (unsigned int n = 0; n < sottopoligoni.size(); ++n) {
            if (sottopoligoni[n].size() > 2) {
                fracturemesh.IdCell2d.push_back(n);
                vector<unsigned int> verticescell2d;

                for (const auto& vertex : sottopoligoni[n]) {
                    auto it = coord_to_id.find(vertex);
                    if (it != coord_to_id.end()) {
                        verticescell2d.push_back(it->second);
                    }
                }

                fracturemesh.IdVerticesCell2d.push_back(verticescell2d);
                vector<unsigned int> edgescell2d;

                for (unsigned int p = 0; p < fracturemesh.IdVerticesCell1d.size(); ++p) {
                    for (unsigned int j = 0; j < verticescell2d.size(); ++j) {
                        unsigned int next_j = (j + 1) % verticescell2d.size(); // next index, wraps around to 0
                        vector<unsigned int> extremesids1 = {verticescell2d[j], verticescell2d[next_j]};
                        vector<unsigned int> extremesids2 = {verticescell2d[next_j], verticescell2d[j]};

                        if (extremesids1 == fracturemesh.IdVerticesCell1d[p] || extremesids2 == fracturemesh.IdVerticesCell1d[p]) {
                            edgescell2d.push_back(fracturemesh.IdCell1d[p]);
                        }
                    }
                }

                fracturemesh.IdEdgesCell2d.push_back(edgescell2d);
            }
        }
        fracturemesh.NumberCell2d = fracturemesh.IdCell2d.size();

        fileout << "Id frattura: " << id << endl;
        fileout << "Numero di tracce: " << idtraces.size() << endl;
        fileout << "Numero di sottopoligoni: " << fracturemesh.NumberCell2d << endl;
        fileout << "Il numero di celle 0d è: " << fracturemesh.NumberCell0d << endl;
        fileout << "IdCell0d: X; Y; Z" << endl;


        for (unsigned int n = 0; n < fracturemesh.NumberCell0d; ++n) {
            fileout << fracturemesh.IdCell0d[n] << ": ";
            for (unsigned int j = 0; j < 3; ++j) {
                fileout << fracturemesh.CoordCell0d[n](j) << "; ";
            }
            fileout << endl;
        }


        fileout << "Il numero di celle 1d è: " << fracturemesh.NumberCell1d << endl;
        fileout << "Id Cell1d: Origin; End" << endl;
        for (unsigned int n = 0; n < fracturemesh.NumberCell1d; ++n) {
            fileout << fracturemesh.IdCell1d[n] << ": ";
            for (unsigned int j = 0; j < 2; ++j) {
                fileout << fracturemesh.IdVerticesCell1d[n][j] << "; ";
            }
            fileout << endl;
        }

        fileout << "Il numero di celle 2d è: " << fracturemesh.NumberCell2d << endl;
        fileout << "IdCells2d: NumVertices; Vertices; NumEdges; Edges" << endl;
        for (unsigned int n = 0; n < fracturemesh.NumberCell2d; ++n) {
            fileout  << fracturemesh.IdCell2d[n] << ": " ;
            fileout << fracturemesh.IdVerticesCell2d[n].size() << "; " ;
            for (unsigned int j = 0; j < fracturemesh.IdVerticesCell2d[n].size(); ++j) {
                fileout << fracturemesh.IdVerticesCell2d[n][j] << "; ";
            }

            fileout << fracturemesh.IdEdgesCell2d[n].size() << "; " ;
            for (unsigned int j = 0; j < fracturemesh.IdEdgesCell2d[n].size(); ++j) {
                fileout << fracturemesh.IdEdgesCell2d[n][j] << "; ";
            }
            fileout << endl;
        }

        fileout << endl;

        mesh.push_back(fracturemesh);

        chrono::steady_clock::time_point t_end = chrono::steady_clock::now();
        double tempoTrascorso = chrono::duration_cast<chrono::microseconds>(t_end-t_begin).count();
        cout << "Tempo mesh: " << tempoTrascorso << "micros. " << endl;
    }
    fileout.close();
    return true;
}
}

namespace Export{
void exportMesh(const vector<MeshLibrary::PolygonalMesh>& mesh){
    unsigned int numCols = mesh[7].NumberCell0d;


    // Creo la matrice di dimensioni 3xNumCols
    MatrixXd matrix(3, numCols);

    for (unsigned int i = 0; i < numCols; ++i) {
        matrix.col(i) = mesh[7].CoordCell0d[i];
    }

    Gedim::UCDUtilities exporter;
    string fileName = "./ExportMeshVertices.inp";
    exporter.ExportPoints( fileName, matrix,{},{});


    unsigned int numCols2 = mesh[7].NumberCell1d;

    // Crea la matrice di dimensioni 2xNumCols
    MatrixXi matrix2(2, numCols2);

    for (unsigned int i = 0; i < numCols2; ++i) {
        matrix2.col(i) = Vector2i(mesh[7].IdVerticesCell1d[i][0],mesh[7].IdVerticesCell1d[i][1]);
    }
    string fileName1 = "./ExportMeshEdges.inp";
    exporter.ExportSegments( fileName1, matrix,matrix2);

}
}
