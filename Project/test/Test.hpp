#pragma once

#include <gtest/gtest.h>
#include "Eigen/Eigen"
#include "Utils.hpp"
#include "Fractures.hpp"
#include <iostream>
#include <cmath>

using namespace Eigen;
using namespace std;
using namespace Polygons;

double tol=10*numeric_limits<double>::epsilon();

TEST(DFNImportTest, CorrectNumberOfFractures) {
    Fractures fractures;
    importdfn("/home/alberto/Documenti/Materiale scuola Alberto/Appunti/Programmazione e calcolo scientifico/Consegne/Progetto_PCS_2024/Project/DFN/FR3_data.txt", fractures);

    ASSERT_EQ(fractures.FracturesNumber, 3); // Il file contiene veramente 3 fratture
}

TEST(DFNImportTest, CorrectFractures) {

    Fractures fractures;
    importdfn("/home/alberto/Documenti/Materiale scuola Alberto/Appunti/Programmazione e calcolo scientifico/Consegne/Progetto_PCS_2024/Project/DFN/FR3_data.txt", fractures);

    ASSERT_EQ(fractures.FracturesId[0], 0); // Controlliamo che gli id siano memorizzati correttamente
    ASSERT_EQ(fractures.FracturesId[1], 1);
    ASSERT_EQ(fractures.FracturesId[2], 2);

    ASSERT_EQ(fractures.CoordVertices[0].size(), 4); // La prima frattura ha 4 vertici
    ASSERT_EQ(fractures.CoordVertices[1].size(), 4); // La seconda frattura ha 4 vertici

    ASSERT_TRUE(abs(fractures.CoordVertices[0][0].x() - 0) < tol); // Controlliamo che alcune coordinate siano memorizzate correttamente
    ASSERT_TRUE(abs(fractures.CoordVertices[0][0].y() - 0) < tol);
    ASSERT_TRUE(abs(fractures.CoordVertices[0][0].z() - 0) < tol);
    ASSERT_TRUE(abs(fractures.CoordVertices[0][1].x() - 1) < tol);
    ASSERT_TRUE(abs(fractures.CoordVertices[0][2].x() - 1) < tol);
    ASSERT_TRUE(abs(fractures.CoordVertices[0][3].x() - 0) < tol);
}

TEST(DFNAnalytics, Spheres){
    // testo la funzione calcspheres su un quadrato
    //Fractures fractures;

    vector<Vector3d> quadrato = {Vector3d(0,1,1),Vector3d(-1,0,0),Vector3d(0,-1,-1),Vector3d(1,0,0)};
    Vector4d sferaquad = Analytics::calcsphere(quadrato);

    ASSERT_TRUE(abs(sferaquad(0)) < tol);
    ASSERT_TRUE(abs(sferaquad(1)) < tol);
    ASSERT_TRUE(abs(sferaquad(2)) < tol);

    ASSERT_NEAR(sferaquad(3), sqrt(2.0), tol);

    vector<Vector3d> quadrato2 = {Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(1,1,0),Vector3d(0,1,0)};
    Vector4d sferaquad2 = Analytics::calcsphere(quadrato2);

    ASSERT_NEAR(sferaquad2(0), 0.5, tol);
    ASSERT_NEAR(sferaquad2(1), 0.5, tol);
    ASSERT_NEAR(sferaquad2(2), 0, tol);
    ASSERT_NEAR(sferaquad2(3), 1/sqrt(2.0), tol);
}

TEST(DFNAnalytics, Distance){
    double dist1 = Analytics::distance(Vector3d(1,1,1), Vector3d(0,0,0));
    ASSERT_NEAR(dist1, sqrt(3.0), tol);
    double dist2 = Analytics::distance(Vector3d(-3,2,4), Vector3d(1,2,3));
    ASSERT_NEAR(dist2, sqrt(17.0), tol);
}

TEST(DFNAnalytics, Normals){
    Vector3d norm1 = Analytics::normal({Vector3d(1,0,0), Vector3d(0,1,0), Vector3d(-1,0,0), Vector3d(0,-1,0)});
    ASSERT_NEAR(norm1(0),0,tol);
    ASSERT_NEAR(norm1(1),0,tol);
    ASSERT_NEAR(norm1(2),1,tol);
}

TEST(DFNAnalytics, Intersectrettaretta){
    Vector3d inter;

    //test  quando le rette si intersecano
    bool int1 = Analytics::intersectrettaretta(Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(0,1,0),Vector3d(1,0,0), inter);
    ASSERT_TRUE(int1);//al bordo del segmento
    ASSERT_NEAR(inter.x(),1,tol);
    ASSERT_NEAR(inter.y(),1,tol);
    ASSERT_NEAR(inter.z(),0,tol);

    bool int4 = Analytics::intersectrettaretta(Vector3d(0.5,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), inter);
    ASSERT_TRUE(int4);
    ASSERT_NEAR(inter.x(),0.5,tol);
    ASSERT_NEAR(inter.y(),2,tol);
    ASSERT_NEAR(inter.z(),0,tol);

    bool int6 = Analytics::intersectrettaretta(Vector3d(0,0,0),Vector3d(1,1,1),Vector3d(1,0,0),Vector3d(0,1,1), inter);
    ASSERT_TRUE(int6); // retta 3D
    ASSERT_NEAR(inter.x(),1,tol);
    ASSERT_NEAR(inter.y(),1,tol);
    ASSERT_NEAR(inter.z(),1,tol);

    //test quando le rette non si intersecano
    bool int2 = Analytics::intersectrettaretta(Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(1,0,0), inter);
    ASSERT_FALSE(int2);
    ASSERT_EQ(inter, Vector3d(1000,1000,1000));


    //le rette si intersecano al di fuori del segmento
    bool int3 = Analytics::intersectrettaretta(Vector3d(2,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), inter);
    ASSERT_FALSE(int3);
    ASSERT_EQ(inter, Vector3d(1000,1000,1000));

    bool int5 = Analytics::intersectrettaretta(Vector3d(0,2,0),Vector3d(1,0,0),Vector3d(0.5,0,0),Vector3d(0,1,0), inter);
    ASSERT_FALSE(int5);
    ASSERT_EQ(inter, Vector3d(1000,1000,1000));
}


TEST(DFNAnalytics, Intersectrettasemiretta){
    Vector3d control;

    //test  quando le rette si intersecano
    bool s_int1 = Analytics::intersectrettasemiretta(Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(0,1,0),Vector3d(1,0,0), control);
    ASSERT_TRUE(s_int1);//al bordo del segmento
    ASSERT_NEAR(control.x(),1,tol);
    ASSERT_NEAR(control.y(),1,tol);
    ASSERT_NEAR(control.z(),0,tol);

    bool s_int4 = Analytics::intersectrettasemiretta(Vector3d(0.5,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), control);
    ASSERT_TRUE(s_int4);
    ASSERT_NEAR(control.x(),0.5,tol);
    ASSERT_NEAR(control.y(),2,tol);
    ASSERT_NEAR(control.z(),0,tol);


    //test quando le rette non si intersecano
    bool s_int2 = Analytics::intersectrettasemiretta(Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(1,0,0), control);
    ASSERT_FALSE(s_int2);
    ASSERT_EQ(control, Vector3d(1000,1000,1000));


    //le rette si intersecano al di fuori del segmento
    bool s_int3 = Analytics::intersectrettasemiretta(Vector3d(2,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), control);
    ASSERT_FALSE(s_int3);
    ASSERT_EQ(control, Vector3d(1000,1000,1000));


    //le rette si intersecano con t negativo
    bool s_int5 = Analytics::intersectrettasemiretta(Vector3d(0.5,0,0),Vector3d(0,1,0),Vector3d(0,-1,0),Vector3d(1,0,0), control);
    ASSERT_FALSE(s_int5);
    ASSERT_EQ(control, Vector3d(1000,1000,1000));
}

TEST(DFNOutputfilesTools, Printtraces){
    string tipsfileout = "traces_3.txt";

    // Legge il contenuto del file riga per riga
    std::ifstream infile(tipsfileout);
    if (!infile.is_open()) {
        cerr << "Failed to open file " << tipsfileout << endl;
    }

    // Vettore contenente l'output atteso riga per riga
    vector<string> expected_lines = {
        "# Number of Traces",
        "2",
        "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2",
        "0; 0; 1; 0.8; 0; 0; 0.8; 1; 0; ",
        "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2",
        "1; 0; 2; 0; 0.5; 0; 0.316184; 0.5; 0; "
    };

    string line;
    size_t line_index = 0;
    while (getline(infile, line)) {
        if (line_index < expected_lines.size()) {
            ASSERT_EQ(line, expected_lines[line_index]) << "Mismatch in line " << line_index + 1;
        }
        else {
            FAIL() << "Unexpected additional line: " << line;
        }
        ++line_index;
    }

    // Controlla che non ci siano meno righe di quelle attese
    ASSERT_EQ(line_index, expected_lines.size()) << "File has fewer lines than expected.";

    infile.close();

    // Rimuove il file temporaneo dal filesystem (remove prende in input un puntatore (const char*),
    // per ottenerlo devo usare il metodo .c_str() sulla stringa)
    remove(tipsfileout.c_str());
}

TEST(PrintToFileTest, Printtips) {
    string tipsfileout = "tips_3.txt";

    // Legge il contenuto del file riga per riga
    std::ifstream infile(tipsfileout);
    if (!infile.is_open()) {
        cerr << "Failed to open file " << tipsfileout << endl;
    }

    // Vettore contenente l'output atteso riga per riga
    vector<string> expected_lines = {
        "# FractureId; NumTraces",
        "0; 2",
        "# TraceId; Tips; Length",
        "0; True; 1",
        "# TraceId; Tips; Length",
        "1; False; 0.316184",
        "",
        "# FractureId; NumTraces",
        "1; 1",
        "# TraceId; Tips; Length",
        "0; True; 1",
        "",
        "# FractureId; NumTraces",
        "2; 1",
        "# TraceId; Tips; Length",
        "1; False; 0.316184",
        ""
    };

    string line;
    size_t line_index = 0;
    while (getline(infile, line)) {
        if (line_index < expected_lines.size()) {
            ASSERT_EQ(line, expected_lines[line_index]) << "Mismatch in line " << line_index + 1;
        }
        else {
            FAIL() << "Unexpected additional line: " << line;
        }
        ++line_index;
    }

    // Controlla che non ci siano meno righe di quelle attese
    ASSERT_EQ(line_index, expected_lines.size()) << "File has fewer lines than expected.";

    infile.close();

    // Rimuove il file temporaneo dal filesystem
    remove(tipsfileout.c_str());
}
