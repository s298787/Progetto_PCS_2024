#pragma once

#include <gtest/gtest.h>
#include "Eigen/Eigen"
#include "Utils.hpp"
#include "Fractures.hpp"
#include <iostream>
#include <sstream>
#include <cmath>

using namespace Eigen;
using namespace std;
using namespace Polygons;

double tol = 100*numeric_limits<double>::epsilon();


TEST(DFNImportTest, CorrectNumberOfFractures) {
    Fractures fractures;
    importdfn("./FR3_data.txt", fractures);

    EXPECT_EQ(fractures.FracturesNumber, 3); // Il file contiene veramente 3 fratture
}

TEST(DFNImportTest, CorrectFractures) {

    Fractures fractures;
    importdfn("./FR3_data.txt", fractures);

    EXPECT_EQ(fractures.FracturesId[0], 0); // Controlliamo che gli id siano memorizzati correttamente
    EXPECT_EQ(fractures.FracturesId[1], 1);
    EXPECT_EQ(fractures.FracturesId[2], 2);

    EXPECT_EQ(fractures.CoordVertices[0].size(), 4); // La prima frattura ha 4 vertici
    EXPECT_EQ(fractures.CoordVertices[1].size(), 4); // La seconda frattura ha 4 vertici

    EXPECT_TRUE(abs(fractures.CoordVertices[0][0].x() - 0) < tol); // Controlliamo che alcune coordinate siano memorizzate correttamente
    EXPECT_TRUE(abs(fractures.CoordVertices[0][0].y() - 0) < tol);
    EXPECT_TRUE(abs(fractures.CoordVertices[0][0].z() - 0) < tol);
    EXPECT_TRUE(abs(fractures.CoordVertices[0][1].x() - 1) < tol);
    EXPECT_TRUE(abs(fractures.CoordVertices[0][2].x() - 1) < tol);
    EXPECT_TRUE(abs(fractures.CoordVertices[0][3].x() - 0) < tol);
}


TEST(DFNAnalytics, BasicSpheres){
    vector<Vector3d> quadrato2 = {Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(1,1,0),Vector3d(0,1,0)};
    Vector4d sferaquad2 = Analytics::calcsphere(quadrato2);

    EXPECT_NEAR(sferaquad2(0), 0.5, tol);
    EXPECT_NEAR(sferaquad2(1), 0.5, tol);
    EXPECT_NEAR(sferaquad2(2), 0, tol);
    EXPECT_NEAR(sferaquad2(3), 1/sqrt(2.0), tol);
}

TEST(DFNAnalytics, 3DSpheres){
    // testo la funzione calcspheres su un quadrato
    //Fractures fractures;

    vector<Vector3d> quadrato = {Vector3d(0,1,1),Vector3d(-1,0,0),Vector3d(0,-1,-1),Vector3d(1,0,0)};
    Vector4d sferaquad = Analytics::calcsphere(quadrato);

    EXPECT_TRUE(abs(sferaquad(0)) < tol);
    EXPECT_TRUE(abs(sferaquad(1)) < tol);
    EXPECT_TRUE(abs(sferaquad(2)) < tol);

    EXPECT_NEAR(sferaquad(3), sqrt(2.0), tol);
}

TEST(DFNAnalytics, BasicDistance){
    double dist1 = Analytics::distance(Vector3d(1,1,1), Vector3d(0,0,0));
    EXPECT_NEAR(dist1, sqrt(3.0), tol);
}
TEST(DFNAnalytics, 3DDistance){
    double dist2 = Analytics::distance(Vector3d(-3,2,4), Vector3d(1,2,3));
    EXPECT_NEAR(dist2, sqrt(17.0), tol);
}

TEST(DFNAnalytics, Normals){
    Vector3d norm1 = Analytics::normal({Vector3d(1,0,0), Vector3d(0,1,0), Vector3d(-1,0,0), Vector3d(0,-1,0)});
    EXPECT_NEAR(norm1(0),0,tol);
    EXPECT_NEAR(norm1(1),0,tol);
    EXPECT_NEAR(norm1(2),1,tol);
}

TEST(DFNAnalytics, RettaRettaIntersAlBordo){
    Vector3d inter;

    //test  quando le rette si intersecano
    bool int1 = Analytics::intersectrettaretta(Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(0,1,0),Vector3d(1,0,0), inter);
    EXPECT_TRUE(int1);//al bordo del segmento
    EXPECT_NEAR(inter.x(),1,tol);
    EXPECT_NEAR(inter.y(),1,tol);
    EXPECT_NEAR(inter.z(),0,tol);
}

TEST(DFNAnalytics, RettaRettaInters){
    //test  quando le rette si intersecano
    Vector3d inter;
    bool int4 = Analytics::intersectrettaretta(Vector3d(0.5,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), inter);
    EXPECT_TRUE(int4);
    EXPECT_NEAR(inter.x(),0.5,tol);
    EXPECT_NEAR(inter.y(),2,tol);
    EXPECT_NEAR(inter.z(),0,tol);
}

TEST(DFNAnalytics, RettaRettaInters3D){
    Vector3d inter;
    bool int6 = Analytics::intersectrettaretta(Vector3d(0,0,0),Vector3d(1,1,1),Vector3d(1,0,0),Vector3d(0,1,1), inter);
    EXPECT_TRUE(int6); // retta 3D
    EXPECT_NEAR(inter.x(),1,tol);
    EXPECT_NEAR(inter.y(),1,tol);
    EXPECT_NEAR(inter.z(),1,tol);
}

TEST(DFNAnalytics, RettaRettaNoInters){
    //test quando le rette non si intersecano
    Vector3d inter;
    bool int2 = Analytics::intersectrettaretta(Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(1,0,0), inter);
    EXPECT_FALSE(int2);
    EXPECT_EQ(inter, Vector3d(1000,1000,1000));
}

TEST(DFNAnalytics, RettaRettaIntersFuori){
    //le rette si intersecano al di fuori del segmento
    Vector3d inter;
    bool int3 = Analytics::intersectrettaretta(Vector3d(2,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), inter);
    EXPECT_FALSE(int3);
    EXPECT_EQ(inter, Vector3d(1000,1000,1000));

    bool int5 = Analytics::intersectrettaretta(Vector3d(0,2,0),Vector3d(1,0,0),Vector3d(0.5,0,0),Vector3d(0,1,0), inter);
    EXPECT_FALSE(int5);
    EXPECT_EQ(inter, Vector3d(1000,1000,1000));
}


TEST(DFNAnalytics, RettaSemirettaIntersAlBordo){
    Vector3d control;

    //test  quando le rette si intersecano
    bool s_int1 = Analytics::intersectrettasemiretta(Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(0,1,0),Vector3d(1,0,0), control);
    EXPECT_TRUE(s_int1);//al bordo del segmento
    EXPECT_NEAR(control.x(),1,tol);
    EXPECT_NEAR(control.y(),1,tol);
    EXPECT_NEAR(control.z(),0,tol);

}
TEST(DFNAnalytics, RettaSemirettaInters){
    Vector3d control;
    //test  quando le rette si intersecano
    bool s_int4 = Analytics::intersectrettasemiretta(Vector3d(0.5,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), control);
    EXPECT_TRUE(s_int4);
    EXPECT_NEAR(control.x(),0.5,tol);
    EXPECT_NEAR(control.y(),2,tol);
    EXPECT_NEAR(control.z(),0,tol);
}
TEST(DFNAnalytics, RettaSemirettaNoInters){
    //test quando le rette non si intersecano
    Vector3d control;
    bool s_int2 = Analytics::intersectrettasemiretta(Vector3d(0,0,0),Vector3d(1,0,0),Vector3d(0,1,0),Vector3d(1,0,0), control);
    EXPECT_FALSE(s_int2);
    EXPECT_EQ(control, Vector3d(1000,1000,1000));
}

TEST(DFNAnalytics, RettaSemirettaIntersFuori){
    //le rette si intersecano al di fuori del segmento
    Vector3d control;
    bool s_int3 = Analytics::intersectrettasemiretta(Vector3d(2,0,0),Vector3d(0,1,0),Vector3d(0,2,0),Vector3d(1,0,0), control);
    EXPECT_FALSE(s_int3);
    EXPECT_EQ(control, Vector3d(1000,1000,1000));
}

TEST(DFNAnalytics, RettaSemirettatNegativo){
    //le rette si intersecano con t negativo
    Vector3d control;
    bool s_int5 = Analytics::intersectrettasemiretta(Vector3d(0.5,0,0),Vector3d(0,1,0),Vector3d(0,-1,0),Vector3d(1,0,0), control);
    EXPECT_FALSE(s_int5);
    EXPECT_EQ(control, Vector3d(1000,1000,1000));
}

TEST(DFNAnalytics, Angolo_xy){
    double ang = Analytics::calcolangolo(Vector3d(1,0,0) , Vector3d(0,1,0), Vector3d(0,0,1));
    EXPECT_NEAR(ang, M_PI/2, tol);
}
TEST(DFNAnalytics, Angolo_xz){
    double ang = Analytics::calcolangolo(Vector3d(1,0,0) , Vector3d(0,0,1), Vector3d(0,1,0));
    EXPECT_NEAR(ang, 3*M_PI/2, tol);
}
TEST(DFNAnalytics, Angolo_yz){
    double ang = Analytics::calcolangolo(Vector3d(0,1,0) , Vector3d(0,0,1), Vector3d(1,0,0));
    EXPECT_NEAR(ang, M_PI/2, tol);
}

TEST(DFNAnalytics, Angolo3D){
    double ang = Analytics::calcolangolo(Vector3d(-1,1,0) , Vector3d(0,-1,1), Vector3d(1,1,1));
    EXPECT_NEAR(ang, 2*M_PI/3, tol);
}

TEST(DFNAnalytics, BasicAntiorario){
    vector<Vector3d> sottopol;
    sottopol.push_back(Vector3d(1,0,0));
    sottopol.push_back(Vector3d(0,1,0));
    sottopol.push_back(Vector3d(0,0,0));
    sottopol.push_back(Vector3d(0.5,2,0));
    sottopol.push_back(Vector3d(1,1,0));
    sottopol.push_back(Vector3d(2,0.5,0));

    Analytics::antiorario(sottopol, Vector3d(0,0,1));

    EXPECT_EQ(sottopol[0], Vector3d(1,1,0));
    EXPECT_EQ(sottopol[1], Vector3d(0.5,2,0));
    EXPECT_EQ(sottopol[2], Vector3d(0,1,0));
    EXPECT_EQ(sottopol[3], Vector3d(0,0,0));
    EXPECT_EQ(sottopol[4], Vector3d(1,0,0));
    EXPECT_EQ(sottopol[5], Vector3d(2,0.5,0));

}
