#include <stdio.h>
#include "Screen.hpp"
#include "Primitives.hpp"
#include "Render.hpp"
#include <iostream>

#define WIDTH  400
#define HEIGHT 300

#define RATIO (double(WIDTH) / HEIGHT)
#define SAMPLES  100

int main(int, char **)
{
    go::Screen s(WIDTH, HEIGHT);

    go::Pixel* px = new go::Pixel(2,2);
    
    uint8_t pp[] = {
        222,222,222,0,222,0,
        0,0,222,222,222,222

    };
    px->assign(pp,sizeof(pp));

    auto sf = s.createSurface();

    auto l1 = std::make_shared<go::Lambertian>(Vector3d(0.1,0.5,0.7));

    auto l2 = std::make_shared<go::Lambertian>(px);

    auto l3 = std::make_shared<go::Lambertian>(new go::TestColor());

    auto m1 = std::make_shared<go::Metal>(Vector3d(0.9,0.9,0.9),0);

    auto d1 = std::make_shared<go::Dielectric>(1.5);

    auto db = std::make_shared<go::NormalColor>();

    auto light1 = std::make_shared<go::Light>(Vector3d(10,10,10));
    
    go::Camera c(Eigen::Vector3d(0, 0, 3.5), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0.3, 0),pi / 4, RATIO,4.9,0 * pi / 100);
    

    Vector4d ambient = Vector4d(0.1,0.1,0.1,1);

    go::Scene sc(1000,ambient);

    go::Sphere* sq = new go::Sphere(Vector3d(-1.5, 0, -2.5), 0.5,l1);

    go::Sphere* sq4 = new go::Sphere(Vector3d(1.3, 0, 0), 0.5,l2);

    go::Sphere* sq2 = new go::Sphere(Vector3d(0., 0,0.), 0.5,m1);

    go::Sphere* sq5 = new go::Sphere(Vector3d(0, 0.7,-0.7), 0.3,light1);

    go::Sphere* sqlt = new go::Sphere(Vector3d(0., 3,0.), 1,light1);

     go::Sphere* sq3 = new go::Sphere(Vector3d(-1.1, 0,0.), 0.5,d1);

    go::Planer* pla = new go::Planer(Vector3d(0,-0.5,0),Vector3d(0,1,0),l3);

    go::Triangle* tra = new go::Triangle(Vector3d(-1,0.2,-1),Vector3d(1,1,-1),Vector3d(-2,0.5,0),db);

    go::Quad* qua = new go::Quad(Vector3d(1,1,-1),Vector3d(-1,1,-1),Vector3d(-1,-1,-1),db);

    go::Quad* qua1 = new go::Quad(Vector3d(-1,1,-1),Vector3d(-1,1,1),Vector3d(-1,-1,1),db);

    go::Quad* qua2 = new go::Quad(Vector3d(1,-1,-1),Vector3d(1,-1,1),Vector3d(1,1,1),db);

    go::Quad* qua3 = new go::Quad(Vector3d(-1,-1,1),Vector3d(1,-1,1),Vector3d(1,-1,-1),db);

    go::Quad* qua4 = new go::Quad(Vector3d(1,1,1),Vector3d(-1,1,1),Vector3d(-1,1,-1),db);

    go::Quad* qua5 = new go::Quad(Vector3d(0.6,0.99,0.6),Vector3d(-0.6,0.99,0.6),Vector3d(-0.6,0.99,-0.6),light1);


    go::Sphere* qsq = new go::Sphere(Vector3d(-0.5, -0.8, 0), 0.2,l2);

    go::Sphere* qsq1 = new go::Sphere(Vector3d(0.0, -0.8, 0.5), 0.2,m1);

    go::Sphere* qsq2 = new go::Sphere(Vector3d(0.5, -0.6,0.5), 0.2,d1);

    // sc.add(sqlt);
    // sc.add(tra);
    // sc.add(sq5);
    // sc.add(sq2);
    // sc.add(sq);
    // sc.add(sq3);
    // sc.add(sq4);
    // sc.add(pla);
    sc.add(qua);
    sc.add(qua2);
    sc.add(qua1);
    sc.add(qua3);
    sc.add(qua4);
    sc.add(qua5);
    sc.add(qsq);
    sc.add(qsq1);
    sc.add(qsq2);
    go::Render render(sf,c,sc);


    
    

    render.draw(SAMPLES);
    // sf.smooth();
    sf.ppm("m.ppm");
    s.showWindow(sf);
    sf.invalid();


    s.wait([](auto &sceen) {

    });
}
