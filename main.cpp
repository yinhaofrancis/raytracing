#include <stdio.h>
#include "Screen.hpp"
#include "Primitives.hpp"
#include "Render.hpp"
#include <iostream>
// #include <immintrin.h>
#define WIDTH  400
#define HEIGHT 300

#define RATIO (double(WIDTH) / HEIGHT)
#define SAMPLES  10


// void scene1(){

// }


int main(int, char **)
{

    
    go::Screen s(WIDTH, HEIGHT);

    go::Pixel* px = new go::Pixel(2,2);
    
    uint8_t pp[] = {
        222,222,222, 0,222,0,
        0,0,222,     222,222,222

    };
    px->assign(pp,sizeof(pp));

    auto sf = s.createSurface();

    auto l1 = std::make_shared<go::Lambertian>(Vector3d(0.8,0.4,0.3));

    auto l2 = std::make_shared<go::Lambertian>(px);

    auto l3 = std::make_shared<go::Lambertian>(new go::TestColor());

    auto m1 = std::make_shared<go::Metal>(Vector3d(0.9,0.9,0.9),0);

    auto d1 = std::make_shared<go::Dielectric>(1.1);

    auto db = std::make_shared<go::Lambertian>(new go::TestColor());

    auto light1 = std::make_shared<go::Light>(Vector3d(9,9,9));

    auto n = std::make_shared<go::NormalColor>();
    
    go::Camera c(Eigen::Vector3d(0,0,8), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0.3, 0),pi / 4, RATIO,4.8,0 * pi / 100);
    

    Vector4d ambient = Vector4d(0,0,0,1);

    go::Scene sc(100,ambient);

    go::Sphere* sq = new go::Sphere(Vector3d(0, -1, 0), 0.5,l1);

    go::Sphere* sq4 = new go::Sphere(Vector3d(1.3, 0, 0), 0.5,l2);

    go::Sphere* sq2 = new go::Sphere(Vector3d(-1.5, 0,-2.5), 0.5,m1);

    go::Sphere* sq5 = new go::Sphere(Vector3d(0, 0.7,-0.7), 0.3,light1);

    go::Sphere* sqlt = new go::Sphere(Vector3d(0., 3,0.), 1,light1);

     go::Sphere* sq3 = new go::Sphere(Vector3d(-1, 0,1), 0.5,d1);

    go::Planer* pla = new go::Planer(Vector3d(0,-0.5,0),Vector3d(0,1,0),l3);

    go::Triangle* tri = new go::Triangle(
        go::Vertex( 1., 1.,-1.,1.,0),
        go::Vertex(-1., 1.,-1.,0,0),
        go::Vertex(-1.,-1.,-1.,0,1.),
        l2
    );
    
    go::Quad* qua = new go::Quad(
        go::Vertex( 1, 1,0,1,0),
        go::Vertex(-1, 1,0,0,0),
        go::Vertex(-1,-1,0,0,1),
        go::Vertex( 1,-1,0,1,1),m1);

    
    qua->transform(rotate(1,0,0,M_PI / 6));
    // qua->transform(translate(Vector3d(1,1,1)));
    // qua->transform(scale(Vector3d(2,2,1)));


    // go::Quad* qua1 = new go::Quad(Vector3d(-1,1,-1),Vector3d(-1,1,1),Vector3d(-1,-1,1),db);

    // go::Quad* qua2 = new go::Quad(Vector3d(1,-1,-1),Vector3d(1,-1,1),Vector3d(1,1,1),db);

    // go::Quad* qua3 = new go::Quad(Vector3d(-1,-1,1),Vector3d(1,-1,1),Vector3d(1,-1,-1),db);

    // go::Quad* qua4 = new go::Quad(Vector3d(1,1,1),Vector3d(-1,1,1),Vector3d(-1,1,-1),db);

    go::Quad* qua5 = new go::Quad(light1);
    qua5->transform(translate(-2,4,-1));
    qua5->transform(scale(1,1,1));
    qua5->double_face() = true;


    auto box = new go::Box(d1);
    box->transform(scale(2,2,0.1));
    box->transform(translate(0,1.5,0));

    auto room = new go::Room(l2);
    room->transform(scale(2,2,2));

    auto light = new go::Quad(light1);
    light->transform(translate(0,1.99,0));
    // light->transform(scale(-1,1,1));
    light->double_face() = true;
    sc.light(light);
    sc.add(room);
    // sc.add(sqlt);
    // sc.add(tra);
    // sc.add(sq5);
    // sc.add(sq2);
    sc.add(sq);
    // sc.add(sq3);
    // sc.add(sq4);
    // sc.add(pla);
    // sc.add(tri);
    // sc.add(qua);
    // sc.add(box);
    // sc.add(qua2);
    // sc.add(qua1);
    // sc.add(qua3);
    // sc.add(qua4);
    // sc.add(qua5);
    // sc.add(qsq);
    // sc.add(qsq1);
    // sc.add(qsq2);
    go::Render render(sf,c,sc);


    
    

    render.draw(SAMPLES,1);
    // sf.smooth();
    sf.ppm("m.ppm");
    s.showWindow(sf);
    sf.invalid();


    s.wait([](auto &sceen) {

    });
}
