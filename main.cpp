#include <stdio.h>
#include "Screen.hpp"
#include "Math.hpp"
#include "Render.hpp"
#include <iostream>

#define WIDTH  300
#define HEIGHT 300

#define RATIO (double(WIDTH) / HEIGHT)
#define SAMPLES  1000

int main(int, char **)
{
    go::Screen s(WIDTH, HEIGHT);

    

    auto sf = s.createSurface();

    

    auto l1 = std::make_shared<go::Lambertian>(Vector4d(0.1,0.5,0.7,1));

    auto l2 = std::make_shared<go::Lambertian>(Vector4d(0.5,0.7,0.1,1));

    auto m1 = std::make_shared<go::Metal>(Vector4d(0.7,0.7,0.7,1),0);

    auto d1 = std::make_shared<go::Dielectric>(1.5);
    
    go::Camera c(Eigen::Vector3d(-1, 1, 2), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0.3, 0),pi / 2, WIDTH / HEIGHT,2.1,pi / 2);

    go::Scene sc(100);

    go::Sphere* sq = new go::Sphere(Vector3d(-1.5, 0, -2.5), 0.5,l1);

    go::Sphere* sq4 = new go::Sphere(Vector3d(-1.5, 3.5, -2.5), 0.5,l1);

    go::Sphere* sq2 = new go::Sphere(Vector3d(0., 0,0.),Vector3d(0., 0,0.), 0.5,m1);

     go::Sphere* sq3 = new go::Sphere(Vector3d(-1., 0,0.), 0.5,d1);

    go::Planer* pla = new go::Planer(Vector3d(0,-0.5,0),Vector3d(0,1,0),l2);

    


    sc.add(sq2);
    sc.add(sq);
    sc.add(sq3);
    sc.add(sq4);
    sc.add(pla);

    go::Render render(sf,c,sc);


    render.draw(SAMPLES);

    s.showWindow(sf);

    s.wait([](auto &sceen) {

    });
}
