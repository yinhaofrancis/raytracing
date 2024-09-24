#include <iostream>
#include <array>
#include "ray/Render.hpp"
int main2();

int main()
{
    const char* s1 = "v:111";
    const char* s2 = "v";


    return main2();
}

int main2()
{
    std::cout << "Hello, from rt!\n";

    vr::Camera cam(800, 450, {0, 2, 9}, {0, 1.5, 0}, {0, 1, 0}, M_PI / 3.);

    vr::Quadrilateral lightp({0, -1, 0}, {0, 0, 1}, {0, 4.99, 0}, 1);

    vr::Sphere sp2({2., 1.0, 1.}, 1.0);
    vr::Sphere sp1({-2., 1.0, 1.}, 1.0);
    vr::Sphere sp({0., 1.0, 0.}, 1.0);

    vr::Quadrilateral left({1, 0, 0}, {0, 1, 0}, {-5, 0, 0}, 5);
    vr::Quadrilateral right({-1, 0, 0}, {0, 1, 0}, {5, 0, 0}, 5);
    vr::Quadrilateral top({0, -1, 0}, {0, 0, 1}, {0, 5, 0}, 5);
    vr::Quadrilateral back({0, 0, 1}, {0, 1, 0}, {0, 0, -5}, 5);

    vr::Quadrilateral vp({0, 1, 0}, {1, 0, 0}, {0, -1, 0}, 100);

    vr::Lambertain color3({1, 1., 1});
    vr::Lambertain color4({0.7, 0.7, 0.7});
    vr::Lambertain color2({0.3, .4, 0.3});
    vr::Lambertain color1({0.5, .1, 0.1});
    vr::EmitLight light({10, 10, 10});
    vr::Metal metal1({0.9, .9, 0.9}, 0.);
    vr::Dielectric glass1(1.2);


    vr::ObjectReader vor("/home/yinhao/projects/raytracing/m.obj");
    vr::Scene sc;
    // sc.add(vor,&color1);
    sc.add(&sp, &color1);
    sc.add(&vp, &color2);
    sc.add(&sp2, &metal1);
    sc.add(&left, &color4);
    sc.add(&right, &color4);
    sc.add(&top, &color4);
    sc.add(&back, &color4);
    sc.add(&lightp, &light);
    sc.add(&sp1, &glass1);
    sc.end();
    vr::Render render(100);
    render.gamma() = 1.0 / 1.9;
    render.render(cam, sc);

    return 0;
}