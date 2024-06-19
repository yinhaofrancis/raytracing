#pragma once 

#include "Screen.hpp"
#include "Math.hpp"
#include <thread>
#include "Raytracing.hpp"



namespace go{

class Render{
public:
    Render(Surface&,Camera&,Scene &);
    ~Render();
    Scene&   scene();
    Camera&  camera();
    void     draw(int,double);
private:
    void renderPixel(int x,int y,int samples,double gamma);
    Surface  m_surface;
    Scene   m_scene;
    Camera  m_camera;
};


}