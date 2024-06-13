#pragma once 

#include "Screen.hpp"
#include "Math.hpp"




namespace go{

class Render{
public:
    Render(Surface&,Camera&,Scene &);
    Scene&   scene();
    Camera&  camera();
    void     draw(int);
private:
    Surface  m_surface;
    Scene   m_scene;
    Camera  m_camera;
};


}