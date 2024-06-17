#include "Render.hpp"
#include <iostream>
#include <atomic>
go::Render::Render(go::Surface &sur, go::Camera &c, go::Scene &s) : m_surface(sur), m_camera(c), m_scene(s)
{
}

go::Render::~Render()
{
}

go::Scene &go::Render::scene()
{
    return m_scene;
}

go::Camera &go::Render::camera()
{
    return m_camera;
}

void go::Render::draw(int samples)
{

    int max = 20;
    std::atomic_int b = 0;
    for (size_t t = 0; t < max; t++)
    {
        int index = t;
        std::thread m([this, samples,index,max,&b](){
            for (size_t i = index; i < m_surface.width(); i+=max)
            {
               
                for (size_t j = 0; j < m_surface.height(); j++)
                {
                    renderPixel(i,j,samples);
                }
                std::cout << "col:" << i << std::endl;
            } 
            b+=1;
        });
        m.detach();
    }
    while (b < max)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

void go::Render::renderPixel(int x, int y, int samples)
{
    Vector4d color(0, 0, 0, 0);
    for (int k = 0; k < samples; k++)
    {
        auto a = m_surface.pixelToScreenMultiSample(x, y);
        auto ray = m_camera.createRay(a);
        Vector4d tcolor = m_scene.hit(ray);
        color += tcolor;
    }
    color = color / samples;
    m_surface.draw(gamma(color, 1.0), x, y);
}
