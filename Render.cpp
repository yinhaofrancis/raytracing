#include "Render.hpp"
#include <iostream>
go::Render::Render(go::Surface& sur,go::Camera& c,go::Scene& s):m_surface(sur),m_camera(c),m_scene(s)
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
    for (size_t i = 0; i < m_surface.width(); i++)
    {
        std::cout << "col:" << i << std::endl;
        
        for (size_t j = 0; j < m_surface.height(); j++)
        {
            Vector4d color(0,0,0,0);
            
            for (int k = 0; k < samples; k++)
            {
                
                auto a = m_surface.pixelToScreenMultiSample(i, j);
                auto ray = m_camera.createRay(a);
                Vector4d tcolor = m_scene.hit(ray);
                color += tcolor;
            }

            color = color / samples;
        
            m_surface.draw(gamma(color,1),i,j);
        }
    }
}
