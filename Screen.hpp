#pragma once

#include <functional>

#include <SDL2/SDL.h>
#include <SDL2/SDL_video.h>
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL_render.h>

#include<eigen3/Eigen/Eigen>

namespace go
{
    class Surface{
    public:
        Surface(int w,int h);
        ~Surface();
        void draw(uint32_t color,int x,int y);

        void draw(double r,double g,double b,double a,int x,int y);

        void draw(Eigen::Vector3d color, int x, int y);

        void draw(Eigen::Vector4d color, int x, int y);


        Eigen::Vector2d pixelToScreen(int x,int y);

        Eigen::Vector2d pixelToScreenMultiSample(int x,int y);

        friend class Screen;

        int width();
        int height();
        void ppm(const char *path);
    private:
        int m_width,m_height;
        uint32_t *buffer = nullptr;
    };
    class Screen{
    public:
        Screen(int w,int h);
        ~Screen();


        void showWindow(Surface &s);
        
        Surface createSurface();

        void wait(std::function<void(Screen&)> drawcall);
    private:
        int m_width,m_height;
        SDL_Window *m_window = nullptr;
        SDL_Renderer *m_renderer = nullptr;
        
    };

} // namespace go