#include "Screen.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include "Math.hpp"

go::Screen::Screen(int w, int h) : m_height(h), m_width(w)
{
}

go::Screen::~Screen()
{
    SDL_DestroyRenderer(m_renderer);
    SDL_DestroyWindow(m_window);
}

void go::Screen::showWindow(go::Surface &s)
{
    m_window = SDL_CreateWindow("title", 0, 0, m_width, m_height, 0);

    m_renderer = SDL_CreateRenderer(m_window, -1, SDL_RENDERER_ACCELERATED);

    SDL_RenderClear(m_renderer);
    for (int x = 0; x < s.m_width; x++)
    {
        for (int y = 0; y < s.m_height; y++)
        {
            uint32_t color = s.buffer[x * s.m_height + y];
            uint8_t r = color >> 24;
            uint8_t g = (color & 0x00ff0000) >> 16;
            uint8_t b = (color & 0x0000ff00) >> 8;
            uint8_t a = (color & 0x000000ff);
            SDL_SetRenderDrawColor(m_renderer, r, g, b, a);
            SDL_RenderDrawPoint(m_renderer, x, m_height - y);
        }
    }
    SDL_RenderPresent(m_renderer);
}

go::Surface go::Screen::createSurface()
{
    return Surface(m_width, m_height);
}

Eigen::Vector2d go::Surface::pixelToScreen(int x, int y)
{
    double sx = x / double(m_width);
    double sy = y / double(m_height);
    return Eigen::Vector2d(sx * 2.0 - 1.0, sy * 2.0 - 1.0);
}

Eigen::Vector2d go::Surface::pixelToScreenMultiSample(int x, int y,int i,int sample)
{
    auto r = random_in_square(i,sample);
    double sx = (x + r.x()) / double(m_width);
    double sy = (y + r.y()) / double(m_height);
    return Eigen::Vector2d(sx * 2.0 - 1.0, sy * 2.0 - 1.0);
}

int go::Surface::width()
{
    return m_width;
}

int go::Surface::height()
{
    return m_height;
}

void go::Surface::ppm(const char *path)
{
    std::ofstream f(path);
    f << "P3" << std::endl
      << m_width << " " << m_height << std::endl
      << "255" << std::endl;

    for (int y = m_height - 1; y >= 0; y--)
    {
        for (int x = 0; x < m_width; x++)
        {
            uint32_t color = buffer[x * m_height + y];
            uint8_t r = color >> 24;
            uint8_t g = (color & 0x00ff0000) >> 16;
            uint8_t b = (color & 0x0000ff00) >> 8;
            uint8_t a = (color & 0x000000ff);
            f << uint16_t(r) << " " << uint16_t(g) << " " << uint16_t(b) << std::endl;
        }
    }
    f.close();
}

void go::Surface::invalid()
{
    delete [] buffer;
    buffer = nullptr;
}

void go::Surface::smooth()
{
    int size = 1;
    uint32_t *n_buffer = new uint32_t[m_height * m_width];

    for (int x = 0; x < m_width; x++)
    {
        for (int y = 0; y < m_height; y++)
        {
            uint32_t sumr = 0;
            uint32_t sumg = 0;
            uint32_t sumb = 0;
            uint32_t suma = 0;
            
            for(int i = -size; i <= size; i++){
                for (int j = -size; j <= size; j++)
                {
                    int index = (x + i)* m_height + (y + j);
                    index = index < 0 ? 0 : index;
                    index = index >= m_height * m_width ? m_height * m_width - 1 : index;
                    uint32_t color = buffer[index];
                    sumr += color >> 24;
                    sumg += (color & 0x00ff0000) >> 16;
                    sumb += (color & 0x0000ff00) >> 8;
                    suma += (color & 0x000000ff);
                }
                
            }
            sumr = sumr / (size * 2 + 1) / (size * 2 + 1);
            sumg = sumg / (size * 2 + 1) / (size * 2 + 1);
            sumb = sumb / (size * 2 + 1) / (size * 2 + 1);
            suma = suma / (size * 2 + 1) / (size * 2 + 1);
            uint32_t sum = ((sumr & 0x000000ff) << 24) + ((sumg & 0x000000ff) << 16) + ((sumb & 0x000000ff) << 8) + (suma & 0x000000ff);
            n_buffer[x * m_height + y] = sum;
        }
    }
    delete [] buffer;
    buffer = n_buffer;
}

void go::Screen::wait(std::function<void(Screen &)> drawcall)
{
    while (true)
    {
        SDL_Event e;
        if (SDL_PollEvent(&e))
        {
            if (e.type == SDL_QUIT)
            {
                break;
            }
            drawcall(*this);
        }
    }
}

go::Surface::Surface(int w, int h) : m_width(w), m_height(h)
{
    buffer = new uint32_t[w * h];
}

go::Surface::~Surface()
{
}

void go::Surface::draw(uint32_t color, int x, int y)
{
    int index = x * m_height + y;
    buffer[index] = color;
}

void go::Surface::draw(double r, double g, double b, double a, int x, int y)
{
    uint8_t cr = SDL_clamp(r, 0, 1) * 255;
    uint8_t cg = SDL_clamp(g, 0, 1) * 255;
    uint8_t cb = SDL_clamp(b, 0, 1) * 255;
    uint8_t ca = SDL_clamp(a, 0, 1) * 255;
    uint32_t v = ca;

    uint32_t color = ((0x00000000 | cr) << 24) + ((0x00000000 | cg) << 16) + ((0x00000000 | cb) << 8) + ca;

    draw(color, x, y);
}

void go::Surface::draw(Eigen::Vector3d color, int x, int y)
{
    draw(color.x(), color.y(), color.z(), 1, x, y);
}

void go::Surface::draw(Eigen::Vector4d color, int x, int y)
{
    draw(color.x(), color.y(), color.z(), color.w(), x, y);
}