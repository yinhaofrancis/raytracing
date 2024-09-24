#include "Window.hpp"
#include <fstream>
vr::Window::Window(Uint width, Uint height) : m_buffer(new uint32_t[width * height]), m_width(width), m_height(height)
{

    m_window = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_INIT_VIDEO);
    m_render = SDL_CreateRenderer(m_window, -1, 0);
    m_texture = SDL_CreateTexture(m_render, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, width, height);
}

uint32_t &vr::Window::pixel(Uint x, Uint y)
{
    return m_buffer[m_width * y + x];
}

void vr::Window::pixel(uint8_t r, uint8_t g, uint8_t b, Uint x, Uint y)
{
    uint32_t pixel = 0xff;
    pixel <<= 8;
    pixel |= r;
    pixel <<= 8;
    pixel |= g;
    pixel <<= 8;
    pixel |= b;
    this->pixel(x, y) = pixel;
}

void vr::Window::pixel(Vec3 pixel, Uint x, Uint y)
{

    auto pix = glm::clamp(pixel, {0, 0, 0}, {1, 1, 1});
    this->pixel(pix.r * 255, pix.g * 255, pix.b * 255, x, y);
}

void vr::Window::update() const
{
    SDL_UpdateTexture(m_texture, nullptr, m_buffer, sizeof(uint32_t) * m_width);
    SDL_RenderClear(m_render);
    SDL_RenderCopy(m_render, m_texture, nullptr, nullptr);
    SDL_RenderPresent(m_render);
}

void vr::Window::wait() const
{
    bool quit = false;
    SDL_Event event;
    while (!quit)
    {
        SDL_WaitEvent(&event);
        if (event.type == SDL_QUIT)
        {
            quit = true;
        }
    }
}

void vr::Window::write(const char *name) const
{
    std::fstream fstream(name, std::ios::out);
    if (fstream.is_open())
    {
        fstream << "P3" << std::endl;
        fstream << m_width << " "<< m_height << std::endl;
        fstream << "255" << std::endl;
        for (size_t i = 0; i < m_width * m_height; i++)
        {
            Uint r = (m_buffer[i] & 0x00ff0000) >> 16;
            Uint g = (m_buffer[i] & 0x0000ff00) >> 8;
            Uint b = (m_buffer[i] & 0x000000ff);
            fstream << r << " " << g << " " << b << std::endl;
        }
        fstream.close();
        
    }
}

vr::Window::~Window()
{
    delete[] m_buffer;
}
