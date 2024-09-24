#include <SDL2/SDL.h>
#include "Ray.hpp"

namespace vr
{
    class Window
    {
    public:
        Window(Uint width,Uint height);
        uint32_t& pixel(Uint x,Uint y);
        void pixel(uint8_t r,uint8_t g,uint8_t b,Uint x,Uint y);
        void pixel(Vec3 pixel,Uint x,Uint y);
        void update() const;
        void wait() const;
        void write(const char* name) const;
        ~Window();
    private:
        uint32_t* m_buffer;
        Uint m_width,m_height;
        SDL_Window *m_window = nullptr;
        SDL_Renderer *m_render = nullptr;
        SDL_Texture *m_texture = nullptr;
    };
} // namespace vr
