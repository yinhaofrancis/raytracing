#include "Raytracing.hpp"

go::Color::Color(Vector3d &albedo) : m_albedo(albedo)
{
}

go::Color::~Color() {}

Vector3d go::Color::color(Vector2d uv, Vector3d &point)
{
    return m_albedo;
}

Vector3d go::TestColor::color(Vector2d uv, Vector3d &point)
{
    return uv.homogeneous();
    // return (point + Vector3d(1,1,1)) / 2.0;
}

go::Pixel::Pixel(uint32_t w, uint32_t h) : m_pixels(new uint8_t[w * h * 3]), m_w(w), m_h(h), m_p(3)
{
}

Vector3d go::Pixel::color(Vector2d uv, Vector3d &point)
{
    Vector2d size = Vector2d(double(m_w - 1), double(m_h - 1));
    Vector2d center = uv;
    Vector2d index = center.cwiseProduct(size);

    if (is_linear)
    {
        Vector2d p00 = Vector2d(floor(index.x()), floor(index.y()));
        Vector2d p10 = Vector2d(ceil(index.x()), floor(index.y()));
        Vector2d p01 = Vector2d(floor(index.x()), ceil(index.y()));
        Vector2d p11 = Vector2d(ceil(index.x()), ceil(index.y()));

        Vector3d p00_color = (*this)(p00.cast<int>());
        Vector3d p10_color = (*this)(p10.cast<int>());
        Vector3d p01_color = (*this)(p01.cast<int>());
        Vector3d p11_color = (*this)(p11.cast<int>());

        double dx = index.x() - floor(index.x());
        double dy = index.y() - floor(index.y());

        return p00_color * dx * dy + p10_color * (1 - dx) * dy + p01_color * dx * (1 - dy) + (1 - dx) * (1 - dy) * p11_color;
    }
    else
    {

        int x = ceil(index.x()) - index.x() > 0.5 ? floor(index.x()) : ceil(index.x());
        int y = ceil(index.y()) - index.y() > 0.5 ? floor(index.y()) : ceil(index.y());
        Vector3d color = (*this)(Vector2i(x, y));
        // std::cout <<uv.x() << "|" << uv.y() << std::endl;
        return color;
    }
}

Vector3d go::Pixel::operator()(Vector2i uv)
{
    if (uv.x() < 0)
    {
        uv[0] = 0;
    }
    if (uv.y() < 0)
    {
        uv[1] = 0;
    }
    if (uv.x() >= m_w)
    {
        uv[0] = m_w - 1;
    }
    if (uv.y() >= m_h)
    {
        uv[1] = m_h - 1;
    }
    int x = (uv.y() * m_w + uv.x()) * m_p;
    u_int8_t *p = m_pixels + x;
    u_int8_t r = *p;
    u_int8_t g = *(p + 1);
    u_int8_t b = *(p + 2);
    return Vector3d(r / 255.0, g / 255.0, b / 255.0);
}

void go::Pixel::assign(uint8_t *buffer, uint32_t size)
{

    size_t copys = std::min(size, m_h * m_p * m_w);

    memcpy(m_pixels, buffer, copys);
}