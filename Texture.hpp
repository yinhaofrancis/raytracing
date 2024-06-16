#pragma once
#include <eigen3/Eigen/Eigen>

#include <eigen3/Eigen/Dense>

#include <random>
#include <memory>
#include <vector>
#include "Math.hpp"

using namespace Eigen;

namespace go
{
    class Texture
    {
    public:
        virtual ~Texture() = default;
        virtual Vector3d color(Vector2d uv, Vector3d &point) = 0;
    };
    class Color : public Texture
    {
    public:
        Color(Vector3d &albedo);
        ~Color();
        Vector3d color(Vector2d uv, Vector3d &point);

    private:
        Vector3d m_albedo;
    };
    class Pixel : public Texture
    {
    public:
        Pixel(uint32_t w, uint32_t h);
        Vector3d color(Vector2d uv, Vector3d &point);
        Vector3d operator()(Vector2i uv);
        void assign(uint8_t *buffer, uint32_t size);

    private:
        uint8_t *m_pixels = nullptr;
        uint32_t m_w;
        uint32_t m_h;
        uint32_t m_p = 3;
        bool is_linear = true;
    };

    class TestColor : public Texture
    {
    public:
        Vector3d color(Vector2d uv, Vector3d &point);
    };
}