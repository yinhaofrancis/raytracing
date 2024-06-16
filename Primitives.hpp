#pragma once
#include "Raytracing.hpp"

namespace go
{
    class Sphere : public Hitable, public Textureable
    {
    public:
        Sphere(Vector3d center, double radius, std::shared_ptr<Material>);
        Sphere(Vector3d center, Vector3d center2, double radius, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        virtual void uv(Vector2d &uv, const Vector3d &point);

    private:
        Vector3d center(Ray &);
        std::shared_ptr<Material> m_mat;
        Vector3d m_center;
        Vector3d m_center2;
        double m_radius;
    };

    class Planer : public Hitable, public Textureable
    {
    public:
        Planer(Vector3d point, Vector3d normal, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        virtual void uv(Vector2d &uv, const Vector3d &point);

    private:
        std::shared_ptr<Material> m_mat;
        Vector3d m_point;
        Vector3d m_normal;
        Matrix3d m_uv_back;
    };

    class Triangle : public Hitable, public Textureable
    {
    public:
        Triangle(Vector3d point1, Vector3d point2, Vector3d point3, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        virtual void uv(Vector2d &uv, const Vector3d &point);
    };
} // namespace go
