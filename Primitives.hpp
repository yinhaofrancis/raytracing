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


    struct Vertex{
        Vector3d point;
        Vector2d uv;
    };
    class Triangle : public Hitable, public Textureable
    {
    public:
        Triangle(Vertex point1, Vertex point2, Vertex point3, std::shared_ptr<Material>);
        Triangle(Vector3d&& point1, Vector3d&& point2, Vector3d&& point3, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        virtual void uv(Vector2d &uv, const Vector3d &point);
        bool contain(Vector3d &p);

    private:

        bool same_size(Vector3d &point1, Vector3d &point2, Vector3d &point3,Vector3d &p);
        std::shared_ptr<Material> m_mat;
        Vertex m_point[3];
        Vector3d m_normal;
    };


} // namespace go
