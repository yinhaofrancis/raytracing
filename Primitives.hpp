#pragma once
#include "Raytracing.hpp"

namespace go
{
    class Sphere : public Hitable
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

    class Planer : public Hitable
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

    struct Vertex
    {
        Vector3d point;
        Vector2d uv;
        Vertex(double x,double y,double z,double u,double v);
        Vertex() {}
    };
    class Triangle : public Hitable
    {
    public:
        Triangle();
        Triangle(Vertex point1, Vertex point2, Vertex point3, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        virtual void uv(Vector2d &uv, const Vector3d &point);
        bool contain(const Vector3d &p);
        void transform(const Matrix4d& transform);
        bool& double_face();
    private:
        bool same_size(const Vector3d &point1, const Vector3d &point2, const Vector3d &point3, const Vector3d &p);
        std::shared_ptr<Material> m_mat;
        Vertex m_point[3];
        bool is_double_face = false;
        Vector3d m_normal;
    };
    class Quad : public Hitable
    {
    public:
    
        Quad() {}
        Quad(std::shared_ptr<Material> m);
        Quad(Vertex point1, Vertex point2, Vertex point3, Vertex point4, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        virtual void uv(Vector2d &uv, const Vector3d &point);
        void transform(const Matrix4d& transform);
        bool& double_face();
    private:
        Triangle m_triangles[2];
        
        bool is_double_face = false;
    };

    class Box : public Hitable
    {
    public:
        Box(std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        virtual void uv(Vector2d &uv, const Vector3d &point);
        void transform(const Matrix4d&  transform);
    private:
        Quad m_quad[6];
    };

} // namespace go
