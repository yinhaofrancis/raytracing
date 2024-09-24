#pragma once
#include "Util.hpp"
#include "Ray.hpp"
#include <list>
#include <memory>
namespace vr
{
    enum HitMode
    {
        front = 1,
        back = 2
    };
    struct HitRecord
    {
        Vec3 position;
        Vec3 normal;
        Scalar t;
        Vec3 uvinterp;
        HitMode hitMode;
        Vec2 uv;
    };

    class Primitives
    {
    public:
        virtual bool hitTest(
            const Ray &ray,
            const Region &region, HitRecord &record) const = 0;

        virtual const AABB &bbox() const = 0;

        inline void setHitMode(HitMode mode) { m_hit_mode = mode; }

        inline HitMode hitMode() const { return m_hit_mode; }

        virtual Scalar area() const = 0;

        virtual Vec3 random() const = 0;

        virtual Vec3 lightNormal(const Vec3 &lightDirection) const = 0;

    protected:
        static Random::Generate s_random;

    private:
        HitMode m_hit_mode = vr::front;
    };

    class Triangle : public Primitives
    {
    public:
        Triangle(std::array<Vec3, 3> &points);
        Triangle(std::array<Vec3, 3> &points, std::array<Vec3, 3> &normal);
        Triangle(std::array<Vec3, 3> &points, std::array<Vec2, 3> &uv, std::array<Vec3, 3> &normal);

        virtual bool hitTest(
            const Ray &ray,
            const Region &region,
            HitRecord &record) const;

        virtual const AABB &bbox() const;

        virtual Scalar area() const;

        virtual Vec3 random() const;

        virtual Vec3 normal(const Vec3 &point) const;

        virtual Vec3 lightNormal(const Vec3 &lightDirection) const;

    private:
        const std::array<Vec3, 3> m_points;
        const std::vector<Vec2> m_uv;
        Vec3 m_face_normal;
        const std::vector<Vec3> m_normal;
        AABB m_bbox;
    };

    class Quadrilateral : public Primitives
    {
    public:
        Quadrilateral(std::array<Vec3, 4> &points);
        Quadrilateral(const Vec3 &normal, const Vec3 &up, const Vec3 &center, const Scalar radius);

        virtual bool hitTest(
            const Ray &ray,
            const Region &region,
            HitRecord &record) const;

        virtual const AABB &bbox() const;

        virtual Scalar area() const;

        virtual Vec3 random() const;

        virtual Vec3 lightNormal(const Vec3 &lightDirection) const;

    private:
        void createTriangles(std::array<vr::Vec3, 4> &points);
        std::vector<Triangle> m_triangles;
        AABB m_aabb;
    };

    class Sphere : public Primitives
    {
    public:
        Sphere(const Vec3 &center, Scalar radius);
        virtual bool hitTest(
            const Ray &ray,
            const Region &region,
            HitRecord &record) const;

        virtual const AABB &bbox() const;

        virtual Vec3 random() const;

        virtual Scalar area() const;

        virtual Vec3 lightNormal(const Vec3 &lightDirection) const;

    private:
        Vec3 m_center;

        Scalar m_radius;

        AABB m_bbox;
    };

    class ObjectReader{
    public:
        ObjectReader(const char *path);
        void addBias(std::size_t &size, std::string &p);
        const std::list<std::shared_ptr<Triangle>>& tiangles() const;

    private:
        std::list<std::shared_ptr<Triangle>> m_tri;
    };
    // class Plant : public Primitives
    // {
    // public:
    //     Plant(Vec3 center, Vec3 normal);
    //     virtual bool hitTest(
    //         const Ray &ray,
    //         const Region &region,
    //         HitRecord &record) const;

    //     virtual const AABB &bbox() const;

    //     virtual Vec3 random() const;

    //     virtual Scalar area() const;

    //     virtual Vec3 lightNormal(const Vec3 &lightDirection) const;

    // private:
    //     Vec3 m_center;
    //     Vec3 m_normal;
    //     AABB m_bbox;
    // };
} // namespace vr
