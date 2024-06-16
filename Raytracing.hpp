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
    class Textureable
    {
    public:
        virtual void uv(Vector2d &uv, const Vector3d &point) = 0;
    };

    class Ray
    {
    public:
        Ray();
        Ray(Vector3d location, Vector3d direction);
        void lambertian(Vector3d normal, Vector3d hitPoint);
        void reflect(Vector3d normal, Vector3d hitPoint, double fuzz);
        void refract(Vector3d normal, Vector3d hitPoint, double fract);
        void hitRandom(Vector3d normal, Vector3d hitPoint);
        Vector3d location() const;
        Vector3d direction() const;
        Vector3d exec(double t);

        int depth();

        double time();

    private:
        int m_max = max_depth;
        Vector3d m_location;
        Vector3d m_direction;
        double m_time;
    };

    class Camera
    {
    public:
        Camera(Vector3d location, Vector3d target, Vector3d up, double fov, double ratio, double length, double defocus);
        Ray createRay(Vector2d &screenUV);

    private:
        Matrix3d m_lookAt;
        Vector3d m_location;
        double m_length;
        double m_ratio;
        double m_fov;
        double m_disk_radius;
    };

    struct HitResult;
    class Material
    {
    public:
        Material() {}
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);
        virtual Vector3d emitted(HitResult &hit);
    };

    struct HitResult
    {
        Vector3d hit, normal;
        Vector2d uv;
        std::shared_ptr<Material> mat;
        double t;
        bool isFront;
    };

    class Hitable
    {
    public:
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        bool &isFrontFace();
        virtual ~Hitable();

    private:
        bool m_front = true;
    };

    class Lambertian : public Material
    {
    public:
        Lambertian(Vector3d albedo);
        Lambertian(Texture *texure);
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);

    private:
        std::shared_ptr<Texture> m_texture;
    };
    class Metal : public Material
    {
    public:
        Metal(Vector3d albedo, double fuzz);
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);

    private:
        std::shared_ptr<Texture> m_texture;
        double m_fuzz;
    };

    class Dielectric : public Material
    {

    public:
        Dielectric(double index);
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);

    private:
        double m_index;
    };

    class Scene
    {
    public:
        Scene(double max_distace);
        ~Scene();
        void add(Hitable *);
        bool hitOnce(Ray &ray, HitResult &out);
        Vector4d hit(Ray &ray);

    private:
        std::vector<Hitable *> items;
        double m_max_distance;
    };

    class Light : public Material
    {
    public:
        Light(Vector3d light);
        virtual Vector3d emitted(HitResult &hit);
    private:
        Vector3d m_light;
    };

    class Sphere : public Hitable, public Textureable
    {
    public:
        Sphere(Vector3d center, double radius, std::shared_ptr<Material>);
        Sphere(Vector3d center, Vector3d center2, double radius, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);
        void uv(Vector2d &uv, const Vector3d &point);

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
        void uv(Vector2d &uv, const Vector3d &point);

    private:
        std::shared_ptr<Material> m_mat;
        Vector3d m_point;
        Vector3d m_normal;
        Matrix3d m_uv_back;
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

} // namespace go
