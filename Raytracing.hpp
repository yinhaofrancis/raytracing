#pragma once
#include <eigen3/Eigen/Eigen>

#include <eigen3/Eigen/Dense>

#include <random>
#include <memory>
#include <vector>
#include "Math.hpp"

#include "Texture.hpp"
using namespace Eigen;

namespace go
{
    class Material;
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

    class Scene
    {
    public:
        Scene(double max_distace,Vector4d& ambient);
        ~Scene();
        void add(Hitable *);
        bool hitOnce(Ray &ray, HitResult &out);
        Vector4d hit(Ray &ray);

    private:
        std::vector<Hitable *> items;
        double m_max_distance;
        Vector4d m_ambient;
    };



} // namespace go


// material
namespace go
{
    class Material
    {
    public:
        Material() {}
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);
        virtual Vector3d emitted(HitResult &hit);
    };

    class NormalColor : public Material{
    public:
        NormalColor();
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);
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
    class Light : public Material
    {
    public:
        Light(Vector3d light);
        virtual Vector3d emitted(HitResult &hit);

    private:
        Vector3d m_light;
    };
} // namespace go
