#pragma once
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>

#include <random>
#include <memory>
#include <vector>

using namespace Eigen;

const double infinity = std::numeric_limits<double>::infinity();

const double pi = 3.1415926535897932385;

const double epsilon = 1e-20;

Vector3d reflect(Vector3d &incident, Vector3d &normal);

double random_double(double min, double max);

Vector3d random_in_unit_sphere();

Vector3d random_in_unit_semisphere(Vector3d &normal);

Vector3d random_in_unit_disk();

const int max_depth = 20;

Vector4d gamma(Vector4d &v, double gamma);

namespace go
{
    class Interval
    {
    public:
        Interval(double min, double max);
        Interval();
        static Interval max(double max);
        double size() const;
        bool contains(double x) const;
        bool surrounds(double x) const;
        Interval expands(double x) const;

    private:
        double m_min;
        double m_max;
    };
    class Ray
    {
    public:
        Ray();
        Ray(Vector3d location, Vector3d direction);
        void lambertian(Vector3d normal, Vector3d hitPoint);
        void reflect(Vector3d normal, Vector3d hitPoint,double fuzz);
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
        Camera(Vector3d location, Vector3d target, Vector3d up, double fov,double ratio,double length,double defocus);
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
    };

    struct HitResult
    {
        Vector3d hit, normal;
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
        Lambertian(Vector4d albedo);
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);

    private:
        Vector4d m_albedo;
    };
    class Metal : public Material
    {
    public:
        Metal(Vector4d albedo,double fuzz);
        virtual bool scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out);
    private:
        Vector4d m_albedo;
        double m_fuzz;
    };


    class Dielectric:public Material{
    
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

    class Sphere : public Hitable
    {
    public:
        Sphere(Vector3d center, double radius, std::shared_ptr<Material>);
        Sphere(Vector3d center,Vector3d center2, double radius, std::shared_ptr<Material>);
        virtual bool hit(Ray &ray, Interval ray_t, HitResult &result);

    private:
        Vector3d center(Ray&);
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

    private:
        std::shared_ptr<Material> m_mat;
        Vector3d m_point;
        Vector3d m_normal;
    };

} // namespace go