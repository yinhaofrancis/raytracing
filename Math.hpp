#pragma once
#include <eigen3/Eigen/Eigen>

#include <eigen3/Eigen/Dense>

#include <random>
#include <memory>
#include <vector>
#include <math.h>

using namespace Eigen;

const double infinity = std::numeric_limits<double>::infinity();

const double pi = 3.1415926535897932385;

const double epsilon = 1e-20;

Vector3d reflect(Vector3d &incident, Vector3d &normal);

Vector3d refract(Vector3d &incident, Vector3d &normal, double eta);

double random_double(double min, double max);

Vector3d mix(Vector3d start, Vector3d end, double s);

Vector3d random_in_unit_sphere();

Vector3d random_in_unit_semisphere(Vector3d &normal);

Vector3d random_in_unit_disk();

Vector3d random_in_triangle(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3);

Vector2d random_in_square(int sampleIndex, int sample);

const int max_depth = 8;

Vector4d gamma(Vector4d &v, double gamma);

Matrix4d translate(double x, double y, double z);

Matrix4d rotate(double x, double y, double z, double angle);

Matrix4d scale(double x, double y, double z);

Vector2d interpolateUV(const Vector3d &p, const Vector3d &v1, const Vector3d &v2, const Vector3d &v3, const Vector2d &uv1, const Vector2d &uv2, const Vector2d &uv3);

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

    class Random
    {
    public:
        Random(const Vector3d &w);
        Vector3d cosine_direction() const;
        Vector3d w() const;
    private:
        Vector3d m_w;
        Matrix3d m_mat;
    };
    class pdf
    {
    public:
        virtual ~pdf() {}
        virtual double value(const Vector3d &direction) const = 0;
        virtual Vector3d generate() const = 0;
    };
    class sphere_pdf:public pdf{
    public:
        ~sphere_pdf() {}
        virtual double value(const Vector3d &direction) const;
        virtual Vector3d generate() const;
    };
    class cosine_pdf:public pdf{
    public:
        ~cosine_pdf() {}
        cosine_pdf(const Vector3d& w);
        virtual double value(const Vector3d &direction) const;
        virtual Vector3d generate() const;
    private:
        Random m_r;

    };


} // namespace go
