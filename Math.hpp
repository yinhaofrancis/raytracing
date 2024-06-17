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

Vector3d refract(Vector3d &incident, Vector3d &normal, double eta);

double random_double(double min, double max);

Vector3d mix(Vector3d start, Vector3d end, double s);

Vector3d random_in_unit_sphere();

Vector3d random_in_unit_semisphere(Vector3d &normal);

Vector3d random_in_unit_disk();

const int max_depth = 20;

Vector4d gamma(Vector4d &v, double gamma);

double crossProductZ(const Vector3d& a, const Vector3d& b, const Vector3d& c);

Vector2d interpolateUV(const Vector3d& p, const Vector3d& v1, const Vector3d& v2, const Vector3d& v3,
                      const Vector2d& uv1, const Vector2d& uv2, const Vector2d& uv3);

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
} // namespace go
