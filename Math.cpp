
#include "Math.hpp"

double random_double(double min, double max)
{
    static std::uniform_real_distribution<double> distribution(min, max);
    static std::mt19937 generator;
    return distribution(generator);
}

Vector3d mix(Vector3d start, Vector3d end, double s)
{
    Vector3d out = (end - start) * s;
    return start + out;
}

Vector3d random_in_unit_sphere()
{
    double x = random_double(-1, 1);
    double y = random_double(-1, 1);
    double z = random_double(-1, 1);
    Vector3d d(x, y, z);
    if (d.size() > 1)
    {
        d.normalize();
    }
    return d;
}

Vector3d random_in_unit_semisphere(Vector3d &normal)
{
    Vector3d v = random_in_unit_sphere();
    return v.dot(normal) > 0 ? v : -v;
}

Vector3d random_in_unit_disk()
{

    Vector3d v(random_double(-1, 1), random_double(1, 1), 0);
    if (v.size() > 1)
    {
        v.normalize();
        return v;
    }
    return v;
}

Vector3d reflect(Vector3d &incident, Vector3d &normal)
{
    return incident - 2 * incident.dot(normal) * normal;
}

Vector3d refract(Vector3d &incident, Vector3d &normal, double eta)
{

    auto cos_theta = std::min(normal.dot(-incident), 1.0);
    Vector3d r_out_perp = (incident + cos_theta * normal) * eta;
    Vector3d r_out_parallel = -std::sqrt(std::fabs(1 - r_out_perp.size())) * normal;
    return r_out_parallel + r_out_perp;
}

Vector4d gamma(Vector4d &v, double gamma)
{
    return Vector4d(std::pow(v.x(), gamma), std::pow(v.y(), gamma), std::pow(v.z(), gamma), std::pow(v.w(), gamma));
}