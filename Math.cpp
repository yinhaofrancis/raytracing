
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

double crossProductZ(const Vector3d &a, const Vector3d &b, const Vector3d &c)
{
    return (b.x() - a.x()) * (c.z() - a.z()) - (b.z() - a.z()) * (c.x() - a.x());
}
Vector2d interpolateUV(const Vector3d &p, const Vector3d &v1, const Vector3d &v2, const Vector3d &v3,
                       const Vector2d &uv1, const Vector2d &uv2, const Vector2d &uv3)
{
    // 判断点P是否在以v1v2为底的三角形内
    double s = crossProductZ(p, v1, v2) / crossProductZ(v3, v1, v2);
    // 判断点P是否在以v2v3为底的三角形内
    double t = crossProductZ(p, v2, v3) / crossProductZ(v1, v2, v3);

    // 确保s和t都在[0, 1]之间，即点P确实在三角形内部
    Vector2d uv;
    uv[0] = uv1.x() * (1 - s) * (1 - t) + uv2.x() * s * (1 - t) + uv3.x() * t;
    uv[1] = uv1.y() * (1 - s) * (1 - t) + uv2.y() * s * (1 - t) + uv3.y() * t;
    return uv;
}

go::Interval::Interval(double min, double max) : m_min(min), m_max(max)
{
}

go::Interval::Interval() : Interval(+infinity, -infinity)
{
}

go::Interval go::Interval::max(double max)
{
    return Interval(0.000001, max);
}

double go::Interval::size() const
{
    return m_max - m_min;
}

bool go::Interval::contains(double x) const
{
    return m_min <= x && m_max >= x;
}

bool go::Interval::surrounds(double x) const
{
    return m_min < x && m_max > x;
    ;
}

go::Interval go::Interval::expands(double x) const
{
    Interval m(*this);
    m.m_min -= x;
    m.m_max += x;
    return m;
}