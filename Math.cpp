
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

Vector2d random_in_square(int i, int sample)
{
    double w = sqrt(sample);
    double y = i / w + random_double(-0.5,0.5);
    double x = fmod(i,w) + random_double(-0.5,0.5);
    Vector2d r(x,y);
    return r / w;
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

Matrix4d translate(double x,double y,double z)
{
    Matrix4d m;
    m << 1,0,0,x,
         0,1,0,y,
         0,0,1,z,
         0,0,0,1;
    return m;
}

Matrix4d rotate(double x,double y,double z,double angle)
{
    AngleAxisd ax(angle,Vector3d(x,y,z).normalized());
    Eigen::Quaterniond m(ax);
    m.normalize();
    Matrix4d m4 = Matrix4d::Identity();
    m4.block<3,3>(0,0) = m.toRotationMatrix();
    return m4;
}

Matrix4d scale(double x,double y,double z)
{
    Matrix4d m4;
    m4 << x,0,0,0,
          0,y,0,0,
          0,0,z,0,
          0,0,0,1;
    return m4;
}

Vector2d interpolateUV(const Vector3d &p, const Vector3d &a, const Vector3d &b, const Vector3d &c,
                       const Vector2d &uv1, const Vector2d &uv2, const Vector2d &uv3)
{
    double at = -(p.x() - b.x()) * (c.y() - b.y()) + (p.y() - b.y()) * (c.x() - b.x());
    double ab = -(a.x() - b.x()) * (c.y() - b.y()) + (a.y() - b.y()) * (c.x() - b.x());
    double bt = -(p.x() - c.x()) * (a.y() - c.y()) + (p.y() - c.y()) * (a.x() - c.x());
    double bb = -(b.x() - c.x()) * (a.y() - c.y()) + (b.y() - c.y()) * (a.x() - c.x());
    double alpha = at / ab;
    double beta = bt / bb;
    double gama = 1 - alpha - beta;
    return uv1 * alpha + uv2 * beta + uv3 * gama;
}

go::Interval::Interval(double min, double max) : m_min(min), m_max(max)
{
}

go::Interval::Interval() : Interval(+infinity, -infinity)
{
}

go::Interval go::Interval::max(double max)
{
    return Interval(0.00001, max);
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

go::Random::Random(const Vector3d &w)
{
    auto s = w.normalized();
    auto a = s.x() > 0.9 ? Vector3d(0,1,0) : Vector3d(1,0,0);
    auto v = w.cross(a);
    auto u = w.cross(v);
    m_mat << u,v,w;
}

Vector3d go::Random::cosine_direction()
{
    auto r1 = random_double(0,1);
    auto r2 = random_double(0,1);
    auto phi = 2 * pi * r1;
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);
    auto z = sqrt(1 - r2);
    return Vector3d(x,y,z);
}
