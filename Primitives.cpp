#include "Primitives.hpp"


bool go::Hitable::hit(Ray &ray, Interval ray_t, HitResult &result)
{
    return false;
}

bool &go::Hitable::isFrontFace()
{
    return m_front;
}

go::Hitable::~Hitable()
{
}

go::Sphere::Sphere(Vector3d center, double radius, std::shared_ptr<Material> mat) : m_center(center), m_radius(radius), m_mat(mat), m_center2(center)
{
}

go::Sphere::Sphere(Vector3d center, Vector3d center2, double radius, std::shared_ptr<Material> mat) : m_center(center), m_radius(radius), m_mat(mat), m_center2(center2)
{
}

bool go::Sphere::hit(Ray &ray, Interval ray_t, HitResult &result)
{
    double A = ray.direction().dot(ray.direction());
    Vector3d delta = this->center(ray) - ray.location();
    double B = delta.dot(ray.direction()) * -2.;
    double C = delta.dot(delta) - m_radius * m_radius;

    double d = B * B - 4 * A * C;
    if (d >= 0)
    {
        double x1 = (-B + sqrt(d)) / (2 * A);
        double x2 = (-B - sqrt(d)) / (2 * A);
        double t1 = std::min(x1, x2);
        if (!ray_t.surrounds(t1))
        {
            t1 = std::max(x1, x2);
            ;
            if (!ray_t.surrounds(t1))
            {
                return false;
            }
        }

        Vector3d h1 = ray.exec(t1);

        Vector3d n1 = (isFrontFace() ? 1 : -1) * (h1 - center(ray)).normalized();
        result.t = t1;
        result.normal = n1;
        result.hit = h1;
        result.mat = m_mat;
        this->uv(result.uv, result.hit);
        result.isFront = n1.dot(ray.direction()) <= 0;
        return true;
    }
    return false;
}

void go::Sphere::uv(Vector2d &uv, const Vector3d &point)
{
    auto theta = acos(-point.y());
    auto phi = atan2(-point.z(), point.x()) + pi;

    double u = phi / (2 * pi);
    double v = theta / pi;
    uv[0] = u;
    uv[1] = v;
}

Vector3d go::Sphere::center(Ray &t)
{
    return mix(m_center, m_center2, t.time());
}

go::Planer::Planer(Vector3d point, Vector3d normal, std::shared_ptr<Material> mat) : m_point(point), m_normal(normal), m_mat(mat)
{
    Vector3d tx = Vector3d(1.0, 0, 0);
    Vector3d y = tx.cross(m_normal);
    Vector3d x = m_normal.cross(y);
    Matrix3d tra;
    tra << x, y, m_normal;
    m_uv_back = tra.inverse();
}

bool go::Planer::hit(Ray &ray, Interval ray_t, HitResult &result)
{
    bool hit = ray.direction().dot(m_normal) < 0;
    if (hit)
    {
        double A = m_normal.dot(ray.direction());
        double B = (m_point - ray.location()).dot(m_normal);
        double t = B / A;
        result.t = t;
        result.hit = ray.exec(t);
        result.normal = m_normal;
        result.mat = m_mat;
        result.isFront = true;
        this->uv(result.uv, result.hit);
        if (ray_t.contains(t))
        {
            return true;
        }
    }
    return false;
}

void go::Planer::uv(Vector2d &uv, const Vector3d &point)
{
    Vector3d uvw = m_uv_back * point;

    uv[0] = abs(std::fmod(uvw[0], 1));
    uv[1] = abs(std::fmod(uvw[1], 1));
}