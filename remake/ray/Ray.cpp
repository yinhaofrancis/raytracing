#include "Ray.hpp"
#include <iostream>
using namespace vr;

std::random_device Random::m_rd{};

const Random Random::shared{};

const Random::Sphere Ray::sphere(Random::shared.sphere());

Random::Generate vr::Random::normal(Scalar mean, Scalar stddev) const
{
    std::mt19937 m_gen(m_rd());
    std::normal_distribution dis{mean, stddev};
    Generate m([dis, m_gen]() mutable
               { return dis(m_gen); });
    return m;
}

vr::Random::Random()
{
}

Random::Generate vr::Random::uniform(Scalar min, Scalar max) const
{
    std::mt19937 m_gen(m_rd());
    std::uniform_real_distribution<Scalar> dis(min, max);
    Generate m([dis, m_gen]() mutable
               { return dis(m_gen); });
    return m;
}

Random::Sphere vr::Random::sphere() const
{
    std::mt19937 m_gen(m_rd());
    std::uniform_real_distribution<Scalar> dis(-1, 1);
    return Random::Sphere([dis, m_gen]() mutable
                          {
        Scalar x = dis(m_gen);
        Scalar y = dis(m_gen);
        Scalar z = dis(m_gen);
        return glm::normalize(Vec3(x,y,z)); });
}

Random::Sphere vr::Random::circle() const
{
    std::mt19937 m_gen(m_rd());
    std::uniform_real_distribution<Scalar> radius(0, 1);
    std::uniform_real_distribution<Scalar> theta(0, M_PI * 2);
    return Random::Sphere([radius, theta, m_gen]() mutable
                          {
        Scalar r = radius(m_gen);
        Scalar angle = theta(m_gen);
        Scalar x = r * cos(angle);
        Scalar y = r * sin(angle);
        return Vec3(x,y,0); });
}

Random::Sphere vr::Random::sphere(Scalar zmin, Scalar zmax)
{
    std::mt19937 m_gen(m_rd());
    std::uniform_real_distribution<Scalar> dis(-1, 1);
    std::uniform_real_distribution<Scalar> zdis(zmin, zmax);
    return Random::Sphere([zdis, dis, m_gen]() mutable
                          {
        Scalar x = dis(m_gen);
        Scalar y = dis(m_gen);
        Scalar z = zdis(m_gen);
        return glm::normalize(Vec3(x,y,z)); });
}

Ray vr::Ray::reflect(const Vec3 &position, const Vec3 &normal) const
{
    auto ndirection = glm::reflect(direction, normal);
    return {position, ndirection};
}

Ray vr::Ray::refract(const Vec3 &position, const Vec3 &normal, const Scalar eta) const
{
    auto ndirection = glm::refract(direction, normal, eta);
    return {position, ndirection};
}

Ray vr::Ray::hemisphere(const Vec3 &position, const Vec3 &normal) const
{
    return Ray{position, Ray::sphere(normal)};
}

Scalar vr::Random::Generate::operator()() const
{
    return m_call();
}

vr::Random::Generate::Generate(std::function<Scalar()> call) : m_call(call)
{
}

Vec3 vr::Random::Sphere::operator()() const
{
    return m_call();
}

Vec3 vr::Random::Sphere::operator()(const Vec3 &normal) const
{
    auto v = m_call();
    return glm::dot(v, normal) > 0 ? v : -v;
}

vr::Random::Sphere::Sphere(std::function<Vec3()> call) : m_call(call) {}

vr::Random::Sphere::Sphere(Sphere &&s) : m_call(s.m_call)
{
}

Vec3 vr::Sensor::operator()(Uint x, Uint y, Uint w, Uint h) const
{
    auto px = x * m_pixel_size + m_pixel_size / 2 + m_random() - Scalar(w) * m_pixel_size / 2;
    auto py = y * m_pixel_size + m_pixel_size / 2 + m_random() - Scalar(h) * m_pixel_size / 2;
    return Vec3(px, -py, 0);
}

vr::Camera::Camera(
    Uint width,
    Uint height,
    const Vec3 &location,
    const Vec3 &lookAt,
    const Vec3 &up,
    Scalar fov)
    : m_width(width),
      m_height(height),
      m_location(location),
      m_look_at(lookAt),
      m_fov(fov),
      m_front(0.001),
      m_up(up)
{
    auto tanv = tan(fov / 2);
    Vec3 z = glm::normalize(m_location - m_look_at);
    Vec3 x = glm::normalize(glm::cross(up, z));
    Vec3 y = glm::normalize(glm::cross(z, x));
    m_lookAt_matrix = Mat3(x, y, z);
    m_length = std::min(width, height) * m_front.pixelSize() / 2.0 / tanv;
}

Ray vr::Camera::ray(Uint x, Uint y) const
{
    Vec3 front_point = m_front(x, y, m_width, m_height);
    Vec3 back_point{0, 0, m_length};
    Vec3 direction = glm::normalize(m_lookAt_matrix * (front_point - back_point));
    return Ray(m_location, direction);
}

Uint vr::Camera::width() const
{
    return m_width;
}

Uint vr::Camera::height() const
{
    return m_height;
}

std::ostream &vr::operator<<(std::ostream &out, const Ray &ray)
{
    return out << "ray: " << glm::to_string(ray.location) << "|" << glm::to_string(ray.direction);
}

bool vr::AABB::hitTest(const Ray &ray, Scalar &t) const
{
    Scalar tminx = -INFINITY;
    auto bminx = this->hitAxisAlign(x.min, ray.location.x, ray.direction.x, tminx);

    Scalar tminy = -INFINITY;
    auto bminy = this->hitAxisAlign(y.min, ray.location.y, ray.direction.y, tminy);

    Scalar tminz = -INFINITY;
    auto bminz = this->hitAxisAlign(z.min, ray.location.z, ray.direction.z, tminz);
    Scalar tmin = std::max(std::max(tminx, tminy), tminz);

    if (bminx && bminy && bminz)
    {
        t = tmin;
        return true;
    }
    Vec3 pmin = ray.eval(tmin);

    if (bminx && y.contain(pmin.y) && z.contain(pmin.z))
    {
        t = tmin;
        return true;
    }
    if (bminy && x.contain(pmin.x) && z.contain(pmin.z))
    {
        t = tmin;
        return true;
    }
    if (bminz && x.contain(pmin.x) && y.contain(pmin.y))
    {
        t = tmin;
        return true;
    }
    return false;
}
