#include "Material.hpp"
#include <iostream>
using namespace vr;

Random::Generate Material::s_random = Random::shared.uniform(0, 1);

vr::Lambertain::Lambertain(const Color &color) : m_color(color)
{
}

Color vr::Lambertain::albedo(const Vec2 &uvInterp) const
{
    return m_color * 1.0 / M_PI;
}

Color vr::Lambertain::normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const
{
    return Color(0, 0, 1);
}

Ray vr::Lambertain::nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const
{
    auto r = ray.hemisphere(position, normal);
    r.direction *= s_random();
    return r;
}

Color vr::Lambertain::lightDirect(
    const Ray &current,
    const Vec2 &uv,
    const Color &emitLight,
    const Vec3 &pointOnLight,
    const Vec3 &hitPosition,
    const Vec3 &hitNormal,
    const Vec3 &lightNormal,
    const Scalar Area) const
{
    Scalar distance = glm::distance(hitPosition, pointOnLight);
    Vec3 lightDirection = glm::normalize(pointOnLight - hitPosition);
    Scalar hitPointCos = std::max(glm::dot(lightDirection, hitNormal), 0.0);
    Scalar lightCos = std::max(glm::dot(lightNormal, -lightDirection), 0.0);
    Scalar lightPdf = 1 / Area;
    auto res = emitLight * albedo(uv) * lightCos * hitPointCos / (distance * distance) / lightPdf;
    return res;
}

Color vr::Lambertain::lightInDirect(const Vec2 &uv, const Vec3 &rayDirection, const Vec3 &hitNormal, const Scalar prr) const
{
    Scalar param = glm::dot(rayDirection, hitNormal) / pdf() / prr;
    return albedo(uv) * param;
}

vr::Metal::Metal(const Color &color, const Scalar blur) : m_color(color), m_blur(blur)
{
}

Color vr::Metal::albedo(const Vec2 &uvInterp) const
{
    return m_color;
}

Color vr::Metal::normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const
{
    return Color(0, 0, 1);
}

Ray vr::Metal::nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const
{
    auto r = ray.reflect(position, normal);
    if (m_blur > 0.00001)
    {
        Vec3 delta = Vec3(s_random() * 2 - 1, s_random() * 2 - 1, s_random() * 2 - 1) * m_blur;
        r.direction += delta;
        r.direction = glm::normalize(r.direction);
    }
    return r;
}

Color vr::Metal::lightInDirect(const Vec2 &uv, const Vec3 &rayDirection, const Vec3 &hitNormal, const Scalar prr) const
{
    return albedo(uv);
}

vr::EmitLight::EmitLight(const Color &color, const Scalar scale) : m_color(color), m_scale(scale)
{
}

Color vr::EmitLight::albedo(const Vec2 &uvInterp) const
{
    return m_color;
}

Color vr::EmitLight::normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const
{
    return {0, 0, 1};
}

Ray vr::EmitLight::nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const
{
    return Ray();
}

Color vr::EmitLight::lightInDirect(const Vec2 &uv, const Vec3 &rayDirection, const Vec3 &hitNormal, const Scalar prr) const
{
    return Vec3{1};
}

vr::Dielectric::Dielectric(Scalar eta):m_eta(eta)
{
}

Color vr::Dielectric::albedo(const Vec2 &uvInterp) const
{
    return Color(1);
}

Color vr::Dielectric::normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const
{
    return Color(0,0,1);
}

Ray vr::Dielectric::nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const
{
    auto direct = glm::dot(normal,ray.direction);
    Scalar eta = direct < 0 ? 1.0 / m_eta : m_eta;
    auto n = direct < 0 ? normal : -normal;
    auto cos = glm::dot(n,-ray.direction);
    auto sin = std::sqrt(1 - cos * cos);
    bool notRefract = eta * sin > 1.0;
    auto r0 = (1 - eta) / (1 + eta);
    r0 = r0*r0;
    auto reflection = r0 + (1-r0)*std::pow((1 - cos),5);
    if(notRefract || reflection > s_random()){
        return ray.reflect(position,n);
    }else{
        return ray.refract(position,n,eta);
    }
    
}

Color vr::Dielectric::lightInDirect(const Vec2 &uv, const Vec3 &rayDirection, const Vec3 &hitNormal, const Scalar prr) const
{
    return Color(1);
}
