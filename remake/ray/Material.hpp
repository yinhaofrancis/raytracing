#pragma once
#include "Util.hpp"
#include "Ray.hpp"

namespace vr
{
    class Material
    {
    public:
        virtual Color albedo(const Vec2 &uvInterp) const = 0;

        virtual bool isLight() const = 0;

        virtual bool hasDirectLight() const = 0;

        virtual Color lightDirect(
            const Ray &current,
            const Vec2 &uv,
            const Color &emitLight,
            const Vec3 &pointOnLight,
            const Vec3 &hitPosition,
            const Vec3 &hitNormal,
            const Vec3 &lightNormal,
            const Scalar Area) const = 0;

        virtual Color lightInDirect(
            const Vec2 &uv,
            const Vec3 &rayDirection,
            const Vec3 &hitNormal,
            const Scalar prr) const = 0;

        virtual Ray nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const = 0;

        virtual Scalar pdf() const
        {
            return 1.0 / 2.0 / M_PI;
        }

        virtual Scalar pdf(Scalar area) const
        {
            return 1.0 / area;
        }

        virtual Scalar scale() const
        {
            return 1;
        }
        virtual Color normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const = 0;

        virtual bool volume() const { return false; }

    protected:
        static Random::Generate s_random;
    };

    class Lambertain : public Material
    {
    public:
        Lambertain(const Color &color);

        virtual Color albedo(const Vec2 &uvInterp) const;

        virtual Color normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const;

        virtual Ray nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const;

        virtual bool isLight() const { return false; }

        virtual bool hasDirectLight() const { return true; }

        virtual Color lightDirect(
            const Ray &current,
            const Vec2 &uv,
            const Color &emitLight,
            const Vec3 &pointOnLight,
            const Vec3 &hitPosition,
            const Vec3 &hitNormal,
            const Vec3 &lightNormal,
            const Scalar Area) const;

        virtual Color lightInDirect(
            const Vec2 &uv,
            const Vec3 &rayDirection,
            const Vec3 &hitNormal,
            const Scalar prr) const;

    private:
        Color m_color;
    };

    class EmitLight : public Material
    {
    public:
        EmitLight(const Color &color, const Scalar scale = 1);

        virtual Color albedo(const Vec2 &uvInterp) const;

        virtual Color normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const;

        virtual Ray nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const;

        virtual bool isLight() const { return true; }

        virtual bool hasDirectLight() const { return false; }

        virtual Color lightInDirect(
            const Vec2 &uv,
            const Vec3 &rayDirection,
            const Vec3 &hitNormal,
            const Scalar prr) const;
        virtual Color lightDirect(
            const Ray &current,
            const Vec2 &uv,
            const Color &emitLight,
            const Vec3 &pointOnLight,
            const Vec3 &hitPosition,
            const Vec3 &hitNormal,
            const Vec3 &lightNormal,
            const Scalar Area) const { return Vec3{0, 0, 0}; }

        Scalar scale() const
        {
            return m_scale;
        }

    private:
        Color m_color;
        Scalar m_scale = 1;
    };

    class Metal : public Material
    {
    public:
        Metal(const Color &color, const Scalar blur);

        virtual Color albedo(const Vec2 &uvInterp) const;

        virtual Color normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const;

        virtual Ray nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const;

        virtual bool isLight() const { return false; }

        virtual bool hasDirectLight() const { return false; }

        virtual Color lightInDirect(
            const Vec2 &uv,
            const Vec3 &rayDirection,
            const Vec3 &hitNormal,
            const Scalar prr) const;

        virtual Color lightDirect(
            const Ray &current,
            const Vec2 &uv,
            const Color &emitLight,
            const Vec3 &pointOnLight,
            const Vec3 &hitPosition,
            const Vec3 &hitNormal,
            const Vec3 &lightNormal,
            const Scalar Area) const
        {
            return {0, 0, 0};
        };

    private:
        Color m_color;
        Scalar m_blur;
    };

    class Dielectric : public Material
    {
    public:
        Dielectric(Scalar eta);

         virtual Color albedo(const Vec2 &uvInterp) const;

        virtual Color normal(const Vec2 &uvInerp, const Vec3 &plantNormal) const;

        virtual Ray nextRay(const Ray &ray, const Vec3 &position, const Vec3 &normal) const;

        virtual bool isLight() const { return false; }

        virtual bool hasDirectLight() const { return false; }

        virtual Color lightInDirect(
            const Vec2 &uv,
            const Vec3 &rayDirection,
            const Vec3 &hitNormal,
            const Scalar prr) const;

        virtual Color lightDirect(
            const Ray &current,
            const Vec2 &uv,
            const Color &emitLight,
            const Vec3 &pointOnLight,
            const Vec3 &hitPosition,
            const Vec3 &hitNormal,
            const Vec3 &lightNormal,
            const Scalar Area) const
        {
            return {0, 0, 0};
        };
    private:
        Scalar m_eta;
    };
} // namespace vr
