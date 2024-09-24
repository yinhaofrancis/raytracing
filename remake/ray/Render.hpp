#pragma once
#include "Util.hpp"
#include "Ray.hpp"
#include "Window.hpp"
#include "Geometry.hpp"
#include "Material.hpp"
#include <memory>
#include <atomic>

namespace vr
{
    
    struct Model
    {
        Primitives *primitives;
        Material *material;
    };
    struct Light
    {
        Vec3 pointOfLight;
        Color emmitLight;
        Vec3 lightNormal;
        Scalar areaOfLight;
        Scalar scalar;
    };
    class Bvh;
    class Scene
    {
    public:
        void add(const Model &model);
        void add(const ObjectReader& r, Material *material);
        void add(Primitives *primitives, Material *material);
        virtual bool hitTest(
            const Ray &ray,
            const bool includeLight,
            const Region &region,
            HitRecord &record,
            Material **material) const;
        void end();

        virtual std::vector<Light> lights(const Vec3 &position, const Region &region) const;

        virtual bool hitLight(const Ray &ray, const Region &region) const;

    private:
        std::vector<Model> m_primitives;
        std::shared_ptr<Bvh> m_bvh;
        std::vector<Model> m_lights;
    };
    class Render
    {
    public:
        Render(Uint count = 10);
        void render(const Camera &cam, const Scene &scene) const;
        void pixelShade(const vr::Camera &cam, vr::Uint j, vr::Uint i, const vr::Scene &scene, vr::Window &win) const;
        Scalar &gamma()
        {
            return m_gamma;
        }

    private:
        Color rayHit(const Ray &ray, const Scene &scene) const;
        Color rayRecursion(const Ray &ray, const HitRecord hr, const Scene &scene, const Material *material) const;
        Color specular(const Ray &ray, const vr::Scene &scene, const vr::HitRecord &hitRecord, const vr::Material *material) const;
        Color diffUse(Ray &ray, const HitRecord &hitRecord, const Scalar prr, const Scene &scene, const Material *material) const;
        const Uint m_count;
        Region m_region = {0.00001, 1000000};
        Scalar m_full = 8;
        Random::Generate m_gen = Random::shared.uniform(0.0, m_full);
        Scalar m_max = 7;
        std::atomic_int m_thread_count = 16;
        Scalar m_gamma = 1.0 / 2.2;
    };

    class Bvh 
    {
    public:
        std::vector<Bvh*> hit(const Ray& ray);
        Bvh(std::vector<Model> m);
        static std::shared_ptr<Bvh> make(std::vector<Model>& m);
        friend class Scene;
        void cut();
    private:
        std::vector<Model> m_leaf;
        std::shared_ptr<Bvh> m_right;
        std::shared_ptr<Bvh> m_left;
        AABB m_box;
        bool m_cut = false;
    };
} // namespace vr
