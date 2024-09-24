#include "Render.hpp"
#include <iostream>
#include <thread>
#include <stack>
using namespace vr;

void vr::Scene::add(const Model &model)
{
    if (model.material->isLight())
    {
        if (model.material->volume())
        {
            auto m = HitMode(vr::front | vr::back);
            auto p = model.primitives;
            p->setHitMode(m);
        }
        m_lights.push_back(model);
        m_primitives.push_back(model);
    }
    else
    {
        m_primitives.push_back(model);
    }
}

void vr::Scene::add(const ObjectReader &r, Material *material)
{
    auto list = r.tiangles();
    for (auto &&i : list)
    {
        add(i.get(), material);
    }
}

void vr::Scene::add(Primitives *primitives, Material *material)
{
    Model m = {primitives, material};
    add(m);
}

bool vr::Scene::hitTest(const Ray &ray, const bool includeLight, const Region &region, HitRecord &record, Material **material) const
{
    auto boxs = m_bvh->hit(ray);
    Region local = region;
    bool hit = false;
    for (auto &&b : boxs)
    {
        for (auto i : b->m_leaf)
        {
            if (i.material->isLight() && !includeLight)
            {
                continue;
            }
            if (i.primitives->hitTest(ray, local, record))
            {
                local.max = record.t;
                hit = true;
                *material = i.material;
            }
        }
    }

    return hit;
}

void vr::Scene::end()
{
    m_bvh = Bvh::make(m_primitives);
    m_primitives.clear();
}

std::vector<Light> vr::Scene::lights(const Vec3 &position, const Region &region) const
{
    std::vector<Light> lights;
    for (auto &&light : m_lights)
    {
        Light l;
        auto pt = light.primitives->random();
        auto nor = light.primitives->lightNormal(position - pt);
        l.areaOfLight = light.primitives->area();
        l.emmitLight = light.material->albedo({0, 0});
        l.lightNormal = nor;
        l.pointOfLight = pt;
        l.scalar = light.material->scale();
        Ray r(position, glm::normalize(pt - position));
        HitRecord hr;
        Material *m;
        if (this->hitTest(r, true, region, hr, &m))
        {
            if (m->isLight())
            {
                lights.push_back(l);
            }
        }
    }
    return lights;
}

bool vr::Scene::hitLight(const Ray &ray, const Region &region) const
{
    HitRecord hr;
    Material *m;
    if (this->hitTest(ray, true, region, hr, &m))
    {
        if (m->isLight())
        {
            return true;
        }
    }
    return false;
}

vr::Render::Render(Uint count) : m_count(count)
{
}

void vr::Render::render(const Camera &cam, const Scene &scene) const
{
    Uint w = cam.width();
    Uint h = cam.height();
    std::atomic_int m = 0;
    Window win(w, h);
    Uint c = ceil(h / m_thread_count) + 1;
    for (Uint s = 0; s < m_thread_count; s++)
    {
        std::thread(
            [s, this, c, &cam, &scene, &win, h, w, &m]()
            {
                for (Uint i = s * c; i < s * c + c; i++)
                {
                    if (i < h)
                    {
                        for (Uint j = 0; j < w; j++)
                        {
                            pixelShade(cam, j, i, scene, win);
                        }
                    }
                }
                m++;
            })
            .detach();
    }
    while (m < m_thread_count)
    {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        win.update();
    }

    win.update();
    win.write("temp.ppm");
    win.wait();
}

void vr::Render::pixelShade(const vr::Camera &cam, vr::Uint j, vr::Uint i, const vr::Scene &scene, vr::Window &win) const
{
    Int count = 0;
    Color sum;
    while (count < m_count)
    {
        Ray r = cam.ray(j, i);
        sum += rayHit(r, scene);
        count++;
    }
    sum /= m_count;
    sum = glm::pow(sum, Vec3(m_gamma));
    win.pixel(sum, j, i);
}

vr::Color vr::Render::rayHit(const Ray &ray, const Scene &scene) const
{
    HitRecord hitRecord;
    Material *material;
    Ray newRay = ray;
    Scalar m_rr = m_max / m_full;
    if (scene.hitTest(newRay, true, m_region, hitRecord, &material))
    {
        if (material->isLight())
        {
            return Vec3(1);
        }
        else
        {
            return rayRecursion(ray, hitRecord, scene, material);
        }
    }
    return Vec3{0};
}

Color vr::Render::rayRecursion(const Ray &ray, const HitRecord hitRecord, const Scene &scene, const Material *material) const
{
    Scalar m_rr = m_max / m_full;

    Color scolor = material->hasDirectLight() ? specular(ray, scene, hitRecord, material) : Vec3(0);
    Ray newRay = ray;
    Color dcolor = diffUse(newRay, hitRecord, m_rr, scene, material);
    if (m_gen() < m_max)
    {
        HitRecord newHitRecord;
        Material *newMaterial;
        if (scene.hitTest(newRay, true, m_region, newHitRecord, &newMaterial))
        {
            if (scene.hitLight(newRay, m_region))
            {
                if (!material->hasDirectLight())
                {
                    return material->albedo({0, 0});
                }
            }
            newRay.count += 1;
            return scolor + dcolor * rayRecursion(newRay, newHitRecord, scene, newMaterial);
        }
        else
        {
            return scolor;
        }
    }
    else
    {
        return scolor;
    }
}

Color vr::Render::specular(const Ray &ray, const vr::Scene &scene, const vr::HitRecord &hitRecord, const vr::Material *material) const
{

    Color specularColor = Vec3{0};

    auto lights = scene.lights(hitRecord.position, m_region);

    for (auto &&li : lights)
    {
        specularColor += material->lightDirect(ray, hitRecord.uv, li.emmitLight, li.pointOfLight, hitRecord.position, hitRecord.normal, li.lightNormal, li.areaOfLight * li.scalar);
    }
    return specularColor;
}

Color vr::Render::diffUse(Ray &ray, const HitRecord &hitRecord, Scalar prr, const Scene &scene, const Material *material) const
{
    ray = material->nextRay(ray, hitRecord.position, hitRecord.normal);
    return material->lightInDirect(hitRecord.uv, ray.direction, hitRecord.normal, prr);
}

std::vector<Bvh *> vr::Bvh::hit(const Ray &ray)
{
    Scalar t;
    if (m_box.hitTest(ray, t))
    {
        if (m_left == nullptr && m_right == nullptr)
        {
            return {this};
        }
        else
        {
            auto a = m_left->hit(ray);
            auto b = m_right->hit(ray);
            a.insert(a.end(), b.begin(), b.end());
            return a;
        }
    }
    return {};
}

vr::Bvh::Bvh(std::vector<Model> m) : m_leaf(m)
{
    for (auto &&i : m)
    {
        m_box.add(i.primitives->bbox());
    }
}
std::shared_ptr<Bvh> vr::Bvh::make(std::vector<Model> &m)
{
    std::stack<std::shared_ptr<Bvh>> stack;
    std::shared_ptr<Bvh> mm(new Bvh(m));
    stack.push(mm);
    while (stack.size() > 0)
    {
        if (stack.top()->m_cut == false || stack.top()->m_leaf.size() > 20)
        {
            stack.top()->cut();
            auto last = stack.top();
            stack.pop();
            if(last->m_left != nullptr){
                stack.push(last->m_left);
            }
            if(last->m_right != nullptr){
                stack.push(last->m_right);
            }
        }
        else
        {
            stack.pop();
        }
    }

    return mm;
}
void vr::Bvh::cut()
{
    AABB left, right;
    std::vector<Model> leftData, rightData;
    int c = 0, b = 0;
    AABB current_box = m_box;
    current_box.cut(left, right);

    for (auto &&i : m_leaf)
    {
        if (left.intersection(i.primitives->bbox()) && c < m_leaf.size() / 2)
        {
            leftData.push_back(i);
            c++;
        }
        else
        {
            rightData.push_back(i);
            b++;
        }
    }
    if (c != 0 && b != 0)
    {
        m_left = std::shared_ptr<Bvh>(new Bvh(leftData));
        m_right = std::shared_ptr<Bvh>(new Bvh(rightData));
        m_leaf.clear();
        m_leaf.resize(0);
    }
    m_cut = true;
}
