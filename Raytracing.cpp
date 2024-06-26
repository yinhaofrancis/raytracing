#include "Raytracing.hpp"
#include <math.h>
#include <algorithm>
#include <iostream>
#include "Math.hpp"

go::Ray::Ray()
{
}
go::Ray::Ray(Vector3d location, Vector3d direction) : m_location(location), m_direction(direction), m_time(random_double(0, 1))
{
}

void go::Ray::reflect(Vector3d normal, Vector3d hitPoint, double fuzz)
{
    m_location = hitPoint;
    if (fuzz < epsilon)
    {
        m_direction = ::reflect(m_direction, normal).normalized();
    }
    else
    {
        m_direction = ::reflect(m_direction, normal).normalized() + random_in_unit_sphere() * fuzz;
    }

    this->m_max -= 1;
}

void go::Ray::refract(Vector3d normal, Vector3d hitPoint, double fract)
{
    m_location = hitPoint;
    m_direction = ::refract(m_direction, normal, fract);
    this->m_max -= 1;
}

void go::Ray::hitRandom(Vector3d normal, Vector3d hitPoint)
{
    m_direction = random_in_unit_semisphere(normal);
    m_location = hitPoint;
    this->m_max -= 1;
}
void go::Ray::fullRandom(Vector3d hitPoint)
{
    m_direction = random_in_unit_sphere();
    m_location = hitPoint;
    this->m_max -= 1;
}
void go::Ray::lambertian(Vector3d normal, Vector3d hitPoint)
{
    
    m_direction = Random(normal).cosine_direction();
    if (m_direction.dot(m_direction) < epsilon)
    {
        m_direction = normal;
    }
    m_location = hitPoint;
    this->m_max -= 1;
}

Vector3d go::Ray::location() const
{
    return m_location;
}

Vector3d go::Ray::direction() const
{
    return m_direction;
}

Vector3d go::Ray::exec(double t)
{
    return m_location + t * m_direction;
}

void go::Ray::setDirection(const Vector3d &d)
{
    m_direction = d;
    m_direction.normalize();
}

int go::Ray::depth()
{
    return m_max;
}
double go::Ray::time()
{
    return m_time;
}
go::Camera::Camera(Vector3d location, Vector3d target, Vector3d up, double fov, double ratio, double length, double defocus)
    : m_location(location),
      m_ratio(ratio),
      m_fov(fov),
      m_length(length)
{
    auto z = (location - target).normalized();
    auto x = up.cross(z).normalized();
    auto y = z.cross(x).normalized();
    m_disk_radius = length * tan(defocus / 2);
    m_lookAt << x, y, z;
}

go::Ray go::Camera::createRay(Vector2d &screenUV)
{
    double length = m_length * tan(m_fov / 2);
    Vector3d offset = random_in_unit_disk() * m_disk_radius;
    Vector3d filmCoo(
        screenUV.x() * m_ratio * length,
        screenUV.y() * length, -m_length);
    Vector3d rd = filmCoo - offset;
    filmCoo = (m_lookAt * rd).normalized();
    offset = (m_lookAt * offset);

    return go::Ray(m_location + offset, filmCoo);
}

bool go::Hitable::hit(Ray &ray, Interval ray_t, HitResult &result)
{
    return false;
}

bool &go::Hitable::isFrontFace()
{
    return m_front;
}

go::LightStatus go::Hitable::get_light_status(const HitResult &result)
{
    return LightStatus();
}

Vector3d go::Hitable::random_serface(const Vector3d &point)
{
    return Vector3d();
}

go::Hitable::~Hitable()
{
}



go::Scene::Scene(double distance,Vector4d& m) : m_max_distance(distance),m_ambient(m)
{
}

go::Scene::~Scene()
{
    items.clear();
    for (auto i : items)
    {
        delete i;
    }
}

void go::Scene::add(Hitable *i)
{
    items.push_back(i);
}

void go::Scene::light(Hitable *i)
{
    lights.push_back(i);
    add(i);
}

bool go::Scene::hitOnce(Ray &ray, HitResult &out,LightStatus& status)
{
    Interval m = Interval::max(m_max_distance);
    bool hit_some = false;
    HitResult result;
    for (auto &&i : items)
    {
        if (i->hit(ray, m, result))
        {
            m = Interval::max(result.t);
            hit_some = true;
            out = result;
        }
    }
    return hit_some;
}

Vector4d go::Scene::hit(Ray &ray)
{
    if (ray.depth() > 0)
    {
        HitResult result;
        LightStatus light;
        if (hitOnce(ray, result,light))
        {
            Ray next;
            Vector4d color;
            auto em = result.mat->emitted(result);
            auto m = result.mat->scatter(ray, color, result, next);
            auto scatter_pdf = result.mat->scatter_pdf(ray,result,next);
            Vector4d c;
            c << em,1.0;
            if(!m){
                return c;
            }
            
            auto status = lights.front()->get_light_status(result);
            status.light_area = status.light_area ;
            bool drawLight = false;
            cosine_pdf surface_pdf(result.normal);
            double pdf_val = surface_pdf.value(surface_pdf.generate());
            if(hasLight()){
                auto toL = status.toLight(result);
                toL.normalize();
                if(result.normal.dot(toL) > 0){
                    next.setDirection(toL);
                    drawLight = true;
                }
            }
            auto pdf = !drawLight ? 1 / (2 * pi) :status.pdf(result) * pdf_val;
            auto nextColor = hit(next);
            return c + (color*scatter_pdf).cwiseProduct(nextColor) / pdf;
        }
        return m_ambient;
    }

    return Vector4d(0, 0, 0, 1);
}

bool go::Scene::hasLight()
{
    return lights.size() > 0;
}

bool go::Material::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    return false;
}

Vector3d go::Material::emitted(HitResult &hit)
{
    return Vector3d(0,0,0);
}

double go::Material::scatter_pdf(const Ray &in, const HitResult &hit, const Ray &out)
{
    return 1.0;
}

go::Lambertian::Lambertian(Vector3d albedo) : m_texture(std::make_shared<Color>(albedo))
{
}

go::Lambertian::Lambertian(Texture *texure) : m_texture(texure)
{
}

bool go::Lambertian::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    // std::cout << "Lambertian" << std::endl;
    out = in;
    auto tColor = m_texture->color(hit.uv, hit.hit);
    color << tColor, 1.0;
    out.lambertian(hit.normal, hit.hit);
    return true;
}

double go::Lambertian::scatter_pdf(const Ray &in, const HitResult &hit, const Ray &out)
{
    return 1 / (2 * pi);
}

go::Metal::Metal(Vector3d albedo, double fuzz) : m_texture(std::make_shared<Color>(albedo)), m_fuzz(fuzz)
{
}

bool go::Metal::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    out = in;

    out.reflect(hit.normal, hit.hit, m_fuzz);
    auto tColor = m_texture->color(hit.uv, hit.hit);
    color << tColor, 1.0;
    return true;
}

go::Dielectric::Dielectric(double index) : m_index(index)
{
}

static double reflectance(double cosine, double refraction_index)
{
    // Use Schlick's approximation for reflectance.
    auto r0 = (1 - refraction_index) / (1 + refraction_index);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}

bool go::Dielectric::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{

    double index = !hit.isFront ? m_index : 1 / m_index;
    Vector3d n = hit.isFront ? hit.normal : -hit.normal;
    double ct = (-in.direction()).dot(n);
    double cos_theta = std::min(ct, 1.0);
    double sin_theta = sqrt(1 - cos_theta * cos_theta);

    out = in;
    bool ro = reflectance(cos_theta, index) > random_double(0, 1);
    if (index * sin_theta > 1 || ro)
    {
        out.reflect(n, hit.hit, 0);
    }
    else
    {
        out.refract(n, hit.hit, index);
    }
    color = Vector4d(1, 1, 1, 1);

    return true;
}

go::Light::Light(Vector3d light):m_light(light)
{
}

Vector3d go::Light::emitted(HitResult &hit)
{
    // std::cout << "Light" << std::endl;
    // if(!hit.isFront){
    //     return Vector3d(0,0,0);
    // }
    return m_light;
}

go::NormalColor::NormalColor()
{
}

bool go::NormalColor::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    out = in;
    out.lambertian(hit.normal, hit.hit);
    Vector4d m = hit.uv.homogeneous().homogeneous();
    m[3] = 1;
    color = m;
    
    return true;
}

go::isotropic::isotropic(Vector3d albedo):m_texture(std::make_shared<Color>(albedo))
{

}

go::isotropic::isotropic(std::shared_ptr<Texture> tex)
{
}

bool go::isotropic::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    Vector3d m = m_texture->color(hit.uv,hit.hit);
    color = m.homogeneous();
    color[3] = 1;
    out.fullRandom(hit.hit);
    return true;
}

double go::isotropic::scatter_pdf(const Ray &in, const HitResult &hit, const Ray &out)
{
    return 1 / (4 * pi);
}

double go::LightStatus::pdf(const go::HitResult &hit)
{
    auto to_light = toLight(hit);
    if (to_light.normalized().dot(hit.normal) < 0){
        return -1;
    }
    auto light_cosine = fabs(to_light.y());
    if (light_cosine < 0.000001)
        return -1;
    auto length_squared = to_light.dot(to_light);
    auto pdf = length_squared / (light_cosine * light_area);
    return pdf;
}

Vector3d go::LightStatus::toLight(const go::HitResult & hit)
{
    return on_light - hit.hit;
}
