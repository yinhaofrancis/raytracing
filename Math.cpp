#include "Math.hpp"
#include <math.h>
#include <iostream>


double random_double(double min, double max)
{
    static std::uniform_real_distribution<double> distribution(min, max);
    static std::mt19937 generator;
    return distribution(generator);
}

Vector3d mix(Vector3d start,Vector3d end,double s){
    Vector3d out = (end - start) * s;
    return start + out;
}

Vector3d random_in_unit_sphere()
{
    double x = random_double(-1, 1);
    double y = random_double(-1, 1);
    double z = random_double(-1, 1);
    Vector3d d(x, y, z);
    if (d.size() > 1){
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
    
    Vector3d v(random_double(-1,1),random_double(1,1),0);
    if (v.size() > 1){
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
go::Ray::Ray()
{
}
go::Ray::Ray(Vector3d location, Vector3d direction) : m_location(location), m_direction(direction),m_time(random_double(0,1))
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

void go::Ray::lambertian(Vector3d normal, Vector3d hitPoint)
{
    m_direction = random_in_unit_sphere() + normal;
    if (m_direction.dot(m_direction) < epsilon)
    {
        m_direction = normal;
    }
    // m_direction.normalize();
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
        screenUV.x() * m_ratio * length ,
        screenUV.y()* length,-m_length);
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

go::Hitable::~Hitable()
{
}

go::Sphere::Sphere(Vector3d center, double radius, std::shared_ptr<Material> mat) : m_center(center), m_radius(radius), m_mat(mat),m_center2(center)
{
}

go::Sphere::Sphere(Vector3d center, Vector3d center2, double radius, std::shared_ptr<Material> mat): m_center(center), m_radius(radius), m_mat(mat),m_center2(center2)
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
        this->uv(result.uv,result.hit);
        result.isFront = n1.dot(ray.direction()) <= 0;
        return true;
    }
    return false;
}

void go::Sphere::uv(Vector2d &uv, const Vector3d &point)
{
    auto theta = acos(-point.y());
    auto phi = atan2(-point.z(), point.x()) + pi;

    double u = phi / (2*pi);        
    double v = theta / pi;
    uv[0] = u;
    uv[1] = v;
}

Vector3d go::Sphere::center(Ray & t)
{
    return mix(m_center,m_center2,t.time());
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

go::Planer::Planer(Vector3d point, Vector3d normal, std::shared_ptr<Material> mat) : m_point(point), m_normal(normal), m_mat(mat)
{
    Vector3d tx = Vector3d(1.0,0,0);
    Vector3d y = tx.cross(m_normal);
    Vector3d x = m_normal.cross(y);
    Matrix3d tra;
    tra << x,y,m_normal;
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
        this->uv(result.uv,result.hit);
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

   uv[0] = abs(std::fmod(uvw[0],1));
   uv[1] = abs(std::fmod(uvw[1],1));
}

go::Scene::Scene(double distance) : m_max_distance(distance)
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

bool go::Scene::hitOnce(Ray &ray, HitResult &out)
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
        if (hitOnce(ray, result))
        {
            Ray next;
            Vector4d color;
            auto m = result.mat->scatter(ray, color, result, next);
            auto nextColor = hit(next);
            return color.cwiseProduct(nextColor);
        }
        return Vector4d(1, 1, 1, 1);
    }

    return Vector4d(0, 0, 0, 1);
}

bool go::Material::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    return false;
}

go::Lambertian::Lambertian(Vector3d albedo) : m_texture(std::make_shared<Color>(albedo))
{
}

go::Lambertian::Lambertian(Texture *texure):m_texture(texure)
{
}

bool go::Lambertian::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    out = in;
    auto tColor = m_texture->color(hit.uv,hit.hit);
    color << tColor ,1.0;
    out.lambertian(hit.normal, hit.hit);
    return true;
}

go::Metal::Metal(Vector3d albedo, double fuzz) : m_texture(std::make_shared<Color>(albedo)), m_fuzz(fuzz)
{
}

bool go::Metal::scatter(const Ray &in, Vector4d &color, HitResult &hit, Ray &out)
{
    out = in;

    out.reflect(hit.normal, hit.hit, m_fuzz);
    auto tColor = m_texture->color(hit.uv,hit.hit);
    color << tColor ,1.0;
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
    bool ro = reflectance(cos_theta,index) > random_double(0,1);
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

go::Color::Color(Vector3d &albedo):m_albedo(albedo)
{

}

go::Color::~Color(){}

Vector3d go::Color::color(Vector2d uv, Vector3d &point)
{
    return m_albedo;
}



Vector3d go::TestColor::color(Vector2d uv, Vector3d &point)
{
    if(uv.x() > 0.5 && uv.y() > 0.5){
        return Vector3d(1,1,1);
    }else if(uv.x() <= 0.5 && uv.y() <= 0.5){
        return Vector3d(1,1,1);
    }else{
        return Vector3d(0.5,0.5,0.5);
    }
}
