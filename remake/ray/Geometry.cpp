#include "Geometry.hpp"
#include <cmath>
#include <fstream>
#include <string>
using namespace vr;

Random::Generate vr::Primitives::s_random = Random::shared.uniform(0, 1);

vr::Triangle::Triangle(std::array<Vec3, 3> &points)
    : m_points(points)
{
    for (auto &&i : points)
    {
        m_bbox.add(i);
    }
    auto v = glm::normalize(glm::cross(m_points[1] - m_points[0], m_points[2] - m_points[0]));
    m_face_normal = v;
}

vr::Triangle::Triangle(std::array<Vec3, 3> &points, std::array<Vec3, 3> &normal)
    : m_points(points), m_normal(normal.begin(), normal.end())
{
    for (auto &&i : points)
    {
        m_bbox.add(i);
    }
}

vr::Triangle::Triangle(std::array<Vec3, 3> &points, std::array<Vec2, 3> &uv, std::array<Vec3, 3> &normal)
    : m_points(points), m_uv(uv.begin(), uv.end()), m_normal(normal.begin(), normal.end())
{
    for (auto &&i : points)
    {
        m_bbox.add(i);
    }
}

bool vr::Triangle::hitTest(
    const Ray &ray,
    const Region &region,
    HitRecord &record) const
{
    auto mode = hitMode();
    if (mode != back | front)
    {
        if (mode == front && glm::dot(ray.direction, m_face_normal) >= 0)
            return false;
        if (mode == back && glm::dot(ray.direction, m_face_normal) <= 0)
            return false;
    }

    auto s = ray.location - m_points[0];
    auto e1 = m_points[1] - m_points[0];
    auto e2 = m_points[2] - m_points[0];
    auto s1 = glm::cross(ray.direction, e2);
    auto s2 = glm::cross(s, e1);

    auto s1e1 = glm::dot(s1, e1);
    auto temp = glm::dot(s2, e2) / s1e1;
    auto b1 = glm::dot(s1, s) / s1e1;
    auto b2 = glm::dot(s2, ray.direction) / s1e1;
    if (temp >= 0 && b1 >= 0 && b2 >= 0 && (1 - b1 - b2) >= 0)
    {
        if (region.contain(temp))
        {
            record.t = temp;
            record.uvinterp.y = b1;
            record.uvinterp.z = b2;
            record.uvinterp.x = 1 - b1 - b2;
            record.position = ray.eval(temp);

            if (m_normal.size() < 3)
            {
                record.normal = m_face_normal;
            }
            else
            {
                record.normal = m_normal[0] * (1 - b1 - b2) + m_normal[1] * b1 + m_normal[2] * b2;
            }
            record.hitMode = glm::dot(ray.direction, record.normal) <= 0 ? front : back;
            if (m_uv.size() >= 3)
            {
                record.uv = m_uv[0] * (1 - b1 - b2) + m_uv[1] * b1 + m_uv[2] * b2;
            }
            return true;
        }
        return false;
    }
    return false;
}

const AABB &vr::Triangle::bbox() const
{
    return m_bbox;
}

Scalar vr::Triangle::area() const
{
    Scalar a = glm::distance(m_points[1], m_points[0]);
    Scalar b = glm::distance(m_points[2], m_points[1]);
    Scalar c = glm::distance(m_points[0], m_points[2]);
    if (a + b > c && a + c > b && b + c > a)
    {
        Scalar s = (a + b + c) / 2;
        Scalar A = std::sqrt(s * (s - a) * (s - b) * (s - c));
        return a;
    }
    return 0;
}

Vec3 vr::Triangle::random() const
{
    Scalar a = s_random();
    Scalar b = s_random();
    Scalar c = s_random();
    Scalar sum = a + b + c;
    Vec3 p = a / sum * m_points[0] + b / sum * m_points[1] + c / sum * m_points[2];
    return p;
}

Vec3 vr::Triangle::normal(const Vec3 &point) const
{
    return m_face_normal;
}

Vec3 vr::Triangle::lightNormal(const Vec3 &lightDirection) const
{
    return m_face_normal;
}

vr::Sphere::Sphere(const Vec3 &center, Scalar radius) : m_center(center), m_radius(radius)
{
    auto p1 = m_center + radius;
    auto p2 = m_center - radius;
    m_bbox.x.min = std::min(p1.x, p2.x);
    m_bbox.y.min = std::min(p1.x, p2.x);
    m_bbox.z.min = std::min(p1.x, p2.x);

    m_bbox.x.max = std::max(p1.x, p2.x);
    m_bbox.y.max = std::max(p1.x, p2.x);
    m_bbox.z.max = std::max(p1.x, p2.x);
}

bool vr::Sphere::hitTest(const Ray &ray, const Region &region, HitRecord &record) const
{
    auto so = ray.location - m_center;
    auto a = glm::dot(ray.direction, ray.direction);
    auto b = 2 * glm::dot(so, ray.direction);
    auto c = glm::dot(so, so) - m_radius * m_radius;
    auto mode = hitMode();
    auto delta = b * b - 4 * a * c;
    if (delta >= 0)
    {
        auto t1 = (-b + sqrt(delta)) / (a * 2);
        auto t2 = (-b - sqrt(delta)) / (a * 2);
        auto tmin = std::min(t1, t2);
        auto tmax = std::max(t2, t1);
        auto t = region.contain(tmin) ? tmin : tmax;
        if (!region.contain(t))
            return false;
        record.t = t;
        record.position = ray.eval(t);
        record.normal = record.position - m_center;
        record.hitMode = glm::dot(ray.direction, record.normal) <= 0 ? front : back;
        record.uv = Vec2(0, 0);
        record.uvinterp = {1, 0, 0};
        return true;
    }
    return false;
}

const AABB &vr::Sphere::bbox() const
{
    return m_bbox;
}

Vec3 vr::Sphere::random() const
{
    Scalar a = s_random() * 2 - 1;
    Scalar b = s_random() * 2 - 1;
    Scalar c = s_random() * 2 - 1;
    Vec3 point = glm::normalize(Vec3{a, b, c}) * m_radius;
    return point * s_random() + m_center;
}

Scalar vr::Sphere::area() const
{
    return M_PI * m_radius * m_radius;
}

Vec3 vr::Sphere::lightNormal(const Vec3 &lightDirection) const
{
    return glm::normalize(lightDirection);
}

// vr::Plant::Plant(Vec3 center, Vec3 normal) : m_center(center), m_normal(normal)
// {
// }

// bool vr::Plant::hitTest(const Ray &ray, const Region &region, HitRecord &record) const
// {
//     auto mode = hitMode();
//     if (glm::dot(ray.direction, m_normal) >= 0 && mode & vr::front)
//     {
//         return false;
//     }

//     if (glm::dot(ray.direction, m_normal) <= 0 && mode & vr::back)
//     {
//         return false;
//     }

//     Scalar a = glm::dot(m_center, m_normal);
//     Scalar b = glm::dot(ray.location, m_normal);
//     Scalar c = glm::dot(ray.direction, m_normal);
//     auto t = (a - b) / c;
//     if (region.contain(t))
//     {
//         record.position = ray.eval(t);
//         record.normal = m_normal;
//         record.hitMode = glm::dot(ray.direction, record.normal) <= 0 ? front : back;
//         record.t = t;
//         record.uvinterp = Vec3(1, 0, 0);
//         record.uv = {0, 0};
//         return true;
//     }
//     return false;
// }

// const AABB &vr::Plant::bbox() const
// {
//     return m_bbox;
// }

// Vec3 vr::Plant::random() const
// {
//     return Vec3(0);
// }

// Scalar vr::Plant::area() const
// {
//     return 0;
// }

// Vec3 vr::Plant::lightNormal(const Vec3 &lightDirection) const
// {
//     return m_normal;
// }

vr::Quadrilateral::Quadrilateral(std::array<Vec3, 4> &points)
{
    createTriangles(points);
}

void vr::Quadrilateral::createTriangles(std::array<vr::Vec3, 4> &points)
{
    std::array<Vec3, 3> a1 = {points[0], points[1], points[2]};
    vr::Triangle tr1(a1);

    std::array<Vec3, 3> a2 = {points[0], points[2], points[3]};
    vr::Triangle tr2(a2);
    m_triangles.push_back(tr1);
    m_triangles.push_back(tr2);
    m_aabb.add(tr1.bbox());
    m_aabb.add(tr2.bbox());
}

vr::Quadrilateral::Quadrilateral(const Vec3 &normal, const Vec3 &up, const Vec3 &center, const Scalar radius)
{
    Vec3 x = glm::normalize(glm::cross(up, normal));
    Vec3 y = glm::normalize(glm::cross(normal, x));
    Vec3 z = glm::normalize(normal);
    Mat3 matx = {x, y, z};
    Vec3 p1 = matx * Vec3{-radius, radius, 0} + center;
    Vec3 p2 = matx * Vec3{-radius, -radius, 0} + center;
    Vec3 p3 = matx * Vec3{radius, -radius, 0} + center;
    Vec3 p4 = matx * Vec3{radius, radius, 0} + center;
    std::array<vr::Vec3, 4> ps = {p1, p2, p3, p4};
    createTriangles(ps);
}

bool vr::Quadrilateral::hitTest(const Ray &ray, const Region &region, HitRecord &record) const
{
    return m_triangles[0].hitTest(ray, region, record) || m_triangles[1].hitTest(ray, region, record);
}

const AABB &vr::Quadrilateral::bbox() const
{
    return m_aabb;
}

Scalar vr::Quadrilateral::area() const
{
    return m_triangles[0].area() + m_triangles[1].area();
}

Vec3 vr::Quadrilateral::random() const
{
    return s_random() > 0.5 ? m_triangles[0].random() : m_triangles[1].random();
}

Vec3 vr::Quadrilateral::lightNormal(const Vec3 &lightDirection) const
{
    return m_triangles[0].lightNormal(lightDirection);
}

vr::ObjectReader::ObjectReader(const char *path)
{
    std::fstream f(path, std::ios::in);
    std::vector<Vec3> temp_vector;
    std::vector<Vec3> temp_normal;
    std::string buffer;
    buffer.resize(100);
    if (f.is_open())
    {
        while (!f.eof())
        {
            f.getline(buffer.data(), buffer.size());
            if (buffer.compare(0, 1, "v") == 0)
            {
                auto p = buffer.substr(2, buffer.size() - 2);
                std::string::size_type size;
                Scalar x = std::stod(p, &size);
                p = p.substr(size);
                Scalar y = std::stod(p, &size);
                p = p.substr(size);
                Scalar z = std::stod(p, &size);
                temp_vector.push_back({x, y, z});
            }
            if (buffer.compare(0, 2, "vn") == 0)
            {
                auto p = buffer.substr(2, buffer.size() - 2);
                std::string::size_type size;
                Scalar x = std::stod(p, &size);
                p = p.substr(size);
                Scalar y = std::stod(p, &size);
                p = p.substr(size);
                Scalar z = std::stod(p, &size);
                temp_normal.push_back({x, y, z});
            }
            if (buffer.compare(0, 1, "f") == 0)
            {
                auto p = buffer.substr(2, buffer.size() - 2);
                std::string::size_type size;
                Uint p1v = std::stol(p, &size);
                addBias(size, p);
                p = p.substr(size);
                Uint p1n = std::stoi(p, &size);
                addBias(size, p);
                p = p.substr(size);
                Uint p2v = std::stol(p, &size);
                addBias(size, p);
                p = p.substr(size);
                Uint p2n = std::stoi(p, &size);
                addBias(size, p);
                p = p.substr(size);
                Uint p3v = std::stol(p, &size);
                addBias(size, p);
                p = p.substr(size);
                Uint p3n = std::stoi(p, &size);
                std::array<Vec3,3> va = {temp_vector[p1v],temp_vector[p2v],temp_vector[p3v]};
                std::array<Vec3,3> na = {temp_normal[p1n],temp_normal[p2n],temp_normal[p3n]};

                m_tri.push_back(std::shared_ptr<Triangle>(new Triangle(va,na)));
            }
        }
    }

    f.close();
}

void vr::ObjectReader::addBias(std::size_t &size, std::string &p)
{
    for (int i = size; i < p.size(); i++)
    {
        if (p[i] == '/')
        {
            size++;
        }
        else{
            return;
        }
    }
}

const std::list<std::shared_ptr<Triangle>> &vr::ObjectReader::tiangles() const
{
    return m_tri;
}
