#include "Primitives.hpp"
#include <iostream>

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

bool go::Triangle::hit(Ray &ray, Interval ray_t, HitResult &result)
{
    if (ray.direction().dot(m_normal) >= 0 && is_double_face == false)
    {
        return false;
    }
    double A = m_normal.dot(ray.direction());
    double B = (m_point[0].point - ray.location()).dot(m_normal);
    double t = B / A;
    if(!ray_t.contains(t)){
        return false;
    }
    Vector3d hit = ray.exec(t);
    if (!contain(hit))
    {
        return false;
    }

    result.t = t;
    result.hit = ray.exec(t);
    result.normal = m_normal;
    result.mat = m_mat;
    result.isFront = A <= 0;
    this->uv(result.uv, result.hit);
    return true;
}

void go::Triangle::uv(Vector2d &uv, const Vector3d &point)
{
    uv = interpolateUV(
        point,
        this->m_point[0].point,
        this->m_point[1].point,
        this->m_point[2].point,
        this->m_point[0].uv,
        this->m_point[1].uv,
        this->m_point[2].uv
        );
}

bool go::Triangle::contain(const Vector3d &p)
{
    return same_size(m_point[0].point, m_point[1].point, m_point[2].point, p) &&
           same_size(m_point[1].point, m_point[2].point, m_point[0].point, p) &&
           same_size(m_point[2].point, m_point[0].point, m_point[1].point, p);
}

void go::Triangle::transform(const Matrix4d& transform)
{
    for (size_t i = 0; i < 3; i++)
    {
        Vector4d temp = (transform * (m_point[i].point.homogeneous()));
        m_point[i].point << temp.head<3>();
    }
    Matrix4d norm = transform.inverse().transpose();
    Vector4d temp = norm * m_normal.homogeneous();
    m_normal = temp.head<3>().normalized();
}

bool &go::Triangle::double_face()
{
    return is_double_face;
}

bool go::Triangle::same_size(const Vector3d &point0, const Vector3d &point1, const Vector3d &point2, const Vector3d &p)
{
    Vector3d ab = point1 - point0;
    Vector3d ac = point2 - point0;
    Vector3d ap = p - point0;
    Vector3d v1 = ab.cross(ac);
    Vector3d v2 = ab.cross(ap);

    return v1.dot(v2) >= 0 && v1.normalized().dot(v2.normalized()) > 0.99 ;
}

go::Triangle::Triangle()
{
}

go::Triangle::Triangle(go::Vertex point1, go::Vertex point2, go::Vertex point3, std::shared_ptr<Material> m) : m_mat(m)
{
    m_point[0] = point1;
    m_point[1] = point2;
    m_point[2] = point3;
    m_normal = (point2.point - point1.point).cross(point3.point - point2.point);
    m_normal.normalize();
}

go::Quad::Quad(std::shared_ptr<Material> m)
{
    Vertex v0( 1,0,-1,1,0);
    Vertex v1(-1,0,-1,0,0);
    Vertex v2(-1,0, 1,0,1);
    Vertex v3( 1,0, 1,1,1);

    m_triangles[0] = go::Triangle(v0,v1,v2,m);
    m_triangles[1] = go::Triangle(v0,v2,v3,m);
}

go::Quad::Quad(Vertex point1, Vertex point2, Vertex point3, Vertex point4, std::shared_ptr<Material> m)
{
    m_triangles[0] = go::Triangle(point1,point2,point3,m);
    m_triangles[1] = go::Triangle(point1,point3,point4,m);
}

bool go::Quad::hit(Ray &ray, Interval ray_t, HitResult &result)
{
    for (size_t i = 0; i < 2; i++)
    {
        m_triangles[i].double_face() = is_double_face;
    }
    HitResult ret1,ret2;
    bool t1 = m_triangles[0].hit(ray,ray_t,ret1);
    bool t2 = m_triangles[1].hit(ray,ray_t,ret2);
    if(!t1 && !t2 ){
        return false;
    }
    if(t1 && t2){
        result = ret1.t > ret2.t ? ret2 : ret1;
        return true;
    }else if(t2){
        result = ret2;
        return true;
    }else{
        result = ret1;
        return true;
    }
}

void go::Quad::uv(Vector2d &uv, const Vector3d &point)
{
    if(m_triangles[0].contain(point)){
        m_triangles[0].uv(uv,point);
    }
    if(m_triangles[1].contain(point)){
        m_triangles[1].uv(uv,point);
    }
}

void go::Quad::transform(const Matrix4d &transform)
{
    for (size_t i = 0; i < 2; i++)
    {
        m_triangles[i].transform(transform);
    }
    
}

bool &go::Quad::double_face()
{
    return is_double_face;
}

go::Box::Box(std::shared_ptr<Material> m)
{
    Vector4d v[5] = {
        Vector4d(0,1,0,pi / 2),
        Vector4d(0,1,0,pi / -2),
        Vector4d(0,1,0,pi),
        Vector4d(1,0,0,pi / 2),
        Vector4d(1,0,0,pi / -2)
    };
    go::Quad q(
            go::Vertex( 1, 1,1,1,0),
            go::Vertex(-1, 1,1,0,0),
            go::Vertex(-1,-1,1,0,1),
            go::Vertex( 1,-1,1,1,1),m);
    q.double_face() = true;
    m_quad[0] = q;
    
    for (size_t i = 0; i < 5; i++)
    {
        go::Quad q(
            go::Vertex( 1, 1,1,1,0),
            go::Vertex(-1, 1,1,0,0),
            go::Vertex(-1,-1,1,0,1),
            go::Vertex( 1,-1,1,1,1),m);
        q.double_face() = true;
        q.transform(rotate(v[i].x(),v[i].y(),v[i].z(),v[i].w()));
        m_quad[i + 1] = q;
    }
    
}

bool go::Box::hit(Ray &ray, Interval ray_t, HitResult &result)
{
    HitResult ret;
    bool has_hit = false;
    ret.t = infinity;
    for (size_t i = 0; i < 6; i++)
    {
        Ray tray = ray;
        HitResult r;
        if (!m_quad[i].hit(tray,ray_t,r)){
            continue;
        }
        
        if(r.t < ret.t && ray_t.contains(r.t)){
            ret = r;
            has_hit = true;
        }
    }
    if(!has_hit){
        return false;
    }
    result = ret;
    return true;
}

void go::Box::uv(Vector2d &uv, const Vector3d &point)
{

}

void go::Box::transform(const Matrix4d &transform)
{
    for (size_t i = 0; i < 6; i++)
    {
        m_quad[i].transform(transform);
    }
}

go::Vertex::Vertex(double x, double y, double z, double u, double v):point(x,y,z),uv(u,v)
{

}
