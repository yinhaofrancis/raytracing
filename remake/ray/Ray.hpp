#pragma once
#include "Util.hpp"

namespace vr
{
    class Ray;

    struct Region
    {
        Scalar min = -INFINITY;
        Scalar max = INFINITY;
        inline Scalar size() const
        {
            return max - min;
        }
        inline bool contain(Scalar t) const
        {
            return t > min && t < max;
        }
        inline bool contain(const Region &t) const
        {
            return this->contain(t.min) && this->contain(t.max);
        }
        inline void addmin(Scalar v)
        {
            min = min > v ? v : min;
        }
        inline void addmax(Scalar v)
        {
            max = max < v ? v : max;
        }
        inline void add(Scalar v){
            addmin(v);
            addmax(v);
        }
        inline void add(Region region)
        {
            addmin(region.min);
            addmax(region.max);
        }
        inline bool hasIntersection(const Region &t) const
        {
            return this->contain(t.min) || this->contain(t.max);
        }
        inline bool isEmpty() const
        {
            return min >= max;
        }
        inline void intersectionSet(Region v)
        {
            if (isEmpty())
            {
            }
            else if (v.isEmpty())
            {
                min = v.min;
                max = v.max;
            }
            else
            {
                min = std::max(min, v.min);
                max = std::min(max, v.max);
            }
        }

        inline void unionSet(Region v)
        {
            if (isEmpty())
            {
                min = v.min;
                max = v.max;
            }
            else if (v.isEmpty())
            {
            }
            else
            {
                min = std::min(min, v.min);
                max = std::max(max, v.max);
            }
        }
    };
    struct AABB
    {
        Region x{1, -1};
        Region y{1, -1};
        Region z{1, -1};

        inline void intersectionSet(AABB b)
        {
            x.intersectionSet(b.x);
            y.intersectionSet(b.y);
            z.intersectionSet(b.z);
        }
        inline void unionSet(AABB b)
        {
            x.unionSet(b.x);
            y.unionSet(b.y);
            z.unionSet(b.z);
        }
        inline void add(const AABB& box)
        {
            x.add(box.x);
            y.add(box.y);
            z.add(box.z);
        }
        inline void add(const Vec3& v){
            x.add(v.x);
            y.add(v.y);
            z.add(v.z);
        }
        inline bool isEmpty() const
        {
            return x.isEmpty() && y.isEmpty() && z.isEmpty();
        }
        inline bool contain(const AABB &box) const
        {
            return x.contain(box.x) && y.contain(box.y) && z.contain(box.z);
        }
        inline bool intersection(const AABB &box) const
        {
            return x.contain(box.x) || y.contain(box.y) || z.contain(box.z);
        }
        void cut(AABB &left, AABB &right) const
        {
            if (x.size() > y.size() && x.size() > z.size())
            {
                Scalar center = (x.min + x.max) / 2;
                left = *this;
                right = *this;
                left.x.max = center;
                right.x.min = center;
            }
            else if (y.size() > z.size() && y.size() > x.size())
            {
                Scalar center = (y.min + y.max) / 2;
                left = *this;
                right = *this;
                left.y.max = center;
                right.y.min = center;
            }
            else
            {
                Scalar center = (z.min + z.max) / 2;
                left = *this;
                right = *this;
                left.z.max = center;
                right.z.min = center;
            }
        }
        inline bool hitAxisAlign(Scalar align, Scalar o, Scalar d, Scalar &t) const
        {
            if (d == 0)
            {
                return false;
            }
            t = (align - o) / d;
            return true;
        }
        bool hitTest(const Ray &ray, Scalar &t) const;
    };

    class Random
    {
    public:
        class Generate
        {
        public:
            Scalar operator()() const;
            Generate(std::function<Scalar()>);
            Generate(Generate &&g) : m_call(g.m_call)
            {
            }

        private:
            std::function<Scalar()> m_call;
        };

        class Sphere
        {
        public:
            Vec3 operator()() const;
            Vec3 operator()(const Vec3 &normal) const;
            Sphere(std::function<Vec3()>);
            Sphere(Sphere &&);

        private:
            std::function<Vec3()> m_call;
        };

        Random();
        Generate normal(Scalar mean, Scalar stddev) const;
        Generate uniform(Scalar min, Scalar max) const;
        Sphere sphere() const;
        Sphere circle() const;
        Sphere sphere(Scalar zmin, Scalar zmax);
        const static Random shared;

    private:
        static std::random_device m_rd;
    };

    class Sensor
    {
    public:
        Sensor(Scalar pixSize = 1) : m_pixel_size(pixSize), m_random(Random::shared.uniform(-pixSize / 2, pixSize / 2)) {};
        Vec3 operator()(Uint x, Uint y, Uint w, Uint h) const;
        Scalar &pixelSize() { return m_pixel_size; }

    private:
        Scalar m_pixel_size;
        Random::Generate m_random;
    };
    class Camera
    {
    public:
        Camera(Uint width, Uint height, const Vec3 &location, const Vec3 &lookAt, const Vec3 &up, Scalar fov);

        Ray ray(Uint x, Uint y) const;

        Uint width() const;

        Uint height() const;

    private:
        Sensor m_front;
        Vec3 m_location, m_look_at, m_up;
        Scalar m_width, m_height, m_fov, m_length;
        Mat3 m_lookAt_matrix;
    };

    struct Ray
    {
        Vec3 location;
        Vec3 direction;
        Uint count = 0;

        Ray reflect(const Vec3 &position, const Vec3 &normal) const;
        Ray refract(const Vec3 &position, const Vec3 &norma, const Scalar eta) const;
        Ray hemisphere(const Vec3 &position, const Vec3 &normal) const;

        Vec3 eval(const Scalar t) const
        {
            return location + direction * t;
        }
        Ray move(Scalar delta_t) const
        {
            return Ray{eval(delta_t), direction};
        }
        friend std::ostream &operator<<(std::ostream &out, const Ray &ray);
        static const Random::Sphere sphere;
    };

    std::ostream &operator<<(std::ostream &out, const Ray &ray);
} // namespace vr
