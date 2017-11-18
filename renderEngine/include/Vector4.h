#pragma once
#include "basics.h"

// even though this has a time component, we don't give it a Lorentzian signature here, that is held separately in the metric g
class Vector4
{
public:
  Vector4(){}
  inline Vector4(double T, double X, double Y, double Z)
  {
    t = T;  x = X; y = Y; z = Z;
  }
  double t, x, y, z;

  inline const double &operator[](int i) const
  {
//    ASSERT(i >= 0 && i < 4);
    const double *n[4] = { &t, &x, &y, &z };
    return *(n[i]);
  }
  // TODO: this is a convenience function for maps
  inline bool operator <(const Vector4 &f) const
  {
    double a = t + 100.1*x + 100.1*100.1*y + 100.1*100.1*100.1*z;
    double b = f.t + 100.1*f.x + 100.1*100.1*f.y + 100.1*100.1*100.1*f.z;
    return a < b;
  }
  inline Vector4 operator *(double f) const
  {
    return Vector4(t * f, x * f, y * f, z * f);
  }
  inline Vector4 operator *(const Vector4& vec) const
  {
    return Vector4(t * vec.t, x * vec.x, y * vec.y, z * vec.z);
  }
  inline bool operator ==(const Vector4 &vec)
  {
    return t == vec.t && x == vec.x && y == vec.y && z == vec.z;
  }
  inline void operator *=(double f)
  {
    t *= f;
    x *= f;
    y *= f;
    z *= f;
  }
  inline Vector4 operator /(double f) const
  {
    double div = 1.0 / f;
    return Vector4(t * div, x * div, y * div, z * div);
  }
  inline Vector4 operator /(const Vector4& div) const
  {
    return Vector4(t / div.t, x / div.x, y / div.y, z / div.z);
  }
  inline void operator /=(double f)
  {
    double div = 1.0 / f;
    t *= div;
    x *= div;
    y *= div;
    z *= div;
  }
  inline Vector4 operator +(const Vector4& v) const
  {
    return Vector4(t + v.t, x + v.x, y + v.y, z + v.z);
  }
  inline void operator +=(const Vector4& v)
  {
    t += v.t;
    x += v.x;
    y += v.y;
    z += v.z;
  }
  inline Vector4 operator -() const
  {
    return Vector4(-t, -x, -y, -z);
  }
  inline Vector4 operator -(const Vector4& v) const
  {
    return Vector4(t - v.t, x - v.x, y - v.y, z - v.z);
  }
  inline void operator -=(const Vector4& v)
  {
    t -= v.t;
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }
  inline void set(const Vector4& v)
  {
    t = v.t;
    x = v.x;
    y = v.y;
    z = v.z;
  }
  inline void set(double _t, double _x, double _y, double _z)
  {
    t = _t;
    x = _x;
    y = _y;
    z = _z;
  }
  inline void setToZero()
  {
    t = x = y = z = 0.0;
  }
  inline double magnitudeSquared() const
  {
    return t*t + x*x + y*y + z*z;
  }
  inline double dot(const Vector4& v) const
  {
    return t * v.t + x * v.x + y * v.y + z * v.z;
  }
};
Vector4 operator *(double f, const Vector4& v);
