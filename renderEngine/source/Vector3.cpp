#include "basics.h"

void Vector3::clamp(float minVal, float maxVal)
{
  x = clamped(x, minVal, maxVal);
  y = clamped(y, minVal, maxVal);
  z = clamped(z, minVal, maxVal);
}

Vector3 operator *(float f, const Vector3& v)
{
  return Vector3(v.x * f, v.y * f, v.z * f);
}

Vector3 Vector3::getRotationVector(const Direction& from, const Direction& to)
{
  Vector3 cross = Vector3::cross(from, to);
  float dot = from.dot(to);
  float mag = cross.normalise();
  Angle angle = atan2f(mag, dot);
  return cross * angle;
}

Vector4 operator *(double f, const Vector4& v)
{
  return Vector4(v.t * f, v.x * f, v.y * f, v.z * f);
}