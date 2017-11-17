#include "basics.h"


Vector4 operator *(double f, const Vector4& v)
{
  return Vector4(v.t * f, v.x * f, v.y * f, v.z * f);
}