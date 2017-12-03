#include "basics.h"


Vector4 operator *(double f, const Vector4& v)
{
  return Vector4(v.t * f, v.x * f, v.y * f, v.z * f);
}
Vector4 operator *(const Matrix<double, 4, 4> &mat, const Vector4& v)
{
  Vector4d vec(v.t, v.x, v.y, v.z);
  Vector4d res = mat * vec;
  return Vector4(res[0], res[1], res[2], res[3]);
}
