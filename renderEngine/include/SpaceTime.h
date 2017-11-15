#pragma once
#include "basics.h"
#include "Simplices.h"


struct Spacetime
{
  map<Vector4, Vertex> points;
  map<Vector4, Line> lines;
  map<Vector4, Triangle> triangles;
  map<Vector4, Tetrahedron> tetrahedrons;
  map<Vector4, Pentachoron> pentachorons;

  bool inBounds(const Vector4 &v);
  void init();
  void update();
  VectorXd getEdgeErrors(SparseMatrix<double> &jacobian) const;
  double getDeficitAngle(const Triangle &bone) const;
  void validate() const;

private:
  build();
  connect();
  setTriangleEdgeMatrices();
  calculateSignatures();
};