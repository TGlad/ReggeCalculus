#pragma once
#include "basics.h"
#include "Simplices.h"

struct SpaceTime
{
  map<Vector4, Vertex> points;
  map<Vector4, Line> lines;
  map<Vector4, Triangle> triangles;
  map<Vector4, Tetrahedron> tetrahedrons;
  map<Vector4, Pentachoron> pentachorons;

  bool inBounds(const Vector4 &v);
  void init();
  void update();
  VectorXd getEdgeErrors(SparseMatrix<double> &jacobian);
  double getDeficitAngle(Triangle &bone);

private:
  void build();
  void connect();
  void setTriangleEdgeMatrices();
  void calculateSignatures();
  void calculateEdgeLengths();
};