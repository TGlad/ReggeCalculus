#pragma once
#include "basics.h"

// simplex structure.
struct Vertex
{
  Vertex(){}
  Vertex(const Vector4 &pos) : pos(pos) {}
  Vector4 pos;

  vector<struct Line *> edges;
  vector<struct Triangle *> triangles;
  vector<struct Tetrahedron *> tetrahedrons;
  vector<struct Pentachoron *> pentachorons;
};
struct Line
{
  Line(){}
  Line(Vertex *end0, Vertex *end1) { ends[0] = end0; ends[1] = end1; }

  Vertex *ends[2];
  double lengthSqr; // only value used in metric
  int index; // used in the vector representation

  vector<struct Triangle *> triangles;
  vector<struct Tetrahedron *> tetrahedrons;
  vector<struct Pentachoron *> pentachorons;

  double signature; // 1 is timelike
};
struct Triangle
{
  Triangle(){}
  Triangle(Vertex *corner0, Vertex *corner1, Vertex *corner2) { corners[0] = corner0; corners[1] = corner1; corners[2] = corner2; }

  Vertex *corners[3];
  Line *edges[3];

  vector<struct Tetrahedron *> tetrahedrons;
  vector<struct Pentachoron *> pentachorons;

  struct EdgeMat
  {
    struct Line *edges[5][5];
  };
  vector<EdgeMat> edgeMatrix;

  double signature; // 1 = timelike, -1 = spacelike, 0 = lightlike

  // cached data
  double deficitAngle;
  double areaSquared;
  SparseVector<double> areaSquaredDot, areaDot;
  SparseVector<double> deficitAngleDot; // deficit angle with respect to all edges in full complex
};
struct Tetrahedron
{
  Tetrahedron(){}
  Tetrahedron(Vertex *corner0, Vertex *corner1, Vertex *corner2, Vertex *corner3) { corners[0] = corner0; corners[1] = corner1; corners[2] = corner2; corners[3] = corner3; }

  Vertex *corners[4];
  Line *edges[6];
  Triangle *faces[4];

  vector<struct Pentachoron *> pentachorons;
};
struct Pentachoron
{
  Pentachoron(){}
  Pentachoron(Vertex *corner0, Vertex *corner1, Vertex *corner2, Vertex *corner3, Vertex *corner4) { corners[0] = corner0; corners[1] = corner1; corners[2] = corner2; corners[3] = corner3; corners[4] = corner4; }

  Vertex *corners[5];
  Line *edges[10];
  Triangle *faces[10];
  Tetrahedron *volumes[5];
};
