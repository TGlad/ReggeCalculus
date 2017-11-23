#include "SpaceTime.h"
#include <ostream>

static const int tMax = 4;// 10;
static const int wMax = 4;// 10;

static const int TMax = 2 * tMax - 1;
static const int WMax = 2 * wMax - 1;

bool SpaceTime::inBounds(const Vector4 &v)
{
  return v.t >= 0 && v.t < TMax &&
    v.x >= 0 && v.x < WMax &&
    v.y >= 0 && v.y < WMax &&
    v.z >= 0 && v.z < WMax;
}

void SpaceTime::build()
{
  // 1. points and vertices
  for (int t = 0; t < TMax; t++) for (int x = 0; x < WMax; x++) for (int y = 0; y < WMax; y++) for (int z = 0; z < WMax; z++)
  {
    bool corner = (t % 2) == 0 && (x % 2) == 0 && (y % 2) == 0 && (z % 2) == 0;
    bool middle = (t % 2) == 1 && (x % 2) == 1 && (y % 2) == 1 && (z % 2) == 1;
    if (corner || middle)
      points[Vector4(t, x, y, z)] = Vertex(Vector4(t, x, y, z));
  }

  // all the rest
  for (int t = 0; t < TMax; t += 2) for (int x = 0; x < WMax; x += 2) for (int y = 0; y < WMax; y += 2) for (int z = 0; z < WMax; z += 2)
  {
    bool cornerPick = !((t / 2) % 2 + (x / 2) % 2 + (y / 2) % 2 + (z / 2) % 2);
    if (!cornerPick)
      continue;

    // do this once in every direction, i.e. 16 directions
    for (int st = -1; st < 2; st += 2) for (int sx = -1; sx < 2; sx += 2) for (int sy = -1; sy < 2; sy += 2) for (int sz = -1; sz < 2; sz += 2)
    {
      Vector4 corners[5] = { Vector4(0, 0, 0, 0), Vector4(2 * st, 0, 0, 0), Vector4(0, 2 * sx, 0, 0), Vector4(0, 0, 2 * sy, 0), Vector4(0, 0, 0, 2 * sz) };
      Vector4 event(t, x, y, z);
      Vector4 coords[5];
      Vertex *vertices[5];
      for (int h = 0; h < 2; h++)
      {
        for (int i = 0; i < 5; i++)
        {
          coords[0] = event + corners[i];
          if (i == 0 && h)
            coords[0] += Vector4(st, sx, sy, sz); // second pentachoron goes from middle rather than corner
          if (inBounds(coords[0]))
          {
            vertices[0] = &points[coords[0]];
            for (int j = i + 1; j < 5; j++)
            {
              coords[1] = event + corners[j];
              if (inBounds(coords[1]))
              {
                vertices[1] = &points[coords[1]];
                if (lines.find(coords[0] + coords[1]) == lines.end())
                  lines[coords[0] + coords[1]] = Line(vertices[0], vertices[1]);

                for (int k = j + 1; k < 5; k++)
                {
                  coords[2] = event + corners[k];
                  if (inBounds(coords[2]))
                  {
                    vertices[2] = &points[coords[2]];
                    if (triangles.find(coords[0] + coords[1] + coords[2]) == triangles.end())
                      triangles[coords[0] + coords[1] + coords[2]] = Triangle(vertices[0], vertices[1], vertices[2]);

                    for (int l = k + 1; l < 5; l++)
                    {
                      coords[3] = event + corners[l];
                      if (inBounds(coords[3]))
                      {
                        vertices[3] = &points[coords[3]];
                        Vector4 ev = coords[0] + coords[1] + coords[2] + coords[3];
                        if (tetrahedrons.find(ev) == tetrahedrons.end())
                          tetrahedrons[ev] = Tetrahedron(vertices[0], vertices[1], vertices[2], vertices[3]);

                        for (int m = l + 1; m < 5; m++)
                        {
                          coords[4] = event + corners[m];
                          if (inBounds(coords[4]))
                          {
                            vertices[4] = &points[coords[4]];
                            Vector4 ev = coords[0] + coords[1] + coords[2] + coords[3] + coords[4];
                            if (pentachorons.find(ev) == pentachorons.end())
                            {
                              double sign = st*sx*sy*sz;
                              if (h)
                              {
                                sign = -sign;
                                if (i != 0)
                                  cout << "bll" << endl;
                              }
                              if (sign > 0)
                                pentachorons[ev] = Pentachoron(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4]);
                              else
                                pentachorons[ev] = Pentachoron(vertices[0], vertices[1], vertices[2], vertices[4], vertices[3]); // reverse orientation, is this right?
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void SpaceTime::connect()
{
  // upwards connections
  for (auto &tripair : triangles)
  {
    Triangle &tri = tripair.second;
    for (int i = 0; i < 3; i++)
      tri.edges[i] = &lines[tri.corners[i]->pos + tri.corners[(i + 1) % 3]->pos];
  }
  for (auto &tetpair : tetrahedrons)
  {
    Tetrahedron &tet = tetpair.second;
    int i = 0;
    for (int a = 0; a < 4; a++)
      for (int b = a + 1; b < 4; b++)
        tet.edges[i++] = &lines[tet.corners[a]->pos + tet.corners[b]->pos];
    i = 0;
    for (int a = 0; a < 4; a++)
      for (int b = a + 1; b < 4; b++)
        for (int c = b + 1; c < 4; c++)
          tet.faces[i++] = &triangles[tet.corners[a]->pos + tet.corners[b]->pos + tet.corners[c]->pos];
  }
  for (auto &pentpair : pentachorons)
  {
    Pentachoron &pent = pentpair.second;
    int i = 0;
    for (int a = 0; a < 5; a++)
      for (int b = a + 1; b < 5; b++)
        pent.edges[i++] = &lines[pent.corners[a]->pos + pent.corners[b]->pos];
    i = 0;
    for (int a = 0; a < 5; a++)
      for (int b = a + 1; b < 5; b++)
        for (int c = b + 1; c < 5; c++)
          pent.faces[i++] = &triangles[pent.corners[a]->pos + pent.corners[b]->pos + pent.corners[c]->pos];
    i = 0;
    for (int a = 0; a < 5; a++)
      for (int b = a + 1; b < 5; b++)
        for (int c = b + 1; c < 5; c++)
          for (int d = c + 1; d < 5; d++)
            pent.volumes[i++] = &tetrahedrons[pent.corners[a]->pos + pent.corners[b]->pos + pent.corners[c]->pos + pent.corners[d]->pos];
  }


  // downwards connections
  for (auto &pentpair : pentachorons)
  {
    Pentachoron &pent = pentpair.second;
    for (int i = 0; i < 5; i++)
      pent.corners[i]->pentachorons.push_back(&pent);
    for (int i = 0; i < 10; i++)
      pent.edges[i]->pentachorons.push_back(&pent);
    for (int i = 0; i < 10; i++)
      pent.faces[i]->pentachorons.push_back(&pent);
    for (int i = 0; i < 5; i++)
      pent.volumes[i]->pentachorons.push_back(&pent);
  }
  for (auto &tetpair : tetrahedrons)
  {
    Tetrahedron &tet = tetpair.second;
    for (int i = 0; i < 4; i++)
      tet.corners[i]->tetrahedrons.push_back(&tet);
    for (int i = 0; i < 6; i++)
      tet.edges[i]->tetrahedrons.push_back(&tet);
    for (int i = 0; i < 4; i++)
      tet.faces[i]->tetrahedrons.push_back(&tet);
  }
  for (auto &tripair : triangles)
  {
    Triangle &tri = tripair.second;
    for (int i = 0; i < 3; i++)
      tri.corners[i]->triangles.push_back(&tri);
    for (int i = 0; i < 3; i++)
      tri.edges[i]->triangles.push_back(&tri);
  }
  int index = 0;
  for (auto &linepair : lines)
  {
    Line &line = linepair.second;
    for (int i = 0; i < 2; i++)
      line.ends[i]->edges.push_back(&line);
    line.index = index++;
  }
}

void SpaceTime::calculateSignatures()
{
  double eps = 1e-5;
  for (auto &edgePair : lines)
  {
    Line &edge = edgePair.second;
    Vector4 l = edge.ends[1]->pos - edge.ends[0]->pos;
    double S = sqrt(sqr(l.x) + sqr(l.y) + sqr(l.z));
    double T = abs(l.t);
    double eps = 1e-10;
    edge.signature = T > S + eps ? -1.0 : (T < S - eps ? 1 : 0);
  }
  for (auto &triPair : triangles)
  {
    Triangle &tri = triPair.second;
    double difs[3];
    for (int i = 0; i < 3; i++)
      difs[i] = abs(tri.corners[i]->pos.t - tri.corners[(i + 1) % 3]->pos.t);
    if (difs[0]>eps && difs[1] > eps && difs[2] > eps)
      tri.signature = -1; // temporal
    else
    {
      int a = difs[0] < difs[1] ? (difs[0] < difs[2] ? 0 : 2) : (difs[1] < difs[2] ? 1 : 2);
      Vector4 pa = tri.corners[a]->pos;
      Vector4 pb = tri.corners[(a + 1) % 3]->pos;
      Vector4 pc = tri.corners[(a + 2) % 3]->pos;
      double t = (pb - pa).dot(pc - pa) / (pb - pa).magnitudeSquared();
      ASSERT(t >= -eps && t <= 1.0 + eps);
      Vector4 pab = pa + (pb - pa) * t;
      Vector4 dir = pab - pc;
      double x = sqrt(sqr(dir.x) + sqr(dir.y) + sqr(dir.z));
      if (x - abs(dir.t) > eps)
        tri.signature = 1; // spacelike
      else if (x - abs(dir.t) < -eps)
        tri.signature = -1; // timelike
      else
        tri.signature = 0; // lightlile
    }
  }
  /*
  for (auto &tetpair : tetrahedrons)
  {
  Tetrahedron &tet = tetpair.second;
  Vector4 m0 = tet.corners[1].pos - tet.corners[0].pos;
  Vector4 m1 = tet.corners[2].pos - tet.corners[0].pos;
  Vector4 m2 = tet.corners[3].pos - tet.corners[0].pos;
  Matrix33 A, B, C, D;
  A  << m0.x, m0.y, m0.z,
  m1.x, m1.y, m1.z,
  m2.x, m2.y, m2.z;
  B  << m0.t, m0.y, m0.z,
  m1.t, m1.y, m1.z,
  m2.t, m2.y, m2.z;
  C  << m0.t, m0.x, m0.z,
  m1.t, m1.x, m1.z,
  m2.t, m2.x, m2.z;
  D  << m0.t, m0.x, m0.y,
  m1.t, m1.x, m1.y,
  m2.t, m2.x, m2.y;
  // this orthogonal vector is for a Euclidean R4, which is sufficient
  Vector4 ortho = Vector4(1, 0, 0, 0) * A.determinant() -
  Vector4(0, 1, 0, 0) * B.determinant() +
  Vector4(0, 0, 1, 0) * C.determinant() -
  Vector4(0, 0, 0, 1) * D.determinant();
  // how does this tell me the signature?
  double x = sqr(ortho.x) + sqr(ortho.y) + sqr(ortho.z);
  double t = abs(ortho.t);
  tet.signature = t > x + eps ? 1 : (t < x-eps ? -1 : 0); // the opposite sign here is because we are using the orthogonal vector
  }*/
}

void SpaceTime::setTriangleEdgeMatrices()
{
  // Next set some links that are specific to Regge algebra
  // each triangle needs to know the two tetrahedra that border each pentachoron
  for (auto &pentaPair : pentachorons)
  {
    Pentachoron &penta = pentaPair.second;
    for (int j = 0; j <= 4; j++)
    {
      for (int k = 1; k <= 4; k++)
      {
        auto &e = lines.find(penta.corners[j]->pos + penta.corners[k]->pos);
        penta.edgeMatrix[j][k] = e == lines.end() ? NULL : &e->second;
        if (j == 0 && k == 4 && (penta.edgeMatrix[j][k] == NULL || abs(penta.edgeMatrix[j][k]->lengthSqr) < 1e-20))
        {
          cout << "bad edge 4 " << endl;
        }
      }
    }
  }
}

void SpaceTime::calculateEdgeLengths()
{
  for (auto &linePair : lines)
  {
    Line &edge = linePair.second;
    Vector4 dif = edge.ends[1]->pos - edge.ends[0]->pos;
    edge.lengthSqr = -sqr(dif.t) + sqr(dif.x) + sqr(dif.y) + sqr(dif.z);
  }
}

void SpaceTime::init()
{
  build();
  calculateEdgeLengths();
  connect();
  setTriangleEdgeMatrices();
  calculateSignatures();
}

