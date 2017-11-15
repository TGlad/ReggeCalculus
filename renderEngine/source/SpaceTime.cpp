#include "SpaceTime.h"

static const int tMax = 10;
static const int wMax = 10;

static const int TMax = 2 * tMax - 1;
static const int WMax = 2 * wMax - 1;

bool SpaceTime::inBounds(const Vector4 &v)
{
  return v.t >= 0 && v.t < TMax &&
    v.x >= 0 && v.x < WMax &&
    v.y >= 0 && v.y < WMax &&
    v.z >= 0 && v.z < WMax;
}

void SpaceTime::init()
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
                              pentachorons[ev] = Pentachoron(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4]);
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

  // Next set some links that are specific to Regge algebra
  // each triangle needs to know the two tetrahedra that border each pentachoron
  for (auto &bonepair : triangles)
  {
    Triangle &bone = bonepair.second;
    bone.edgeMatrix.resize(bone.pentachorons.size());

    bone.tetraPairs[0].resize(bone.pentachorons.size());
    bone.tetraPairs[1].resize(bone.pentachorons.size());
    for (int i = 0; i < (int)bone.pentachorons.size(); i++)
    {
      Pentachoron *pent = bone.pentachorons[i];
#if 1
      int tetraIndex = 0;
      for (auto &tetra : pent->volumes)
      {
        for (auto &tr : tetra->faces)
        {
          if (tr == &bone)
          {
            ASSERT(tetraIndex < 2);
            bone.tetraPairs[tetraIndex++][i] = tetra;
          }
        }
      }
      ASSERT(tetraIndex == 2);
#else
      Vertex *vs[5] = { bone.corners[0], bone.corners[1], bone.corners[2], 0, 0 }; // vs[0] will always be even in all components or odd in all
      int iv = 3;
      int cornerIndex = 0;
      for (auto &cornerPair : pent->corners)
      {
        Vertex &corner = cornerPair.second;
        if (corner != bone.corners[0] && corner != bone.corners[1] && corner != bone.corners[2])
          vs[iv++] = corner;
      }
      for (int j = 0; j <= 4; j++)
      {
        for (int k = 1; k <= 4; k++)
        {
          auto &e = edges(vs[j].pos + vs[k].pos);
          bone.edgeMatrix[i].edges[j][k] = e == edges.end() ? NULL : &e.second;
        }
      }
#endif
    }
  }

  // Next: calculate signatures:
  double eps = 1e-5;
  for (auto &edgePair : lines)
  {
    Edge &edge = edgePair.second;
    Vector4 l = edge.corners[1].pos - edge.corners[0].pos;
    double S = sqrt(sqr(l.x) + sqr(l.y) + sqr(l.z));
    double T = abs(l.t);
    double eps = 1e-10;
    edge.signature = T > S+eps ? -1.0 : (T < S-eps ? 1 : 0);
  }
  for (auto &triPair : triangles)
  {
    Triangle &tri = triPair.second;
    tri.jacobian.resize(edges.size(), edges.size());
    double difs[3];
    for (int i = 0; i < 3; i++)
      difs[i] = abs(tri.corners[i].t - tri.corners[(i + 1) % 3].t);
    if (difs[0]>eps && difs[1] > eps && difs[2] > eps)
      tri.signature = -1; // temporal
    else
    {
      int a = difs[0] < difs[1] ? (difs[0] < difs[2] ? 0 : 2) : (difs[1] < difs[2] ? 1 : 2);
      Vector4 pa = tri.corners[a].pos;
      Vector4 pb = tri.corners[(a + 1) % 3].pos;
      Vector4 pc = tri.corners[(a + 2) % 3].pos;
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
  }
}

// TODO: I think the best approach is to switch to absolute sparse matrices throughout.

VectorXd SpaceTime::getEdgeErrors(SparseMatrix<double> &jacobian)
{
  for (auto &triPair : triangles)
  {
    Triangle &tri = triPair.second;
    tri.deficitAngle = getDeficitAngle(tri);

    // one facing edge
    int ia = tri->edges[0] == edge ? 0 : (tri->edges[1] == edge ? 1 : 2);
    Edge *e0i = tri->edges[(ia + 1) % 3];
    Edge *e0j = tri->edges[(ia + 2) % 3];
    Edge *eij = tri->edges[ia % 3];
    double sij = eij->lengthSqr;
    double s0i = e0i->lengthSqr;
    double s0j = e0j->lengthSqr;

    // below: Hartle '84 eq 3.7
    tri.areaSquared = (s0i*s0j - 0.25*sqr(s0i + s0j - sij)) / 4.0; // correct
    // TODO: I need to make these partial derivatives absolute by using a sparse vector
    // Hartle '84 eq 3.13. I think this is right:
    tri.areaSquaredDot(e0i->index) = (1.0 / 8.0) * (s0j + sij - s0i);
    tri.areaSquaredDot(e0j->index) = (1.0 / 8.0) * (s0i + sij - s0j);
    tri.areaSquaredDot(eij->index) = (1.0 / 8.0) * (s0i + s0j - sij); 
 //   tri.areaSquaredDot = (s0i + s0j - sij) / 8.0;  correct 
    tri.areaDot = tri.areaSquaredDot / (2.0 * sqrt(tri.areaSquared)); // correct
    SparseMatrix<double, edges.size(), edges.size()> areaSquaredDotDot;

    areaSquaredDotDot(e0i->index, e0j->index) = 1; areaSquaredDotDot(e0i->index, eij->index) = 1; areaSquaredDotDot(e0i->index, e0i->index) = -1;
    areaSquaredDotDot(e0j->index, e0i->index) = 1; areaSquaredDotDot(e0j->index, eij->index) = 1; areaSquaredDotDot(e0j->index, e0j->index) = -1;
    areaSquaredDotDot(eij->index, e0i->index) = 1; areaSquaredDotDot(eij->index, e0j->index) = 1; areaSquaredDotDot(eij->index, eij->index) = -1;
    areaSquaredDotDot *= 1.0 / 8.0;
    // TODO: is this transpose correct?
    tri.areaDotDot = (areaSquaredDotDot * sqrt(tri.areaSquared) - tri.areaSquaredDot * tri.areaDot.transpose()) / (2.0*tri.areaSquared);
  }

  VectorXd edgeErrors(edges.size());
  int index = 0;
  for (auto &edgePair : edges)
  {
    Edge &edge = edgePair.second;
    double sum = 0;
    SparseMatrix<double, edges.size(), edges.size()> jacobianSum;
    jacobianSum.setZero;
    double scale = 0.5 * sqrt(edge->lengthSqr);
    for (auto &tri : edge->triangles)
    {
      sum += tri.deficitAngle * tri.areaDot;
      jacobianSum += tri.deficitAngle * tri.areaDotDot + tri.deficitAngleDot * tri.areaDot.transpose();
    }
    ASSERT(index++ == edge.index);
    edgeErrors[edge.index] = sum; // could subtract some value for a simple mass component
    jacobian += jacobianSum;
  }
  return edgeErrors;
}

double SpaceTime::getDeficitAngle(const Triangle &bone) const
{
  double sum = 0;
  bone.deficitAngleDot.setZero();
  for (int p = 0; p < (int)bone.tetraPairs.size(); p++)
  {
    // in this paper the bone triangle has indices 0,1,2 and the two remaining vertices of this pentachoron attached to the bone are 3,4
    // this means I need edges 0 to i in 1..4, and i to j in 1..4, so I don't need tetrapairs, I need the two other points
    Vertex *vs[5] = { bone.corners[0], bone.corners[1], bone.corners[2], bone.cornerPairs[p][0], bone.cornerPairs[p][1] };
    Tetrahedron &tetra = tetraPairs[p][0];

    double s[5][5]; // square edge lengths
    for (int i = 0; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        s[i][j] = bone.edgeMatrix[p].edges[i][j] ? bone.edgeMatrix[p].edges[i][j]->lengthSqr : 0;
    double g[5][5]; // we don't use the 0 index, to match Brewin
    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        g[i][j] = 0.5*(s[0][i] + s[0][j] - s[i][j]);
    Vector4 n3 = -sign(g[4][4]) * Vector4(g[4][1], g[4][2], g[4][3], g[4][4]) / sqrt(abs(g[3][3]));
    double h[4][4];

    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        h[i][j] = g[i][j] - (g[4][i] * g[4][j]) / g[4][4];
    Vector4 m3 = -sign(h[3][3]) * Vector4(h[4][1], h[4][2], h[4][3], h[4][4]) / sqrt(abs(h[3][3]));
    // m4 just symmetric to the m3 block above
    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        h[i][j] = g[i][j] - (g[3][i] * g[3][j]) / g[3][3];
    Vector4 m4 = -sign(h[4][4]) * Vector4(h[3][1], h[3][2], h[3][3], h[3][4]) / sqrt(abs(h[4][4]));

    double phi34;
    double m3m4 = m3.dot(m4);
    if (bone.signature == -1)
      phi34 = acos(m3m4);
    else
    {
      double n3m4 = n3.dot(m4);;
      double rho12 = abs(m3m4) < abs(n3m4) ? sign(n3m4)*m3m4 : sign(m3m4)*n3m4;
      phi34 = sign(m3.dot(m3)) * asinh(rho12);
    }
    sum += phi34;

    // Jacobian part
    SparseVector<double> mn(edges.size());
    for (int u = 0; u < 4; u++)
    {
      for (int v = 0; v < 4; v++)
      {
        double scale = (m3[u] * n3[v] + m4[u] * n4[v]);
        Edge *edge0 = bone.edgeMatrix[p].edges[0][u + 1];
        if (edge0)
          mn(edge0->index) += 0.5*scale; // dguv/ds = 0.5
        Edge *edge1 = bone.edgeMatrix[p].edges[0][v + 1];
        if (edge1)
          mn(edge1->index) += 0.5*scale; // dguv/ds = 0.5
        Edge *edge2 = bone.edgeMatrix[p].edges[u+1][v + 1];
        if (edge2)
          mn(edge2->index) -= 0.5*scale; // dguv/ds = -0.5
      }
    }
    mn *= 0.5*tetra->signature * bone->signature;
    
    bone.deficitAngleDot += mn;
  }

  return bone.signature == -1 ? 2.0*pi - sum : -sum;
}

void SpaceTime::update()
{
  SparseMatrix<double> jacobian(edges.size(), edges.size());
  for (int i = 0; i < 10; i++) // multiple Newton iterations
  {
    VectorXd error = getEdgeErrors(jacobian);
    SparseLU<SparseMatrix<double> > solver;
    solver.compute(Jacobian);
    VectorXd deltaLengthSqrs = solver.solve(-error);
    int index;
    for (auto &edgePair : lines)
      edgePair.second->lengthSqr += deltaLengthSqrs[index++];
  }
}

// correct processing should maintain certain inequalities. See Khavari ch 5.
void SpaceTime::validate()
{
#if 0
  for (auto &triPair : triangles)
  {
    Triangle &tri = triPair.second;
    if (TTT || SSS)
    {
      double maxS = max(tri.edges[0].s, max(tri.edges[1].s, tri.edges[2].s));
      double sum = 0;
      for (auto &edge: tri.edges)
        if (edge.s != maxS)
          sum += sqrt(edge.s);
      ASSERT(sqrt(maxS) > sum);
    }
    if (SST)
    {
      if (pure timelike twin)
        ASSERT(b > c + a);
      else
        ASSERT(c < b + a);
    }
  }
#endif
}