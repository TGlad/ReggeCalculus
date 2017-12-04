#include "SpaceTime.h"
#include "/Code/Eigen/SparseLU"


#define TEST
VectorXd SpaceTime::getEdgeErrors(SparseMatrix<double> &jacobian)
{
  for (auto &triPair : triangles)
  {
    Triangle &tri = triPair.second;
    tri.deficitAngle = tri.getDeficitAngle(lines.size());
    if (!(tri.deficitAngle == tri.deficitAngle) || (tri.deficitAngle && abs(tri.deficitAngle) < 1e-20))
      cout << "blah" << endl;

    tri.calcAreaDot(lines.size());
    tri.calcAreaDotDot(lines.size());
  }
#if defined(TEST)
  // we know change in angle with change in lengths,
  // so randomly change lengths, calculate suposed change and compare to actual
  SparseVector<double> change(lines.size());
  const double eps = 0.01;
  for (int i = 0; i < change.size(); i++)
    change.coeffRef(i) = random(-eps, eps);
  int i = 0;
  for (auto &edgePair : lines)
  {
    Line &edge = edgePair.second;
    edge.lengthSqr += change.coeff(i);
    i++;
  }
  double totalError = 0;
  double totalValue = 0;
  for (auto &triPair : triangles)
  {
    Triangle &tri = triPair.second;
#if defined(TEST_DEFICIT)
    double expectedAngleChange = tri.deficitAngleDot.dot(change);
    double angleChange = tri.getDeficitAngle(lines.size()) - tri.deficitAngle;
    double error = abs(angleChange - expectedAngleChange);
    if (error > 0.1*eps)
      cout << "bad deficit change prediction for triangle " << i << ", expected: " << expectedAngleChange << ", actual: " << angleChange << endl;
    totalValue += abs(angleChange);
#else // test angle dot
    SparseVector<double> expectedAreaDotChange = tri.areaDotDot * change;
    SparseVector<double> areaDot = tri.areaDot;
    tri.calcAreaDot(lines.size());
    SparseVector<double> areaDotChange = tri.areaDot - areaDot;
    double error = (areaDotChange - expectedAreaDotChange).norm();
    if (error > 0.1*eps)
      cout << "bad area dot prediction for triangle " << i << ", expected: " << expectedAreaDotChange << ", actual: " << areaDotChange << endl;
    totalValue += areaDot.norm();
#endif
    totalError += error;
  }
  if (totalError > totalValue * 0.001)
    cout << "bad overall prediction, error " << totalError << ", total value: " << totalValue << endl;
#endif

  VectorXd edgeErrors(lines.size());
  for (auto &edgePair : lines)
  {
    Line &edge = edgePair.second;
    double sum = 0;
    SparseMatrix<double> jacobianSum(lines.size(), lines.size());
    double scale = 0.5 * sqrt(edge.lengthSqr);
    for (Triangle *tri : edge.triangles)
    {
      sum += tri->deficitAngle * tri->areaDot.coeff(edge.index);
      jacobianSum += tri->deficitAngle * tri->areaDotDot + tri->deficitAngleDot * tri->areaDot.transpose();
    }
    edgeErrors[edge.index] = sum; // could subtract some value for a simple mass component
    jacobian += jacobianSum;
  }
  return edgeErrors;
}

void Triangle::calcAreaDot(int numLines)
{
  areaDot.resize(numLines);
  areaSquaredDot.resize(numLines);
  Line *e0i = edges[0];
  Line *e0j = edges[1];
  Line *eij = edges[2];
  double sij = eij->lengthSqr;
  double s0i = e0i->lengthSqr;
  double s0j = e0j->lengthSqr;

  // below: Hartle '84 eq 3.7
  areaSquared = (s0i*s0j - 0.25*sqr(s0i + s0j - sij)) / 4.0; // correct
  if (abs(areaSquared) < 1e-20)
    cout << "blah" << endl;
  // TODO: I need to make these partial derivatives absolute by using a sparse vector
  // Hartle '84 eq 3.13. I think this is right:
  areaSquaredDot.coeffRef(e0i->index) = (1.0 / 8.0) * (s0j + sij - s0i);
  areaSquaredDot.coeffRef(e0j->index) = (1.0 / 8.0) * (s0i + sij - s0j);
  areaSquaredDot.coeffRef(eij->index) = (1.0 / 8.0) * (s0i + s0j - sij);
  
  areaDot = areaSquaredDot / (2.0 * sqrt(areaSquared)); // correct
}

void Triangle::calcAreaDotDot(int numLines)
{
  Line *e0i = edges[0];
  Line *e0j = edges[1];
  Line *eij = edges[2];
  SparseMatrix<double> areaSquaredDotDot(numLines, numLines);
  areaSquaredDotDot.coeffRef(e0i->index, e0j->index) = 1; areaSquaredDotDot.coeffRef(e0i->index, eij->index) = 1; areaSquaredDotDot.coeffRef(e0i->index, e0i->index) = -1;
  areaSquaredDotDot.coeffRef(e0j->index, e0i->index) = 1; areaSquaredDotDot.coeffRef(e0j->index, eij->index) = 1; areaSquaredDotDot.coeffRef(e0j->index, e0j->index) = -1;
  areaSquaredDotDot.coeffRef(eij->index, e0i->index) = 1; areaSquaredDotDot.coeffRef(eij->index, e0j->index) = 1; areaSquaredDotDot.coeffRef(eij->index, eij->index) = -1;
  areaSquaredDotDot *= 1.0 / 8.0;
  // TODO: is this transpose correct?
  areaDotDot = (areaSquaredDotDot * sqrt(areaSquared) - areaSquaredDot * areaDot.transpose()) / (2.0*areaSquared);
}

double Triangle::getDeficitAngle(int numLines)
{
  double sum = 0;
  deficitAngleDot.resize(numLines);
  static int xx = 0;
  for (int p = 0; p < (int)pentachorons.size(); p++)
  {
    Pentachoron *penta = pentachorons[p];
    double s[5][5]; // square edge lengths
    for (int i = 0; i <= 4; i++)
      for (int j = 0; j <= 4; j++)
        s[i][j] = penta->edgeMatrix[i][j] ? penta->edgeMatrix[i][j]->lengthSqr : 0; // error is in this edgeMatrix
//#if defined(BREWIN2011)
    Matrix<double, 4, 4> gdown;
    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        gdown(i-1,j-1) = 0.5*(s[0][i] + s[0][j] - s[i][j]); // g[4][4] will be zero whenever edge [0][4] is timelike...  so h will all be undef
    Matrix<double, 4, 4> gd = penta->M * gdown * penta->M.transpose();
    
    Matrix<double, 4, 4> gup = gd.inverse();
    double gu[5][5];
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++)
      gu[i + 1][j + 1] = gup(i, j);


    if (abs(gu[4][4]) < 1e-20 || abs(gu[3][3]) < 1e-20)
    {
      cout << "bad g metric" << endl;
      ASSERT(false);
    }
    Vector4 n3u = -sign(gu[4][4]) * Vector4(gu[4][1], gu[4][2], gu[4][3], gu[4][4]) / sqrt(abs(gu[3][3]));
    Vector4 n4u = -sign(gu[3][3]) * Vector4(gu[3][1], gu[3][2], gu[3][3], gu[3][4]) / sqrt(abs(gu[4][4])); // TODO: is this correct?

    double hu[5][5];
    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        hu[i][j] = gu[i][j] - (gu[4][i] * gu[4][j]) / gu[4][4];
    Vector4 m3u = sign(hu[3][3]) * Vector4(hu[3][1], hu[3][2], hu[3][3], hu[3][4]) / sqrt(abs(hu[3][3]) + 1e-10);
    // m4 just symmetric to the m3 block above
    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        hu[i][j] = gu[i][j] - (gu[3][i] * gu[3][j]) / gu[3][3];
    Vector4 m4u = sign(hu[4][4]) * Vector4(hu[4][1], hu[4][2], hu[4][3], hu[4][4]) / sqrt(abs(hu[4][4]) + 1e-10); // TODO: is this correct?
    
    Vector4 m4d = gd * m4u;
    Vector4 m3d = gd * m3u;

    double m3m4 = m3u.dot(m4d); // TODO: try the Hartle eq 3.9 version of this m3m4 value, and see what result we get
    double n3m4 = n3u.dot(m4d);
//    if (((int)penta->corners[0]->pos.t) % 2)
//      cout << "middle" << endl;
#if 1//else // HARTLE85
    Matrix<double, 3, 3> WaWb, Va, Vb;
    int as[] = { 1, 2, 3 };
    int bs[] = { 1, 2, 4 };
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        WaWb(i, j) = 0.5*(s[0][as[i]] + s[0][bs[j]] - s[as[i]][bs[j]]);
        Va(i, j) = 0.5*(s[0][as[i]] + s[0][as[j]] - s[as[i]][as[j]]);
        Vb(i, j) = 0.5*(s[0][bs[i]] + s[0][bs[j]] - s[bs[i]][bs[j]]);
      }
    }
    double wawb = WaWb.determinant() / 36.0;
    double va = Va.determinant() / 36.0;
    double vb = Vb.determinant() / 36.0;
    double m3m42 = wawb / sqrt(va*vb + 1e-10);   // cosAngle // seems more likely to be correct here...
//    if (m3m42 != m3m4)
//      cout << "mismatch" << endl;
    double n3m42 = 0; // TODO: how do we calculate this??
#define calc_sine_angle
#if defined(calc_sine_angle)
    Matrix<double,2,2> Vbone;
    Matrix<double, 4, 4> Vpent;
    for (int i = 1; i <= 2; i++)
      for (int j = 1; j <= 2; j++)
        Vbone(i-1, j-1) = 0.5*(s[0][i] + s[0][j] - s[i][j]);
    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
        Vpent(i-1, j-1) = 0.5*(s[0][i] + s[0][j] - s[i][j]);
    double vbone = Vbone.determinant() / 4.0;
    double vpent = Vpent.determinant() / sqr(24);
    double sineAngle = (4.0 / 3.0) * sqrt(vbone * vpent) / sqrt(va*vb + 1e-10); // Note that this and cosAngle equate (so are correct) for random values
#endif
#endif
    double phi34;
    if (signature == -1)
      phi34 = acos(m3m4);
    else
    {
      double rho12 = abs(m3m4) < abs(n3m4) ? sign(n3m4)*m3m4 : sign(m3m4)*n3m4;
      phi34 = sign(m3u.dot(m3d)) * asinh(rho12);
    }
    if (!(phi34 == phi34) || (phi34 && abs(phi34) < 1e-20))
      cout << "blah" << endl;
    sum += phi34;
    xx++;

    // Jacobian part
    SparseVector<double> mn(numLines);
    for (int u = 0; u < 4; u++)
    {
      for (int v = 0; v < 4; v++)
      {
        double scale = (m3u[u] * n3u[v] + m4u[u] * n4u[v]);
        Line *edge0 = penta->edgeMatrix[0][u + 1];
        if (edge0)
          mn.coeffRef(edge0->index) += 0.5*scale; // dguv/ds = 0.5
        Line *edge1 = penta->edgeMatrix[0][v + 1];
        if (edge1)
          mn.coeffRef(edge1->index) += 0.5*scale; // dguv/ds = 0.5
        Line *edge2 = penta->edgeMatrix[u+1][v + 1];
        if (edge2)
          mn.coeffRef(edge2->index) -= 0.5*scale; // dguv/ds = -0.5
      }
    }
    mn *= 0.5 * signature;
    
    deficitAngleDot += mn;
  }

  return signature == -1 ? 2.0*pi - sum : -sum;
}

void SpaceTime::update()
{
  SparseMatrix<double> jacobian(lines.size(), lines.size());
  for (int i = 0; i < 10; i++) // multiple Newton iterations
  {
    VectorXd error = getEdgeErrors(jacobian);
    SparseLU<SparseMatrix<double> > solver;
    solver.compute(jacobian);
    VectorXd deltaLengthSqrs = solver.solve(-error);
    int index = 0;
    for (auto &edgePair : lines)
      edgePair.second.lengthSqr += deltaLengthSqrs[index++];
  }
}

