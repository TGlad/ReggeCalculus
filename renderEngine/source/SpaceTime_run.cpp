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

