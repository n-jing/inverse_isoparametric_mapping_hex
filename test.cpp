#include <iostream>
#include <Eigen/Core>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Geometry> 

#include "tet_hex_interpolation.h"

using namespace std;
using namespace Eigen;

using FLOAT = double;

using namespace Jing;


int main (int argc, char *argv[])
{
#ifdef DEBUG_BUILD
#error
#endif
  Matrix<FLOAT, -1, -1> verts(8, 3);
  verts << -1, -1, -1,
    1, -1, -1,
    1, 1, -1,
    -1, 1, -1,
    -1, -1, 1,
    1, -1, 1,
    1, 1, 1,
    -1, 1, 1;

  verts += Matrix<FLOAT, -1, -1>::Constant(8, 3, 3);
  verts *= 2;
  verts += Matrix<FLOAT, -1, -1>::Random(8, 3);
  default_random_engine e;
  uniform_real_distribution<double> u(2, 10);
  Matrix<FLOAT, -1, 1> global(3);
  for (int i = 0; i < 200; ++i)
  {
    global << u(e), u(e), u(e);
    Matrix<FLOAT, -1, 1> l = GetHexInterpolateCoeff(global, verts);
    cerr << global.transpose() << endl;
    cerr << l.transpose() << endl;
    cerr << (verts.transpose() * l).transpose() << endl;
    cerr << (global - verts.transpose()*l).norm() << endl;
    assert((global - verts.transpose()*l).norm() < 1e-6);
  }


  Matrix<FLOAT, -1, -1> tet(4, 3);
  tet << -1, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1;
  Matrix<FLOAT, -1, 1> q(3);
  for (int i = 0; i < 300; ++i)
  {
    q << u(e), u(e), u(e);
    Matrix<FLOAT, -1, 1> p = GetTetInterpolateCoeff(q, tet);
    assert((q-tet.transpose() * p).norm() < 1e-8);
    cerr << q.transpose() << endl;
    cerr << (tet.transpose() * p).transpose() << endl;
  }

  return 0;
}
