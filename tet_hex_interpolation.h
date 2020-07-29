#ifndef TET_HEX_INTERPOLATION_JJ_H
#define TET_HEX_INTERPOLATION_JJ_H


#include <Eigen/Dense>
#include <Eigen/Geometry> 
#include <Eigen/Core>
#include <iostream>

using namespace std;

namespace Jing
{
  //! \brief calculate global coordinates responding to local according to FEM.
  //! \param[in] verts is in dimension 8*3
  //! \param[in] local parameter domain coordinates
  //! \param[out] p coordinates
  template<typename T>
  Eigen::Matrix<T, -1, 1> GetGlobalCoordinatesInHex(const Eigen::Matrix<T, -1, 1> &local,
                                                    const Eigen::Matrix<T, -1, -1> &verts)
  {
    assert(verts.cols() == 3 && verts.rows() == 8);
    T epsilon[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    T eta[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    T phi[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

    T n[8];
    for (int i = 0; i < 8; ++i)
      n[i] = (1+local[0]*epsilon[i]) * (1+local[1]*eta[i]) * (1+local[2]*phi[i]) / 8.0;

    Eigen::Matrix<T, -1, 1> p = Eigen::Matrix<T, -1, 1>::Zero(3);
    for (int i = 0; i < 8; ++i)
      p += n[i] * verts.row(i);

    return p;
  }

  //! \brief calcualte parameter domain coordinates according to Inverse Isoparametric(Efficient Inverse Isoparametric Mapping Algorithm for Whole-Body Computed Tomography Registration Using Deformations Predicted by Nonlinear Finite Element Modeling)
  //! \para[in] verts is in dimension 8*3
  //! \param[out] parameter domain coordinates
  template<typename T>
  std::tuple<bool, Eigen::Matrix<T, -1, 1>> GetParameterCoordinatesInHex(
    const Eigen::Matrix<T, -1, 1> &global,
    const Eigen::Matrix<T, -1, -1> &verts)
  {
    assert(verts.cols() == 3 && verts.rows() == 8);
    for (int a = 0; a < 3; ++a)
    {
      T m = verts.col(a).minCoeff();
      T n = verts.col(a).maxCoeff();
      if (global[a] < 1.5*m-0.5*n || global[a] > 1.5*n-0.5*m)
      {
        printf("[ \033[1;35mwarning\033[0m ] interpolation points are outside the hex bounding box!\n");
        return std::make_tuple(false, Eigen::Matrix<T, -1, 1>::Zero(3));
      } 
    }

    Eigen::Matrix<T, -1, 1> local(3);
    local << 0, 0, 0;
    Eigen::Matrix<T, -1, 1> p = GetGlobalCoordinatesInHex(local, verts);
    Eigen::Matrix<T, -1, -1> coeff(3, 8);

    T epsilon[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    T eta[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    T phi[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    int count = 0;

    T leps;
    if (std::is_same<T, float>::value)
      leps = 1e-4;
    else
      leps = 1e-8;
  
    while ((p-global).norm() > leps)
    {
      for (int i = 0; i < 8; ++i)
        coeff(0, i) = (1+local[1]*eta[i]) * (1+local[2]*phi[i]) * epsilon[i] / 8.0;
      for (int i = 0; i < 8; ++i)
        coeff(1, i) = (1+local[0]*epsilon[i]) * (1+local[2]*phi[i]) * eta[i] / 8.0;
      for (int i = 0; i < 8; ++i)
        coeff(2, i) = (1+local[0]*epsilon[i]) * (1+local[1]*eta[i]) * phi[i] / 8.0;

      Eigen::Matrix<T, -1, -1> Jaccobi = coeff * verts;
      Eigen::FullPivLU<Eigen::Matrix<T, -1, -1>> lu(Jaccobi);
      local += lu.solve(global - p);
      p = GetGlobalCoordinatesInHex(local, verts);
      if (count > 100)
      {
        printf("[ \033[1;35mwarning\033[0m ] warning in local coordinates calculation!\n");
        return std::make_tuple(false, Eigen::Matrix<T, -1, 1>::Zero(3));
      }
      ++count;
    }

    printf("[   \033[1;37mlog\033[0m   ] parameter coodinates convergence in %d iteration!\n", count);
    return std::make_tuple(true, local);
  }

  //! \brief get interpolate coefficient in hex. inverse parameter will be returned if global is inside hex, otherwise least norm interpolate coefficient will be return.
  //! \param[in] verts in dimension 8*3
  template<typename T>
  Eigen::Matrix<T, -1, 1> GetHexInterpolateCoeff(const Eigen::Matrix<T, -1, 1> &global,
                                                 const Eigen::Matrix<T, -1, -1> &verts)

  {
    assert(verts.cols() == 3 && verts.rows() == 8);
    bool is_parameter = false;
    Eigen::Matrix<T, -1, 1> local(3);
    std::tie(is_parameter, local) = GetParameterCoordinatesInHex(global, verts);
    if (is_parameter)
    {
      T epsilon[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
      T eta[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
      T phi[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

      Eigen::Matrix<T, -1, 1> n(8);
      for (int i = 0; i < 8; ++i)
        n[i] = (1+local[0]*epsilon[i]) * (1+local[1]*eta[i]) * (1+local[2]*phi[i]) / 8.0;

      return n;
    }
    else
    {
      printf("[ \033[1;35mwarning\033[0m ] extrapolation in hex!\n");
//! \note RowMajor must be added because verts is ColMajor
      Eigen::Map<const Eigen::Matrix<T, -1, -1, Eigen::RowMajor>>
        p(verts.data(), verts.cols(), verts.rows());
      Eigen::Matrix<T, -1, -1> AA = p * p.transpose();
#ifdef DEBUG
      Eigen::FullPivLU<Eigen::Matrix<T, -1, -1>> lu(AA);
      assert(lu.rank() == 3);
#endif
      Eigen::Matrix<T, -1, -1> AAI = AA.inverse();
      Eigen::Matrix<T, -1, 1> n = p.transpose() * AAI * (global);
      return n;
    }

    return Eigen::Matrix<T, -1, 1>::Zero(8);
  }

  //! \param[in] verts in dimension 4*3
  template<typename T>
  Eigen::Matrix<T, -1, 1> GetTetInterpolateCoeff(const Eigen::Matrix<T, -1, 1> &global,
                                                 const Eigen::Matrix<T, -1, -1> &verts)
  {
    assert(verts.cols() == 3 && verts.rows() == 4);
    Eigen::Map<const Eigen::Matrix<T, -1, -1, Eigen::RowMajor>>
      p(verts.data(), verts.cols(), verts.rows()); //! \note RowMajor must be added because verts is ColMajor
    Eigen::Matrix<T, -1, -1> D = Eigen::Matrix<T, -1, -1>::Ones(4, 4);
    D.block(0, 0, 3, 4) = p;
    Eigen::Matrix<T, -1, 1> b = Eigen::Matrix<T, -1, 1>::Ones(4); 
    b.block(0, 0, 3, 1) = global;

    Eigen::Matrix<T, -1, 1> n = D.inverse() * b;
    return n;
  }
  
}



#endif // TET_HEX_INTERPOLATION_JJ_H
