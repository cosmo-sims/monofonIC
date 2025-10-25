// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <vector>
#include <cassert>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>


/// @brief 1D interpolation class
/// @tparam logx static flag to indicate logarithmic interpolation in x
/// @tparam logy static flag to indicate logarithmic interpolation in y
/// @tparam periodic static flag to indicate periodic interpolation in x
template <bool logx, bool logy, bool periodic>
class interpolated_function_1d
{

private:

  double min_val_ = 0. ; ///< minimum value of y (subtracted and then readded to avoid numerical issues) : TOMA
  bool isinit_; ///< flag to indicate whether the interpolation has been initialized
  std::vector<double> data_x_, data_y_; ///< data vectors
  gsl_interp_accel *gsl_ia_; ///< GSL interpolation accelerator
  gsl_spline *gsl_sp_; ///< GSL spline object

  /// @brief deallocate GSL objects
  void deallocate()
  {
    gsl_spline_free(gsl_sp_);
    gsl_interp_accel_free(gsl_ia_);
  }

public:

  /// @brief default copy constructor (deleted)
  interpolated_function_1d(const interpolated_function_1d &) = delete;

  /// @brief empty constructor (without data)
  interpolated_function_1d() : isinit_(false){}

  /// @brief constructor with data
  /// @param data_x x data vector
  /// @param data_y y data vector
  interpolated_function_1d(const std::vector<double> &data_x, const std::vector<double> &data_y)
  : isinit_(false)
  {
    static_assert(!(logx & periodic),"Class \'interpolated_function_1d\' cannot both be periodic and logarithmic in x!");
    this->set_data( data_x, data_y );
  }

  /// @brief destructor
  ~interpolated_function_1d()
  {
    if (isinit_) this->deallocate();
  }

  /// @brief set data
  /// @param data_x x data vector
  /// @param data_y y data vector
  void set_data(const std::vector<double> &data_x, const std::vector<double> &data_y)
  {
    // TOMA
    min_val_ = *std::min_element(data_y.begin(), data_y.end()) -1.0;

    data_x_ = data_x;
    data_y_ = data_y;

    for (size_t i=0; i< data_y.size(); i++) //TOMA
    {
        data_y_[i] = data_y_[i] - min_val_;
    }

    
    assert(data_x_.size() == data_y_.size());
    assert(data_x_.size() > 5);

    if (logx) for (auto &d : data_x_) d = std::log(d);
    if (logy) for (auto &d : data_y_) d = std::log(d);

    if (isinit_) this->deallocate();

    gsl_ia_ = gsl_interp_accel_alloc();
    gsl_sp_ = gsl_spline_alloc(periodic ? gsl_interp_cspline_periodic : gsl_interp_cspline, data_x_.size());
    gsl_spline_init(gsl_sp_, &data_x_[0], &data_y_[0], data_x_.size());

    isinit_ = true;
  }

  /// @brief evaluate the interpolation
  /// @param x x value
  /// @return y value
  double operator()(double x) const noexcept
  {
    assert( isinit_ && !(logx&&x<=0.0) );
    const double xa = logx ? std::log(x) : x;
    const double y(gsl_spline_eval(gsl_sp_, xa, gsl_ia_));
    return logy ? std::exp(y) + min_val_ : y + min_val_;
    //return logy ? std::exp(y) : y;
  }
};
