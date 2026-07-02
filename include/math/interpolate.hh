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

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <logger.hh>


/// @brief 1D interpolation class
/// @tparam logx static flag to indicate logarithmic interpolation in x
/// @tparam logy static flag to indicate logarithmic interpolation in y
/// @tparam periodic static flag to indicate periodic interpolation in x
template <bool logx, bool logy, bool periodic>
class interpolated_function_1d
{

private:

  double y_sign_ = 1.0; ///< common sign restored after logarithmic-y interpolation
  bool use_log_y_ = logy; ///< whether logarithmic-y interpolation is active for the current data
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
    data_x_ = data_x;
    data_y_ = data_y;

    assert(data_x_.size() == data_y_.size());
    assert(data_x_.size() > 5);

    if (logx) {
      for (auto &d : data_x_) {
        if (!std::isfinite(d) || d <= 0.0) {
          music::elog << "Log-x interpolation requires finite, positive data." << std::endl;
          throw std::invalid_argument("Log-x interpolation requires finite, positive data");
        }
        d = std::log(d);
      }
    }

    if (logy) {
      bool has_positive = false;
      bool has_negative = false;
      bool has_zero = false;

      for (const auto d : data_y_) {
        if (!std::isfinite(d)) {
          music::elog << "Log-y interpolation requires finite data." << std::endl;
          throw std::invalid_argument("Log-y interpolation requires finite data");
        }
        has_positive |= d > 0.0;
        has_negative |= d < 0.0;
        has_zero |= d == 0.0;
      }

      use_log_y_ = !has_zero && !(has_positive && has_negative);
      if (use_log_y_) {
        y_sign_ = has_negative ? -1.0 : 1.0;
        for (auto &d : data_y_)
          d = std::log(std::abs(d));
      } else {
        y_sign_ = 1.0;
        music::wlog
          << "Log-y interpolation received data that is not strictly single-signed; "
             "falling back to linear-y interpolation."
          << std::endl;
      }
    }

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
    return use_log_y_ ? y_sign_ * std::exp(y) : y;
  }
};
