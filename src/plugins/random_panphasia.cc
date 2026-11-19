// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2021 by Oliver Hahn and Adrian Jenkins (this file)
// but see distinct licensing for PANPHASIA below
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
//
// IMPORTANT NOTICE:
// Note that PANPHASIA itself is not released under the GPL. Make sure
// to read and agree to its distinct licensing before you use or modify
// the code below or in the /external/panphasia directory which can be
// found here: http://icc.dur.ac.uk/Panphasia.php
// NOTE THAT PANPHASIA REQUIRES REGISTRATION ON THIS WEBSITE PRIOR TO USE

#if defined(USE_PANPHASIA) && defined (USE_MPI)

#include <general.hh>
#include <random_plugin.hh>
#include <config_file.hh>

#include <vector>
#include <cmath>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <grid_fft.hh>


namespace PANPHASIA2
{
  extern "C"
  {
#include <panphasia_ho/panphasia_functions.h>
    extern size_t descriptor_base_size;
  }
}

#include "PANPHASIA.hh"

class RNG_panphasia : public RNG_plugin
{
private:
protected:
  std::string descriptor_string_;
  int num_threads_;
  int panphasia_mode_;
  size_t grid_res_;
  real_t boxlength_;

  PANPHASIA1::RNG *ppan1_rng_;

public:
  explicit RNG_panphasia(config_file &cf) : RNG_plugin(cf)
  {

#ifdef _OPENMP
    num_threads_ = omp_get_max_threads();
#else
    num_threads_ = 1;
#endif

    descriptor_string_ = pcf_->get_value<std::string>("random", "descriptor");
    grid_res_ = pcf_->get_value<size_t>("setup", "GridRes");
    boxlength_ = pcf_->get_value<real_t>("setup", "BoxLength");

    if( pcf_->get_value_safe<bool>("setup", "DoFixing", false) ){
      music::flog << "Fixing all the modes to the mean power negates any advantage of using the Panphasia field.\n";
      music::flog << "With the new panphasia ho it is possible by choosing the descriptor to fix the largest modes without losing the ability to resimulate to much higher resolution.\n";
      throw std::runtime_error("PANPHASIA: incompatible parameter.");
    }

    panphasia_mode_ = 0;
    PANPHASIA2::parse_and_validate_descriptor_(descriptor_string_.c_str(), &panphasia_mode_);

    if (panphasia_mode_ == 0)
    {
      ppan1_rng_ = new PANPHASIA1::RNG(&cf);
    }
  }

  ~RNG_panphasia()
  {
    if (panphasia_mode_ == 0)
    {
      delete ppan1_rng_;
    }
  }

  bool isMultiscale() const { return true; }

  void Run_Panphasia_Highorder(Grid_FFT<real_t> &g);

  void Fill_Grid(Grid_FFT<real_t> &g)
  {
    switch (panphasia_mode_)
    {

    case 0: // old mode
      music::ilog << "PANPHASIA: Old descriptor" << std::endl;
      ppan1_rng_->Fill(g);
      break;

    case 1: // PANPHASIA HO descriptor
      music::ilog << "PANPHASIA: New descriptor" << std::endl;
      this->Run_Panphasia_Highorder(g);
      break;

    default: // unknown PANPHASIA mode
      music::elog << "PANPHASIA: Something went wrong with descriptor" << std::endl;
      abort();
      break;
    }
  }
};

void RNG_panphasia::Run_Panphasia_Highorder(Grid_FFT<real_t> &g)
{
  int verbose = 0;
  int error;
  size_t x0 = 0, y0 = 0, z0 = 0;
  size_t rel_level;
  int fdim = 2; //Option to scale Fourier grid dimension relative to Panphasia coefficient grid

  //char descriptor[300] = "[Panph6,L20,(424060,82570,148256),S1,KK0,CH-999,Auriga_100_vol2]";

  PANPHASIA2::PANPHASIA_init_descriptor_(descriptor_string_.c_str(), &verbose);

  printf("Descriptor %s\n ngrid_load %lu\n", descriptor_string_.c_str(), grid_res_);

  // Choose smallest value of level to equal of exceed grid_res_)

  for (rel_level = 0; fdim * (PANPHASIA2::descriptor_base_size <<rel_level) <  grid_res_; rel_level++);

  printf("Setting relative level = %lu\n", rel_level);

  if ((error = PANPHASIA2::PANPHASIA_init_level_(&rel_level, &x0, &y0, &z0, &verbose)))
  {
    printf("Abort: PANPHASIA_init_level_ :error code %d\n", error);
    abort();
  };

  //======================= FFTW ==============================

  ptrdiff_t alloc_local, local_n0, local_0_start;

  size_t N0 = fdim * (PANPHASIA2::descriptor_base_size << rel_level);

  alloc_local = FFTW_MPI_LOCAL_SIZE_3D(N0, N0, N0 + 2, MPI_COMM_WORLD, &local_n0, &local_0_start);

  Grid_FFT<real_t> pan_grid({{N0, N0, N0}}, {{boxlength_, boxlength_, boxlength_}});

  _unused(alloc_local);

  assert(pan_grid.n_[0] == N0);
  assert(pan_grid.n_[1] == N0);
  assert(pan_grid.n_[2] == N0);
  assert(pan_grid.local_0_start_ == local_0_start);
  assert(pan_grid.local_0_size_ == local_n0);
  assert(pan_grid.ntot_ == size_t(alloc_local));

  pan_grid.FourierTransformForward(false);

  FFTW_COMPLEX *Panphasia_White_Noise_Field = reinterpret_cast<FFTW_COMPLEX *>(&pan_grid.data_[0]);

  // Panphasia_White_Noise_Field = FFTW_ALLOC_COMPLEX(alloc_local);

  if ((error = PANPHASIA2::PANPHASIA_compute_kspace_field_(rel_level, N0, local_n0, local_0_start, Panphasia_White_Noise_Field)))
  {
    music::elog << "Error code from PANPHASIA_compute ...  (ErrCode = " << error << ")" << std::endl;
  };

  pan_grid.FourierInterpolateCopyTo(g);
}
namespace
{
  RNG_plugin_creator_concrete<RNG_panphasia> creator("PANPHASIA");
}
#endif // defined(USE_PANPHASIA)
