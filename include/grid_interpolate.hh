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

#include <array>
#include <vector>

#include <general.hh>

#include <math/vec3.hh>

template <int interp_order, typename grid_t>
struct grid_interpolate
{
  using data_t = typename grid_t::data_t;
  using vec3 = std::array<real_t, 3>;

  static constexpr bool is_distributed_trait = grid_t::is_distributed_trait;
  static constexpr int interpolation_order = interp_order;

  std::vector<data_t> boundary_;
  std::vector<int> local0starts_;
  const grid_t &gridref;
  size_t nx_, ny_, nz_;

  explicit grid_interpolate(const grid_t &g)
      : gridref(g), nx_(g.n_[0]), ny_(g.n_[1]), nz_(g.n_[2])
  {
    static_assert(interpolation_order >= 0 && interpolation_order <= 2, "Interpolation order needs to be 0 (NGP), 1 (CIC), or 2 (TSC).");

    if (is_distributed_trait)
    {
      update_ghosts( g );
    }
  }

  void update_ghosts( const grid_t &g )
  {
  #if defined(USE_MPI)

    int local_0_start = int(gridref.local_0_start_);
    local0starts_.assign(MPI::get_size(), 0);

    MPI_Allgather(&local_0_start, 1, MPI_INT, &local0starts_[0], 1, MPI_INT, MPI_COMM_WORLD);

    //... exchange boundary
    size_t nx = interpolation_order + 1;
    size_t ny = g.n_[1];
    size_t nz = g.n_[2];

    boundary_.assign(nx * ny * nz, data_t{0.0});

    for (size_t i = 0; i < nx; ++i)
    {
      for (size_t j = 0; j < ny; ++j)
      {
        for (size_t k = 0; k < nz; ++k)
        {
          boundary_[(i * ny + j) * nz + k] = g.relem(i, j, k);
        }
      }
    }

    int sendto = (MPI::get_rank() + MPI::get_size() - 1) % MPI::get_size();
    int recvfrom = (MPI::get_rank() + MPI::get_size() + 1) % MPI::get_size();

    MPI_Status status;
    status.MPI_ERROR = MPI_SUCCESS;

    int err = MPI_Sendrecv_replace(&boundary_[0], nx * ny * nz, MPI::get_datatype<data_t>(), sendto,
                          MPI::get_rank() + 1000, recvfrom, recvfrom + 1000, MPI_COMM_WORLD, &status);

    if( err != MPI_SUCCESS ){
      char errstr[256]; int errlen=256;
      MPI_Error_string(err, errstr, &errlen ); 
      music::elog << "MPI_ERROR #" << err << " : " << errstr << std::endl;
    }
#endif
  }

  data_t get_ngp_at(const std::array<real_t, 3> &pos, std::vector<data_t> &val) const noexcept
  {
    size_t ix = static_cast<size_t>(pos[0]);
    size_t iy = static_cast<size_t>(pos[1]);
    size_t iz = static_cast<size_t>(pos[2]);
    return gridref.relem(ix - gridref.local_0_start_, iy, iz);
  }

  data_t get_cic_at(const std::array<real_t, 3> &pos) const noexcept
  {
    size_t ix = static_cast<size_t>(pos[0]);
    size_t iy = static_cast<size_t>(pos[1]);
    size_t iz = static_cast<size_t>(pos[2]);
    real_t dx = pos[0] - real_t(ix), tx = 1.0 - dx;
    real_t dy = pos[1] - real_t(iy), ty = 1.0 - dy;
    real_t dz = pos[2] - real_t(iz), tz = 1.0 - dz;
    size_t iy1 = (iy + 1) % ny_;
    size_t iz1 = (iz + 1) % nz_;

    data_t val{0.0};
    
    if( is_distributed_trait ){
      ptrdiff_t localix = ix-gridref.local_0_start_;
      val += gridref.relem(localix, iy, iz) * tx * ty * tz;
      val += gridref.relem(localix, iy, iz1) * tx * ty * dz;
      val += gridref.relem(localix, iy1, iz) * tx * dy * tz;
      val += gridref.relem(localix, iy1, iz1) * tx * dy * dz;

      if( localix+1 >= gridref.local_0_size_ ){
        size_t localix1 = localix+1 - gridref.local_0_size_;
        val += boundary_[(localix1*ny_+iy)*nz_+iz] * dx * ty * tz;
        val += boundary_[(localix1*ny_+iy)*nz_+iz1] * dx * ty * dz;
        val += boundary_[(localix1*ny_+iy1)*nz_+iz] * dx * dy * tz;
        val += boundary_[(localix1*ny_+iy1)*nz_+iz1] * dx * dy * dz;
      }else{
        size_t localix1 = localix+1;
        val += gridref.relem(localix1, iy, iz) * dx * ty * tz;
        val += gridref.relem(localix1, iy, iz1) * dx * ty * dz;
        val += gridref.relem(localix1, iy1, iz) * dx * dy * tz;
        val += gridref.relem(localix1, iy1, iz1) * dx * dy * dz;
      }
    }else{
      size_t ix1 = (ix + 1) % nx_;
      val += gridref.relem(ix, iy, iz) * tx * ty * tz;
      val += gridref.relem(ix, iy, iz1) * tx * ty * dz;
      val += gridref.relem(ix, iy1, iz) * tx * dy * tz;
      val += gridref.relem(ix, iy1, iz1) * tx * dy * dz;
      val += gridref.relem(ix1, iy, iz) * dx * ty * tz;
      val += gridref.relem(ix1, iy, iz1) * dx * ty * dz;
      val += gridref.relem(ix1, iy1, iz) * dx * dy * tz;
      val += gridref.relem(ix1, iy1, iz1) * dx * dy * dz;
    }
    return val;
  }

  // data_t get_tsc_at(const std::array<real_t, 3> &pos, std::vector<data_t> &val) const
  // {
  // }

  int get_task(const vec3 &x) const noexcept
  {
    const auto it = std::upper_bound(local0starts_.begin(), local0starts_.end(), int(x[0]));
    return std::distance(local0starts_.begin(), it)-1;
  }

  void domain_decompose_pos(std::vector<vec3> &pos) const noexcept
  {
    if (is_distributed_trait)
    {
#if defined(USE_MPI)
      std::sort(pos.begin(), pos.end(), [&](auto x1, auto x2) { return this->get_task(x1) < this->get_task(x2); });
      std::vector<int> sendcounts(MPI::get_size(), 0), sendoffsets(MPI::get_size(), 0);
      std::vector<int> recvcounts(MPI::get_size(), 0), recvoffsets(MPI::get_size(), 0);
      for (auto x : pos)
      {
        sendcounts[this->get_task(x)] += 3;
      }

      MPI_Alltoall(&sendcounts[0], 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

      size_t tot_receive = recvcounts[0];
//      size_t tot_send = sendcounts[0];
      for (int i = 1; i < MPI::get_size(); ++i)
      {
        sendoffsets[i] = sendcounts[i - 1] + sendoffsets[i - 1];
        recvoffsets[i] = recvcounts[i - 1] + recvoffsets[i - 1];
        tot_receive += recvcounts[i];
//        tot_send += sendcounts[i];
      }

      std::vector<vec3> recvbuf(tot_receive/3,{0.,0.,0.});

      MPI_Alltoallv(&pos[0], &sendcounts[0], &sendoffsets[0], MPI::get_datatype<real_t>(),
                    &recvbuf[0], &recvcounts[0], &recvoffsets[0], MPI::get_datatype<real_t>(), MPI_COMM_WORLD);

      pos.swap( recvbuf );
#endif
    }
  }

  ccomplex_t compensation_kernel( const vec3_t<real_t>& k ) const noexcept
  {
    auto sinc = []( real_t x ){ return (std::fabs(x)>1e-10)? std::sin(x)/x : 1.0; };
    real_t dfx = sinc(0.5*M_PI*k[0]/gridref.kny_[0]);
    real_t dfy = sinc(0.5*M_PI*k[1]/gridref.kny_[1]);
    real_t dfz = sinc(0.5*M_PI*k[2]/gridref.kny_[2]);
    real_t del = std::pow(dfx*dfy*dfz,1+interpolation_order);

    real_t shift = 0.5 * k[0] * gridref.get_dx()[0] + 0.5 * k[1] * gridref.get_dx()[1] + 0.5 * k[2] * gridref.get_dx()[2];

    return std::exp(ccomplex_t(0.0, shift)) / del;
  }

};
