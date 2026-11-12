// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2022 by Oliver Hahn
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
#include <numeric>

#include <general.hh>
#include <math/vec3.hh>

/// @brief implements a wrapper class for grids with ghost zones for MPI communication (slab decomposition)
/// @tparam numghosts number of ghost zones on each side
/// @tparam haveleft flag whether to have left ghost zone
/// @tparam haveright flag whether to have right ghost zone
/// @tparam grid_t grid type to wrap
template <int numghosts, bool haveleft, bool haveright, typename grid_t>
struct grid_with_ghosts
{
  using data_t = typename grid_t::data_t;
  using vec3 = std::array<real_t, 3>;

  static constexpr bool is_distributed_trait = grid_t::is_distributed_trait;
  static constexpr int num_ghosts = numghosts;
  static constexpr bool have_left = haveleft, have_right = haveright;

  std::vector<data_t> boundary_left_, boundary_right_;
  std::vector<int> local0starts_;
  const grid_t &gridref;
  size_t nx_, ny_, nz_, nzp_;


  //... determine communication offsets
  std::vector<ptrdiff_t> offsets_, sizes_;

  /// @brief get task index for a given index
  /// @param index index
  /// @return task index
  int get_task(ptrdiff_t index) const
  {
    int itask = 0;
#if defined(USE_MPI)
    while (itask < MPI::get_size() - 1 && offsets_[itask + 1] <= index)
        ++itask;
#endif
    return itask;
  }

  /// @brief constructor for grid with ghosts
  /// @param g grid to wrap
  explicit grid_with_ghosts(const grid_t &g)
  : gridref(g), nx_(g.n_[0]), ny_(g.n_[1]), nz_(g.n_[2]), nzp_(g.n_[2]+2)
  {
#if defined(USE_MPI)
    if (is_distributed_trait)
    {
      int ntasks(MPI::get_size());

      offsets_.assign(ntasks+1, 0);
      sizes_.assign(ntasks, 0);
      
      MPI_Allgather(&g.local_0_size_, 1, MPI_LONG_LONG, &sizes_[0], 1,
                      MPI_LONG_LONG, MPI_COMM_WORLD);
      MPI_Allgather(&g.local_0_start_, 1, MPI_LONG_LONG, &offsets_[0], 1,
                      MPI_LONG_LONG, MPI_COMM_WORLD);
      
      for( int i=0; i< CONFIG::MPI_task_size; i++ ){
          if( offsets_[i+1] < offsets_[i] + sizes_[i] ) offsets_[i+1] = offsets_[i] + sizes_[i];
      }

      update_ghosts_allow_multiple( g );
    }
#endif
  }

  /// @brief update ghost zones via MPI communication
  /// @param g grid to wrap
  void update_ghosts_allow_multiple( const grid_t &g )
  {
  #if defined(USE_MPI)
    //... exchange boundary
    if( have_left  ) boundary_left_.assign(num_ghosts * ny_ * nzp_, data_t{0.0});
    if( have_right ) boundary_right_.assign(num_ghosts * ny_ * nzp_, data_t{0.0});

    size_t slicesz = ny_ * nzp_;

    MPI_Status status;
    std::vector<MPI_Request> req;
    MPI_Request temp_req;

    if( have_right ){
      for( int itask=0; itask<CONFIG::MPI_task_size; ++itask ){
        for( size_t i=0; i<num_ghosts; ++i ){
          ptrdiff_t iglobal_request = (offsets_[itask] + sizes_[itask] + i) % g.n_[0];
          if( iglobal_request >= g.local_0_start_ && iglobal_request < g.local_0_start_ + g.local_0_size_ ){
            size_t ii = iglobal_request - g.local_0_start_;
            MPI_Isend( &g.relem(ii*slicesz), slicesz, MPI::get_datatype<data_t>(), itask, iglobal_request, MPI_COMM_WORLD, &temp_req);
            req.push_back(temp_req);
          }
        }
      }

      //--- receive data ------------------------------------------------------------
      #pragma omp parallel if(CONFIG::MPI_threads_ok)
      {
          MPI_Status status;

          #pragma omp for 
          for( size_t i=0; i<num_ghosts; ++i ){
            ptrdiff_t iglobal_request = (g.local_0_start_ + g.local_0_size_ + i) % g.n_[0];
          
            int recvfrom = get_task(iglobal_request);

            //#pragma omp critical // need critical region here if we do "MPI_THREAD_FUNNELED", 
            {
                // receive data slice and check for MPI errors when in debug mode
                status.MPI_ERROR = MPI_SUCCESS;
                MPI_Recv(&boundary_right_[i*slicesz], (int)slicesz, MPI::get_datatype<data_t>(), recvfrom, (int)iglobal_request, MPI_COMM_WORLD, &status);
                assert(status.MPI_ERROR == MPI_SUCCESS);
            }
          }
      }
    }

    MPI_Barrier( MPI_COMM_WORLD );

    if( have_left ){
      for( int itask=0; itask<CONFIG::MPI_task_size; ++itask ){
        for( size_t i=0; i<num_ghosts; ++i ){
          ptrdiff_t iglobal_request = (offsets_[itask] + g.n_[0] - num_ghosts + i) % g.n_[0];
          if( iglobal_request >= g.local_0_start_ && iglobal_request < g.local_0_start_ + g.local_0_size_ ){
            size_t ii = iglobal_request - g.local_0_start_;
            MPI_Isend( &g.relem(ii*slicesz), slicesz, MPI::get_datatype<data_t>(), itask, iglobal_request, MPI_COMM_WORLD, &temp_req);
            req.push_back(temp_req);
          }
        }
      }

      //--- receive data ------------------------------------------------------------
      #pragma omp parallel if(CONFIG::MPI_threads_ok)
      {
          MPI_Status status;

          #pragma omp for 
          for( size_t i=0; i<num_ghosts; ++i ){
            ptrdiff_t iglobal_request = (g.local_0_start_ + g.n_[0] - num_ghosts + i) % g.n_[0];
          
            int recvfrom = get_task(iglobal_request);

            //#pragma omp critical // need critical region here if we do "MPI_THREAD_FUNNELED", 
            {
                // receive data slice and check for MPI errors when in debug mode
                status.MPI_ERROR = MPI_SUCCESS;
                MPI_Recv(&boundary_left_[i*slicesz], (int)slicesz, MPI::get_datatype<data_t>(), recvfrom, (int)iglobal_request, MPI_COMM_WORLD, &status);
                assert(status.MPI_ERROR == MPI_SUCCESS);
            }
          }
      }
    }

    MPI_Barrier( MPI_COMM_WORLD );
    
    for (size_t i = 0; i < req.size(); ++i)
    {
        // need to set status as wait does not necessarily modify it
        // c.f. http://www.open-mpi.org/community/lists/devel/2007/04/1402.php
        status.MPI_ERROR = MPI_SUCCESS;
        // std::cout << "task " << CONFIG::MPI_task_rank << " : checking request No" << i << std::endl;
        int flag(1);
        MPI_Test(&req[i], &flag, &status);
        if( !flag ){
            std::cout << "task " << CONFIG::MPI_task_rank << " : request No" << i << " unsuccessful" << std::endl;
        }

        MPI_Wait(&req[i], &status);
        // std::cout << "---> ok!" << std::endl;
        assert(status.MPI_ERROR == MPI_SUCCESS);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    #endif
  }

  /// @brief return the element at position (i,j,k) in the grid
  /// @param i index in x direction
  /// @param j index in y direction
  /// @param k index in z direction
  /// @return grid element at position (i,j,k)
  data_t relem(const ptrdiff_t& i, const ptrdiff_t& j, const ptrdiff_t&k ) const noexcept
  {
    return this->relem({i,j,k});
  }

  /// @brief return the element at position (i,j,k) in the grid
  /// @param pos position in the grid (array of size 3: {i,j,k})
  /// @return grid element at position (i,j,k)
  data_t relem(const std::array<ptrdiff_t, 3> &pos) const noexcept
  {
    const ptrdiff_t ix = pos[0];
    const ptrdiff_t iy = (pos[1]+gridref.n_[1])%gridref.n_[1];
    const ptrdiff_t iz = (pos[2]+gridref.n_[2])%gridref.n_[2];

    if( is_distributed_trait ){
      const ptrdiff_t localix = ix;
      if( localix < 0 ){
        return boundary_left_[((localix+num_ghosts)*ny_+iy)*nzp_+iz];
      }else if( localix >= gridref.local_0_size_ ){
        return boundary_right_[((localix-gridref.local_0_size_)*ny_+iy)*nzp_+iz];
      }else{
        return gridref.relem(localix, iy, iz);
      }
    }

    return gridref.relem((ix+gridref.n_[0])%gridref.n_[0], iy, iz);
  }
};

