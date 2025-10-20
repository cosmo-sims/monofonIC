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
#include <algorithm>

#include <general.hh>
#include <config_file.hh>
#include <grid_interpolate.hh>

#if defined(USE_HDF5)
#include "HDF_IO.hh"
#endif

namespace particle
{
    using vec3 = std::array<real_t,3>;

    template <typename field_t>
    struct glass_loader
    {
        using data_t = typename field_t::data_t;
        size_t num_p, off_p;
        grid_interpolate<1, field_t> interp_;
        std::vector<vec3> glass_posr;

        glass_loader( config_file& cf, const field_t &field )
        : num_p(0), off_p(0), interp_( field )
        {
            real_t lglassbox = 1.0;
            std::string glass_fname = cf.get_value<std::string>("setup", "GlassFileName");
            size_t ntiles = cf.get_value<size_t>("setup", "GlassTiles");

#if defined(USE_HDF5)
            // Read box size attribute (all ranks)
            HDFReadGroupAttribute(glass_fname, "Header", "BoxSize", lglassbox);

            // Get dataset size without reading full data
            std::vector<int> glass_extent;
            HDFGetDatasetExtent(glass_fname, "/PartType1/Coordinates", glass_extent);
            size_t np_in_file = glass_extent[0];
#else
            throw std::runtime_error("Glass lattice requires HDF5 support. Enable and recompile.");
#endif

#if defined(USE_MPI)
            // Broadcast number of particles in file to all ranks
            MPI_Bcast(&np_in_file, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

            const int nranks = MPI::get_size();
            const int myrank = MPI::get_rank();

            const size_t ntiles3 = ntiles * ntiles * ntiles;
            const size_t total_glass_particles = np_in_file * ntiles3;
            const size_t fft_grid_size = field.n_[0] * field.n_[1] * field.n_[2];
            const real_t fft_oversampling = real_t(total_glass_particles) / real_t(fft_grid_size);

            music::ilog << "Glass file contains " << np_in_file << " particles." << std::endl;
            music::ilog << "Tiling " << ntiles << "^3 = " << ntiles3 << "x results in " << total_glass_particles << " particles total." << std::endl;
            music::ilog << "FFT grid oversampling factor: " << fft_oversampling << "x (vs. " << field.n_[0] << "^3 grid)" << std::endl;

            // Calculate chunk sizes for each rank (handle remainder)
            const size_t base_chunk = np_in_file / static_cast<size_t>(nranks);
            const size_t remainder = np_in_file % static_cast<size_t>(nranks);

            // Determine how many base glass particles this rank will receive
            const size_t my_base_count = base_chunk + (static_cast<size_t>(myrank) < remainder ? 1 : 0);

            // Allocate space for local base glass particles
            std::vector<real_t> local_base_glass(my_base_count * 3);

            // Chunked read and scatter loop
            for (int ichunk = 0; ichunk < nranks; ++ichunk)
            {
                const size_t chunk_count = base_chunk + (static_cast<size_t>(ichunk) < remainder ? 1 : 0);
                const size_t chunk_offset = static_cast<size_t>(ichunk) * base_chunk + std::min(static_cast<size_t>(ichunk), remainder);

                std::vector<real_t> chunk_data;

                // Rank 0 reads the chunk from HDF5 using existing HDFReadDatasetSlab function
                if (myrank == 0)
                {
                    HDFReadDatasetSlab(glass_fname, "/PartType1/Coordinates",
                                      chunk_offset, chunk_count, chunk_data);
                }

                // Scatter chunk to the appropriate rank
                if (ichunk == myrank)
                {
                    if (myrank == 0)
                    {
                        local_base_glass = std::move(chunk_data);
                    }
                    else
                    {
                        MPI_Recv(local_base_glass.data(), chunk_count * 3, MPI::get_datatype<real_t>(),
                                 0, ichunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                else if (myrank == 0)
                {
                    MPI_Send(chunk_data.data(), chunk_count * 3, MPI::get_datatype<real_t>(),
                             ichunk, ichunk, MPI_COMM_WORLD);
                }
            }

            // Now tile the local base glass particles
            glass_posr.assign(my_base_count * ntiles3, {0.0, 0.0, 0.0});

            std::array<real_t, 3> ng({real_t(field.n_[0]), real_t(field.n_[1]), real_t(field.n_[2])});

            #pragma omp parallel for
            for (size_t ibase = 0; ibase < my_base_count; ++ibase)
            {
                for (size_t itile = 0; itile < ntiles3; ++itile)
                {
                    size_t tile_x = itile / (ntiles * ntiles);
                    size_t tile_y = (itile / ntiles) % ntiles;
                    size_t tile_z = itile % ntiles;

                    size_t idx = ibase * ntiles3 + itile;
                    glass_posr[idx][0] = std::fmod((local_base_glass[3 * ibase + 0] / lglassbox + real_t(tile_x)) / ntiles * ng[0] + ng[0], ng[0]);
                    glass_posr[idx][1] = std::fmod((local_base_glass[3 * ibase + 1] / lglassbox + real_t(tile_y)) / ntiles * ng[1] + ng[1], ng[1]);
                    glass_posr[idx][2] = std::fmod((local_base_glass[3 * ibase + 2] / lglassbox + real_t(tile_z)) / ntiles * ng[2] + ng[2], ng[2]);
                }
            }

            // Apply spatial domain decomposition
            interp_.domain_decompose_pos(glass_posr);

            num_p = glass_posr.size();
            std::vector<size_t> all_num_p( MPI::get_size(), 0 );
            MPI_Allgather( &num_p, 1, MPI_UNSIGNED_LONG_LONG, &all_num_p[0], 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD );
            off_p = 0;
            for( int itask=0; itask<myrank; ++itask ){
                off_p += all_num_p[itask];
            }
#else
            const size_t ntiles3 = ntiles * ntiles * ntiles;
            const size_t total_glass_particles = np_in_file * ntiles3;
            const size_t fft_grid_size = field.n_[0] * field.n_[1] * field.n_[2];
            const real_t glass_oversampling = real_t(total_glass_particles) / real_t(np_in_file);
            const real_t fft_oversampling = real_t(total_glass_particles) / real_t(fft_grid_size);

            music::ilog << "Glass file contains " << np_in_file << " particles." << std::endl;
            music::ilog << "Tiling " << ntiles << "^3 = " << ntiles3 << "x results in " << total_glass_particles << " particles total." << std::endl;
            music::ilog << "Glass oversampling factor: " << glass_oversampling << "x (vs. base glass)" << std::endl;
            music::ilog << "FFT grid oversampling factor: " << fft_oversampling << "x (vs. " << field.n_[0] << "^3 grid)" << std::endl;

            // Non-MPI version: read full glass file
            std::vector<real_t> glass_pos;
            HDFReadDataset(glass_fname, "/PartType1/Coordinates", glass_pos);

            num_p = total_glass_particles;
            off_p = 0;

            glass_posr.assign(num_p, {0.0, 0.0, 0.0});
            std::array<real_t, 3> ng({real_t(field.n_[0]), real_t(field.n_[1]), real_t(field.n_[2])});

            #pragma omp parallel for
            for (size_t i = 0; i < num_p; ++i)
            {
                size_t idx_in_glass = i % np_in_file;
                size_t idxtile = i / np_in_file;
                size_t tile_x = idxtile / (ntiles * ntiles);
                size_t tile_y = (idxtile / ntiles) % ntiles;
                size_t tile_z = idxtile % ntiles;
                glass_posr[i][0] = std::fmod((glass_pos[3 * idx_in_glass + 0] / lglassbox + real_t(tile_x)) / ntiles * ng[0] + ng[0], ng[0]);
                glass_posr[i][1] = std::fmod((glass_pos[3 * idx_in_glass + 1] / lglassbox + real_t(tile_y)) / ntiles * ng[1] + ng[1], ng[1]);
                glass_posr[i][2] = std::fmod((glass_pos[3 * idx_in_glass + 2] / lglassbox + real_t(tile_z)) / ntiles * ng[2] + ng[2], ng[2]);
            }
#endif
        }

        void update_ghosts( const field_t &field )
        {
            interp_.update_ghosts( field );
        }

        data_t get_at( const vec3& x ) const noexcept
        {
            return interp_.get_cic_at( x );
        }

        size_t size() const noexcept
        {
            return num_p;
        }

        size_t offset() const noexcept
        {
            return off_p;
        }
    };

} // namespace particle
