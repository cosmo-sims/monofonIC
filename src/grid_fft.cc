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

#include <general.hh>
#include <grid_fft.hh>
#include <thread>

#include "memory_stat.hh"

void memory_report(void) 
{
    //... report memory usage
    size_t curr_mem_high_mark = 0;
    local_mem_high_mark = memory::getCurrentRSS();
#if defined(USE_MPI)
    MPI_Allreduce(&local_mem_high_mark, &curr_mem_high_mark, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
    curr_mem_high_mark = local_mem_high_mark;
#endif
    if( curr_mem_high_mark > 1.1*global_mem_high_mark ){
        music::ilog << std::setw(57) << std::setfill(' ') << std::right << "mem high: " << std::setw(8) << curr_mem_high_mark/(1ull<<20) << " MBytes / task" << std::endl;
        global_mem_high_mark = curr_mem_high_mark;
    }
}

template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::allocate(void)
{
    if (!bdistributed)
    {
        ntot_ = (n_[2] + 2) * n_[1] * n_[0];

        music::dlog.Print("[FFT] Setting up a shared memory field %lux%lux%lu\n", n_[0], n_[1], n_[2]);
        if (typeid(data_t) == typeid(real_t))
        {
            data_ = reinterpret_cast<data_t *>(FFTW_API(malloc)(ntot_ * sizeof(real_t)));
            cdata_ = reinterpret_cast<ccomplex_t *>(data_);

            plan_ = FFTW_API(plan_dft_r2c_3d)(n_[0], n_[1], n_[2], (real_t *)data_, (complex_t *)data_, FFTW_RUNMODE);
            iplan_ = FFTW_API(plan_dft_c2r_3d)(n_[0], n_[1], n_[2], (complex_t *)data_, (real_t *)data_, FFTW_RUNMODE);
        }
        else if (typeid(data_t) == typeid(ccomplex_t))
        {
            data_ = reinterpret_cast<data_t *>(FFTW_API(malloc)(ntot_ * sizeof(ccomplex_t)));
            cdata_ = reinterpret_cast<ccomplex_t *>(data_);

            plan_ = FFTW_API(plan_dft_3d)(n_[0], n_[1], n_[2], (complex_t *)data_, (complex_t *)data_, FFTW_FORWARD, FFTW_RUNMODE);
            iplan_ = FFTW_API(plan_dft_3d)(n_[0], n_[1], n_[2], (complex_t *)data_, (complex_t *)data_, FFTW_BACKWARD, FFTW_RUNMODE);
        }
        else
        {
            music::elog.Print("invalid data type in Grid_FFT<data_t>::setup_fft_interface\n");
        }

        fft_norm_fac_ = 1.0 / std::sqrt((real_t)((size_t)n_[0] * (real_t)n_[1] * (real_t)n_[2]));

        if (typeid(data_t) == typeid(real_t))
        {
            npr_ = n_[2] + 2;
            npc_ = n_[2] / 2 + 1;
        }
        else
        {
            npr_ = n_[2];
            npc_ = n_[2];
        }

        for (int i = 0; i < 3; ++i)
        {
            nhalf_[i] = n_[i] / 2;
            kfac_[i] = 2.0 * M_PI / length_[i];
            kny_[i] = kfac_[i] * n_[i]/2;
            dx_[i] = length_[i] / n_[i];

            global_range_.x1_[i] = 0;
            global_range_.x2_[i] = n_[i];
        }

        local_0_size_ = n_[0];
        local_1_size_ = n_[1];
        local_0_start_ = 0;
        local_1_start_ = 0;

        if (space_ == rspace_id)
        {
            sizes_[0] = n_[0];
            sizes_[1] = n_[1];
            sizes_[2] = n_[2];
            sizes_[3] = npr_;
        }
        else
        {
            sizes_[0] = n_[1];
            sizes_[1] = n_[0];
            sizes_[2] = npc_;
            sizes_[3] = npc_;
        }
    }
    else
    {
#ifdef USE_MPI //// i.e. ifdef USE_MPI ////////////////////////////////////////////////////////////////////////////////////
        size_t cmplxsz;

        if (typeid(data_t) == typeid(real_t))
        {
            cmplxsz = FFTW_API(mpi_local_size_3d_transposed)(n_[0], n_[1], n_[2], MPI_COMM_WORLD,
                                                             &local_0_size_, &local_0_start_, &local_1_size_, &local_1_start_);
            ntot_ = local_0_size_ * n_[1] * (n_[2]+2);
            data_ = (data_t *)FFTW_API(malloc)(ntot_ * sizeof(real_t));
            cdata_ = reinterpret_cast<ccomplex_t *>(data_);
            plan_ = FFTW_API(mpi_plan_dft_r2c_3d)(n_[0], n_[1], n_[2], (real_t *)data_, (complex_t *)data_,
                                                  MPI_COMM_WORLD, FFTW_RUNMODE | FFTW_MPI_TRANSPOSED_OUT);
            iplan_ = FFTW_API(mpi_plan_dft_c2r_3d)(n_[0], n_[1], n_[2], (complex_t *)data_, (real_t *)data_,
                                                   MPI_COMM_WORLD, FFTW_RUNMODE | FFTW_MPI_TRANSPOSED_IN);
        }
        else if (typeid(data_t) == typeid(ccomplex_t))
        {
            cmplxsz = FFTW_API(mpi_local_size_3d_transposed)(n_[0], n_[1], n_[2], MPI_COMM_WORLD,
                                                             &local_0_size_, &local_0_start_, &local_1_size_, &local_1_start_);
            ntot_ = cmplxsz;
            data_ = (data_t *)FFTW_API(malloc)(ntot_ * sizeof(ccomplex_t));
            cdata_ = reinterpret_cast<ccomplex_t *>(data_);
            plan_ = FFTW_API(mpi_plan_dft_3d)(n_[0], n_[1], n_[2], (complex_t *)data_, (complex_t *)data_,
                                              MPI_COMM_WORLD, FFTW_FORWARD, FFTW_RUNMODE | FFTW_MPI_TRANSPOSED_OUT);
            iplan_ = FFTW_API(mpi_plan_dft_3d)(n_[0], n_[1], n_[2], (complex_t *)data_, (complex_t *)data_,
                                               MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_RUNMODE | FFTW_MPI_TRANSPOSED_IN);
        }
        else
        {
            music::elog.Print("unknown data type in Grid_FFT<data_t>::setup_fft_interface\n");
            abort();
        }

        music::dlog.Print("[FFT] Setting up a distributed memory field %lux%lux%lu\n", n_[0], n_[1], n_[2]);

        fft_norm_fac_ = 1.0 / sqrt((real_t)n_[0] * (real_t)n_[1] * (real_t)n_[2]);

        if (typeid(data_t) == typeid(real_t))
        {
            npr_ = n_[2] + 2;
            npc_ = n_[2] / 2 + 1;
        }
        else
        {
            npr_ = n_[2];
            npc_ = n_[2];
        }

        for (int i = 0; i < 3; ++i)
        {
            nhalf_[i] = n_[i] / 2;
            kfac_[i] = 2.0 * M_PI / length_[i];
            kny_[i] = kfac_[i] * n_[i]/2;
            dx_[i] = length_[i] / n_[i];

            global_range_.x1_[i] = 0;
            global_range_.x2_[i] = n_[i];
        }
        global_range_.x1_[0] = (int)local_0_start_;
        global_range_.x2_[0] = (int)(local_0_start_ + local_0_size_);

        if (space_ == rspace_id)
        {
            sizes_[0] = (int)local_0_size_;
            sizes_[1] = n_[1];
            sizes_[2] = n_[2];
            sizes_[3] = npr_; // holds the physical memory size along the 3rd dimension
        }
        else
        {
            sizes_[0] = (int)local_1_size_;
            sizes_[1] = n_[0];
            sizes_[2] = npc_;
            sizes_[3] = npc_; // holds the physical memory size along the 3rd dimension
        }
#else
        music::flog << "MPI is required for distributed FFT arrays!" << std::endl;
        throw std::runtime_error("MPI is required for distributed FFT arrays!");
#endif //// of #ifdef #else USE_MPI ////////////////////////////////////////////////////////////////////////////////////
    }
    ballocated_ = true;
    memory_report();
}

template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::ApplyNorm(void)
{
#pragma omp parallel for
    for (size_t i = 0; i < ntot_; ++i)
        data_[i] *= fft_norm_fac_;
}

//! Perform a backwards Fourier transformation
template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::FourierTransformForward(bool do_transform)
{
#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (space_ != kspace_id)
    {
        //.............................
        if (do_transform)
        {
            double wtime = get_wtime();
            music::dlog.Print("[FFT] Calling Grid_FFT::to_kspace (%lux%lux%lu)", sizes_[0], sizes_[1], sizes_[2]);
            FFTW_API(execute)
            (plan_);
            this->ApplyNorm();

            wtime = get_wtime() - wtime;
            music::dlog.Print("[FFT] Completed Grid_FFT::to_kspace (%lux%lux%lu), took %f s", sizes_[0], sizes_[1], sizes_[2], wtime);
        }

        sizes_[0] = local_1_size_;
        sizes_[1] = n_[0];
        sizes_[2] = (int)npc_;
        sizes_[3] = npc_;

        space_ = kspace_id;
        //.............................
    }
}

//! Perform a forward Fourier transformation
template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::FourierTransformBackward(bool do_transform)
{
#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (space_ != rspace_id)
    {
        //.............................
        if (do_transform)
        {
            music::dlog.Print("[FFT] Calling Grid_FFT::to_rspace (%dx%dx%d)\n", sizes_[0], sizes_[1], sizes_[2]);
            double wtime = get_wtime();

            FFTW_API(execute)(iplan_);
            this->ApplyNorm();

            wtime = get_wtime() - wtime;
            music::dlog.Print("[FFT] Completed Grid_FFT::to_rspace (%dx%dx%d), took %f s\n", sizes_[0], sizes_[1], sizes_[2], wtime);
        }
        sizes_[0] = local_0_size_;
        sizes_[1] = n_[1];
        sizes_[2] = n_[2];
        sizes_[3] = npr_;

        space_ = rspace_id;
        //.............................
    }
}

//! Perform a copy to another field, not necessarily the same size, using Fourier interpolation
template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::FourierInterpolateCopyTo( grid_fft_t &grid_to )
{
    grid_fft_t &grid_from = *this;

    //... transform both grids to Fourier space
    grid_from.FourierTransformForward(true);
    grid_to.FourierTransformForward(false);
    
    // if grids are same size, we directly copy without the need for interpolation
    if( grid_from.n_[0] == grid_to.n_[0] 
        && grid_from.n_[1] == grid_to.n_[1] 
        && grid_from.n_[2] == grid_to.n_[2] )
    {
        grid_to.copy_from( grid_from );
        return;
    }

    // set to zero so that regions without data are zeroed
    grid_to.zero();

    // if not running MPI, then can do without sending and receiving
#if !defined(USE_MPI)

    // determine effective Nyquist modes representable by both fields and their locations in array
    size_t fny0_left  = std::min(grid_from.n_[0] / 2, grid_to.n_[0] / 2);
    size_t fny0_right = std::max(grid_from.n_[0] - grid_to.n_[0] / 2, grid_from.n_[0] / 2);
    size_t fny1_left  = std::min(grid_from.n_[1] / 2, grid_to.n_[1] / 2);
    size_t fny1_right = std::max(grid_from.n_[1] - grid_to.n_[1] / 2, grid_from.n_[1] / 2);
    size_t fny2_left  = std::min(grid_from.n_[2] / 2, grid_to.n_[2] / 2);
    // size_t fny2_right = std::max(grid_from.n_[2] - grid_to.n_[2] / 2, grid_from.n_[2] / 2);

    const size_t fny0_left_recv  = fny0_left;
    const size_t fny0_right_recv = (fny0_right + grid_to.n_[0]) - grid_from.n_[0];
    const size_t fny1_left_recv  = fny1_left;
    const size_t fny1_right_recv = (fny1_right + grid_to.n_[1]) - grid_from.n_[1];
    const size_t fny2_left_recv  = fny2_left;
    // const size_t fny2_right_recv = (fny2_right + grid_to.n_[2]) - grid_from.n_[2];

    #pragma omp parallel for
    for( size_t i=0; i<grid_to.size(0); ++i )
    {
        if (i < fny0_left_recv || i > fny0_right_recv)
        {
            size_t isend = (i < fny0_left_recv)? i : (i + grid_from.n_[0]) - grid_to.n_[0];

            // copy data slice into new grid, avoiding modes that do not exist in both
            for( size_t j=0; j<grid_to.size(1); ++j )
            {
                if( j < fny1_left_recv || j > fny1_right_recv )
                {
                    const size_t jsend = (j < fny1_left_recv)? j : (j + grid_from.n_[1]) - grid_to.n_[1];

                    for( size_t k=0; k<fny2_left_recv; ++k )
                    {
                        grid_to.kelem(i,j,k) = grid_from.kelem(isend,jsend,k);
                    }
                }
            }
        }
    }

#else
    // if they are not the same size, we use Fourier interpolation to upscale/downscale
    double tstart = get_wtime();
    music::dlog << "[MPI] Started scatter for Fourier interpolation/copy" << std::endl;

    //... determine communication offsets
    std::vector<ptrdiff_t> offsets_send, offsets_recv, sizes_send, sizes_recv;

    // this should use bisection, not linear search
    auto get_task = [](ptrdiff_t index, const std::vector<ptrdiff_t> &offsets, const int ntasks) -> int
    {
        int itask = 0;
        while (itask < ntasks - 1 && offsets[itask + 1] <= index)
            ++itask;
        return itask;
    };

    int ntasks(MPI::get_size());

    offsets_send.assign(ntasks+1, 0);
    sizes_send.assign(ntasks, 0);
    offsets_recv.assign(ntasks+1, 0);
    sizes_recv.assign(ntasks, 0);

    MPI_Allgather(&grid_from.local_1_size_, 1, MPI_LONG_LONG, &sizes_send[0], 1,
                    MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&grid_to.local_1_size_, 1, MPI_LONG_LONG, &sizes_recv[0], 1,
                    MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&grid_from.local_1_start_, 1, MPI_LONG_LONG, &offsets_send[0], 1,
                    MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&grid_to.local_1_start_, 1, MPI_LONG_LONG, &offsets_recv[0], 1,
                    MPI_LONG_LONG, MPI_COMM_WORLD);

    for( int i=0; i< CONFIG::MPI_task_size; i++ ){
        if( offsets_send[i+1] < offsets_send[i] + sizes_send[i] ) offsets_send[i+1] = offsets_send[i] + sizes_send[i];
        if( offsets_recv[i+1] < offsets_recv[i] + sizes_recv[i] ) offsets_recv[i+1] = offsets_recv[i] + sizes_recv[i];
    }

    const MPI_Datatype datatype =
        (typeid(data_t) == typeid(float)) ? MPI_C_FLOAT_COMPLEX 
        : (typeid(data_t) == typeid(double)) ? MPI_C_DOUBLE_COMPLEX 
        : (typeid(data_t) == typeid(long double)) ? MPI_C_LONG_DOUBLE_COMPLEX
        : MPI_BYTE;

    const size_t slicesz = grid_from.size(1) * grid_from.size(3);

    // determine effective Nyquist modes representable by both fields and their locations in array
    const size_t fny0_left  = std::min(grid_from.n_[1] / 2, grid_to.n_[1] / 2);
    const size_t fny0_right = std::max(grid_from.n_[1] - grid_to.n_[1] / 2, grid_from.n_[1] / 2);
    const size_t fny1_left  = std::min(grid_from.n_[0] / 2, grid_to.n_[0] / 2);
    const size_t fny1_right = std::max(grid_from.n_[0] - grid_to.n_[0] / 2, grid_from.n_[0] / 2);
    const size_t fny2_left  = std::min(grid_from.n_[2] / 2, grid_to.n_[2] / 2);
    // size_t fny2_right = std::max(grid_from.n_[2] - grid_to.n_[2] / 2, grid_from.n_[2] / 2);

    const size_t fny0_left_recv  = fny0_left;
    const size_t fny0_right_recv = (fny0_right + grid_to.n_[1]) - grid_from.n_[1];
    const size_t fny1_left_recv  = fny1_left;
    const size_t fny1_right_recv = (fny1_right + grid_to.n_[0]) - grid_from.n_[0];
    const size_t fny2_left_recv  = fny2_left;
    // const size_t fny2_right_recv = (fny2_right + grid_to.n_[2]) - grid_from.n_[2];

    //--- send data from buffer ---------------------------------------------------
    
    std::vector<MPI_Request> req;
    MPI_Request temp_req;

    for (size_t i = 0; i < grid_from.size(0); ++i)
    {
        size_t iglobal_send = i + offsets_send[CONFIG::MPI_task_rank];
        
        if (iglobal_send < fny0_left) 
        {
            size_t iglobal_recv = iglobal_send;
            int sendto = get_task(iglobal_recv, offsets_recv, CONFIG::MPI_task_size);
            MPI_Isend(&grid_from.kelem(i * slicesz), (int)slicesz, datatype, sendto,
                        (int)iglobal_recv, MPI_COMM_WORLD, &temp_req);
            req.push_back(temp_req);
            // std::cout << "task " << CONFIG::MPI_task_rank << " : added request No" << req.size()-1 << ": Isend #" << iglobal_send << " to task " << sendto << ", size = " << slicesz << std::endl;
        }
        if (iglobal_send > fny0_right) 
        {
            size_t iglobal_recv = (iglobal_send + grid_to.n_[1]) - grid_from.n_[1]; 
            int sendto = get_task(iglobal_recv, offsets_recv, CONFIG::MPI_task_size);
            MPI_Isend(&grid_from.kelem(i * slicesz), (int)slicesz, datatype, sendto,
                        (int)iglobal_recv, MPI_COMM_WORLD, &temp_req);
            req.push_back(temp_req);
            // std::cout << "task  " << CONFIG::MPI_task_rank << " : added request No" << req.size()-1 << ": Isend #" << iglobal_send << " to task " << sendto << ", size = " << slicesz<< std::endl;
        }
    }

    //--- receive data ------------------------------------------------------------
    #pragma omp parallel if(CONFIG::MPI_threads_ok)
    {
        MPI_Status status;
        ccomplex_t  * recvbuf = new ccomplex_t[ slicesz ]; 

        #pragma omp for schedule(dynamic) 
        for( size_t i=0; i<grid_to.size(0); ++i )
        {   
            size_t iglobal_recv = i + offsets_recv[CONFIG::MPI_task_rank];

            if (iglobal_recv < fny0_left_recv || iglobal_recv > fny0_right_recv)
            {
                size_t iglobal_send = (iglobal_recv < fny0_left_recv)? iglobal_recv : (iglobal_recv + grid_from.n_[1]) - grid_to.n_[1];

                int recvfrom = get_task(iglobal_send, offsets_send, CONFIG::MPI_task_size);
                
                //#pragma omp critical // need critical region here if we do "MPI_THREAD_FUNNELED", 
                {
                    // receive data slice and check for MPI errors when in debug mode
                    status.MPI_ERROR = MPI_SUCCESS;
                    MPI_Recv(&recvbuf[0], (int)slicesz, datatype, recvfrom, (int)iglobal_recv, MPI_COMM_WORLD, &status);
                    assert(status.MPI_ERROR == MPI_SUCCESS);
                }

                // copy data slice into new grid, avoiding modes that do not exist in both
                for( size_t j=0; j<grid_to.size(1); ++j )
                {
                    if( j < fny1_left_recv || j > fny1_right_recv )
                    {
                        const size_t jsend = (j < fny1_left_recv)? j : (j + grid_from.n_[0]) - grid_to.n_[0];

                        for( size_t k=0; k<fny2_left_recv; ++k )
                        {
                            grid_to.kelem(i,j,k) = recvbuf[jsend * grid_from.sizes_[3] + k];
                        }
                    }
                }
            }
        }
        delete[] recvbuf;
    }
    MPI_Barrier( MPI_COMM_WORLD );
    
    MPI_Status status;
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

    music::dlog.Print("[MPI] Completed scatter for Fourier interpolation/copy, took %fs\n",
                        get_wtime() - tstart);  
#endif //defined(USE_MPI)      
}


#define H5_USE_16_API
#include <hdf5.h>

bool file_exists(std::string Filename)
{
    bool flag = false;
    std::fstream fin(Filename.c_str(), std::ios::in | std::ios::binary);
    if (fin.is_open())
        flag = true;
    fin.close();
    return flag;
}

void create_hdf5(std::string Filename)
{
    hid_t HDF_FileID;
    HDF_FileID =
        H5Fcreate(Filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(HDF_FileID);
}

template <typename T>
hid_t hdf5_get_data_type(void)
{
    if (typeid(T) == typeid(int))
        return H5T_NATIVE_INT;

    if (typeid(T) == typeid(unsigned))
        return H5T_NATIVE_UINT;

    if (typeid(T) == typeid(float))
        return H5T_NATIVE_FLOAT;

    if (typeid(T) == typeid(double))
        return H5T_NATIVE_DOUBLE;
    
    if (typeid(T) == typeid(long double))
        return H5T_NATIVE_LDOUBLE;

    if (typeid(T) == typeid(long long))
        return H5T_NATIVE_LLONG;

    if (typeid(T) == typeid(unsigned long long))
        return H5T_NATIVE_ULLONG;

    if (typeid(T) == typeid(size_t))
        return H5T_NATIVE_ULLONG;

    music::elog << "[HDF_IO] trying to evaluate unsupported type in GetDataType";
    return -1;
}

template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::Read_from_HDF5(const std::string Filename, const std::string ObjName)
{
    if (bdistributed)
    {
        music::elog << "Attempt to read from HDF5 into MPI-distributed array. This is not supported yet!" << std::endl;
        abort();
    }

    hid_t HDF_Type = hdf5_get_data_type<data_t>();

    hid_t HDF_FileID = H5Fopen(Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    //... save old error handler
    herr_t (*old_func)(void *);
    void *old_client_data;

    H5Eget_auto(&old_func, &old_client_data);

    //... turn off error handling by hdf5 library
    H5Eset_auto(NULL, NULL);

    //... probe dataset opening
    hid_t HDF_DatasetID = H5Dopen(HDF_FileID, ObjName.c_str());

    //... restore previous error handler
    H5Eset_auto(old_func, old_client_data);

    //... dataset did not exist or was empty
    if (HDF_DatasetID < 0)
    {
        music::elog << "Dataset \'" << ObjName.c_str() << "\' does not exist or is empty." << std::endl;
        H5Fclose(HDF_FileID);
        abort();
    }

    //... get space associated with dataset and its extensions
    hid_t HDF_DataspaceID = H5Dget_space(HDF_DatasetID);

    int ndims = H5Sget_simple_extent_ndims(HDF_DataspaceID);

    hsize_t dimsize[3];

    H5Sget_simple_extent_dims(HDF_DataspaceID, dimsize, NULL);

    hsize_t HDF_StorageSize = 1;
    for (int i = 0; i < ndims; ++i)
        HDF_StorageSize *= dimsize[i];

    //... adjust the array size to hold the data
    std::vector<data_t> Data;
    Data.reserve(HDF_StorageSize);
    Data.assign(HDF_StorageSize, (data_t)0);

    if (Data.capacity() < HDF_StorageSize)
    {
        music::elog << "Not enough memory to store all data in HDFReadDataset!" << std::endl;
        H5Sclose(HDF_DataspaceID);
        H5Dclose(HDF_DatasetID);
        H5Fclose(HDF_FileID);
        abort();
    }

    //... read the dataset
    H5Dread(HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Data[0]);

    if (Data.size() != HDF_StorageSize)
    {
        music::elog << "Something went wrong while reading!" << std::endl;
        H5Sclose(HDF_DataspaceID);
        H5Dclose(HDF_DatasetID);
        H5Fclose(HDF_FileID);
        abort();
    }

    H5Sclose(HDF_DataspaceID);
    H5Dclose(HDF_DatasetID);
    H5Fclose(HDF_FileID);

    assert(dimsize[0] == dimsize[1] && dimsize[0] == dimsize[2]);
    music::ilog << "Read external constraint data of dimensions " << dimsize[0] << "**3." << std::endl;

    for (size_t i = 0; i < 3; ++i)
        this->n_[i] = dimsize[i];
    this->space_ = rspace_id;

    this->reset();
    this->allocate();

    //... copy data to internal array ...
    real_t sum1{0.0}, sum2{0.0};
    #pragma omp parallel for reduction(+ : sum1, sum2)
    for (size_t i = 0; i < size(0); ++i)
    {
        for (size_t j = 0; j < size(1); ++j)
        {
            for (size_t k = 0; k < size(2); ++k)
            {
                this->relem(i, j, k) = Data[(i * size(1) + j) * size(2) + k];
                sum2 += std::real(this->relem(i, j, k) * this->relem(i, j, k));
                sum1 += std::real(this->relem(i, j, k));
            }
        }
    }
    sum1 /= Data.size();
    sum2 /= Data.size();
    auto stdw = std::sqrt(sum2 - sum1 * sum1);
    music::ilog << "Constraint field has <W>=" << sum1 << ", <W^2>-<W>^2=" << stdw << std::endl;

    #pragma omp parallel for reduction(+ : sum1, sum2)
    for (size_t i = 0; i < size(0); ++i)
    {
        for (size_t j = 0; j < size(1); ++j)
        {
            for (size_t k = 0; k < size(2); ++k)
            {
                this->relem(i, j, k) /= stdw;
            }
        }
    }
}

template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::Write_to_HDF5(std::string fname, std::string datasetname, const unsigned compression_level) const
{
    // FIXME: cleanup duplicate code in this function!
    if (!bdistributed && CONFIG::MPI_task_rank == 0)
    {

        hid_t file_id, dset_id;    /* file and dataset identifiers */
        hid_t filespace, memspace; /* file and memory dataspace identifiers */
        hsize_t offset[3], count[3];
        hid_t dtype_id = H5T_NATIVE_FLOAT;
        hid_t plist_id = H5P_DEFAULT;

        if (!file_exists(fname))
            create_hdf5(fname);

        file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, plist_id);

        for (int i = 0; i < 3; ++i)
            count[i] = size(i);

        if (typeid(data_t) == typeid(float))
            dtype_id = H5T_NATIVE_FLOAT;
        else if (typeid(data_t) == typeid(double))
            dtype_id = H5T_NATIVE_DOUBLE;
        else if (typeid(data_t) == typeid(long double))
            dtype_id = H5T_NATIVE_LDOUBLE;    
        else if (typeid(data_t) == typeid(std::complex<float>))
            dtype_id = H5T_NATIVE_FLOAT;
        else if (typeid(data_t) == typeid(std::complex<double>))
            dtype_id = H5T_NATIVE_DOUBLE;
        else if (typeid(data_t) == typeid(std::complex<long double>))
            dtype_id = H5T_NATIVE_LDOUBLE;

        filespace = H5Screate_simple(3, count, NULL);
        dset_id = H5Dcreate2(file_id, datasetname.c_str(), dtype_id, filespace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        hsize_t slice_sz = size(1) * size(2);

        real_t *buf = new real_t[slice_sz];

        count[0] = 1;
        count[1] = size(1);
        count[2] = size(2);

        offset[1] = 0;
        offset[2] = 0;

        memspace = H5Screate_simple(3, count, NULL);
        filespace = H5Dget_space(dset_id);

        for (size_t i = 0; i < size(0); ++i)
        {
            offset[0] = i;
            for (size_t j = 0; j < size(1); ++j)
            {
                for (size_t k = 0; k < size(2); ++k)
                {
                    if (this->space_ == rspace_id)
                        buf[j * size(2) + k] = std::real(relem(i, j, k));
                    else
                        buf[j * size(2) + k] = std::real(kelem(i, j, k));
                }
            }

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Dwrite(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, buf);
        }

        H5Sclose(filespace);
        H5Sclose(memspace);

        // H5Sclose(filespace);
        H5Dclose(dset_id);

        if (typeid(data_t) == typeid(std::complex<float>) ||
            typeid(data_t) == typeid(std::complex<double>) ||
            typeid(data_t) == typeid(std::complex<long double>) ||
            this->space_ == kspace_id)
        {
            datasetname += std::string(".im");

            for (int i = 0; i < 3; ++i)
                count[i] = size(i);

            filespace = H5Screate_simple(3, count, NULL);
            dset_id = H5Dcreate2(file_id, datasetname.c_str(), dtype_id, filespace,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(filespace);

            count[0] = 1;

            for (size_t i = 0; i < size(0); ++i)
            {
                offset[0] = i;

                for (size_t j = 0; j < size(1); ++j)
                    for (size_t k = 0; k < size(2); ++k)
                    {
                        if (this->space_ == rspace_id)
                            buf[j * size(2) + k] = std::imag(relem(i, j, k));
                        else
                            buf[j * size(2) + k] = std::imag(kelem(i, j, k));
                    }

                memspace = H5Screate_simple(3, count, NULL);
                filespace = H5Dget_space(dset_id);

                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count,
                                    NULL);

                H5Dwrite(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, buf);

                H5Sclose(memspace);
                H5Sclose(filespace);
            }

            H5Dclose(dset_id);

            delete[] buf;
        }

        H5Fclose(file_id);
        return;
    }

    if (!bdistributed && CONFIG::MPI_task_rank != 0)
        return;

    hid_t file_id, dset_id;    /* file and dataset identifiers */
    hid_t filespace, memspace; /* file and memory dataspace identifiers */
    hsize_t offset[3], count[3];
    hid_t dtype_id = H5T_NATIVE_FLOAT;
    hid_t plist_id;

#if defined(USE_MPI)

    int mpi_size, mpi_rank;

    mpi_size = MPI::get_size();
    mpi_rank = MPI::get_rank();

    if (!file_exists(fname) && mpi_rank == 0)
        create_hdf5(fname);
    MPI_Barrier(MPI_COMM_WORLD);

#if defined(USE_MPI_IO)
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    plist_id = H5P_DEFAULT;
#endif

#else

    if (!file_exists(fname))
        create_hdf5(fname);

    plist_id = H5P_DEFAULT;

#endif

#if defined(USE_MPI)
    std::vector<size_t> sizes0(MPI::get_size(), 0);
    std::vector<size_t> offsets0(MPI::get_size()+1, 0);

    MPI_Allgather((this->space_==kspace_id)? &this->local_1_start_ : &this->local_0_start_, 1, 
        MPI_UNSIGNED_LONG_LONG, &offsets0[0], 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    MPI_Allgather((this->space_==kspace_id)? &this->local_1_size_ : &this->local_0_size_, 1, 
        MPI_UNSIGNED_LONG_LONG, &sizes0[0], 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    for( int i=0; i< CONFIG::MPI_task_size; i++ ){
        if( offsets0[i+1] < offsets0[i] + sizes0[i] ) offsets0[i+1] = offsets0[i] + sizes0[i];
    }
    
#endif

#if defined(USE_MPI)
        auto loc_count = size(0), glob_count = size(0);
        MPI_Allreduce( &loc_count, &glob_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
#endif

#if defined(USE_MPI) && !defined(USE_MPI_IO)

    for (int itask = 0; itask < mpi_size; ++itask)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (itask != mpi_rank)
            continue;

#endif

        // file_id = H5Fcreate( fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT,
        // H5P_DEFAULT );
        file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, plist_id);

        for (int i = 0; i < 3; ++i)
            count[i] = size(i);

#if defined(USE_MPI)
        count[0] = glob_count;
#endif

        if (typeid(data_t) == typeid(float))
            dtype_id = H5T_NATIVE_FLOAT;
        else if (typeid(data_t) == typeid(double))
            dtype_id = H5T_NATIVE_DOUBLE;
        else if (typeid(data_t) == typeid(long double))
            dtype_id = H5T_NATIVE_LDOUBLE;
        else if (typeid(data_t) == typeid(std::complex<float>))
            dtype_id = H5T_NATIVE_FLOAT;
        else if (typeid(data_t) == typeid(std::complex<double>))
            dtype_id = H5T_NATIVE_DOUBLE;
        else if (typeid(data_t) == typeid(std::complex<long double>))
            dtype_id = H5T_NATIVE_LDOUBLE;

        hid_t dcpl = H5P_DEFAULT;
        if (size(1) >= 64) {
            dcpl = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk_dims[3] = {64, 64, 64};
            H5Pset_chunk(dcpl, 3, chunk_dims);
            H5Pset_shuffle(dcpl);
            H5Pset_deflate(dcpl, compression_level);
            H5Pset_fletcher32(dcpl);
        }

#if defined(USE_MPI) && !defined(USE_MPI_IO)
        if (itask == 0)
        {
            filespace = H5Screate_simple(3, count, NULL);
            dset_id = H5Dcreate2(file_id, datasetname.c_str(), dtype_id, filespace,
                                 H5P_DEFAULT, dcpl, H5P_DEFAULT);
            H5Sclose(filespace);
        }
        else
        {
            dset_id = H5Dopen1(file_id, datasetname.c_str());
        }
#else
    filespace = H5Screate_simple(3, count, NULL);
    dset_id = H5Dcreate2(file_id, datasetname.c_str(), dtype_id, filespace,
                         H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Sclose(filespace);
#endif

        hsize_t slice_sz = size(1) * size(2);

        real_t *buf = new real_t[slice_sz];

        count[0] = 1;
        count[1] = size(1);
        count[2] = size(2);

        offset[1] = 0;
        offset[2] = 0;

#if defined(USE_MPI) && defined(USE_MPI_IO)
        H5Pclose(plist_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#else
        plist_id = H5P_DEFAULT;
#endif

        memspace = H5Screate_simple(3, count, NULL);
        filespace = H5Dget_space(dset_id);

        for (size_t i = 0; i < size(0); ++i)
        {
#if defined(USE_MPI)
            offset[0] = offsets0[mpi_rank] + i;
#else
            offset[0] = i;
#endif

            for (size_t j = 0; j < size(1); ++j)
            {
                for (size_t k = 0; k < size(2); ++k)
                {
                    if (this->space_ == rspace_id)
                        buf[j * size(2) + k] = std::real(relem(i, j, k));
                    else
                        buf[j * size(2) + k] = std::real(kelem(i, j, k));
                }
            }

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Dwrite(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, buf);
        }

        H5Sclose(filespace);
        H5Sclose(memspace);

#if defined(USE_MPI) && defined(USE_MPI_IO)
        H5Pclose(plist_id);
#endif

        // H5Sclose(filespace);
        H5Dclose(dset_id);

        if (typeid(data_t) == typeid(std::complex<float>) ||
            typeid(data_t) == typeid(std::complex<double>) ||
            typeid(data_t) == typeid(std::complex<long double>) ||
            this->space_ == kspace_id)
        {
            datasetname += std::string(".im");

            for (int i = 0; i < 3; ++i)
                count[i] = size(i);
#if defined(USE_MPI)
            count[0] = glob_count;
#endif

#if defined(USE_MPI) && !defined(USE_MPI_IO)
            if (itask == 0)
            {
                filespace = H5Screate_simple(3, count, NULL);
                dset_id = H5Dcreate2(file_id, datasetname.c_str(), dtype_id, filespace,
                                     H5P_DEFAULT, dcpl, H5P_DEFAULT);
                H5Sclose(filespace);
            }
            else
            {
                dset_id = H5Dopen1(file_id, datasetname.c_str());
            }
#else
        filespace = H5Screate_simple(3, count, NULL);
        dset_id = H5Dcreate2(file_id, datasetname.c_str(), dtype_id, filespace,
                             H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Sclose(filespace);
#endif

#if defined(USE_MPI) && defined(USE_MPI_IO)
            // H5Pclose( plist_id );
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#else
        plist_id = H5P_DEFAULT;
#endif

            count[0] = 1;

            for (size_t i = 0; i < size(0); ++i)
            {
#if defined(USE_MPI)
                offset[0] = offsets0[mpi_rank] + i;//mpi_rank * size(0) + i;
#else
            offset[0] = i;
#endif

                for (size_t j = 0; j < size(1); ++j)
                    for (size_t k = 0; k < size(2); ++k)
                    {
                        if (this->space_ == rspace_id)
                            buf[j * size(2) + k] = std::imag(relem(i, j, k));
                        else
                            buf[j * size(2) + k] = std::imag(kelem(i, j, k));
                    }

                memspace = H5Screate_simple(3, count, NULL);
                filespace = H5Dget_space(dset_id);

                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count,
                                    NULL);

                H5Dwrite(dset_id, dtype_id, memspace, filespace, H5P_DEFAULT, buf);

                H5Sclose(memspace);
                H5Sclose(filespace);
            }

#if defined(USE_MPI) && defined(USE_MPI_IO)
            H5Pclose(plist_id);
#endif

            H5Dclose(dset_id);

            delete[] buf;
        }

        H5Fclose(file_id);
        H5Pclose(dcpl);

#if defined(USE_MPI) && !defined(USE_MPI_IO)
    }
#endif
}

#include <iomanip>

template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::Write_PDF(std::string ofname, int nbins, double scale, double vmin, double vmax)
{
    double logvmin = std::log10(vmin);
    double logvmax = std::log10(vmax);
    double idv = double(nbins) / (logvmax - logvmin);

    std::vector<double> count(nbins, 0.0), scount(nbins, 0.0);

    for (size_t ix = 0; ix < size(0); ix++)
        for (size_t iy = 0; iy < size(1); iy++)
            for (size_t iz = 0; iz < size(2); iz++)
            {
                auto v = this->relem(ix, iy, iz);
                int ibin = int((std::log10(std::fabs(v)) - logvmin) * idv);
                if (ibin >= 0 && ibin < nbins)
                {
                    count[ibin] += 1.0;
                }
                ibin = int(((std::log10((std::fabs(v) - 1.0) * scale + 1.0)) - logvmin) * idv);
                if (ibin >= 0 && ibin < nbins)
                {
                    scount[ibin] += 1.0;
                }
            }

#if defined(USE_MPI)
    if (CONFIG::MPI_task_rank == 0)
    {
#endif
        std::ofstream ofs(ofname.c_str());
        std::size_t numcells = size(0) * size(1) * size(2);

        //ofs << "# a = " << aexp << std::endl;
        ofs << "# " << std::setw(14) << "rho" << std::setw(16) << "d rho / dV" << std::setw(16) << "d (rho/D+) / dV"
            << "\n";

        for (int ibin = 0; ibin < nbins; ++ibin)
        {
            double vmean = std::pow(10.0, logvmin + (double(ibin) + 0.5) / idv);
            double dv = std::pow(10.0, logvmin + (double(ibin) + 1.0) / idv) - std::pow(10.0, logvmin + (double(ibin)) / idv);

            ofs << std::setw(16) << vmean
                << std::setw(16) << count[ibin] / dv / numcells
                << std::setw(16) << scount[ibin] / dv / numcells
                << std::endl;
        }
#if defined(USE_MPI)
    }
#endif
}

template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::Write_PowerSpectrum(std::string ofname)
{
    std::vector<double> bin_k, bin_P, bin_eP;
    std::vector<size_t> bin_count;
    this->Compute_PowerSpectrum(bin_k, bin_P, bin_eP, bin_count);
#if defined(USE_MPI)
    if (CONFIG::MPI_task_rank == 0)
    {
#endif
        std::ofstream ofs(ofname.c_str());

        ofs << "# " << std::setw(14) << "k" << std::setw(16) << "P(k)" << std::setw(16) << "err. P(k)"
            << std::setw(16) << "#modes"
            << "\n";

        for (size_t ibin = 0; ibin < bin_k.size(); ++ibin)
        {
            if (bin_count[ibin] > 0)
                ofs << std::setw(16) << bin_k[ibin]
                    << std::setw(16) << bin_P[ibin]
                    << std::setw(16) << bin_eP[ibin]
                    << std::setw(16) << bin_count[ibin]
                    << std::endl;
        }
#if defined(USE_MPI)
    }
#endif
}

template <typename data_t, bool bdistributed>
void Grid_FFT<data_t, bdistributed>::Compute_PowerSpectrum(std::vector<double> &bin_k, std::vector<double> &bin_P, std::vector<double> &bin_eP, std::vector<size_t> &bin_count)
{
    this->FourierTransformForward();

    real_t kmax = std::max(std::max(kfac_[0] * nhalf_[0], kfac_[1] * nhalf_[1]),
                           kfac_[2] * nhalf_[2]),
           kmin = std::min(std::min(kfac_[0], kfac_[1]), kfac_[2]),
           dk = kmin;

    const int nbins = kmax / kmin;

    bin_count.assign(nbins, 0);
    bin_k.assign(nbins, 0);
    bin_P.assign(nbins, 0);
    bin_eP.assign(nbins, 0);

    for (size_t ix = 0; ix < size(0); ix++)
        for (size_t iy = 0; iy < size(1); iy++)
            for (size_t iz = 0; iz < size(2); iz++)
            {
                vec3_t<double> k3 = get_k<double>(ix, iy, iz);
                double k = k3.norm();
                int idx2 = k / dk; //int((1.0f / dklog * std::log10(k / kmin)));
                auto z = this->kelem(ix, iy, iz);
                double vabs = z.real() * z.real() + z.imag() * z.imag();

                if (k >= kmin && k < kmax)
                {
                    if (iz == 0)
                    {
                        bin_P[idx2] += vabs;
                        bin_eP[idx2] += vabs * vabs;
                        bin_k[idx2] += k;
                        bin_count[idx2]++;
                    }
                    else
                    {
                        bin_P[idx2] += 2.0 * vabs;
                        bin_eP[idx2] += 2.0 * vabs * vabs;
                        bin_k[idx2] += 2.0 * k;
                        bin_count[idx2] += 2;
                    }
                }
            }

#if defined(USE_MPI)
    std::vector<double> tempv(nbins, 0);
    std::vector<size_t> tempvi(nbins, 0);

    MPI_Allreduce(reinterpret_cast<void *>(&bin_k[0]), reinterpret_cast<void *>(&tempv[0]),
                  nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    bin_k.swap(tempv);
    MPI_Allreduce(reinterpret_cast<void *>(&bin_P[0]), reinterpret_cast<void *>(&tempv[0]),
                  nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    bin_P.swap(tempv);
    MPI_Allreduce(reinterpret_cast<void *>(&bin_eP[0]), reinterpret_cast<void *>(&tempv[0]),
                  nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    bin_eP.swap(tempv);
    MPI_Allreduce(reinterpret_cast<void *>(&bin_count[0]), reinterpret_cast<void *>(&tempvi[0]),
                  nbins, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    bin_count.swap(tempvi);

#endif

    const real_t volfac(length_[0] * length_[1] * length_[2] / std::pow(2.0 * M_PI, 3.0));
    const real_t fftfac(this->fft_norm_fac_ * this->fft_norm_fac_);

    for (int i = 0; i < nbins; ++i)
    {
        if (bin_count[i] > 0)
        {
            bin_k[i] /= bin_count[i];
            bin_P[i] = bin_P[i] / bin_count[i] * volfac * fftfac;
            bin_eP[i] = std::sqrt(bin_eP[i] / bin_count[i] - bin_P[i] * bin_P[i]) / std::sqrt(bin_count[i]) * volfac * fftfac;
        }
    }
}

/********************************************************************************************/

template class Grid_FFT<real_t, true>;
template class Grid_FFT<real_t, false>;
template class Grid_FFT<ccomplex_t, true>;
template class Grid_FFT<ccomplex_t, false>;
