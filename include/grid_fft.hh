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

#include <cmath>
#include <array>
#include <vector>

#include <math/vec3.hh>
#include <general.hh>
#include <bounding_box.hh>
#include <typeinfo>

/// @brief enum to indicate whether a grid is currently in real or k-space
enum space_t { kspace_id, rspace_id };

#ifdef USE_MPI
#define GRID_FFT_DISTRIBUTED true
#else
#define GRID_FFT_DISTRIBUTED false
#endif

/// @brief class for FFTable grids
/// @tparam data_t_ data type
/// @tparam bdistributed flag to indicate whether this grid is distributed in memory
template <typename data_t_, bool bdistributed=GRID_FFT_DISTRIBUTED>
class Grid_FFT
{
public:
    using data_t = data_t_; ///< data type
    static constexpr bool is_distributed_trait{bdistributed}; ///< flag to indicate whether this grid is distributed in memory

protected:
    using grid_fft_t = Grid_FFT<data_t,bdistributed>; ///< type of this grid
    
public:
    std::array<size_t, 3> n_, nhalf_;
    std::array<size_t, 4> sizes_;
    size_t npr_, npc_;
    size_t ntot_;
    std::array<real_t, 3> length_, kfac_, kny_, dx_;

    space_t space_;
    data_t *data_;
    ccomplex_t *cdata_;

    bounding_box<size_t> global_range_;

    fftw_plan_t plan_, iplan_;

    real_t fft_norm_fac_;

    bool ballocated_;

    ptrdiff_t local_0_start_, local_1_start_;
    ptrdiff_t local_0_size_, local_1_size_;

    /// @brief constructor for FTable grid object
    /// @param N number of grid points in each dimension
    /// @param L physical size of the grid in each dimension
    /// @param allocate flag to indicate whether to allocate memory for the grid
    /// @param initialspace flag to indicate whether the grid is initially in real or k-space
    Grid_FFT(const std::array<size_t, 3> &N, const std::array<real_t, 3> &L, bool allocate = true, space_t initialspace = rspace_id)
        : n_(N), length_(L), space_(initialspace), data_(nullptr), cdata_(nullptr), plan_(nullptr), iplan_(nullptr), ballocated_( false )
    {
        if( allocate ){
            this->allocate();
        }
    }

    /// @brief copy constructor [deleted] -- to avoid implicit copying of data
    Grid_FFT(const grid_fft_t &g) = delete;

    /// @brief assignment operator [deleted] -- to avoid implicit copying of data
    grid_fft_t &operator=(const grid_fft_t &g) = delete;

    /// @brief destructor
    ~Grid_FFT() { reset(); }

    /// @brief reset grid object (free memory, etc.)
    void reset()
    {
        if (data_ != nullptr)  { FFTW_API(free)(data_); data_ = nullptr; }
        if (plan_ != nullptr)  { FFTW_API(destroy_plan)(plan_); plan_ = nullptr; }
        if (iplan_ != nullptr) { FFTW_API(destroy_plan)(iplan_); iplan_ = nullptr; }
        ballocated_ = false;
        sizes_ = {0, 0, 0, 0};
        n_ = {0, 0, 0};
        nhalf_ = {0, 0, 0};
        ntot_ = 0;
    }

    /// @brief return the grid object for a given refinement level [dummy implementation for backward compatibility with MUSIC1]
    const grid_fft_t *get_grid(size_t ilevel) const { return this; }

    /// @brief return if grid object is distributed in memory
    /// @return true if grid object is distributed in memory
    bool is_distributed( void ) const noexcept { return bdistributed; }

    /// @brief allocate memory for grid object
    void allocate();

    /// @brief return if grid object is allocated
    /// @return true if grid object is allocated
    bool is_allocated( void ) const noexcept { return ballocated_; }

    //! return the number of data_t elements that we store in the container
    size_t memsize( void ) const noexcept { return ntot_; }

    //! return the (local) size of dimension i
    size_t size(size_t i) const noexcept { assert(i<4); return sizes_[i]; }

    //! return locally stored number of elements of field
    size_t local_size(void) const noexcept { return local_0_size_ * n_[1] * n_[2]; }

    //! return globally stored number of elements of field
    size_t global_size(void) const noexcept { return n_[0] * n_[1] * n_[2]; }

    //! return the (global) size of dimension i
    size_t global_size(size_t i) const noexcept { assert(i<3); return n_[i]; }

    size_t rsize( size_t i ) const noexcept { return (i==0)? local_0_size_ : n_[i]; }

    //! return a bounding box of the global extent of the field
    const bounding_box<size_t> &get_global_range(void) const noexcept
    {
        return global_range_;
    }

    bool is_nyquist_mode( size_t i, size_t j, size_t k ) const
    {
        assert( this->space_ == kspace_id );
        bool bres = (i+local_1_start_ == n_[1]/2);
        bres |= (j == n_[0]/2);
        bres |= (k == n_[2]/2);
        return bres;
    }

    //! set all field elements to zero
    void zero() noexcept
    {
        #pragma omp parallel for
        for (size_t i = 0; i < ntot_; ++i)
            data_[i] = 0.0;
    }

    void copy_from(const grid_fft_t &g)
    {
        // make sure the two fields are in the same space
        if (g.space_ != this->space_)
        {
            if (this->space_ == kspace_id)
                this->FourierTransformBackward(false);
            else
                this->FourierTransformForward(false);
        }

        // make sure the two fields have the same dimensions
        assert(this->n_[0] == g.n_[0]);
        assert(this->n_[1] == g.n_[1]);
        assert(this->n_[2] == g.n_[2]);

        // now we can copy all the data over
        #pragma omp parallel for
        for (size_t i = 0; i < ntot_; ++i)
            data_[i] = g.data_[i];
    }

    data_t &operator[](size_t i) noexcept
    {
        return data_[i];
    }

    data_t &relem(size_t i, size_t j, size_t k) noexcept 
    {
        size_t idx = (i * sizes_[1] + j) * sizes_[3] + k;
        return data_[idx];
    }

    const data_t &relem(size_t i, size_t j, size_t k) const noexcept
    {
        size_t idx = (i * sizes_[1] + j) * sizes_[3] + k;
        return data_[idx];
    }

    ccomplex_t &kelem(size_t i, size_t j, size_t k) noexcept
    {
        size_t idx = (i * sizes_[1] + j) * sizes_[3] + k;
        return cdata_[idx];
    }

    const ccomplex_t &kelem(size_t i, size_t j, size_t k) const noexcept
    {
        size_t idx = (i * sizes_[1] + j) * sizes_[3] + k;
        return cdata_[idx];
    }

    ccomplex_t &kelem(size_t idx) noexcept { return cdata_[idx]; }
    const ccomplex_t &kelem(size_t idx) const noexcept { return cdata_[idx]; }
    data_t &relem(size_t idx) noexcept { return data_[idx]; }
    const data_t &relem(size_t idx) const noexcept { return data_[idx]; }

    size_t get_idx(size_t i, size_t j, size_t k) const noexcept
    {
        return (i * sizes_[1] + j) * sizes_[3] + k;
    }

    template <typename ft>
    vec3_t<ft> get_r(const size_t i, const size_t j, const size_t k) const noexcept
    {
        vec3_t<ft> rr;

        rr[0] = real_t(i + local_0_start_) * dx_[0];
        rr[1] = real_t(j) * dx_[1];
        rr[2] = real_t(k) * dx_[2];

        return rr;
    }

    template <typename ft>
    vec3_t<ft> get_unit_r(const size_t i, const size_t j, const size_t k) const noexcept
    {
        vec3_t<ft> rr;

        rr[0] = real_t(i + local_0_start_) / real_t(n_[0]);
        rr[1] = real_t(j) / real_t(n_[1]);
        rr[2] = real_t(k) / real_t(n_[2]);

        return rr;
    }

    template <typename ft>
    vec3_t<ft> get_unit_r_shifted(const size_t i, const size_t j, const size_t k, const vec3_t<real_t> s) const noexcept
    {
        vec3_t<ft> rr;

        rr[0] = (real_t(i + local_0_start_) + s.x) / real_t(n_[0]);
        rr[1] = (real_t(j) + s.y) / real_t(n_[1]);
        rr[2] = (real_t(k) + s.z) / real_t(n_[2]);

        return rr;
    }

    vec3_t<size_t> get_cell_idx_3d(const size_t i, const size_t j, const size_t k) const noexcept
    {
        return vec3_t<size_t>({i + local_0_start_, j, k});
    }

    size_t get_cell_idx_1d(const size_t i, const size_t j, const size_t k) const noexcept
    {
        return ((i + local_0_start_) * n_[1] + j) * n_[2] + k;
    }

    //! deprecated function, was needed for old output plugin
    size_t count_leaf_cells(int, int) const noexcept
    {
        return n_[0] * n_[1] * n_[2];
    }

    real_t get_dx(int idim) const noexcept
    {
        assert(idim<3&&idim>=0);
        return dx_[idim];
    }

    const std::array<real_t, 3> &get_dx(void) const noexcept
    {
        return dx_;
    }

    template <typename ft>
    vec3_t<ft> get_k(const size_t i, const size_t j, const size_t k) const noexcept
    {
        vec3_t<ft> kk;
        if( bdistributed ){
            auto ip = i + local_1_start_;
            kk[0] = (real_t(j) - real_t(j > nhalf_[0]) * n_[0]) * kfac_[0];
            kk[1] = (real_t(ip) - real_t(ip > nhalf_[1]) * n_[1]) * kfac_[1];
        }else{
            kk[0] = (real_t(i) - real_t(i > nhalf_[0]) * n_[0]) * kfac_[0];
            kk[1] = (real_t(j) - real_t(j > nhalf_[1]) * n_[1]) * kfac_[1];
        }
        kk[2] = (real_t(k) - real_t(k > nhalf_[2]) * n_[2]) * kfac_[2];

        return kk;
    }

    template <typename ft>
    vec3_t<ft> get_k(const real_t i, const real_t j, const real_t k) const noexcept
    {
        vec3_t<ft> kk;
        if( bdistributed ){
            auto ip = i + real_t(local_1_start_);
            kk[0] = (j - real_t(j > real_t(nhalf_[0])) * n_[0]) * kfac_[0];
            kk[1] = (ip - real_t(ip > real_t(nhalf_[1])) * n_[1]) * kfac_[1];
        }else{
            kk[0] = (real_t(i) - real_t(i > real_t(nhalf_[0])) * n_[0]) * kfac_[0];
            kk[1] = (real_t(j) - real_t(j > real_t(nhalf_[1])) * n_[1]) * kfac_[1];
        }
        kk[2] = (real_t(k) - real_t(k > real_t(nhalf_[2])) * n_[2]) * kfac_[2];

        return kk;
    }

    std::array<size_t,3> get_k3(const size_t i, const size_t j, const size_t k) const noexcept
    {
        return bdistributed? std::array<size_t,3>({j,i+local_1_start_,k}) : std::array<size_t,3>({i,j,k});
    }

    data_t get_cic( const vec3_t<real_t>& v ) const noexcept
    {
        // warning! this doesn't work with MPI
        vec3_t<real_t> x({real_t(std::fmod(v.x/length_[0]+1.0,1.0)*n_[0]),
                          real_t(std::fmod(v.y/length_[1]+1.0,1.0)*n_[1]),
                          real_t(std::fmod(v.z/length_[2]+1.0,1.0)*n_[2]) });
        size_t ix = static_cast<size_t>(x.x);
        size_t iy = static_cast<size_t>(x.y);
        size_t iz = static_cast<size_t>(x.z);
        real_t dx = x.x-real_t(ix), tx = 1.0-dx;
        real_t dy = x.y-real_t(iy), ty = 1.0-dy;
        real_t dz = x.z-real_t(iz), tz = 1.0-dz;
        size_t ix1 = (ix+1)%n_[0];
        size_t iy1 = (iy+1)%n_[1];
        size_t iz1 = (iz+1)%n_[2];
        data_t val = 0.0;
        val += this->relem(ix ,iy ,iz ) * tx * ty * tz;
        val += this->relem(ix ,iy ,iz1) * tx * ty * dz;
        val += this->relem(ix ,iy1,iz ) * tx * dy * tz;
        val += this->relem(ix ,iy1,iz1) * tx * dy * dz;
        val += this->relem(ix1,iy ,iz ) * dx * ty * tz;
        val += this->relem(ix1,iy ,iz1) * dx * ty * dz;
        val += this->relem(ix1,iy1,iz ) * dx * dy * tz;
        val += this->relem(ix1,iy1,iz1) * dx * dy * dz;
        return val;
    }

    ccomplex_t get_cic_kspace( const vec3_t<real_t> x ) const noexcept
    {
        // warning! this doesn't work with MPI
        int ix = static_cast<int>(std::floor(x.x));
        int iy = static_cast<int>(std::floor(x.y));
        int iz = static_cast<int>(std::floor(x.z));
        real_t dx = x.x-real_t(ix), tx = 1.0-dx;
        real_t dy = x.y-real_t(iy), ty = 1.0-dy;
        real_t dz = x.z-real_t(iz), tz = 1.0-dz;
        size_t ix1 = (ix+1)%size(0);
        size_t iy1 = (iy+1)%size(1);
        size_t iz1 = std::min((iz+1),int(size(2))-1);
        ccomplex_t val = 0.0;
        val += this->kelem(ix ,iy ,iz ) * tx * ty * tz;
        val += this->kelem(ix ,iy ,iz1) * tx * ty * dz;
        val += this->kelem(ix ,iy1,iz ) * tx * dy * tz;
        val += this->kelem(ix ,iy1,iz1) * tx * dy * dz;
        val += this->kelem(ix1,iy ,iz ) * dx * ty * tz;
        val += this->kelem(ix1,iy ,iz1) * dx * ty * dz;
        val += this->kelem(ix1,iy1,iz ) * dx * dy * tz;
        val += this->kelem(ix1,iy1,iz1) * dx * dy * dz;
        // if( val != val ){
           //auto k = this->get_k<real_t>(ix,iy,iz);
           //std::cerr << ix << " " << iy << " " << iz << " " << val << " " <<  this->gradient(0,{ix,iy,iz}) << " " <<  this->gradient(1,{ix,iy,iz}) << " " <<  this->gradient(2,{ix,iy,iz}) << std::endl;
        // }
        return val;
    }

    
    inline ccomplex_t gradient( const int idim, std::array<size_t,3> ijk ) const
    {
        if( bdistributed ){
            ijk[0] += local_1_start_;
            std::swap(ijk[0],ijk[1]);
        }
        real_t rgrad = 
            (ijk[idim]!=nhalf_[idim])? (real_t(ijk[idim]) - real_t(ijk[idim] > nhalf_[idim]) * n_[idim]) * kfac_[idim] : 0.0; 
        return ccomplex_t(0.0,rgrad);
    }

    inline real_t laplacian( const std::array<size_t,3>& ijk ) const noexcept
    {
        return -this->get_k<real_t>(ijk[0],ijk[1],ijk[2]).norm_squared();
    }

    grid_fft_t &operator*=(data_t x)
    {
        if (space_ == kspace_id)
        {
            this->apply_function_k([&](ccomplex_t &f) { return f * x; });
        }
        else
        {
            this->apply_function_r([&](data_t &f) { return f * x; });
        }
        return *this;
    }

    grid_fft_t &operator/=(data_t x)
    {
        if (space_ == kspace_id)
        {
            this->apply_function_k([&](ccomplex_t &f) { return f / x; });
        }
        else
        {
            this->apply_function_r([&](data_t &f) { return f / x; });
        }
        return *this;
    }

    grid_fft_t &apply_Laplacian(void)
    {
        this->FourierTransformForward();
        this->apply_function_k_dep([&](auto x, auto k) {
            real_t kmod2 = k.norm_squared();
            return -x * kmod2;
        });
        this->zero_DC_mode();
        return *this;
    }

    grid_fft_t &apply_negative_Laplacian(void)
    {
        this->FourierTransformForward();
        this->apply_function_k_dep([&](auto x, auto k) {
            real_t kmod2 = k.norm_squared();
            return x * kmod2;
        });
        this->zero_DC_mode();
        return *this;
    }

    grid_fft_t &apply_InverseLaplacian(void)
    {
        this->FourierTransformForward();
        this->apply_function_k_dep([&](auto x, auto k) {
            real_t kmod2 = k.norm_squared();
            return -x / kmod2;
        });
        this->zero_DC_mode();
        return *this;
    }

    template <typename functional>
    void apply_function_k(const functional &f)
    {
#pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    auto &elem = this->kelem(i, j, k);
                    elem = f(elem);
                }
            }
        }
    }

    template <typename functional>
    void apply_function_r(const functional &f)
    {
#pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    auto &elem = this->relem(i, j, k);
                    elem = f(elem);
                }
            }
        }
    }

    real_t compute_2norm(void) const
    {
        real_t sum1{0.0};
        #pragma omp parallel for reduction(+ : sum1)
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    const auto re = std::real(this->relem(i, j, k));
                    const auto im = std::imag(this->relem(i, j, k));
                    sum1 += re * re + im * im;
                }
            }
        }

        sum1 /= sizes_[0] * sizes_[1] * sizes_[2];

        return sum1;
    }

    real_t std(void) const
    {
        double sum1{0.0}, sum2{0.0};
        size_t count{0};

        #pragma omp parallel for reduction(+ : sum1, sum2)
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    const auto elem = (space_==kspace_id)? this->kelem(i, j, k) : this->relem(i, j, k);
                    sum1 += std::real(elem);
                    sum2 += std::norm(elem);// * elem;
                }
            }
        }
        count = sizes_[0] * sizes_[1] * sizes_[2];

#ifdef USE_MPI
        if( bdistributed ){
            double globsum1{0.0}, globsum2{0.0};
            size_t globcount{0};

            MPI_Allreduce(reinterpret_cast<const void *>(&sum1),
                        reinterpret_cast<void *>(&globsum1),
                        1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            MPI_Allreduce(reinterpret_cast<const void *>(&sum2),
                        reinterpret_cast<void *>(&globsum2),
                        1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            MPI_Allreduce(reinterpret_cast<const void *>(&count),
                        reinterpret_cast<void *>(&globcount),
                        1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

            sum1 = globsum1;
            sum2 = globsum2;
            count = globcount;
        }
#endif
        sum1 /= count;
        sum2 /= count;

        return real_t(std::sqrt(sum2 - sum1 * sum1));
    }

    real_t mean(void) const
    {
        double sum1{0.0};
        size_t count{0};

        #pragma omp parallel for reduction(+ : sum1)
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    const auto elem = std::real(this->relem(i, j, k));
                    sum1 += elem;
                }
            }
        }
        count = sizes_[0] * sizes_[1] * sizes_[2];

#ifdef USE_MPI
        if( bdistributed ){
            double globsum1{0.0};
            size_t globcount{0};

            MPI_Allreduce(reinterpret_cast<const void *>(&sum1),
                        reinterpret_cast<void *>(&globsum1),
                        1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            MPI_Allreduce(reinterpret_cast<const void *>(&count),
                        reinterpret_cast<void *>(&globcount),
                        1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

            sum1 = globsum1;
            count = globcount;
        }
#endif

        sum1 /= count;

        return real_t(sum1);
    }
/*
    real_t absmax(void) const
    {
        double locmax{-1e30};

        #pragma omp parallel for reduction(max : locmax)
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    const auto elem = std::abs(this->relem(i, j, k));
                    locmax = (elem>locmax)? elem : locmax;
                }
            }
        }

#ifdef USE_MPI
        if( bdistributed ){
            double globmax{locmax};
            

            MPI_Allreduce(reinterpret_cast<const void *>(&locmax),
                        reinterpret_cast<void *>(&globmax),
                        1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            locmax  = globmax;
        }
#endif

        return real_t(locmax);
    }

    real_t max(void) const
    {
        double locmax{-1e30};

        #pragma omp parallel for reduction(max : locmax)
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    const auto elem = std::real(this->relem(i, j, k));
                    locmax = (elem>locmax)? elem : locmax;
                }
            }
        }

#ifdef USE_MPI
        if( bdistributed ){
            double globmax{locmax};
            

            MPI_Allreduce(reinterpret_cast<const void *>(&locmax),
                        reinterpret_cast<void *>(&globmax),
                        1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            locmax  = globmax;
        }
#endif

        return real_t(locmax);
    }

    real_t min(void) const
    {
        double locmin{+1e30};

        #pragma omp parallel for reduction(min : locmin)
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    const auto elem = std::real(this->relem(i, j, k));
                    locmin = (elem<locmin)? elem : locmin;
                }
            }
        }

#ifdef USE_MPI
        if( bdistributed ){
            double globmin{locmin};
            

            MPI_Allreduce(reinterpret_cast<const void *>(&locmin),
                        reinterpret_cast<void *>(&globmin),
                        1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

            locmin  = globmin;
        }
#endif

        return real_t(locmin);
    }
*/
    
    //! In real space, assigns the value of a functional of arbitrarily many grids, i.e. f(x) = f(g1(x),g2(x),...)
    template <typename functional, typename... Grids>
    void assign_function_of_grids_r(const functional &f, Grids&... grids)
    {
        list_assert_all( { ((grids.size(0)==this->size(0))&&(grids.size(1)==this->size(1))&&(grids.size(2)==this->size(2)))... } );

        #pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    this->relem(i, j, k) = f((grids.relem(i, j, k))...);
                }
            }
        }
    }

    //! In Fourier space, assigns the value of a functional of arbitrarily many grids, i.e. f(k) = f(g1(k),g2(k),...)
    template <typename functional, typename... Grids>
    void assign_function_of_grids_k(const functional &f, Grids&... grids)
    {
        list_assert_all( { ((grids.size(0)==this->size(0))&&(grids.size(1)==this->size(1))&&(grids.size(2)==this->size(2)))... } );

        #pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    this->kelem(i, j, k) = f((grids.kelem(i, j, k))...);
                }
            }
        }
    }

    //! In Fourier space, assigns the value of a functional of arbitrarily many grids where first argument is the 3d array index
    //! i.e. f[ijk] = f({ijk}, g1[ijk], g2[ijk], ...)
    template <typename functional, typename... Grids>
    void assign_function_of_grids_ijk(const functional &f, Grids&... grids)
    {
        list_assert_all( { ((grids.size(0)==this->size(0))&&(grids.size(1)==this->size(1))&&(grids.size(2)==this->size(2)))... } );

        #pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    this->kelem(i, j, k) = f({i, j, k}, (grids.kelem(i, j, k))...);
                }
            }
        }
    }

    //! In Fourier space, assigns the value of a functional of arbitrarily many grids where first argument is the k vector
    //! i.e. f(k) = f(k, g1(k), g2(k), ...)
    template <typename functional, typename... Grids>
    void assign_function_of_grids_kdep(const functional &f, Grids&... grids)
    {
        // check that all grids are same size
        list_assert_all( { ((grids.size(0)==this->size(0))&&(grids.size(1)==this->size(1))&&(grids.size(2)==this->size(2)))... } );

        #pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    this->kelem(i, j, k) = f(this->get_k<real_t>(i, j, k), (grids.kelem(i, j, k))...);
                }
            }
        }
    }

    template <typename functional>
    void apply_function_k_dep(const functional &f)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    auto &elem = this->kelem(i, j, k);
                    elem = f(elem, this->get_k<real_t>(i, j, k));
                }
            }
        }
    }

    template <typename functional>
    void apply_function_r_dep(const functional &f)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < sizes_[0]; ++i)
        {
            for (size_t j = 0; j < sizes_[1]; ++j)
            {
                for (size_t k = 0; k < sizes_[2]; ++k)
                {
                    auto &elem = this->relem(i, j, k);
                    elem = f(elem, this->get_r<real_t>(i, j, k));
                }
            }
        }
    }

    //! perform a backwards Fourier transform
    void FourierTransformBackward(bool do_transform = true);

    //! perform a forwards Fourier transform
    void FourierTransformForward(bool do_transform = true);

    //! perform a copy operation between to FFT grids that might not be of the same size
    void FourierInterpolateCopyTo( grid_fft_t &grid_to );

    //! normalise field
    void ApplyNorm(void);

    void Write_to_HDF5(std::string fname, std::string datasetname) const;

    void Read_from_HDF5( std::string fname, std::string datasetname );

    void Write_PowerSpectrum(std::string ofname);

    void Compute_PowerSpectrum(std::vector<double> &bin_k, std::vector<double> &bin_P, std::vector<double> &bin_eP, std::vector<size_t> &bin_count);

    void Write_PDF(std::string ofname, int nbins = 1000, double scale = 1.0, double rhomin = 1e-3, double rhomax = 1e3);

    void shift_field( const vec3_t<real_t>& s, bool transform_back=true )
    {
        FourierTransformForward();
        apply_function_k_dep([&](auto x, auto k) -> ccomplex_t {
            real_t shift = s.x * k[0] * this->get_dx()[0] + s.y * k[1] * this->get_dx()[1] + s.z * k[2] * this->get_dx()[2];
            return x * std::exp(ccomplex_t(0.0, shift));
        });
        if( transform_back ){
            FourierTransformBackward();
        }
    }

    void zero_DC_mode(void)
    {
        if (space_ == kspace_id)
        {
            if (CONFIG::MPI_task_rank == 0 || !bdistributed )
                cdata_[0] = (data_t)0.0;
        }
        else
        {
            data_t sum = 0.0;
            // #pragma omp parallel for reduction(+:sum)
            for (size_t i = 0; i < sizes_[0]; ++i)
            {
                for (size_t j = 0; j < sizes_[1]; ++j)
                {
                    for (size_t k = 0; k < sizes_[2]; ++k)
                    {
                        sum += this->relem(i, j, k);
                    }
                }
            }
            if( bdistributed ){
#if defined(USE_MPI)
                data_t glob_sum = 0.0;
                MPI_Allreduce(reinterpret_cast<void *>(&sum), reinterpret_cast<void *>(&glob_sum),
                            1, MPI::get_datatype<data_t>(), MPI_SUM, MPI_COMM_WORLD);
                sum = glob_sum;
#endif
            }
            sum /= sizes_[0] * sizes_[1] * sizes_[2];

#pragma omp parallel for
            for (size_t i = 0; i < sizes_[0]; ++i)
            {
                for (size_t j = 0; j < sizes_[1]; ++j)
                {
                    for (size_t k = 0; k < sizes_[2]; ++k)
                    {
                        this->relem(i, j, k) -= sum;
                    }
                }
            }
        }
    }

    void dealias(void)
    {
        static const real_t kmax[3] =
            {(real_t)(2.0 / 3.0 * (real_t)n_[0] / 2 * kfac_[0]),
             (real_t)(2.0 / 3.0 * (real_t)n_[1] / 2 * kfac_[1]),
             (real_t)(2.0 / 3.0 * (real_t)n_[2] / 2 * kfac_[2])};

        //static const real_t kmax2 = kmax*kmax;

        for (size_t i = 0; i < this->size(0); ++i)
        {
            for (size_t j = 0; j < this->size(1); ++j)
            {
                // size_t idx = (i * this->size(1) + j) * this->size(3);
                for (size_t k = 0; k < this->size(2); ++k)
                {
                    auto kk = get_k<real_t>(i, j, k);
                    //if (std::abs(kk[0]) > kmax[0] || std::abs(kk[1]) > kmax[1] || std::abs(kk[2]) > kmax[2])
                    if( kk.norm() > kmax[0] )
                        this->kelem(i,j,k) = 0.0;
                    // ++idx;
                }
            }
        }
    }
};
