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
#include <complex>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <random_plugin.hh>
#include "random_music_wnoise_generator.hh"

inline double Meyer_scaling_function( double k, double kmax )
{
  constexpr double twopithirds{2.0*M_PI/3.0};
  constexpr double fourpithirds{4.0*M_PI/3.0};
  auto nu = []( double x ){ return x<0.0?0.0:(x<1.0?x:1.0); };

  k = std::abs(k)/kmax * 2 * M_PI;

  if( k < twopithirds ) return 1.0;
  else if( k< fourpithirds ){
    return std::cos( 0.5*M_PI * nu(3*k/(2*M_PI)-1.0) );
  }
  return 0.0;
}

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx)
    : res_(res), cubesize_(cubesize), ncubes_(1), baseseed_(baseseed)
{
  music::ilog.Print("Generating random numbers (1) with seed %ld", baseseed);

  initialize();
  fill_subvolume(x0, lx);
}

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(unsigned res, unsigned cubesize, long baseseed, bool zeromean)
    : res_(res), cubesize_(cubesize), ncubes_(1), baseseed_(baseseed)
{
  music::ilog.Print("Generating random numbers (2) with seed %ld", baseseed);

  double mean = 0.0;
  size_t res_l = res;

  initialize();
  mean = fill_all();

  if (zeromean)
  {
    mean = 0.0;

#pragma omp parallel for reduction(+ \
                                   : mean)
    for (int i = 0; i < (int)res_; ++i)
      for (unsigned j = 0; j < res_; ++j)
        for (unsigned k = 0; k < res_; ++k)
          mean += (*this)(i, j, k);

    mean *= 1.0 / (double)(res_l * res_l * res_l);

#pragma omp parallel for
    for (int i = 0; i < (int)res_; ++i)
      for (unsigned j = 0; j < res_; ++j)
        for (unsigned k = 0; k < res_; ++k)
          (*this)(i, j, k) = (*this)(i, j, k) - mean;
  }
}

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(unsigned res, std::string randfname, bool randsign)
    : res_(res), cubesize_(res), ncubes_(1)
{
  rnums_.push_back(new Meshvar<T>(res, 0, 0, 0));
  cubemap_[0] = 0; // create dummy map index

  std::ifstream ifs(randfname.c_str(), std::ios::binary);
  if (!ifs)
  {
    music::elog.Print("Could not open random number file \'%s\'!", randfname.c_str());
    throw std::runtime_error(std::string("Could not open random number file \'") + randfname + std::string("\'!"));
  }

  unsigned vartype;
  unsigned nx, ny, nz, blksz32;
  size_t blksz64;
  int iseed;
  //long seed;

  float sign4 = -1.0f;
  double sign8 = -1.0;

  int addrtype = 32;

  if (randsign) // use grafic2 sign convention
  {
    sign4 = 1.0f;
    sign8 = 1.0;
  }

  //... read header and check if 32bit or 64bit block size .../
  ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
  ifs.read(reinterpret_cast<char *>(&nx), sizeof(unsigned));
  if (blksz32 != 4 * sizeof(int) || nx != res_)
  {
    addrtype = 64;

    ifs.seekg(0);
    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&nx), sizeof(unsigned));

    if (blksz64 != 4 * sizeof(int) || nx != res_)
      addrtype = -1;
  }
  ifs.seekg(0);

  if (addrtype < 0)
    throw std::runtime_error("corrupt random number file");

  if (addrtype == 32)
    ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
  else
    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));

  ifs.read(reinterpret_cast<char *>(&nx), sizeof(unsigned));
  ifs.read(reinterpret_cast<char *>(&ny), sizeof(unsigned));
  ifs.read(reinterpret_cast<char *>(&nz), sizeof(unsigned));
  ifs.read(reinterpret_cast<char *>(&iseed), sizeof(int));
  //seed = (long)iseed;

  if (nx != res_ || ny != res_ || nz != res_)
  {
    char errmsg[128];
    sprintf(errmsg, "White noise file dimensions do not match level dimensions: %ux%ux%u vs. %u**3", nx, ny, nz, res_);
    throw std::runtime_error(errmsg);
  }

  if (addrtype == 32)
    ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
  else
    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));

  //... read data ...//
  //check whether random numbers are single or double precision numbers
  if (addrtype == 32)
  {
    ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
    if (blksz32 == nx * ny * sizeof(float))
      vartype = 4;
    else if (blksz32 == nx * ny * sizeof(double))
      vartype = 8;
    else
      throw std::runtime_error("corrupt random number file");
  }
  else
  {

    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
    if (blksz64 == nx * ny * sizeof(float))
      vartype = 4;
    else if (blksz64 == nx * ny * sizeof(double))
      vartype = 8;
    else
      throw std::runtime_error("corrupt random number file");
  }

  //rewind to beginning of block
  if (addrtype == 32)
    ifs.seekg(-sizeof(int), std::ios::cur);
  else
    ifs.seekg(-sizeof(size_t), std::ios::cur);

  std::vector<float> in_float;
  std::vector<double> in_double;

  music::ilog.Print("Random number file \'%s\'\n   contains %ld numbers. Reading...", randfname.c_str(), nx * ny * nz);

  long double sum = 0.0, sum2 = 0.0;
  size_t count = 0;

  //perform actual reading
  if (vartype == 4)
  {
    for (int ii = 0; ii < (int)nz; ++ii)
    {

      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }

      in_float.assign(nx * ny, 0.0f);
      ifs.read((char *)&in_float[0], nx * ny * sizeof(float));

      for (int jj = 0, q = 0; jj < (int)ny; ++jj)
        for (int kk = 0; kk < (int)nx; ++kk)
        {
          sum += in_float[q];
          sum2 += in_float[q] * in_float[q];
          ++count;

          (*rnums_[0])(kk, jj, ii) = sign4 * in_float[q++];
        }

      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }
    }
  }
  else if (vartype == 8)
  {
    for (int ii = 0; ii < (int)nz; ++ii)
    {
      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }

      in_double.assign(nx * ny, 0.0f);
      ifs.read((char *)&in_double[0], nx * ny * sizeof(double));

      for (int jj = 0, q = 0; jj < (int)ny; ++jj)
        for (int kk = 0; kk < (int)nx; ++kk)
        {
          sum += in_double[q];
          sum2 += in_double[q] * in_double[q];
          ++count;
          (*rnums_[0])(kk, jj, ii) = sign8 * in_double[q++];
        }

      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }
    }
  }

  double mean, var;
  mean = sum / count;
  var = sum2 / count - mean * mean;

  music::ilog.Print("Random numbers in file have \n     mean = %f and var = %f", mean, var);
}

//... copy construct by averaging down
template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(/*const*/ music_wnoise_generator<T> &rc, bool kdegrade)
{
  //if( res > rc.m_res || (res/rc.m_res)%2 != 0 )
  //			throw std::runtime_error("Invalid restriction in random number container copy constructor.");

  long double sum = 0.0, sum2 = 0.0;
  size_t count = 0;

  music::ilog.Print("Generating a coarse white noise field by k-space degrading");
  //... initialize properties of container
  res_ = rc.res_ / 2;
  cubesize_ = res_;
  ncubes_ = 1;
  baseseed_ = -2;

  if (sizeof(real_t) != sizeof(T))
  {
    music::elog.Print("type mismatch with real_t in k-space averaging");
    throw std::runtime_error("type mismatch with real_t in k-space averaging");
  }

  real_t
      *rfine = new real_t[(size_t)rc.res_ * (size_t)rc.res_ * 2 * ((size_t)rc.res_ / 2 + 1)],
      *rcoarse = new real_t[(size_t)res_ * (size_t)res_ * 2 * ((size_t)res_ / 2 + 1)];

  complex_t
      *ccoarse = reinterpret_cast<complex_t *>(rcoarse),
      *cfine = reinterpret_cast<complex_t *>(rfine);

  int nx(rc.res_), ny(rc.res_), nz(rc.res_), nxc(res_), nyc(res_), nzc(res_);

  FFTW_API(plan)
      pf = FFTW_API(plan_dft_r2c_3d)(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
      ipc = FFTW_API(plan_dft_c2r_3d)(nxc, nyc, nzc, ccoarse, rcoarse, FFTW_ESTIMATE);

#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
      {
        size_t q = ((size_t)i * ny + (size_t)j) * (nz + 2) + (size_t)k;
        rfine[q] = rc(i, j, k);
      }

  FFTW_API(execute)(pf);

  double fftnorm = 1.0 / ((double)nxc * (double)nyc * (double)nzc);

#pragma omp parallel for
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc / 2 + 1; k++)
      {
        int ii(i), jj(j), kk(k);

        if (i > nxc / 2)
          ii += nx / 2;
        if (j > nyc / 2)
          jj += ny / 2;

        size_t qc, qf;

        double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
        double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
        double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

        qc = ((size_t)i * nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
        qf = ((size_t)ii * ny + (size_t)jj) * (nz / 2 + 1) + (size_t)kk;

        std::complex<double> val_fine(cfine[qf][0],cfine[qf][1]);
        double phase = (kx / nxc + ky / nyc + kz / nzc) * 0.5 * M_PI;
        std::complex<double> val_phas(cos(phase), sin(phase));

        val_fine *= val_phas * fftnorm / sqrt(8.0);

        ccoarse[qc][0] = val_fine.real();
        ccoarse[qc][1] = val_fine.imag();
      }

  delete[] rfine;

  FFTW_API(execute)(ipc);

  rnums_.push_back(new Meshvar<T>(res_, 0, 0, 0));
  cubemap_[0] = 0; // map all to single array

#pragma omp parallel for reduction(+ \
                                  : sum, sum2, count)
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
      {
        size_t q = ((size_t)i * nyc + (size_t)j) * (nzc + 2) + (size_t)k;
        (*rnums_[0])(i, j, k) = rcoarse[q];
        sum += (*rnums_[0])(i, j, k);
        sum2 += (*rnums_[0])(i, j, k) * (*rnums_[0])(i, j, k);
        ++count;
      }

  delete[] rcoarse;

  FFTW_API(destroy_plan)(pf);
  FFTW_API(destroy_plan)(ipc);


  double rmean, rvar;
  rmean = sum / count;
  rvar = sum2 / count - rmean * rmean;

  music::ilog.Print("Restricted random numbers have\n       mean = %f, var = %f", rmean, rvar);
}

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(music_wnoise_generator<T> &rc, unsigned cubesize, long baseseed,
                                                  bool music2_rng, bool kspace, bool isolated, int *x0_, int *lx_, bool zeromean)
    : res_(2 * rc.res_), cubesize_(cubesize), ncubes_(1), baseseed_(baseseed)
{
  initialize();

  int x0[3], lx[3];
  if (x0_ == NULL || lx_ == NULL)
  {
    for (int i = 0; i < 3; ++i)
    {
      x0[i] = 0;
      lx[i] = res_;
    }
    fill_all();
  }
  else
  {
    for (int i = 0; i < 3; ++i)
    {
      x0[i] = x0_[i];
      lx[i] = lx_[i];
    }
    fill_subvolume(x0, lx);
  }

  if (kspace)
  {

    music::ilog.Print("Generating a constrained random number set with seed %ld\n    using coarse mode replacement...", baseseed);
    assert(lx[0] % 2 == 0 && lx[1] % 2 == 0 && lx[2] % 2 == 0);
    size_t nx = lx[0], ny = lx[1], nz = lx[2],
           nxc = lx[0] / 2, nyc = lx[1] / 2, nzc = lx[2] / 2;

    real_t *rfine = new real_t[nx * ny * (nz + 2l)];
    complex_t *cfine = reinterpret_cast<complex_t *>(rfine);

    FFTW_API(plan)
        pf = FFTW_API(plan_dft_r2c_3d)(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
        ipf = FFTW_API(plan_dft_c2r_3d)(nx, ny, nz, cfine, rfine, FFTW_ESTIMATE);

#pragma omp parallel for
    for (int i = 0; i < (int)nx; i++)
      for (int j = 0; j < (int)ny; j++)
        for (int k = 0; k < (int)nz; k++)
        {
          size_t q = ((size_t)i * (size_t)ny + (size_t)j) * (size_t)(nz + 2) + (size_t)k;
          rfine[q] = (*this)(x0[0] + i, x0[1] + j, x0[2] + k);
        }
    //this->free_all_mem();	// temporarily free memory, allocate again later

    real_t *rcoarse = new real_t[nxc * nyc * (nzc + 2)];
    complex_t *ccoarse = reinterpret_cast<complex_t *>(rcoarse);


    FFTW_API(plan) pc = FFTW_API(plan_dft_r2c_3d)(nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE);

#pragma omp parallel for
    for (int i = 0; i < (int)nxc; i++)
      for (int j = 0; j < (int)nyc; j++)
        for (int k = 0; k < (int)nzc; k++)
        {
          size_t q = ((size_t)i * (size_t)nyc + (size_t)j) * (size_t)(nzc + 2) + (size_t)k;
          rcoarse[q] = rc(x0[0] / 2 + i, x0[1] / 2 + j, x0[2] / 2 + k);
        }

    FFTW_API(execute)(pc);
    FFTW_API(execute)(pf);

    double fftnorm = 1.0 / ((double)nx * (double)ny * (double)nz);
    double sqrt8 = sqrt(8.0);
    double phasefac = -0.5; //-1.0;//-0.125;

    //if( isolated ) phasefac *= 1.5;

    // embedding of coarse white noise by fourier interpolation

#pragma omp parallel for
    for (int i = 0; i < (int)nxc; i++)
      for (int j = 0; j < (int)nyc; j++)
        for (int k = 0; k < (int)nzc / 2 + 1; k++)
        {
          int ii(i), jj(j), kk(k);

          //if( i==(int)nxc/2 ) continue;
          //if( j==(int)nyc/2 ) continue;

          if (i > (int)nxc / 2)
            ii += (int)nx / 2;
          if (j > (int)nyc / 2)
            jj += (int)ny / 2;

          size_t qc, qf;

          double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
          double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
          double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

          qc = ((size_t)i * nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
          qf = ((size_t)ii * ny + (size_t)jj) * (nz / 2 + 1) + (size_t)kk;

          std::complex<double> val(ccoarse[qc][0], ccoarse[qc][1]);
          double phase = (kx / nxc + ky / nyc + kz / nzc) * phasefac * M_PI;

          std::complex<double> val_phas(cos(phase), sin(phase));

          val *= val_phas * sqrt8;

          if(i != (int)nxc / 2 && j != (int)nyc / 2 && k != (int)nzc / 2){
            double blend_coarse_x = Meyer_scaling_function(kx, nxc / 2);
            double blend_coarse_y = Meyer_scaling_function(ky, nyc / 2);
            double blend_coarse_z = Meyer_scaling_function(kz, nzc / 2);

            double blend_coarse = blend_coarse_x*blend_coarse_y*blend_coarse_z;
            double blend_fine = std::sqrt(1.0-blend_coarse*blend_coarse);

            cfine[qf][0] = blend_fine * cfine[qf][0] + blend_coarse * val.real();
            cfine[qf][1] = blend_fine * cfine[qf][1] + blend_coarse * val.imag();
          }
        }

    delete[] rcoarse;

#pragma omp parallel for
    for (int i = 0; i < (int)nx; i++)
      for (int j = 0; j < (int)ny; j++)
        for (int k = 0; k < (int)nz / 2 + 1; k++)
        {
          size_t q = ((size_t)i * ny + (size_t)j) * (nz / 2 + 1) + (size_t)k;

          cfine[q][0] *= fftnorm;
          cfine[q][1] *= fftnorm;
        }

    FFTW_API(execute)(ipf);

#pragma omp parallel for
    for (int i = 0; i < (int)nx; i++)
      for (int j = 0; j < (int)ny; j++)
        for (int k = 0; k < (int)nz; k++)
        {
          size_t q = ((size_t)i * ny + (size_t)j) * (nz + 2) + (size_t)k;
          (*this)(x0[0] + i, x0[1] + j, x0[2] + k, false) = rfine[q];
        }

    delete[] rfine;

    FFTW_API(destroy_plan)(pf);
    FFTW_API(destroy_plan)(pc);
    FFTW_API(destroy_plan)(ipf);
  }
  else
  {
    music::ilog.Print("Generating a constrained random number set with seed %ld\n    using Hoffman-Ribak constraints...", baseseed);

    double fac = 1.0 / sqrt(8.0); //1./sqrt(8.0);

    for (int i = x0[0], ii = x0[0] / 2; i < x0[0] + lx[0]; i += 2, ++ii)
      for (int j = x0[1], jj = x0[1] / 2; j < x0[1] + lx[1]; j += 2, ++jj)
        for (int k = x0[2], kk = x0[2] / 2; k < x0[2] + lx[2]; k += 2, ++kk)
        {
          double topval = rc(ii, jj, kk);
          double locmean = 0.125 * ((*this)(i, j, k) + (*this)(i + 1, j, k) + (*this)(i, j + 1, k) + (*this)(i, j, k + 1) +
                                    (*this)(i + 1, j + 1, k) + (*this)(i + 1, j, k + 1) + (*this)(i, j + 1, k + 1) + (*this)(i + 1, j + 1, k + 1));
          double dif = fac * topval - locmean;

          (*this)(i, j, k) += dif;
          (*this)(i + 1, j, k) += dif;
          (*this)(i, j + 1, k) += dif;
          (*this)(i, j, k + 1) += dif;
          (*this)(i + 1, j + 1, k) += dif;
          (*this)(i + 1, j, k + 1) += dif;
          (*this)(i, j + 1, k + 1) += dif;
          (*this)(i + 1, j + 1, k + 1) += dif;
        }
  }
}

template <typename T>
void music_wnoise_generator<T>::register_cube(int i, int j, int k)
{
  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;
  size_t icube = ((size_t)i * ncubes_ + (size_t)j) * ncubes_ + (size_t)k;

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    rnums_.push_back(NULL);
    cubemap_[icube] = rnums_.size() - 1;
#ifdef DEBUG
    music::dlog.Print("registering new cube %d,%d,%d . ID = %ld, memloc = %ld", i, j, k, icube, cubemap_[icube]);
#endif
  }
}

template <typename T>
double music_wnoise_generator<T>::fill_cube(int i, int j, int k)
{

  gsl_rng *RNG = gsl_rng_alloc(gsl_rng_mt19937);

  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;

  size_t icube = ((size_t)i * ncubes_ + (size_t)j) * ncubes_ + (size_t)k;
  long cubeseed = baseseed_ + icube; //... each cube gets its unique seed

  gsl_rng_set(RNG, cubeseed);

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    music::elog.Print("Attempt to access non-registered random number cube!");
    throw std::runtime_error("Attempt to access non-registered random number cube!");
  }

  size_t cubeidx = it->second;

  if (rnums_[cubeidx] == NULL)
    rnums_[cubeidx] = new Meshvar<T>(cubesize_, 0, 0, 0);

  double mean = 0.0;

  for (int ii = 0; ii < (int)cubesize_; ++ii)
    for (int jj = 0; jj < (int)cubesize_; ++jj)
      for (int kk = 0; kk < (int)cubesize_; ++kk)
      {
        (*rnums_[cubeidx])(ii, jj, kk) = gsl_ran_ugaussian_ratio_method(RNG);
        mean += (*rnums_[cubeidx])(ii, jj, kk);
      }

  gsl_rng_free(RNG);

  return mean / (cubesize_ * cubesize_ * cubesize_);
}

template <typename T>
void music_wnoise_generator<T>::subtract_from_cube(int i, int j, int k, double val)
{
  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;

  size_t icube = ((size_t)i * ncubes_ + (size_t)j) * ncubes_ + (size_t)k;

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    music::elog.Print("Attempt to access unallocated RND cube %d,%d,%d in music_wnoise_generator::subtract_from_cube", i, j, k);
    throw std::runtime_error("Attempt to access unallocated RND cube in music_wnoise_generator::subtract_from_cube");
  }

  size_t cubeidx = it->second;

  for (int ii = 0; ii < (int)cubesize_; ++ii)
    for (int jj = 0; jj < (int)cubesize_; ++jj)
      for (int kk = 0; kk < (int)cubesize_; ++kk)
        (*rnums_[cubeidx])(ii, jj, kk) -= val;
}

template <typename T>
void music_wnoise_generator<T>::free_cube(int i, int j, int k)
{

  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;

  size_t icube = ((size_t)i * (size_t)ncubes_ + (size_t)j) * (size_t)ncubes_ + (size_t)k;

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    music::elog.Print("Attempt to access unallocated RND cube %d,%d,%d in music_wnoise_generator::free_cube", i, j, k);
    throw std::runtime_error("Attempt to access unallocated RND cube in music_wnoise_generator::free_cube");
  }

  size_t cubeidx = it->second;

  if (rnums_[cubeidx] != NULL)
  {
    delete rnums_[cubeidx];
    rnums_[cubeidx] = NULL;
  }
}

template <typename T>
void music_wnoise_generator<T>::initialize(void)
{

  ncubes_ = std::max((int)((double)res_ / cubesize_), 1);
  if (res_ < cubesize_)
  {
    ncubes_ = 1;
    cubesize_ = res_;
  }

  music::ilog.Print("Generating random numbers w/ sample cube size of %d", cubesize_);
}

template <typename T>
double music_wnoise_generator<T>::fill_subvolume(int *i0, int *n)
{
  int i0cube[3], ncube[3];

  i0cube[0] = (int)((double)(res_ + i0[0]) / cubesize_);
  i0cube[1] = (int)((double)(res_ + i0[1]) / cubesize_);
  i0cube[2] = (int)((double)(res_ + i0[2]) / cubesize_);

  ncube[0] = (int)(n[0] / cubesize_) + 2;
  ncube[1] = (int)(n[1] / cubesize_) + 2;
  ncube[2] = (int)(n[2] / cubesize_) + 2;

#ifdef DEBUG
  music::dlog.Print("random numbers needed for region %d,%d,%d ..+ %d,%d,%d", i0[0], i0[1], i0[2], n[0], n[1], n[2]);
  music::dlog.Print("filling cubes %d,%d,%d ..+ %d,%d,%d", i0cube[0], i0cube[1], i0cube[2], ncube[0], ncube[1], ncube[2]);
#endif

  double mean = 0.0;

  for (int i = i0cube[0]; i < i0cube[0] + ncube[0]; ++i)
    for (int j = i0cube[1]; j < i0cube[1] + ncube[1]; ++j)
      for (int k = i0cube[2]; k < i0cube[2] + ncube[2]; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        register_cube(ii, jj, kk);
      }

#pragma omp parallel for reduction(+ \
                                   : mean)
  for (int i = i0cube[0]; i < i0cube[0] + ncube[0]; ++i)
    for (int j = i0cube[1]; j < i0cube[1] + ncube[1]; ++j)
      for (int k = i0cube[2]; k < i0cube[2] + ncube[2]; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        mean += fill_cube(ii, jj, kk);
      }
  return mean / (ncube[0] * ncube[1] * ncube[2]);
}

template <typename T>
double music_wnoise_generator<T>::fill_all(void)
{
  double sum = 0.0;

  for (int i = 0; i < (int)ncubes_; ++i)
    for (int j = 0; j < (int)ncubes_; ++j)
      for (int k = 0; k < (int)ncubes_; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        register_cube(ii, jj, kk);
      }

#pragma omp parallel for reduction(+ \
                                   : sum)
  for (int i = 0; i < (int)ncubes_; ++i)
    for (int j = 0; j < (int)ncubes_; ++j)
      for (int k = 0; k < (int)ncubes_; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        sum += fill_cube(ii, jj, kk);
      }

//... subtract mean
#pragma omp parallel for reduction(+ \
                                   : sum)
  for (int i = 0; i < (int)ncubes_; ++i)
    for (int j = 0; j < (int)ncubes_; ++j)
      for (int k = 0; k < (int)ncubes_; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;
        subtract_from_cube(ii, jj, kk, sum / (ncubes_ * ncubes_ * ncubes_));
      }

  return sum / (ncubes_ * ncubes_ * ncubes_);
}

template <typename T>
void music_wnoise_generator<T>::print_allocated(void)
{
  unsigned ncount = 0, ntot = rnums_.size();
  for (size_t i = 0; i < rnums_.size(); ++i)
    if (rnums_[i] != NULL)
      ncount++;

  music::ilog.Print(" -> %d of %d random number cubes currently allocated", ncount, ntot);
}

#if defined(USE_PRECISION_FLOAT)
template class music_wnoise_generator<float>;
#elif defined(USE_PRECISION_DOUBLE)
template class music_wnoise_generator<double>;
#elif defined(USE_PRECISION_LONGDOUBLE)
template class music_wnoise_generator<long double>;
#endif
