#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <fftw3-mpi.h>

#include "PAN_FFTW3.h"
#include "panphasia_functions.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

extern const int Nbasis;
extern const int irank_p[3][84];

extern size_t descriptor_order;
extern size_t descriptor_kk_limit;
extern size_t descriptor_base_size;

int PANPHASIA_compute_kspace_field_(size_t relative_level, ptrdiff_t N0_fourier_grid,
                                    ptrdiff_t local_n0_fourier_return, ptrdiff_t local_0_start_fourier_return,
                                    FFTW_COMPLEX *return_field)
{

  size_t copy_list[Nbasis];

  int pmax = 6;

  int nsubdivide = 21; //(pmax%2==0)?pmax+1:pmax+2;

  size_t ncopy = (pmax + 1) * (pmax + 2) * (pmax + 3) / 6;

  if (ncopy % nsubdivide != 0)
    return (100010);
  int nchunk = ncopy / nsubdivide;

  int verbose = 1;
  int flag_output_mode = 2;
  int error;
  ptrdiff_t size_to_alloc_fourier;
  ptrdiff_t size_to_alloc_pan;
  ptrdiff_t local_fourier_x_start, local_fourier_x_end;
  FFTW_PLAN output_coeff_forward_plan;

  ptrdiff_t N0_pan_grid = descriptor_base_size << relative_level;

  if (N0_fourier_grid % N0_pan_grid != 0)
    return (100015);

  int fdim = N0_fourier_grid / N0_pan_grid;
  size_t nfft_dim = N0_fourier_grid;
  size_t npan_dim = N0_pan_grid;

  int SHARED_FOUR_PAN_SPACE = (nsubdivide == 1) && (fdim == 1) && (sizeof(PAN_REAL) == sizeof(FFTW_REAL));

  ////////////////////////////////////////////////////////////////////////////////////

  if (pmax > descriptor_order)
    return (100020);

  for (size_t i = 0; i < Nbasis; i++)
    copy_list[i] = i;

  //printf("Dimensions of FT   (%td,%td,%td)\n",N0_fourier_grid,N0_fourier_grid,N0_fourier_grid);
  //printf("Dimensions of PG   (%td,%td,%td)\n",N0_pan_grid,N0_pan_grid,N0_pan_grid);
  //printf("local_no %td local_0_start_fourier %td\n",local_n0_fourier_return, local_0_start_fourier_return);

  //  Compute 1-D Spherical Bessel coefficients for each order //////////////////
  //  These are needed for the convolutions below              //////////////////
  size_t n4dimen;

  n4dimen = (nfft_dim % 4 == 0) ? 4 * (nfft_dim / 4) + 4 : 4 * (nfft_dim / 4) + 5;

  double complex *sph_bessel_coeff = FFTW_MALLOC(sizeof(double complex) * n4dimen * (pmax + 1));

  if (sph_bessel_coeff == NULL)
    return (100030);

  compute_sph_bessel_coeffs(nfft_dim, pmax, n4dimen, fdim, sph_bessel_coeff);

  //printf("Reached here! ndimen4 %ld\n",n4dimen);
  ///////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  // Determine sizes of Fourier and Panphasia coefficient 3-D arrays //

  ptrdiff_t local_n0_pan;
  ptrdiff_t local_0_start_pan;
  {

    int rank = 3;

    ptrdiff_t local_n0_fourier;
    ptrdiff_t local_0_start_fourier;
    const ptrdiff_t ndimens_alloc_fourier[] = {N0_fourier_grid, N0_fourier_grid, N0_fourier_grid + 2}; // Allocated for r2c
    ptrdiff_t howmany = nchunk;
    size_to_alloc_fourier = FFTW_MPI_LOCAL_SIZE_MANY(rank, ndimens_alloc_fourier, howmany,
                                                     FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD,
                                                     &local_n0_fourier, &local_0_start_fourier);

    if (local_0_start_fourier != local_0_start_fourier_return)
    {
      printf("Error local_0_start_fourier!=local_0_start_fourier_return\n");
      return (100032);
    };

    if (local_n0_fourier != local_n0_fourier_return)
    {
      printf("Error local_n0_fourier!=local_n0_fourier_return\n");
      return (100033);
    };

    ptrdiff_t abs_fourier_x_start, abs_fourier_x_end;

    if (local_0_start_fourier_return % fdim == 0)
    {
      abs_fourier_x_start = local_0_start_fourier_return;
    }
    else
    {
      abs_fourier_x_start = local_0_start_fourier_return + fdim - (local_0_start_fourier_return % fdim);
    };

    if ((local_0_start_fourier_return + local_n0_fourier_return - 1) % fdim == 0)
    {
      abs_fourier_x_end = local_0_start_fourier_return + local_n0_fourier_return - 1;
    }
    else
    {
      abs_fourier_x_end = (local_0_start_fourier_return + local_n0_fourier_return - 1) -
                          ((local_0_start_fourier_return + local_n0_fourier_return - 1) % fdim);
    };

    if ((abs_fourier_x_end - abs_fourier_x_start) % fdim != 0)
      return (100036);

    local_0_start_pan = abs_fourier_x_start / fdim;
    local_n0_pan = 1 + (abs_fourier_x_end - abs_fourier_x_start) / fdim;

    local_fourier_x_start = abs_fourier_x_start - local_0_start_fourier_return;
    local_fourier_x_end = abs_fourier_x_end - local_0_start_fourier_return;

    // const ptrdiff_t ndimens_alloc_pan[] = {N0_pan_grid, N0_pan_grid, N0_pan_grid+2};
    howmany = ncopy;

    size_to_alloc_pan = howmany * local_n0_pan * N0_pan_grid * (N0_pan_grid + 2);

    //printf("size_to_alloc_fourier = %td\n",size_to_alloc_fourier);
    //printf("size_to_alloc_pan     = %td\n",size_to_alloc_pan);
  };
  /////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Allocate arrays to store Panphasia coefficients and Fourier information
  // If nsubdivide==1 then use the same structure to store both.

  void *panphasia_coefficients = FFTW_MALLOC(sizeof(PAN_REAL) * size_to_alloc_pan);
  void *fourier_grids;

  if (panphasia_coefficients == NULL)
    return (100040);

  void *mode_weightings = FFTW_MALLOC(sizeof(FFTW_REAL) * size_to_alloc_fourier / nchunk);

  FFTW_REAL *ptr_mode_weightings;
  ptr_mode_weightings = mode_weightings;
  memset(ptr_mode_weightings, 0, sizeof(FFTW_REAL) * size_to_alloc_fourier / nchunk);

  FFTW_REAL *ptr_real_fourier_grid;
  FFTW_REAL *ptr_panphasia_coefficients = panphasia_coefficients;

  FFTW_COMPLEX *ptr_cmplx_fourier_grid;

  if (SHARED_FOUR_PAN_SPACE)
  {
    ptr_real_fourier_grid = panphasia_coefficients;
    ptr_cmplx_fourier_grid = panphasia_coefficients;
  }
  else
  {

    fourier_grids = FFTW_MALLOC(sizeof(FFTW_REAL) * size_to_alloc_fourier);

    if (fourier_grids == NULL)
      return (100041);

    ptr_real_fourier_grid = fourier_grids;
    ptr_cmplx_fourier_grid = fourier_grids;
  };
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////// Compute the Panphasia coefficients   ////////////////////////////////

  {

    size_t xorigin = local_0_start_pan, yorigin = 0, zorigin = 0;

    size_t xextent = local_n0_pan, yextent = N0_pan_grid, zextent = N0_pan_grid;

    if ((error = PANPHASIA_compute_coefficients_(&xorigin, &yorigin, &zorigin, &xextent, &yextent,
                                                 &zextent, copy_list, &ncopy, panphasia_coefficients,
                                                 &flag_output_mode, &verbose)))
      return (100100 + error);
  };
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////

  //////////  Output diagnostics for small grids only

  {

    if (N0_pan_grid < -1)
    {

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      char filename[100];
      sprintf(filename, "output_real_space_field.%d", rank);

      FILE *fp;
      printf("local_n0_pan %ld npan_dim %ld N0_pan_grid %ld\n", local_n0_pan,
             npan_dim, N0_pan_grid);

      fp = fopen(filename, "w");

      for (int ix = 0; ix < local_n0_pan; ix++)
        for (int iy = 0; iy < npan_dim; iy++)
          for (int iz = 0; iz < npan_dim; iz++)
          {
            int index = ix * N0_pan_grid * (N0_pan_grid + 2) + iy * (N0_pan_grid + 2) + iz;
            fprintf(fp, "%6lu%6d%6d %14.8lf %d\n", ix + local_0_start_pan, iy, iz,
                    ptr_panphasia_coefficients[index], index);
          };

      fclose(fp);
    };
  };

  //----------------------------------------------------------------------------------
  ////////////// Set up FTTW plan //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  // printf("Making the plan ... \n");

  //////////////////// Make plan for ncopy interleaved FTs ///////////////////////////
  //////////////////////////////////////////////////////////////////////////

  {
    int rank = 3;
    const ptrdiff_t ndimens[3] = {N0_fourier_grid, N0_fourier_grid, N0_fourier_grid};
    ptrdiff_t howmany = nchunk;
    ptrdiff_t block = FFTW_MPI_DEFAULT_BLOCK;
    ptrdiff_t tblock = FFTW_MPI_DEFAULT_BLOCK;
    unsigned flags = FFTW_ESTIMATE;

    output_coeff_forward_plan = FFTW_MPI_PLAN_MANY_DTF_R2C(rank, ndimens,
                                                           howmany, block, tblock,
                                                           ptr_real_fourier_grid, ptr_cmplx_fourier_grid,
                                                           MPI_COMM_WORLD, flags);
    if (output_coeff_forward_plan == NULL)
    {
      printf("Null plan\n");
      return (100051);
    };
  };

  //printf("Plan completed ... \n");

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------------------------------

  memset(return_field, 0, local_n0_fourier_return * N0_fourier_grid * (N0_fourier_grid + 2) * sizeof(FFTW_REAL));

  for (int iter = 0; iter < nsubdivide; iter++)
  {

    int moffset = iter * nchunk;

    if (!SHARED_FOUR_PAN_SPACE)
    {

      memset(ptr_real_fourier_grid, 0, sizeof(FFTW_REAL) * size_to_alloc_fourier);

      // Copy Panphasia coefficients to Fourier grid with appropriate stride - fdim

      size_t index_p, index_f;
      int m;

      for (int ix_f = local_fourier_x_start, ix_p = 0; ix_f <= local_fourier_x_end; ix_f += fdim, ix_p++)
#ifdef USE_OPENMP
#pragma omp parallel for collapse(2) private(index_p, index_f, m)
#endif

        for (int iy_f = 0, iy_p = 0; iy_f < nfft_dim; iy_f += fdim, iy_p++)
          for (int iz_f = 0, iz_p = 0; iz_f < nfft_dim; iz_f += fdim, iz_p++)
          {
            index_p = ix_p * N0_pan_grid * (N0_pan_grid + 2) + iy_p * (N0_pan_grid + 2) + iz_p;
            index_f = ix_f * N0_fourier_grid * (N0_fourier_grid + 2) + iy_f * (N0_fourier_grid + 2) + iz_f;

            for (m = 0; m < nchunk; m++)
            {
              ptr_real_fourier_grid[nchunk * index_f + m] =
                  ptr_panphasia_coefficients[ncopy * index_p + moffset + m];
            };
          };
    };

    FFTW_MPI_EXECUTE_DFT_R2C(output_coeff_forward_plan,
                             ptr_real_fourier_grid, ptr_cmplx_fourier_grid);

    {
      // Convolve and combine the FT of the Panphasia coefficient field

      size_t index1, index2;
      complex weight;
      //int m;

#ifdef USE_OPENMP
#pragma omp parallel for collapse(3) private(index1, index2, weight, m)
#endif
      for (int ix = 0; ix < local_n0_fourier_return; ix++)
        for (int iy = 0; iy < nfft_dim; iy++)
          for (int iz = 0; iz <= nfft_dim / 2; iz++)
          {
            index1 = ix * N0_fourier_grid * (N0_fourier_grid / 2 + 1) + iy * (N0_fourier_grid / 2 + 1) + iz;
            for (int m = 0; m < nchunk; m++)
            {
              index2 = nchunk * index1 + m;

              weight = sph_bessel_coeff[n4dimen * irank_p[0][m + moffset] + ix + local_0_start_fourier_return] *
                       sph_bessel_coeff[n4dimen * irank_p[1][m + moffset] + iy] *
                       sph_bessel_coeff[n4dimen * irank_p[2][m + moffset] + iz];
              return_field[index1] += weight * ptr_cmplx_fourier_grid[index2];
              ptr_mode_weightings[index1] += cabs(weight) * cabs(weight);
            };
          };
    };

  }; // End loop over iter

  // printf("Reached here 10!\n");

  //Add phase shift and normalise field
  {

    double complex phase_shift_and_scale;
    int kx, ky, kz;
    const double pi = 4.0 * atan(1.0);
    size_t index1;

    double min_weight = 100.0;

#ifdef USE_OPENMP
#pragma omp parallel for collapse(3) private(index1, kx, ky, kz, phase_shift_and_scale)
#endif
    for (int ix = 0; ix < local_n0_fourier_return; ix++)
      for (int iy = 0; iy < nfft_dim; iy++)
        for (int iz = 0; iz <= nfft_dim / 2; iz++)
        {
          index1 = ix * N0_fourier_grid * (N0_fourier_grid / 2 + 1) + iy * (N0_fourier_grid / 2 + 1) + iz;
          kx = (ix + local_0_start_fourier_return > nfft_dim / 2) ? ix + local_0_start_fourier_return - nfft_dim : ix + local_0_start_fourier_return;
          ky = (iy > nfft_dim / 2) ? iy - nfft_dim : iy;
          kz = (iz > nfft_dim / 2) ? iz - nfft_dim : iz;

          if ((kx == nfft_dim / 2) || (ky == nfft_dim / 2) || (kz == nfft_dim / 2))
          {
            // Set Nyquist modes to zero - not used by IC_Gen anyway.
            phase_shift_and_scale = 0.0;       //1.0/pow((double)nfft_dim,1.5);  // No phase shift
            ptr_mode_weightings[index1] = 0.0; // Set squared weight to zero
          }
          else
          {
            phase_shift_and_scale = sqrt((double)(fdim * fdim * fdim)) *
                                    cexp((double)fdim * (-I) * pi * (double)(kx + ky + kz) / sqrt(ptr_mode_weightings[index1]) /
                                         (double)nfft_dim) /
                                    pow((double)nfft_dim, 1.5);
          };

          return_field[index1] *= phase_shift_and_scale;
          if (ptr_mode_weightings[index1] < min_weight)
            min_weight = ptr_mode_weightings[index1];
        };

    //  printf("Minimum weight %lf\n",sqrt(min_weight));
  };

  // printf("Reached here 11!\n");

  // Rescale selected Fourier modes to unit amplitude.
  // By default this part is not executed.

  if (descriptor_kk_limit > 0)
  {
    size_t index1;
    complex weight;
    size_t ksquared;
    int kx, ky, kz;
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3) private(index1, kx, ky, kz, ksquared, weight)
#endif
    for (int ix = 0; ix < local_n0_fourier_return; ix++)
      for (int iy = 0; iy < nfft_dim; iy++)
        for (int iz = 0; iz <= nfft_dim / 2; iz++)
        {
          kx = (ix + local_0_start_fourier_return > nfft_dim / 2) ? ix + local_0_start_fourier_return - nfft_dim : ix + local_0_start_fourier_return;
          ky = (iy > nfft_dim / 2) ? iy - nfft_dim : iy;
          kz = (iz > nfft_dim / 2) ? iz - nfft_dim : iz;
          ksquared = kx * kx + ky * ky + kz * kz;
          if ((kx != nfft_dim / 2) && (ky != nfft_dim / 2) && (kz != nfft_dim / 2))
          { //Omit Nyquist modes

            if ((ksquared <= descriptor_kk_limit) && (ksquared != 0))
            {
              index1 = ix * N0_fourier_grid * (N0_fourier_grid / 2 + 1) + iy * (N0_fourier_grid / 2 + 1) + iz;
              weight = cabs(return_field[index1]);
              return_field[index1] /= weight;
            };
          };
        };
  };

  //printf("Reached here 12!\n");

  // Diagnostic text dump of the Fourier field, only built in Debug builds
  // (PANPHASIA_DEBUG_OUTPUT, set in CMakeLists.txt). The guard below was meant
  // to disable the full dump, but nfft_dim is a size_t so (nfft_dim < -1) is
  // always true and every run wrote the whole field.
#ifdef PANPHASIA_DEBUG_OUTPUT
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char filename[100];
  sprintf(filename, "output_k_space_alt.%d", rank);
  FILE *fp;

  if (nfft_dim < -1)
  {

    FILE *fp;

    fp = fopen(filename, "w");

    for (int ix = 0; ix < local_n0_fourier_return; ix++)
      for (int iy = 0; iy < nfft_dim; iy++)
        for (int iz = 0; iz <= nfft_dim / 2; iz++)
        {

          int index = ix * N0_fourier_grid * (N0_fourier_grid / 2 + 1) + iy * (N0_fourier_grid / 2 + 1) + iz;
          fprintf(fp, "%6lu%6d%6d %14.8lf %14.8lf %14.8lf \n", ix + local_0_start_fourier_return, iy, iz,
                  creal(return_field[index]), cimag(return_field[index]), sqrt(ptr_mode_weightings[index]));
          //		ptr_mode_weightings[index]);
        };
    fclose(fp);
  }
  else
  {

    fp = fopen(filename, "w");

    for (int ix = 0; ix < local_n0_fourier_return; ix++)
      for (int iy = 0; iy < nfft_dim; iy++)
        for (int iz = 0; iz <= nfft_dim / 2; iz++)
        {
          if (ix + iy + iz + local_0_start_fourier_return < 100)
          {
            int index = ix * N0_fourier_grid * (N0_fourier_grid / 2 + 1) + iy * (N0_fourier_grid / 2 + 1) + iz;
            fprintf(fp, "%6lu%6d%6d %14.8lf %14.8lf %14.8lf \n", ix + local_0_start_fourier_return, iy, iz,
                    creal(return_field[index]), cimag(return_field[index]), sqrt(ptr_mode_weightings[index]));
            //   ptr_mode_weightings[index]);
          };
        };
    fclose(fp);
  };
#endif // PANPHASIA_DEBUG_OUTPUT

  // Transpose output field

  {

    FFTW_PLAN transpose_output_plan;
    unsigned flags = FFTW_ESTIMATE;

    void *ptr_inter = return_field;

    FFTW_REAL *ptr_return_as_real_field;
    ptr_return_as_real_field = ptr_inter;

    // int rank = 2;

    //  ptrdiff_t size_to_transpose;
    ptrdiff_t howmany = N0_fourier_grid + 2;
    // ptrdiff_t local_n0, local_0_start;
    // ptrdiff_t local_n1, local_1_start;

    // const ptrdiff_t ndimens[] = {N0_fourier_grid, N0_fourier_grid};

    //  size_to_transpose = FFTW_MPI_LOCAL_SIZE_MANY_TRANSPOSED(rank, ndimens, howmany,
    //                  FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
    //                      MPI_COMM_WORLD,&local_n0, &local_0_start,&local_n1, &local_1_start);

    // printf("size_to_transpose = %td\n",size_to_transpose);

    transpose_output_plan = FFTW_MPI_PLAN_MANY_TRANSPOSE(N0_fourier_grid, N0_fourier_grid,
                                                         howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                                         ptr_return_as_real_field, ptr_return_as_real_field,
                                                         MPI_COMM_WORLD, flags);

    //printf("Transpose plan completed.\n");

    FFTW_MPI_EXECUTE_R2R(transpose_output_plan, ptr_return_as_real_field, ptr_return_as_real_field);

    //printf("Transpose completed.\n");
  };

  // Free all memory assigned by FFTW_MALLOC
  FFTW_FREE(mode_weightings);
  FFTW_FREE(sph_bessel_coeff);
  FFTW_FREE(panphasia_coefficients);

  if (!SHARED_FOUR_PAN_SPACE)
    FFTW_FREE(fourier_grids);

  FFTW_DESTROY_PLAN(output_coeff_forward_plan);

  //printf("Reached end of PANPHASIA_compute_kspace_field_\n");
  return (0);
}

//==========================================================================================
//==========================================================================================
