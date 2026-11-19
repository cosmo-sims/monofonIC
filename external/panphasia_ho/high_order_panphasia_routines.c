#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <fftw3-mpi.h>
#include "PAN_FFTW3.h"
#include "panphasia_functions.h"
#include <gsl/gsl_sf_bessel.h>
#include <time.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

int verbose_warnings_only = 0;
static int start_panph_method = 0;
static int panphasia_rel_origin_set = 0;

// Record descriptor parameters //

size_t descriptor_order;
size_t descriptor_base_level;
size_t descriptor_xorigin, descriptor_yorigin, descriptor_zorigin;
size_t descriptor_base_size;
size_t descriptor_kk_limit;
long long int descriptor_check_digit;
char descriptor_name[100];
char full_descriptor[300];

size_t descriptor_read_in;

// Record relative coordinates for a particular descriptor

size_t rel_level;
size_t rel_origin_x, rel_origin_y, rel_origin_z;
size_t rel_coord_max;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// Matrix operations to solve for individual octree cells
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#ifdef MATRICES_ORDER_1
#include "pan_matrices_order1.h"
#elif MATRICES_ORDER_2
#include "pan_matrices_order2.h"
#elif MATRICES_ORDER_3
#include "pan_matrices_order3.h"
#elif MATRICES_ORDER_4
#include "pan_matrices_order4.h"
#elif MATRICES_ORDER_5
#include "pan_matrices_order5.h"
#elif MATRICES_ORDER_0
#include "pan_matrices_order0.h"
#else
#include "pan_matrices_order6.h"
#endif

void solve_panphasia_cell_(PAN_REAL *input_vec_parent, PAN_REAL *input_vec_children, PAN_REAL *output_vec_children, int control_flag)
{

  int iparity, iconst, irow, l, i, j;
  PAN_REAL element;
  PAN_REAL const norm = sqrt(0.125);

  PAN_REAL parent_constraint[Nbasis]; //__assume_aligned(&parent_constraint, 64);
  PAN_REAL proj_constraint[Nbasis];   //__assume_aligned(&proj_constraint, 64);
  PAN_REAL work_vec1[8 * Nbasis];     //__assume_aligned(&work_vec1, 64);
  PAN_REAL work_vec2[8 * Nbasis];     //__assume_aligned(&work_vec2, 64);

  //===========================================================================
  // Copy inputs and rearrange parent constraints in parity order and set proj_constraint to zero

  for (i = 0; i < 8 * Nbasis; i++)
    work_vec1[i] = input_vec_children[i];

  for (i = 0; i < Nbasis; i++)
  {
    parent_constraint[i] = input_vec_parent[list_by_parity[i]];
    proj_constraint[i] = 0;
  };

  if (control_flag != -999)
  { // Special value of -999 turns off linear parental constraint

    //(1)=============================================================================
    //  Compute projection of the constraints on the free child cells by parity

    for (iparity = 0, irow = 0; iparity < 8; iparity++)
    {
      for (iconst = 0; iconst < num_constr_by_parity[iparity]; iconst++, irow++)
      {
        element = 0;
        for (i = 0; i < Nbasis; i++)
        {
          element += constraint_matrices[irow][i] * work_vec1[i + Nbasis * iparity];
        };
        proj_constraint[irow] = element;
        //  printf("iparity %d irow %d element %f  \n",iparity,irow,element);
      };
    };

    //(2)==============================================================================
    // Apply parent constraint on children by parity

    for (iparity = 0, irow = 0; iparity < 8; iparity++)
    {
      for (iconst = 0; iconst < num_constr_by_parity[iparity]; iconst++, irow++)
      {
        for (i = 0; i < Nbasis; i++)
        {
          work_vec1[i + Nbasis * iparity] = work_vec1[i + Nbasis * iparity] +
                                            constraint_matrices[irow][i] * (parent_constraint[irow] - proj_constraint[irow]);
        };
      };
    };

  }; // End control_flag

  // Two versions for the 8x8 parity transformation - the latter appears faster

  /*
//(3)==== Now apply  Nbasis  8x8 parity tranformations ===========

for (l=0; l<Nbasis; l++){
 for (i=0; i<8; i++){
   element = 0;
   for (j=0; j<8; j++) element+=oct_cell_matrices[l][i][j]*work_vec1[Nbasis*j+l]; 
   work_vec2[Nbasis*i+l] = norm * element;
  }
 }

//(4) Copy back=====================================================================

for (i=0; i<8*Nbasis; i++) output_vec_children[i] = work_vec2[i]; 

  */

  // Alternate version - is it faster? Reorder data to make matrix operation faster
  // (3) - reorder the data and then evaluate the matrix multiplications

  for (i = 0, l = 0; i < Nbasis; i++)
    for (j = 0; j < 8; j++, l++)
      work_vec2[l] = work_vec1[j * Nbasis + i];

  for (l = 0; l < Nbasis; l++)
  {
    for (i = 0; i < 8; i++)
    {
      element = 0;
      for (j = 0; j < 8; j++)
        element += oct_cell_matrices[l][i][j] * work_vec2[j + 8 * l];
      work_vec1[i + 8 * l] = norm * element;
    }
  }

  for (i = 0, l = 0; i < Nbasis; i++)
    for (j = 0; j < 8; j++, l++)
      output_vec_children[j * Nbasis + i] = work_vec1[l];

  return;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// Box-Mueller transformations
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void box_muller_(PAN_REAL *unif_rand, PAN_REAL *gvar)
{
  int i, count;

  const PAN_REAL pi = 4.0 * atan(1.0);
  const PAN_REAL two_pi = 2.0 * pi;
  PAN_REAL unif_real[8 * Nbasis];  //__assume_aligned(&unif_real, 64);
  PAN_REAL gauss_real[8 * Nbasis]; //__assume_aligned(&gauss_real, 64);

  // Copy input uniform random numbers to unif_real

  for (i = 0; i < 8 * Nbasis; i++)
    unif_real[i] = unif_rand[i];

  count = 8 * Nbasis;

  //__assume(count % 8 == 0);

  for (i = 0; i < count; i += 2)
  {

#ifndef PAN_DOUBLE_PRECISION
    PAN_REAL radius = sqrtf(-2.f * logf(unif_real[i]));
    PAN_REAL angle = two_pi * unif_real[i + 1];

    gauss_real[i] = radius * cosf(angle);
    gauss_real[i + 1] = radius * sinf(angle);
#else
    PAN_REAL radius = sqrt(-2.0 * log(unif_real[i]));
    PAN_REAL angle = two_pi * unif_real[i + 1];

    gauss_real[i] = radius * cos(angle);
    gauss_real[i + 1] = radius * sin(angle);

#endif
  }

  // Copy gauss_real to output

  for (i = 0; i < 8 * Nbasis; i++)
    gvar[i] = gauss_real[i];

  return;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// Misc routines
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#define NUM_CELLS 1000000

void speed_test2_()
{
  size_t lev, j1, j2, j3;
  size_t N_cells = 1e6;

  PAN_REAL parent[Nbasis];
  // PAN_REAL child[8 * Nbasis];
  PAN_REAL output[8 * Nbasis];

  //ticks tick_start = getticks();
  //ticks tic_total;

  for (int i = 0; i < Nbasis; i++)
    parent[i] = 1.0;

  threefry4x64_test_(1);
  set_panphasia_key_(1);

  j1 = 3530439;
  j2 = 4530202;
  lev = 47;

  for (j3 = 0; j3 < N_cells; j3++)
    compute_all_properties_of_a_panphasia_cell_(&lev, &j1, &j2, &j3, parent, output);

  //tic_total = getticks() - tick_start;
  //printf("Computed %ld cells in time %.3f %s\n",N_cells,clocks_from_ticks(tic_total),clocks_getunit());
}

//===================================================================================
void compute_all_properties_of_a_panphasia_cell_(size_t *level, size_t *j1, size_t *j2, size_t *j3,
                                                 PAN_REAL *gauss_rand_parent, PAN_REAL *legendre_rand)
{

  //   void compute_all_properties_of_a_panphasia_cell_(size_t *level, size_t *j3, size_t *j1, size_t *j2//,
  //  					   PAN_REAL *gauss_rand_parent, PAN_REAL *legendre_rand) {
  PAN_REAL unif_randoms[8 * Nbasis];
  PAN_REAL gauss_rand_children[8 * Nbasis];
  size_t seed_value = 0;
  size_t allow_non_zero_seed_safety_catch = 0;

  return_uniform_pseudo_rands_threefry4x64_(*level, *j1, *j2, *j3, unif_randoms, seed_value, allow_non_zero_seed_safety_catch);

  box_muller_(unif_randoms, gauss_rand_children);

  // if (*level>17) {        //ARJ for testing purposes only
  //       for(int i=0; i<8*Nbasis;i++)gauss_rand_children[i]=0; //ARJ for testing purposes only
  //        for(int i=56; i<Nbasis;i++)gauss_rand_parent[i]=0; //ARJ for testing purposes only
  //  };

  // if (*level<23) for(int i=0; i<Nbasis;i++)gauss_rand_parent[i]=0; //ARJ for testing purposes only

  solve_panphasia_cell_(gauss_rand_parent, gauss_rand_children, legendre_rand, 0);

  // if ((*level==62)&&(*j1==0)&&(*j2==0)&&(*j3==2299)){
  // if ((*level==61)&&(*j1==0)&&(*j2==0)&&(*j3==1149)){

  //    printf("Inside compute_all_properties ...\n");
  //   printf("Level %llu Cell:(%llu,%llu,%llu)\n",*level,*j1,*j2,*j3);
  //  printf("Start of parent info %f %f %f\n",gauss_rand_parent[0],gauss_rand_parent[1],
  //	  gauss_rand_parent[2]);
  //  printf("Gauss vals:");for(int i=0;i<8;i++)printf("%f ",gauss_rand_children[Nbasis*i]);
  // printf("\n");
  //  printf("Legendre vals:");for(int i=0;i<8;i++)printf("%f ",legendre_rand[Nbasis*i]);
  //  printf("\n");
  //
  //
  // };
}

//==================================================================================
void test_random_dist_(size_t ishift)
{
  const size_t j2 = 988676, l = 62;

  size_t j1 = 409 + ishift;

  PAN_REAL unif_randoms[8 * Nbasis];

  PAN_REAL gauss_rand_children[8 * Nbasis];
  PAN_REAL legendre_rand[8 * Nbasis];
  PAN_REAL gauss_rand_parent[Nbasis];

  size_t seed_value = 0;
  size_t allow_non_zero_seed_saftey_catch = 0;

  const size_t NC = 1000000;

  double sum_squares = 0;
  double rms_value;
  long long int nrand = 0;

  long long int gauss_dist[100] = {0};
  long long int log_uniform_dist[100] = {0};
  const int array_offset = 50;
  PAN_REAL g_expected;

  printf("Generating random numbers...\n");
  threefry4x64_test_(1);
  set_panphasia_key_(1);
  check_panphasia_key_(0);

  printf("Nbasis: %d \n", Nbasis);

  for (int i = 0; i < Nbasis; i++)
    gauss_rand_parent[i] = 0;
  for (int i = 0; i < 8 * Nbasis; i++)
    gauss_rand_children[i] = 0;

  char str1[100], str2[100];

  sprintf(str1, "Gaussian_random_distribution_%lu.dat", ishift);
  sprintf(str2, "Log_uniform_random_distribution_%lu.dat", ishift);

  FILE *file = fopen(str1, "w");
  FILE *file2 = fopen(str2, "w");

  for (size_t j3 = 0; j3 < NC; j3++)
  {
    if (j3 % 10000000 == 0)
      printf("Looped over %ld\n", j3);

    return_uniform_pseudo_rands_threefry4x64_(l, j1, j2, j3, unif_randoms, seed_value, allow_non_zero_seed_saftey_catch);

    for (int i = 0; i < 8 * Nbasis; i++)
    {
      int j = -logf(unif_randoms[i]);
      log_uniform_dist[j]++;
    };

    box_muller_(unif_randoms, gauss_rand_children);
    solve_panphasia_cell_(gauss_rand_parent, gauss_rand_children, legendre_rand, -999);

    for (int i = 0; i < 8 * Nbasis; i++, nrand++)
    {

      sum_squares += pow(legendre_rand[i], 2);
      int j = (legendre_rand[i] * 5.0 + 0.5) + array_offset;
      gauss_dist[j]++;
    };
  }

  rms_value = sqrt(sum_squares / (double)nrand);

  printf("Number of rands %lld  RMS = %12.10lg   Deviation %lg \n",
         nrand, rms_value, (rms_value - 1.0) * sqrt((double)nrand));
  fprintf(file, "Number of rands %lld  RMS = %12.10lg   Deviation %lg \n",
          nrand, rms_value, (rms_value - 1.0) * sqrt((double)nrand));

  for (int i = 0; i < 100; i++)
  {

    if (gauss_dist[i] != 0)
    {
      g_expected = 0.5 * (erf(0.2 * sqrt(0.5) * ((PAN_REAL)(i - array_offset) + 0.5)) - erf(0.2 * sqrt(0.5) * ((PAN_REAL)(i - array_offset) - 0.5))) * (PAN_REAL)nrand;

      printf("%d  %lld  %f %f \n", i - array_offset, gauss_dist[i], g_expected, (gauss_dist[i] - g_expected) / sqrt(gauss_dist[i]));

      fprintf(file, "%d  %lld  %f %f \n", i - array_offset, gauss_dist[i], g_expected, (gauss_dist[i] - g_expected) / sqrt(gauss_dist[i]));
    };
    if (log_uniform_dist[i] != 0)
    {
      double expected = ((double)nrand) * (1.0 - exp(-1.0)) * exp(-(double)i);
      fprintf(file2, "%d %llu %lf %lf \n", i, log_uniform_dist[i], expected, ((double)log_uniform_dist[i] - expected) / sqrt(expected));
    };
  };

  //  fprintf(file, "%f %f %f \n", unif_randoms[i],gauss_rand_children[i],legendre_rand[i]);

  fclose(file);
  fclose(file2);
}

//===================================================================================
void speed_test_()
{
  const size_t l = 61;

  PAN_REAL unif_randoms[8 * Nbasis];

  PAN_REAL gauss_rand_children[8 * Nbasis], legendre_rand[8 * Nbasis];
  PAN_REAL gauss_rand_parent[Nbasis];

  size_t seed_value = 0;
  size_t allow_non_zero_seed_safety_catch = 0;

  printf("Generating random numbers...\n");
  threefry4x64_test_(1);
  set_panphasia_key_(1);
  check_panphasia_key_(0);

  printf("Nbasis: %d \n", Nbasis);

  for (int i = 0; i < Nbasis; i++)
    gauss_rand_parent[i] = 0; /* Set parent info to zero */

  //ticks tic;

  //ticks tic_total = getticks();

  //ticks uniform_total = 0;

  //ticks box_muller_total = 0;

  //ticks matrix_total = 0;

  size_t j1, j2, j3;

  // Time each individual call and sum

  for (j1 = 609090235558, j2 = 9090000544443, j3 = 0; j3 < NUM_CELLS; j3++, j1++, j2++)
  {

    //tic = getticks();

    return_uniform_pseudo_rands_threefry4x64_(l, j1, j2, j3, unif_randoms, seed_value, allow_non_zero_seed_safety_catch);

    //uniform_total += getticks()-tic;

    //tic = getticks();

    box_muller_(unif_randoms, gauss_rand_children);

    //box_muller_total += getticks() - tic;

    //tic = getticks();

    solve_panphasia_cell_(gauss_rand_parent, gauss_rand_children, legendre_rand, 0);

    //matrix_total += getticks() - tic;
  }

  //  printf("-------------------------------------------------\n");
  // printf(" Timings for %d cells \n", (int)NUM_CELLS);
  //printf("-------------------------------------------------\n");
  //
  //printf("Gen uniform rands:      %.3f %s \n",
  //        clocks_from_ticks(uniform_total), clocks_getunit());
  //printf("Box-Muller time:        %.3f %s.\n",clocks_from_ticks(box_muller_total), clocks_getunit());
  //
  //printf("Matrix time:            %.3f %s.\n",
  //        clocks_from_ticks(matrix_total), clocks_getunit());
  //printf("Total time:             %.3f %s.\n",
  //        clocks_from_ticks(getticks() - tic_total), clocks_getunit());
  //printf("-------------------------------------------------\n");
  //
  //printf("CPU clock frequency %llu\n",clocks_get_cpufreq());
  //
  //printf("-------------------------------------------------\n");

  FILE *file = fopen("ic_rand.dat", "w");

  /* fprintf(file, "Gaussiany Random Numbers %d\n",8*Nbasis);*/

  if (sizeof(PAN_REAL) == 4)
    for (int i = 0; i < 8 * Nbasis; i++)
      fprintf(file, "%15.12f\n", gauss_rand_children[i]);
  if (sizeof(PAN_REAL) == 8)
    for (int i = 0; i < 8 * Nbasis; i++)
      fprintf(file, "%15.12lf\n", gauss_rand_children[i]);

  /*  fprintf(file, "\nLegendre Random Numbers\n");
  fprintf(file, "-----------------------\n");
  
  for(int i=0; i<8*Nbasis; i++) fprintf(file, "%f\n", legendre_rand[i]);*/

  fclose(file);
}

//======================================================================================

void check_randoms_()
{
  const size_t j1 = 4, j2 = 9, l = 34;

  PAN_REAL unif_randoms[8 * Nbasis];

  PAN_REAL gauss_rand_children[8 * Nbasis], legendre_rand[8 * Nbasis];
  PAN_REAL gauss_rand_parent[Nbasis];

  size_t seed_value = 0;
  size_t allow_non_zero_seed_saftey_catch = 0;

  const int NC = 500; // Not too big as results output to a file ...

  printf("Generating random numbers...\n");

  threefry4x64_test_(0);
  set_panphasia_key_(0);
  check_panphasia_key_(0);

  printf("Nbasis: %d \n", Nbasis);

  // ticks tic_total = getticks();
  //ticks rng_total = 0;
  //ticks matrix_total = 0;

  FILE *file = fopen("ic_sample.dat", "w");

  for (int i = 0; i < Nbasis; i++)
    gauss_rand_parent[i] = 0.0; /* Set random parent info */
  for (int i = 0; i < 8 * Nbasis; i++)
    gauss_rand_children[i] = 0.0; /* Zero child info */
  for (int i = 0; i < 8 * Nbasis; i++)
    unif_randoms[i] = 0; /* Zero child info */

  for (size_t j3 = 0; j3 < NC; j3++)
  {
    return_uniform_pseudo_rands_threefry4x64_(l, j1, j2, j3, unif_randoms, seed_value, allow_non_zero_seed_saftey_catch);
    box_muller_(unif_randoms, gauss_rand_children);
    solve_panphasia_cell_(gauss_rand_parent, gauss_rand_children, legendre_rand, 0);

    for (int i = 0; i < 8 * Nbasis; i++)
      fprintf(file, "%14.9f %14.9f %14.9f \n", unif_randoms[i], gauss_rand_children[i], legendre_rand[i]);
  }

  fclose(file);

  check_panphasia_key_(0); // checks key is unchanged.
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// Further misc routines
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#include <gsl/gsl_sort_ulong.h>

int demo_descriptor_()
{

  char str[200] = "[Panph6,L11,(2043,2045,2046),S5,CH-999,Testing_only]"; // xyz

  //char str[200] = "[Panph6,L3,(2,3,4),S8,CH-999,Testing_only]"; // xyz

  //char str[200] = "[Panph6,L3,(4,2,3),S8,CH-999,Testing_only]"; // xyz
  //char str[200] = "[Panph6,L3,(3,4,5),S8,CH-999,Testing_only]"; // xyz

  //    char str[200] = "[Panph6,L56,(0,0,31),S5,CH-999,Testing_only]";
  //   char str[200] = "[Panph6,L21,(1136930,890765,1847934),S3,CH2414478110,Auriga_volume2]";
  //  char str[200] = "[Panph6,L21,(1136930,890765,1847934),S3,CH-999,Auriga_volume2]";

  // char copy[200];
  // const char s[20] = "[,L,(),S,CH,]";
  // char *token;

  // size_t desc_level, desc_x, desc_y, desc_z, desc_size;
  // long long int desc_ch;
  // char desc_name[100];
  // char desc_iden[8];
  int error_code;
  int pan_mode;

  descriptor_read_in = 0;

  if ((error_code = parse_and_validate_descriptor_(str, &pan_mode)))
  {

    printf("Invalid descriptor %s\n", str);
    printf("Descriptor error code %d\n", error_code);
  }
  else
  {
    printf("Valid descriptor parsed %s\n", str);
  };

  if (descriptor_read_in)
  {
    printf("-----------------------------------------\n");
    printf("Descriptor order:      %lu\n", descriptor_order);
    printf("Descriptor base level: %lu\n", descriptor_base_level);
    printf("Descriptor x-origin:   %lu\n", descriptor_xorigin);
    printf("Descriptor y-origin:   %lu\n", descriptor_yorigin);
    printf("Descriptor z-origin:   %lu\n", descriptor_zorigin);
    printf("Descriptor base size:  %lu\n", descriptor_base_size);
    printf("Descriptor check digit:%lld\n", descriptor_check_digit);
    printf("Descriptor name        %s\n", descriptor_name);
    printf("-----------------------------------------\n");
    printf("Descriptor %s\n", full_descriptor);
    printf("-----------------------------------------\n");

    printf("Check digit %lld\n", compute_check_digit_());
  };

  int verbose = 0;
  int flag_output_mode = 0;
  PANPHASIA_init_descriptor_(str, &verbose);

  size_t rel_lev = 3;

  size_t rel_orig_x = 33; //xyz
  size_t rel_orig_y = 17;
  size_t rel_orig_z = 9;

  //size_t rel_orig_x = 9;  //zxy
  //size_t rel_orig_y = 33;
  //size_t rel_orig_z = 17;

  // size_t rel_orig_x = 0;
  // size_t rel_orig_y = 0;
  // size_t rel_orig_z = 0;

  verbose = 0;

  if ((error_code = PANPHASIA_init_level_(&rel_lev,
                                          &rel_orig_x, &rel_orig_y, &rel_orig_z, &verbose)))
  {
    printf("Error %d in initialing PANPHASIA_init_level_\n",
           error_code);
    return (error_code);
  };

  size_t xstart = 3, ystart = 5, zstart = 4;
  size_t xextent = 27, yextent = 29, zextent = 40; // xyz

  //   size_t xstart = 4, ystart = 3, zstart = 5;
  //  size_t xextent = 40,  yextent = 27, zextent=29;

  //     size_t xstart = 0, ystart = 0, zstart = 0;
  //  size_t xextent = 4,  yextent = 4, zextent=4;

  size_t copy_list[Nbasis];
  size_t ncopy = 28;

  PAN_REAL *output_values = malloc(sizeof(PAN_REAL) * ncopy * xextent * yextent * zextent);
  if (output_values == NULL)
  {
    printf("Unable to allocate output_values \n");
    abort();
  };

  for (int i = 0; i < Nbasis / 3; i++)
    copy_list[i] = 3 * i;

  if ((error_code = PANPHASIA_compute_coefficients_(&xstart, &ystart, &zstart,
                                                    &xextent, &yextent, &zextent, copy_list, &ncopy,
                                                    output_values, &flag_output_mode, &verbose)))
  {

    printf("Error %d in PANPHASIA_compute_coefficients \n", error_code);
    return (error_code);
  };

  if (xextent * yextent * zextent < 2097153)
  {

    FILE *file = fopen("Panphasia_sample.tex", "w");

    for (size_t xco = 0; xco < xextent; xco++)
      for (size_t yco = 0; yco < yextent; yco++)
        for (size_t zco = 0; zco < zextent; zco++)
          fprintf(file, "%lu %lu %lu %f\n", xco, yco, zco, output_values[ncopy * (xco * yextent * zextent + yco * zextent + zco)]);

    fclose(file);
  };

  return (0);
}

int PANPHASIA_init_descriptor_(const char *descriptor, int *verbose)
{
  int error;
  int verb;

  PANPHASIA_init_descriptor_checks();

  if (*verbose != 2)
  {
    verb = 0;
  }
  else
  {
    verb = 1;
  };

  if (start_panph_method != 1)
  {
    printf("=========================================\n");
    printf("Unable to start_panphasia_method until\nstart_panphasia_method has been called.\n");
    printf("=========================================\n");
    abort();
  };

  threefry4x64_test_(verb);
  set_panphasia_key_(verb);
  check_panphasia_key_(verb);

  int pan_mode;
  if ((error = parse_and_validate_descriptor_(descriptor, &pan_mode)))
  {
    printf("-----------------------------------------\n");
    printf("Error initating start-up Panphasia routines \n");
    printf("Error code %d\n", error);
    printf("pan_mode   %d\n", pan_mode);
    printf("-----------------------------------------\n");
    abort();
  };

  if (*verbose)
    printf("Sucessfully started Panphasia with the descriptor:\n%s\n", descriptor);
  return (0);
}

/////////////////////////////////////////////////////////////////////////////////

void PANPHASIA_init_descriptor_checks()
{

  int verb = 0;

  if (start_panph_method == 1)
    return; //Only run these tests once

  start_panph_method = 1;

  // If any of these tests fail it is not safe to proceed to make initial conditions

  // Check that Threefry4x64 reproduces correct results for 20 rounds.
  // Takes 3 (counter,key) combinations and checks against the expected
  // answer for 20 rounds.

  threefry4x64_test_(verb);

  // For additional tests below, first check the inverse function for the
  // Threefry4x64 function works. The inverse only works up to and including
  // 20 rounds.

  inverse_threefry4x64_test_(verb);

  // Complex test - testing two completely different aspects at the same time.
  // The routine generates a series of random cells in the Panphasia
  // field. For each random cell it then creates two random descriptors
  // for this one cell, and computes the relative coordinates, and
  // and relative levels for the cell.  For one descriptor the cell
  // itself is chosen, while for the second descriptor the children
  // of the cell as chosen.   The test is to compute the Nbasis
  // Legendre coefficients for the random cell by directly evaluating
  // the Panphasia field and using Gaussian quadrature to evaluate the
  // integrals. The test is passed if the cell and its combined eight
  // child cells yield the same Nbasis Legendre coefficients to single/double
  // precision.   This test both that the descriptor/relative coordinates
  // do correctly point to the same cell, and that parent cell information is
  // being accurately propograted to the child cells.
  //
  // The random descriptors are chosen with a minimum side length of
  // 1 cell, up to the entire dimension of Panphasia at that level.
  //
  // The argument determines the number of tests:
  //
  //        0   - default fast test - choose the random cell at
  //              levels 0,5,10,...60. One test per level.
  //              Run time about 1.5 seconds.
  //
  //       N>0  - Do N iterations of the test with 1.
  //              In May 2020 ran with N=8000 - all tested passed.
  //              This provides a good test that the doubly periodic
  //              boundaries (of Panphasia itself, and the region
  //              covered by the descriptor) are working correctly.

  test_propogation_of_moments_(0);

  printf("===================================================\n");
  printf("Test of Threefry4x64 generator function  -  PASSED\n");
  printf("Test of inverse Threefry4x64 function    -  PASSED\n");
  printf("Test of propogation of moments           -  PASSED\n");
  printf("===================================================\n");

  panphasia_rel_origin_set = 0; // Force user to set rel origin themselves.
}

int PANPHASIA_init_level_(size_t *rel_lev,
                          size_t *rel_orig_x, size_t *rel_orig_y,
                          size_t *rel_orig_z, int *verbose)
{

  if (*rel_lev > 63)
    return (101);
  if (descriptor_base_level + *rel_lev > 63)
    return (102);

  if (*rel_orig_x >= (descriptor_base_size << *rel_lev))
    return (103);
  if (*rel_orig_y >= (descriptor_base_size << *rel_lev))
    return (104);
  if (*rel_orig_z >= (descriptor_base_size << *rel_lev))
    return (105);

  // Copy to global set of relative coordinates

  rel_level = *rel_lev;
  rel_origin_x = *rel_orig_x;
  rel_origin_y = *rel_orig_y;
  rel_origin_z = *rel_orig_z;
  rel_coord_max = descriptor_base_size << *rel_lev;

  if (*verbose)
  {
    printf("-----------------------------------------------------------------\n");
    printf("Initialising a Panphasia subgrid\n");
    printf("Relative level %lu\n", rel_level);
    printf("Relative origin (%lu,%lu,%lu)\n", rel_origin_x, rel_origin_y, rel_origin_z);
    printf("The maximum possible extent of this subgrid is %lu cells\n", rel_coord_max);
    printf("-----------------------------------------------------------------\n");
  };

  panphasia_rel_origin_set = 1;

  return (0);
}

//======================================================================================
//======================================================================================
//======================================================================================
//======================================================================================
//======================================================================================
int PANPHASIA_compute_coefficients_(size_t *xstart, size_t *ystart, size_t *zstart,
                                    size_t *xextent, size_t *yextent, size_t *zextent,
                                    size_t *copy_list,
                                    size_t *ncopy, void *output_values, int *flag_output_mode, int *verbose)
{

  size_t cumulative_cell_index[Nbasis + 1];
  size_t level_max = descriptor_base_level + rel_level;
  size_t cell_memory_to_allocate;

  //ticks tic_tot;

  //ticks tic_start = getticks();

  //==================  Basic error checking of input parameters ==========

  if (panphasia_rel_origin_set != 1)
    return (200);
  if (*xstart >= rel_coord_max)
    return (201);
  if (*ystart >= rel_coord_max)
    return (202);
  if (*zstart >= rel_coord_max)
    return (203);

  if (*xextent > rel_coord_max)
    return (204);
  if (*yextent > rel_coord_max)
    return (205);
  if (*zextent > rel_coord_max)
    return (206);

  if ((*ncopy < 0) || (*ncopy > Nbasis))
    return (207);

  if ((copy_list[0] < 0) || (copy_list[*ncopy - 1] >= Nbasis))
    return (208);

  if ((*xextent == 0) || (*yextent == 0) || (*zextent == 0))
    return (0);

  for (int i = 1; i < *ncopy; i++)
    if (copy_list[i] <= copy_list[i - 1])
      return (209);

  //=======================================================================
  // Allocate storage for one dimensional x,y,z cell coordinate lists
  //=======================================================================

  size_t nreturn_x = 2 * (*xextent) + 200;
  size_t nreturn_y = 2 * (*yextent) + 200;
  size_t nreturn_z = 2 * (*zextent) + 200;

  size_t *ret_x_list_coords = malloc(sizeof(size_t) * nreturn_x);
  if (ret_x_list_coords == NULL)
    return (220);
  size_t *ret_y_list_coords = malloc(sizeof(size_t) * nreturn_y);
  if (ret_y_list_coords == NULL)
    return (221);
  size_t *ret_z_list_coords = malloc(sizeof(size_t) * nreturn_z);
  if (ret_z_list_coords == NULL)
    return (222);

  long long int *child_pointer_x = malloc(sizeof(size_t) * 2 * nreturn_x);
  if (child_pointer_x == NULL)
    return (223);
  long long int *child_pointer_y = malloc(sizeof(size_t) * 2 * nreturn_y);
  if (child_pointer_x == NULL)
    return (224);
  long long int *child_pointer_z = malloc(sizeof(size_t) * 2 * nreturn_z);
  if (child_pointer_z == NULL)
    return (225);

  size_t level_begin_x[64], level_count_x[64];
  size_t level_begin_y[64], level_count_y[64];
  size_t level_begin_z[64], level_count_z[64];

  size_t *index_perm_x = malloc(sizeof(size_t) * nreturn_x);
  if (index_perm_x == NULL)
    return (226);
  size_t *index_perm_y = malloc(sizeof(size_t) * nreturn_y);
  if (index_perm_y == NULL)
    return (226);
  size_t *index_perm_z = malloc(sizeof(size_t) * nreturn_z);
  if (index_perm_z == NULL)
    return (226);

  size_t *list_cell_x_coord = malloc(sizeof(size_t) * (*xextent));
  if (list_cell_x_coord == NULL)
    return (227);
  size_t *list_cell_y_coord = malloc(sizeof(size_t) * (*yextent));
  if (list_cell_y_coord == NULL)
    return (228);
  size_t *list_cell_z_coord = malloc(sizeof(size_t) * (*zextent));
  if (list_cell_z_coord == NULL)
    return (229);

  //================================================================
  // Make x,y,z lists of cell coordinates //
  //================================================================
  {
    for (size_t i = 0; i < *xextent; i++)
    {
      size_t xabs, yabs, zabs;
      calc_absolute_coordinates(*xstart + i, *ystart, *zstart, &xabs, &yabs, &zabs);
      list_cell_x_coord[i] = xabs;
    };

    for (size_t i = 0; i < *yextent; i++)
    {
      size_t xabs, yabs, zabs;
      calc_absolute_coordinates(*xstart, *ystart + i, *zstart, &xabs, &yabs, &zabs);
      list_cell_y_coord[i] = yabs;
    };

    for (size_t i = 0; i < *zextent; i++)
    {
      size_t xabs, yabs, zabs;
      calc_absolute_coordinates(*xstart, *ystart, *zstart + i, &xabs, &yabs, &zabs);
      list_cell_z_coord[i] = zabs;
    };
  };
  //================================================================
  // Generate 1-D binary trees for each of the x,y,z cuboid dimensions
  //================================================================
  {
    int error_code;

    if ((error_code = return_binary_tree_cell_lists(level_max, list_cell_x_coord,
                                                    *xextent, ret_x_list_coords, nreturn_x, child_pointer_x,
                                                    level_count_x, level_begin_x, index_perm_x)))
      return (error_code);
    if ((error_code = return_binary_tree_cell_lists(level_max, list_cell_y_coord,
                                                    *yextent, ret_y_list_coords, nreturn_y, child_pointer_y,
                                                    level_count_y, level_begin_y, index_perm_y)))
      return (error_code);
    if ((error_code = return_binary_tree_cell_lists(level_max, list_cell_z_coord,
                                                    *zextent, ret_z_list_coords, nreturn_z, child_pointer_z,
                                                    level_count_z, level_begin_z, index_perm_z)))
      return (error_code);
  };
  //===================================================================
  // Allocate memory to store cell properties for coarse levels (< level_max)
  //===================================================================
  {
    size_t number_of_cells = 0;

    for (int i = (int)level_max - 1; i >= 0; i--)
    {
      cumulative_cell_index[i] = number_of_cells;
      number_of_cells += level_count_x[i] * level_count_y[i] * level_count_z[i];
    };

    if (*verbose)
      printf("Total number cells: %lu \n", number_of_cells);

    cell_memory_to_allocate = sizeof(PAN_REAL) * number_of_cells * Nbasis;
  };

  PAN_REAL *working_space = NULL;
  if (cell_memory_to_allocate > 0)
  {
    working_space = malloc(cell_memory_to_allocate);
    if (working_space == NULL)
      return (210);
  }

  //========================================================================================
  // Loop over octree starting at the root, for all relevant cells at each level
  //========================================================================================
  // size_t total_number_cells = 0;
  // size_t num_cell_compute = 0;
  // size_t num_level_max_cells = 0;
  // size_t total_num_children = 0;
  {
    size_t cell_index, j1, j2, j3;
    // size_t child_cells[8];
    size_t xoffset, yoffset, zoffset;
    size_t ix, iy, iz;
    size_t xco, yco, zco;
    size_t child_index, work_index, selected_child_index;
    size_t i;

    PAN_REAL parent[Nbasis];
    PAN_REAL children[8 * Nbasis];

    PAN_REAL *ptr_real = output_values;
    PAN_COMPLEX *ptr_cmplx = output_values;
    size_t zdimension = (*flag_output_mode == 2) ? *zextent + 2 : *zextent; // For R2C pad by two in z-dimension

    if (level_max == 0)
    {
      return_root_legendre_coefficients_(parent); // Return root cell coefficients
      size_t xloc = index_perm_x[0], yloc = index_perm_y[0], zloc = index_perm_z[0];
      size_t out_v_index = *ncopy * (xloc * (*yextent) * zdimension + yloc * zdimension + zloc);
      if (*flag_output_mode == 1)
      {
        for (i = 0; i < *ncopy; i++)
          ptr_cmplx[out_v_index + i] = (PAN_COMPLEX)parent[copy_list[i]];
      }
      else
      {
        for (i = 0; i < *ncopy; i++)
          ptr_real[out_v_index + i] = parent[copy_list[i]];
      };
    }

#ifdef USE_OPENMP
    double start, end;
    start = omp_get_wtime();
    if (*verbose)
      printf("Start ...\n");
#endif

    for (size_t level = 0; level < level_max; level++)
    {
#ifdef USE_OPENMP
#pragma omp parallel for collapse(3) private(cell_index, xoffset, yoffset, zoffset, j1, j2, j3, ix, iy, iz,   \
                                             xco, yco, zco, child_index, work_index, selected_child_index, i, \
                                             parent, children)
#endif
      for (int cell_x = 0; cell_x < level_count_x[level]; cell_x++)
        for (int cell_y = 0; cell_y < level_count_y[level]; cell_y++)
          for (int cell_z = 0; cell_z < level_count_z[level]; cell_z++)
          {

            cell_index = cumulative_cell_index[level] + cell_x * level_count_y[level] * level_count_z[level] +
                         cell_y * level_count_z[level] + cell_z;

            xoffset = level_begin_x[level] + cell_x;
            yoffset = level_begin_y[level] + cell_y;
            zoffset = level_begin_z[level] + cell_z;

            j1 = ret_x_list_coords[xoffset];
            j2 = ret_y_list_coords[yoffset];
            j3 = ret_z_list_coords[zoffset];

            if (level == 0)
            {
              return_root_legendre_coefficients_(parent); // Root cell parent information
            }
            else
            {
              for (i = 0; i < Nbasis; i++)
                parent[i] = working_space[(Nbasis * cell_index) + i]; // Copy parent information
            };

            //===================================================================================================
            compute_all_properties_of_a_panphasia_cell_(&level, &j1, &j2, &j3, parent, children);
            //===================================================================================================

            // Determine which child information needs to be stored

            for (ix = 0; ix < 2; ix++)
              for (iy = 0; iy < 2; iy++)
                for (iz = 0; iz < 2; iz++)
                {

                  if ((child_pointer_x[2 * xoffset + ix] != -1) && (child_pointer_y[2 * yoffset + iy] != -1) && (child_pointer_z[2 * zoffset + iz] != -1))
                  {

                    xco = child_pointer_x[2 * xoffset + ix] - level_begin_x[level + 1];
                    yco = child_pointer_y[2 * yoffset + iy] - level_begin_y[level + 1];
                    zco = child_pointer_z[2 * zoffset + iz] - level_begin_z[level + 1];

                    selected_child_index = Nbasis * (4 * ix + 2 * iy + iz);

                    if (level + 1 == level_max)
                    {
                      size_t xloc = index_perm_x[xco], yloc = index_perm_y[yco], zloc = index_perm_z[zco];
                      size_t out_v_index = *ncopy * (xloc * (*yextent) * zdimension + yloc * zdimension + zloc);

                      if (*flag_output_mode == 1)
                      {
                        for (i = 0; i < *ncopy; i++)
                          ptr_cmplx[out_v_index + i] = (PAN_COMPLEX)children[selected_child_index + copy_list[i]];
                      }
                      else
                      {
                        for (i = 0; i < *ncopy; i++)
                          ptr_real[out_v_index + i] = children[selected_child_index + copy_list[i]];
                      };
                    }
                    else
                    {
                      child_index = cumulative_cell_index[level + 1] + xco * level_count_y[level + 1] * level_count_z[level + 1] +
                                    yco * level_count_z[level + 1] + zco;

                      work_index = Nbasis * child_index;
                      for (i = 0; i < Nbasis; i++)
                        working_space[work_index + i] = children[selected_child_index + i];
                    };
                  };
                }; // end loop over possible children

            if (*verbose > 1)
              printf("Cell: L%lu %lu %lu %lu\n", level, j1, j2, j3);

          }; // z/y/x-coordinate/level

      //     if (flag_nochildren!=0) return(211); //All cells should have at least one child
    }; // End loop over levels

#ifdef USE_OPENMP

    end = omp_get_wtime();

    if (*verbose)
      printf("End ...\n");

    double cpu_time_used = ((double)(end - start));

    if (*verbose)
      printf("Time in OMP Section = %lf seconds \n", cpu_time_used);

#endif
  };

  //===========================================(==============================================
  // Free all memory (in order of calls to malloc above)
  //=========================================================================================
  free(ret_x_list_coords);
  free(ret_y_list_coords);
  free(ret_z_list_coords);

  free(child_pointer_x);
  free(child_pointer_y);
  free(child_pointer_z);

  free(index_perm_x);
  free(index_perm_y);
  free(index_perm_z);

  free(list_cell_x_coord);
  free(list_cell_y_coord);
  free(list_cell_z_coord);

  free(working_space);

  //tic_tot = getticks()-tic_start;

  // if (*verbose) printf("Total child cells at deepest level %llu \n",num_level_max_cells);
  //if (*verbose) printf("Total number of cells computed %llu \n",num_cell_compute);
  //if (*verbose) printf("Total number of child cells  %llu \n",total_num_children);
  //if (*verbose) printf("Time to compute %llu cells at level %llu: %.3f %s \n",num_level_max_cells,
  //		       level_max,  clocks_from_ticks(tic_tot), clocks_getunit());

  //=======================================================================================
  return (0);
  //=======================================================================================
}
//======================================================================================
//======================================================================================
//======================================================================================
//======================================================================================
//======================================================================================

int parse_and_validate_descriptor_(const char *descriptor, int *pan_mode)
{

  char *token;
  const char split[20] = "[,()]";
  char copy[300];
  size_t desc_order, desc_level, desc_x, desc_y, desc_z, desc_size;
  char desc_name[100];
  size_t desc_kk_limit = 0;
  long long int desc_ch, comp_ch;
  int kk_limit_set = 0;
  int nelement = 0;
  char descriptor_as_read[300];

  strcpy(copy, descriptor);

  token = strtok(copy, split);

  while (token != NULL)
  {
    nelement++;

    // Read in compulsory elements

    switch (nelement)
    {
    case 1:
      if (sscanf(token, "Panph%lu", &desc_order) != 1)
        return (440001);
      break;
    case 2:
      if (sscanf(token, "L%lu", &desc_level) != 1)
        return 440002;
      break;
    case 3:
      if (sscanf(token, "%lu", &desc_x) != 1)
        return 440003;
      break;
    case 4:
      if (sscanf(token, "%lu", &desc_y) != 1)
        return 440004;
      break;
    case 5:
      if (sscanf(token, "%lu", &desc_z) != 1)
        return 440005;
      break;
    case 6:
      if (sscanf(token, "S%lu", &desc_size) != 1)
        return 440005;
      break;
    case 7:
      if (sscanf(token, "KK%lu", &desc_kk_limit) == 1)
      {
        kk_limit_set = 1;
        token = strtok(NULL, split);
      }
      if (sscanf(token, "CH%lld", &desc_ch) != 1)
        return 440006;
      break;
    case 8:
      if (sscanf(token, "%199s", desc_name) != 1)
        return 440007;
      break;
    }
    token = strtok(NULL, split);
  }

  if (kk_limit_set == 0)
  {
    sprintf(descriptor_as_read, "[Panph%lu,L%lu,(%lu,%lu,%lu),S%lu,CH%lld,%s]",
            desc_order, desc_level, desc_x, desc_y, desc_z, desc_size, desc_ch, desc_name);
  }
  else
  {
    sprintf(descriptor_as_read, "[Panph%lu,L%lu,(%lu,%lu,%lu),S%lu,KK%ld,CH%lld,%s]",
            desc_order, desc_level, desc_x, desc_y, desc_z, desc_size, desc_kk_limit, desc_ch, desc_name);
  }

  if (strcmp(descriptor, descriptor_as_read))
  {
    printf("Error - descriptor mismatch\n");
    printf("As read in: %s\n", descriptor_as_read);
    printf("            %s\n", descriptor);
  }

  // Valid format descriptor has been passed - store values

  descriptor_order = desc_order;
  descriptor_base_level = desc_level;
  descriptor_xorigin = desc_x;
  descriptor_yorigin = desc_y;
  descriptor_zorigin = desc_z;
  descriptor_base_size = desc_size;
  descriptor_kk_limit = desc_kk_limit;
  descriptor_check_digit = desc_ch;
  strcpy(descriptor_name, desc_name);
  strcpy(full_descriptor, descriptor);
  descriptor_read_in = 1;

  *pan_mode = (desc_order == 1) ? 0 : 1; // 0 - Old descriptor: 1 HO descriptor

  comp_ch = compute_check_digit_(); // check the check digit

  if ((desc_ch != -999) && (desc_ch != comp_ch))
  {
    descriptor_read_in = 0;
    printf("Check digit read in %llu\n Check digit expected %llu\n", desc_ch, comp_ch);
    return (44008);
  }

  return (0);
}

void calc_absolute_coordinates(size_t xrel, size_t yrel, size_t zrel, size_t *xabs, size_t *yabs, size_t *zabs)
{

  *xabs = ((descriptor_xorigin << rel_level) + ((rel_origin_x + xrel) % rel_coord_max)) % ((size_t)1 << (descriptor_base_level + rel_level));

  *yabs = ((descriptor_yorigin << rel_level) + ((rel_origin_y + yrel) % rel_coord_max)) % ((size_t)1 << (descriptor_base_level + rel_level));

  *zabs = ((descriptor_zorigin << rel_level) + ((rel_origin_z + zrel) % rel_coord_max)) % ((size_t)1 << (descriptor_base_level + rel_level));

  // printf("descriptor_zorigin %llu rel_level %llu zrel %llu rel_origin_z %llu rel_coord_max %llu \n descriptor_base_level %llu, zabs %llu\n",
  //y	descriptor_zorigin,rel_level,zrel,rel_origin_z,rel_coord_max,descriptor_base_level,*zabs);
}

int cell_information(size_t cell_id, size_t *cumulative_cell_index, size_t *cuboid_x_dimen,
                     size_t *cuboid_y_dimen, size_t *cuboid_z_dimen, size_t *cell_lev,
                     size_t *cell_x, size_t *cell_y, size_t *cell_z, size_t number_children,
                     size_t *child_cell_indices)
{

  if (cell_id >= cumulative_cell_index[descriptor_base_level + rel_level + 1])
    return (301);

  size_t cell_level;
  for (cell_level = descriptor_base_level + rel_level;
       cell_id < cumulative_cell_index[cell_level]; cell_level--)
    ;

  size_t local_id = cell_id - cumulative_cell_index[cell_level];

  *cell_x = local_id / (cuboid_y_dimen[cell_level] * cuboid_z_dimen[cell_level]);
  *cell_y = (local_id - *cell_x * cuboid_y_dimen[cell_level] * cuboid_z_dimen[cell_level]) / cuboid_z_dimen[cell_level];
  *cell_z = local_id % cuboid_z_dimen[cell_level];

  //printf("Cell level %llu x %llu y %llu z %llu\n",cell_level,*cell_x,*cell_y,*cell_z);

  return (0);
}

int return_binary_tree_cell_lists(size_t level_max, size_t *list_cell_coordinates,
                                  size_t extent, size_t *return_tree_list_coordinates, size_t nreturn,
                                  long long int *child_pointer, size_t *level_count,
                                  size_t *level_begin, size_t *index_perm)
{

  if (extent == 0)
    return (401);
  if (nreturn < 2 * extent + 192)
    return (402);

  for (size_t i = 0; i < 2 * nreturn; i++)
    child_pointer[i] = -1;

  {
    size_t stride = 1;
    for (size_t i = 0; i < extent; i++)
    {
      index_perm[i] = i;
      return_tree_list_coordinates[i] = list_cell_coordinates[i];
    };
    gsl_sort2_ulong(return_tree_list_coordinates, stride, index_perm, stride, extent);
  }

  //----------------------------------------------------------------------------
  level_begin[level_max] = 0;
  level_count[level_max] = extent;

  size_t offset, counter;
  size_t abs_coord;

  for (size_t level = level_max; level > 0; level--)
  {

    offset = level_begin[level] + level_count[level];
    counter = 0;

    abs_coord = return_tree_list_coordinates[level_begin[level]];

    return_tree_list_coordinates[offset] = abs_coord >> 1;
    child_pointer[2 * offset + abs_coord % 2] = level_begin[level];

    for (size_t cell = 1; cell < level_count[level]; cell++)
    {

      abs_coord = return_tree_list_coordinates[level_begin[level] + cell];

      if (abs_coord >> 1 == return_tree_list_coordinates[offset + counter])
      {
        child_pointer[2 * offset + 2 * counter + abs_coord % 2] = level_begin[level] + cell;
      }
      else
      {
        counter++;
        return_tree_list_coordinates[offset + counter] = abs_coord >> 1;
        child_pointer[2 * offset + 2 * counter + abs_coord % 2] = level_begin[level] + cell;
      }

    } //cell loop

    level_count[level - 1] = ++counter;
    level_begin[level - 1] = level_begin[level] + level_count[level];
  }; // level loop

  return (0);
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// Test code for checking the appropriate moments are preserved
// between levels in Panphasia
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#include <gsl/gsl_sf_legendre.h>

void integrate_cell(int, int, int, size_t, size_t, size_t, FFTW_REAL *, double *);
int compute_panphasia_(double, double, double, size_t, size_t, size_t, FFTW_REAL *, double *);

void test_cell_moments(const char *, size_t, size_t, size_t, size_t, size_t, double *);

////////////////////////////////////////////////////////////////////////////////
void test_moments_()
{

  // int lev = 10;
  char descriptor_demo[300] = "Hello!";
  printf("Demo string %s\n", descriptor_demo);

  //  descriptor_pair_generate_();//, descriptor_demo);
  printf("Parameters: %s\n", descriptor_demo);

  // size_t nlevel = 1;
  double coefficients1[Nbasis];
  double coefficients2[Nbasis];

  double max_diff2 = 0.0;
  double rms_diff2 = 0.0;

  char descriptor[200];

  size_t const xco_full = 0x7504f333f9de6497;
  size_t const yco_full = 0x67ea73c992a3355c;
  size_t const zco_full = 0x5ab50a5892e98768;

  size_t xco = 0;
  size_t yco = 0;
  size_t zco = 0;

  verbose_warnings_only = 1; // Minimize output to screen.

  for (size_t level = 0; level < 63; level++)
  {

    xco = (xco_full) >> (63 - level);
    yco = (yco_full) >> (63 - level);
    zco = (zco_full) >> (63 - level);

    sprintf(descriptor, "[Panph6,L%ld,(%lu,%lu,%lu),S1,CH-999,test]", level, xco, yco, zco);
    // printf("%s\n",descriptor);

    test_cell_moments(descriptor, 0, 0, 0, 0, 1, coefficients1);

    test_cell_moments(descriptor, 1, 0, 0, 0, 2, coefficients2);

    for (int i = 0; i < Nbasis; i++)
    {
      double diff2 = pow(coefficients2[i] - coefficients1[i], 2);
      if (diff2 > max_diff2)
        max_diff2 = diff2;
      rms_diff2 += diff2;
    }

    rms_diff2 /= (double)Nbasis;

    //  for (int i=0; i<Nbasis; i++) printf("X: %d %16.10f %16.10f %e \n",i,
    //   coefficients1[i],coefficients2[i],coefficients2[i]-coefficients1[i]);

    if (sizeof(PAN_REAL) == 4)
    {

      if ((max_diff2 > 1.e-12) || (rms_diff2 > 1.e-12))
      {
        printf("Moments not accurately recovered at single precision\n");
        abort();
      }
    }
    else
    {

      if ((max_diff2 > 1.e-24) || (rms_diff2 > 1.e-24))
      {
        printf("Moments not accurately recovered at double precision\n");
        abort();
      }
    }

    //printf("Acceptable differences:  %e   RMS difference %e\n",sqrt(max_diff2),sqrt(rms_diff2));
  }

  printf("Completed moment test successfully.\n");
}

void test_cell_moments(const char *root_descriptor, size_t rel_lev, size_t rel_orig_x,
                       size_t rel_orig_y, size_t rel_orig_z, size_t extent, double *coeff)
{

  int error_code;
  int verbose = 0;
  int flag_output_mode = 0;

  PANPHASIA_init_descriptor_(root_descriptor, &verbose);

  verbose = 0;

  if ((error_code = PANPHASIA_init_level_(&rel_lev,
                                          &rel_orig_x, &rel_orig_y, &rel_orig_z, &verbose)))
  {
    printf("Error %d in initialing PANPHASIA_init_level_\n",
           error_code);
  }

  size_t xstart = 0, ystart = 0, zstart = 0;

  size_t xextent, yextent, zextent;

  xextent = extent;
  yextent = extent;
  zextent = extent;
  size_t copy_list[Nbasis];
  for (int i = 0; i < Nbasis; i++)
    copy_list[i] = i;
  size_t ncopy = Nbasis;

  FFTW_REAL *output_values = malloc(sizeof(FFTW_REAL) * ncopy * xextent * yextent * zextent);
  if (output_values == NULL)
  {
    printf("Unable to allocate output_values \n");
    abort();
  }

  if ((error_code = PANPHASIA_compute_coefficients_(&xstart, &ystart, &zstart,
                                                    &xextent, &yextent, &zextent, copy_list, &ncopy,
                                                    output_values, &flag_output_mode, &verbose)))
  {

    printf("Error %d in PANPHASIA_compute_coefficients_ \n", error_code);
  }

  /*  // Plot Panphasia field ...

   double panphasia_value;
     FILE *file = fopen("Panphasia_function.tex","w");
     for (int i=0; i<1000; i++){
       double x = 0.0009*(double)(xextent*(i+1)); double y = 0.5; double z=0.5;
       compute_panphasia_(x,y,z,xextent,yextent,zextent,output_values,&panphasia_value);
       fprintf(file,"%12.6f %14.8f\n",x,panphasia_value);
     }
     fclose(file);

     printf("Finished writing file ...\n");
   
   */

  double sum_coefficients[Nbasis];
  for (size_t i = 0; i < Nbasis; ++i)
    sum_coefficients[i] = 0.0;
  double results[Nbasis];
  size_t num_cells = 0;

  for (int ix = 0; ix < extent; ix++)
    for (int iy = 0; iy < extent; iy++)
      for (int iz = 0; iz < extent; iz++)
      {

        integrate_cell(ix, iy, iz, xextent, yextent, zextent, output_values, results);

        for (int i = 0; i < Nbasis; i++)
          sum_coefficients[i] += results[i];
        num_cells++;
      };

  for (int i = 0; i < Nbasis; i++)
    coeff[i] = sum_coefficients[i] / sqrt(xextent * yextent * zextent);

  free(output_values);
}

void integrate_cell(int ix, int iy, int iz, size_t xextent, size_t yextent, size_t zextent, FFTW_REAL *output_values, double *results)
{

  /////////////////////////////////////////////////////////////////////////////
  //
  // This function computes the integral over a cell of the product of the
  // Panphasia field with an 'analysing' Legendre polynomial. As the
  // integrand is a polynomial, Gaussian quadrature can be used for
  // integration as it is exact up to rounding error provide p_order
  // is less than 10.
  //
  /////////////////////////////////////////////////////////////////////////////

  const double GQ_weights[5] = {0.2955242247147529, 0.2692667193099963,
                                0.2190863625159820, 0.1494513491505806,
                                0.0666713443086881};

  const double GQ_abscissa[5] = {0.1488743389816312, 0.4333953941292472,
                                 0.6794095682990244, 0.8650633666889845,
                                 0.9739065285171717};

  double weights[10];
  double abscissa[10];

  for (int i = 0; i < 5; i++)
  {
    weights[i] = GQ_weights[4 - i];
    weights[i + 5] = GQ_weights[i];
    abscissa[i] = -GQ_abscissa[4 - i];
    abscissa[i + 5] = GQ_abscissa[i];
  }

  if (yextent != zextent)
  {
    printf("Non-cubic domain in integrate cell\n");
    abort();
  }

  if (p_order >= 10)
  {
    printf("Higher order Gaussian Quadrature needed!\n");
    abort();
  }

  double a = 0.0;
  double b = 1.0;

  double middle = 0.5 * (b + a);
  double range = 0.5 * (b - a);

  double sum[Nbasis];
  for (size_t i = 0; i < Nbasis; ++i)
    sum[i] = 0.0;

  // double test_sum = 0.0;

  for (int i = 0; i < 10; i++)
  {
    double xp = middle + range * abscissa[i];
    for (int j = 0; j < 10; j++)
    {
      double yp = middle + range * abscissa[j];
      for (int k = 0; k < 10; k++)
      {
        double zp = middle + range * abscissa[k];
        ////////////////////////////////////////////////////////////////////////////////////////

        double panphasia_value;
        double xv = (double)ix + xp;
        double yv = (double)iy + yp;
        double zv = (double)iz + zp;

        if (compute_panphasia_(xv, yv, zv, xextent, yextent, zextent, output_values,
                               &panphasia_value) == 1)
        {
          printf("Call to compute_panphasia_ out of range \n");
          abort();
        }

        double uq, vq, wq;

        uq = 2.0 * (xv / (double)yextent) - 1.0;
        vq = 2.0 * (yv / (double)yextent) - 1.0;
        wq = 2.0 * (zv / (double)yextent) - 1.0;

        int p = p_order;

        double lgp_uq[p_order + 1];
        double lgp_vq[p_order + 1];
        double lgp_wq[p_order + 1];

        gsl_sf_legendre_Pl_array(p, uq, lgp_uq);
        gsl_sf_legendre_Pl_array(p, vq, lgp_vq);
        gsl_sf_legendre_Pl_array(p, wq, lgp_wq);

        for (int ii = 0; ii < p + 1; ii++)
        {
          lgp_uq[ii] *= sqrt(2 * ii + 1);
          lgp_vq[ii] *= sqrt(2 * ii + 1);
          lgp_wq[ii] *= sqrt(2 * ii + 1);
        }

        ////////////////////////////////////////////////////////////////////////////////////////
        double combined_weight = pow(range, 3) * weights[i] * weights[j] * weights[k];

        for (int jj = 0; jj < Nbasis; jj++)
        {
          int lq = irank_p[0][jj];
          int mq = irank_p[1][jj];
          int nq = irank_p[2][jj];
          double af = combined_weight * lgp_uq[lq] * lgp_vq[mq] * lgp_wq[nq]; // Analysing function

          sum[jj] += af * panphasia_value;
          ////////////////////////////////////////////////////////////////////////////////////////
        }
      }
    }
  }

  for (int i = 0; i < Nbasis; i++)
    results[i] = sum[i];
}

int compute_panphasia_(double x, double y, double z, size_t xextent, size_t yextent,
                       size_t zextent, FFTW_REAL *output_values, double *panphasia_value)
{

  if ((x < 0) || (x >= (double)xextent))
    return (1);
  if ((y < 0) || (y >= (double)yextent))
    return (1);
  if ((z < 0) || (z >= (double)zextent))
    return (1);

  int ix = (int)x;
  int iy = (int)y;
  int iz = (int)z;

  double up = 2.0 * (x - ix) - 1.0;
  double vp = 2.0 * (y - iy) - 1.0;
  double wp = 2.0 * (z - iz) - 1.0;

  double lgp_up[p_order + 1];
  double lgp_vp[p_order + 1];
  double lgp_wp[p_order + 1];

  int p = p_order;

  gsl_sf_legendre_Pl_array(p, up, lgp_up);
  gsl_sf_legendre_Pl_array(p, vp, lgp_vp);
  gsl_sf_legendre_Pl_array(p, wp, lgp_wp);

  for (int i = 0; i < p + 1; i++)
  {
    lgp_up[i] *= sqrt(2 * i + 1);
    lgp_vp[i] *= sqrt(2 * i + 1);
    lgp_wp[i] *= sqrt(2 * i + 1);
  };

  *panphasia_value = 0;
  for (int i = 0; i < Nbasis; i++)
  {
    int lp = irank_p[0][i];
    int mp = irank_p[1][i];
    int np = irank_p[2][i];
    *panphasia_value += output_values[Nbasis * (ix * yextent * zextent + iy * zextent + iz) + i] *
                        lgp_up[lp] * lgp_vp[mp] * lgp_wp[np];
  }
  return (0);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// Fortran wrapper for gsl function to evaluate Spherical Bessel functions
void spherical_bessel_(int *l, double *x, double *result)
{

  *result = gsl_sf_bessel_jl(*l, *x);
}

void compute_sph_bessel_coeffs(int nfft, int pmax, int n4dimen, int fdim, double complex *sph_bessel_coeff)
{
  const double pi = 4.0 * atan(1.0);
  for (int l = 0; l <= pmax; l++)
  {
    double norm = sqrt((double)(2 * l + 1));
    double complex phase_shift = cpow(-I, l);
    for (int i = 0; i < nfft; i++)
    {
      int j = (i <= nfft / 2) ? i : i - nfft;
      int k = abs(j);
      double sign = (j < 0) ? pow(-1.0, l) : 1.0;
      double x = pi * (double)fdim * (double)k / (double)nfft;
      double result;
      spherical_bessel_(&l, &x, &result);

      if (i != nfft / 2)
      {
        sph_bessel_coeff[l * n4dimen + i] = norm * phase_shift * sign * result;
      }
      else
      {
        sph_bessel_coeff[l * n4dimen + i] = (l == 0) ? 1.0 : 0.0; //Treat Nyquist mode specially
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//
// End
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
