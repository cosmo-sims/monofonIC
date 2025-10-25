#pragma once

namespace PANPHASIA1
{
  const int maxdim = 60, maxlev = 50, maxpow = 3 * maxdim;
  typedef int rand_offset_[5];
  typedef struct
  {
    int state[133]; // Nstore = Nstate (=5) + Nbatch (=128)
    int need_fill;
    int pos;
  } rand_state_;

  typedef struct
  {
    int base_state[5], base_lev_start[5][maxdim + 1];
    rand_offset_ poweroffset[maxpow + 1], superjump;
    rand_state_ current_state[maxpow + 2];

    int layer_min, layer_max, indep_field;

    long long xorigin_store[2][2][2], yorigin_store[2][2][2], zorigin_store[2][2][2];
    int lev_common, layer_min_store, layer_max_store;
    long long ix_abs_store, iy_abs_store, iz_abs_store, ix_per_store, iy_per_store, iz_per_store, ix_rel_store,
        iy_rel_store, iz_rel_store;
    double exp_coeffs[8][8][maxdim + 2];
    long long xcursor[maxdim + 1], ycursor[maxdim + 1], zcursor[maxdim + 1];
    int ixshift[2][2][2], iyshift[2][2][2], izshift[2][2][2];

    double cell_data[9][8];
    int ixh_last, iyh_last, izh_last;
    int init;

    int init_cell_props;
    int init_lecuyer_state;
    long long p_xcursor[62], p_ycursor[62], p_zcursor[62];

  } pan_state_;

  extern "C"
  {
    void start_panphasia_(pan_state_ *lstate, const char *descriptor, int *ngrid, int *bverbose);

    void parse_descriptor_(const char *descriptor, int16_t *l, int32_t *ix, int32_t *iy, int32_t *iz, int16_t *side1,
                           int16_t *side2, int16_t *side3, int32_t *check_int, char *name);

    void panphasia_cell_properties_(pan_state_ *lstate, int *ixcell, int *iycell, int *izcell, double *cell_prop);

    void adv_panphasia_cell_properties_(pan_state_ *lstate, int *ixcell, int *iycell, int *izcell, int *layer_min,
                                        int *layer_max, int *indep_field, double *cell_prop);

    void set_phases_and_rel_origin_(pan_state_ *lstate, const char *descriptor, int *lev, long long *ix_rel,
                                    long long *iy_rel, long long *iz_rel, int *VERBOSE);
  }

  struct RNG
  {

    config_file *pcf_;

    struct panphasia_descriptor
    {
      int16_t wn_level_base;
      int32_t i_xorigin_base, i_yorigin_base, i_zorigin_base;
      int16_t i_base, i_base_y, i_base_z;
      int32_t check_rand;
      std::string name;

      explicit panphasia_descriptor(std::string dstring)
      {
        char tmp[100];
        std::memset(tmp, ' ', 100);
        parse_descriptor_(dstring.c_str(), &wn_level_base, &i_xorigin_base, &i_yorigin_base, &i_zorigin_base, &i_base,
                          &i_base_y, &i_base_z, &check_rand, tmp);
        for (int i = 0; i < 100; i++)
          if (tmp[i] == ' ')
          {
            tmp[i] = '\0';
            break;
          }
        name = tmp;
        name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
      }
    };

    // greatest common divisor
    int gcd(int a, int b)
    {
      if (b == 0)
        return a;
      return gcd(b, a % b);
    }

    // least common multiple
    int lcm(int a, int b) { return abs(a * b) / gcd(a, b); }

    // Two or largest power of 2 less than the argument
    int largest_power_two_lte(int b)
    {
      int a = 1;
      if (b <= a)
        return a;
      while (2 * a < b)
        a = 2 * a;
      return a;
    }

    std::string descriptor_string_;
    int num_threads_;
    int levelmin_, levelmin_final_, levelmax_, ngrid_, ngrid_panphasia_;
    bool incongruent_fields_;
    double inter_grid_phase_adjustment_;
    // double translation_phase_;
    pan_state_ *lstate;
    int grid_p_;
    int coordinate_system_shift_[3];
    int ix_abs_[3], ix_per_[3], ix_rel_[3], level_p_, lextra_;

    void clear_panphasia_thread_states(void)
    {
      for (int i = 0; i < num_threads_; ++i)
      {
        lstate[i].init = 0;
        lstate[i].init_cell_props = 0;
        lstate[i].init_lecuyer_state = 0;
      }
    }

    void initialize_for_grid_structure(void)
    {
      // if ngrid is not a multiple of i_base, then we need to enlarge and then sample down
      ngrid_ = pcf_->get_value<size_t>("setup", "GridRes");
      int ngridminsize_panphasia = pcf_->get_value_safe<int>("random", "PanphasiaMinRootResolution", 512);

      grid_p_ = pdescriptor_->i_base;

      lextra_ = (log10((double)ngrid_ / (double)grid_p_) + 0.001) / log10(2.0);
      // lmin

      ngrid_panphasia_ = (1 << lextra_) * grid_p_;

      while (ngrid_panphasia_ < ngridminsize_panphasia)
      {
        lextra_++;
        ngrid_panphasia_ *= 2;
      }
      assert(ngrid_panphasia_ >= ngridminsize_panphasia);

      clear_panphasia_thread_states();

      music::ilog.Print("PANPHASIA: using grid size %d (level=%d)", ngrid_panphasia_, lextra_);
      if (ngridminsize_panphasia < 512)
        music::ilog.Print("PANPHASIA WARNING: PanphasiaMinRootResolution = %d below minimum recommended of 512", ngridminsize_panphasia);
      music::ilog.Print("PANPHASIA: running with %d threads", num_threads_);

      coordinate_system_shift_[0] = -pcf_->get_value_safe<int>("setup", "shift_x", 0);
      coordinate_system_shift_[1] = -pcf_->get_value_safe<int>("setup", "shift_y", 0);
      coordinate_system_shift_[2] = -pcf_->get_value_safe<int>("setup", "shift_z", 0);
    }

    std::unique_ptr<panphasia_descriptor> pdescriptor_;

    RNG(config_file *pcf)
        : pcf_(pcf)
    {
      descriptor_string_ = pcf_->get_value<std::string>("random", "descriptor");

#ifdef _OPENMP
      num_threads_ = omp_get_max_threads();
#else
      num_threads_ = 1;
#endif

      // create independent state descriptions for each thread
      lstate = new pan_state_[num_threads_];

      // parse the descriptor for its properties
      pdescriptor_ = std::make_unique<panphasia_descriptor>(descriptor_string_);

      music::ilog.Print("PANPHASIA: descriptor \'%s\' is base %d,", pdescriptor_->name.c_str(), pdescriptor_->i_base);

      // write panphasia base size into config file for the grid construction
      // as the gridding unit we use the least common multiple of 2 and i_base
      std::stringstream ss;
      //ARJ  ss << lcm(2, pdescriptor_->i_base);
      //ss <<  two_or_largest_power_two_less_than(pdescriptor_->i_base);//ARJ
      ss << 2; //ARJ - set gridding unit to two
      pcf_->insert_value("setup", "gridding_unit", ss.str());
      ss.str(std::string());
      ss << pdescriptor_->i_base;
      pcf_->insert_value("random", "base_unit", ss.str());

      this->initialize_for_grid_structure();
    }

    ~RNG() { delete[] lstate; }

    void Fill(Grid_FFT<real_t> &g)
    {
      auto sinc = [](real_t x)
      { return (std::fabs(x) > 1e-16) ? std::sin(x) / x : 1.0; };
      auto dsinc = [](real_t x)
      { return (std::fabs(x) > 1e-16) ? (x * std::cos(x) - std::sin(x)) / (x * x) : 0.0; };
      const real_t sqrt3{std::sqrt(3.0)}, sqrt27{std::sqrt(27.0)};

      // we will overwrite 'g', we can deallocate it while we prepare the panphasia field
      g.reset();

      clear_panphasia_thread_states();

      // temporaries
      Grid_FFT<real_t> g0({size_t(ngrid_panphasia_), size_t(ngrid_panphasia_), size_t(ngrid_panphasia_)}, g.length_);
      Grid_FFT<real_t> g1({size_t(ngrid_panphasia_), size_t(ngrid_panphasia_), size_t(ngrid_panphasia_)}, g.length_);
      Grid_FFT<real_t> g2({size_t(ngrid_panphasia_), size_t(ngrid_panphasia_), size_t(ngrid_panphasia_)}, g.length_);
      Grid_FFT<real_t> g3({size_t(ngrid_panphasia_), size_t(ngrid_panphasia_), size_t(ngrid_panphasia_)}, g.length_);
      Grid_FFT<real_t> g4({size_t(ngrid_panphasia_), size_t(ngrid_panphasia_), size_t(ngrid_panphasia_)}, g.length_);

      double t1 = get_wtime();
      // double tp = t1;

#pragma omp parallel
      {
#ifdef _OPENMP
        const int mythread = omp_get_thread_num();
#else
        const int mythread = 0;
#endif

        //int odd_x, odd_y, odd_z;
        //int ng_level = ngrid_ * (1 << (level - levelmin_)); // full resolution of current level

        int verbosity = (mythread == 0);
        char descriptor[100];
        std::memset(descriptor, 0, 100);
        std::memcpy(descriptor, descriptor_string_.c_str(), descriptor_string_.size());

        start_panphasia_(&lstate[mythread], descriptor, &ngrid_panphasia_, &verbosity);

        {
          panphasia_descriptor d(descriptor_string_);

          int level_p = d.wn_level_base + lextra_;

          lstate[mythread].layer_min = 0;
          lstate[mythread].layer_max = level_p;
          lstate[mythread].indep_field = 1;

          long long ix_rel[3];
          ix_rel[0] = 0; //ileft_corner_p[0];
          ix_rel[1] = 0; //ileft_corner_p[1];
          ix_rel[2] = 0; //ileft_corner_p[2];

          set_phases_and_rel_origin_(&lstate[mythread], descriptor, &level_p, &ix_rel[0], &ix_rel[1], &ix_rel[2],
                                     &verbosity);
        }

        if (verbosity)
          t1 = get_wtime();

        std::array<double, 9> cell_prop;
        pan_state_ *ps = &lstate[mythread];

#pragma omp for //nowait
        for (size_t i = 0; i < g0.size(0); i += 2)
        {
          const int ixmax(std::min<int>(2, g0.size(0) - i));
          for (size_t j = 0; j < g0.size(1); j += 2)
          {
            const int iymax(std::min<int>(2, g0.size(1) - j));
            for (size_t k = 0; k < g0.size(2); k += 2)
            {
              const int izmax(std::min<int>(2, g0.size(2) - k));

              // ARJ - added inner set of loops to speed up evaluation of Panphasia
              for (int ix = 0; ix < ixmax; ++ix)
              {
                for (int iy = 0; iy < iymax; ++iy)
                {
                  for (int iz = 0; iz < izmax; ++iz)
                  {
                    int ilocal = i + ix;
                    int jlocal = j + iy;
                    int klocal = k + iz;

                    int iglobal = ilocal + g0.local_0_start_;
                    int jglobal = jlocal;
                    int kglobal = klocal;

                    adv_panphasia_cell_properties_(ps, &iglobal, &jglobal, &kglobal, &ps->layer_min,
                                                   &ps->layer_max, &ps->indep_field, &cell_prop[0]);

                    g0.relem(ilocal, jlocal, klocal) = cell_prop[0];
                    g1.relem(ilocal, jlocal, klocal) = cell_prop[4];
                    g2.relem(ilocal, jlocal, klocal) = cell_prop[2];
                    g3.relem(ilocal, jlocal, klocal) = cell_prop[1];
                    g4.relem(ilocal, jlocal, klocal) = cell_prop[8];
                  }
                }
              }
            }
          }
        }
      } // end omp parallel region

      g0.FourierTransformForward();
      g1.FourierTransformForward();
      g2.FourierTransformForward();
      g3.FourierTransformForward();
      g4.FourierTransformForward();

#pragma omp parallel for
      for (size_t i = 0; i < g0.size(0); i++)
      {
        for (size_t j = 0; j < g0.size(1); j++)
        {
          for (size_t k = 0; k < g0.size(2); k++)
          {
            if (!g0.is_nyquist_mode(i, j, k))
            {
              auto kvec = g0.get_k<real_t>(i, j, k);

              auto argx = 0.5 * M_PI * kvec[0] / g0.kny_[0];
              auto argy = 0.5 * M_PI * kvec[1] / g0.kny_[1];
              auto argz = 0.5 * M_PI * kvec[2] / g0.kny_[2];

              auto fx = real_t(sinc(argx));
              auto gx = ccomplex_t(0.0, dsinc(argx));
              auto fy = real_t(sinc(argy));
              auto gy = ccomplex_t(0.0, dsinc(argy));
              auto fz = real_t(sinc(argz));
              auto gz = ccomplex_t(0.0, dsinc(argz));

              auto temp = (fx + sqrt3 * gx) * (fy + sqrt3 * gy) * (fz + sqrt3 * gz);
              auto magnitude = real_t(std::sqrt(1.0 - std::fabs(temp * temp)));

              auto y0(g0.kelem(i, j, k)), y1(g1.kelem(i, j, k)), y2(g2.kelem(i, j, k)), y3(g3.kelem(i, j, k)), y4(g4.kelem(i, j, k));

              g0.kelem(i, j, k) = y0 * fx * fy * fz + sqrt3 * (y1 * gx * fy * fz + y2 * fx * gy * fz + y3 * fx * fy * gz) + y4 * magnitude;
            }
            else
            {
              g0.kelem(i, j, k) = 0.0;
            }
          }
        }
      }

      // music::ilog.Print("\033[31mtiming [build panphasia field]: %f s\033[0m", get_wtime() - tp);
      // tp = get_wtime();

      g1.FourierTransformBackward(false);
      g2.FourierTransformBackward(false);
      g3.FourierTransformBackward(false);
      g4.FourierTransformBackward(false);

#pragma omp parallel
      {
#ifdef _OPENMP
        const int mythread = omp_get_thread_num();
#else
        const int mythread = 0;
#endif

        // int odd_x, odd_y, odd_z;
        int verbosity = (mythread == 0);
        char descriptor[100];
        std::memset(descriptor, 0, 100);
        std::memcpy(descriptor, descriptor_string_.c_str(), descriptor_string_.size());

        start_panphasia_(&lstate[mythread], descriptor, &ngrid_panphasia_, &verbosity);

        {
          panphasia_descriptor d(descriptor_string_);

          int level_p = d.wn_level_base + lextra_;

          lstate[mythread].layer_min = 0;
          lstate[mythread].layer_max = level_p;
          lstate[mythread].indep_field = 1;

          long long ix_rel[3];
          ix_rel[0] = 0; //ileft_corner_p[0];
          ix_rel[1] = 0; //ileft_corner_p[1];
          ix_rel[2] = 0; //ileft_corner_p[2];

          set_phases_and_rel_origin_(&lstate[mythread], descriptor, &level_p, &ix_rel[0], &ix_rel[1], &ix_rel[2],
                                     &verbosity);
        }

        if (verbosity)
          t1 = get_wtime();

        //***************************************************************
        // Process Panphasia values: p110, p011, p101, p111
        //****************************************************************
        std::array<double, 9> cell_prop;
        pan_state_ *ps = &lstate[mythread];

#pragma omp for //nowait
        for (size_t i = 0; i < g1.size(0); i += 2)
        {
          const int ixmax(std::min<int>(2, g1.size(0) - i));
          for (size_t j = 0; j < g1.size(1); j += 2)
          {
            const int iymax(std::min<int>(2, g1.size(1) - j));
            for (size_t k = 0; k < g1.size(2); k += 2)
            {
              const int izmax(std::min<int>(2, g1.size(2) - k));

              // ARJ - added inner set of loops to speed up evaluation of Panphasia
              for (int ix = 0; ix < ixmax; ++ix)
              {
                for (int iy = 0; iy < iymax; ++iy)
                {
                  for (int iz = 0; iz < izmax; ++iz)
                  {
                    int ilocal = i + ix;
                    int jlocal = j + iy;
                    int klocal = k + iz;

                    int iglobal = ilocal + g1.local_0_start_;
                    int jglobal = jlocal;
                    int kglobal = klocal;

                    adv_panphasia_cell_properties_(ps, &iglobal, &jglobal, &kglobal, &ps->layer_min,
                                                   &ps->layer_max, &ps->indep_field, &cell_prop[0]);

                    g1.relem(ilocal, jlocal, klocal) = cell_prop[6];
                    g2.relem(ilocal, jlocal, klocal) = cell_prop[3];
                    g3.relem(ilocal, jlocal, klocal) = cell_prop[5];
                    g4.relem(ilocal, jlocal, klocal) = cell_prop[7];
                  }
                }
              }
            }
          }
        }
      } // end omp parallel region

      // music::ilog.Print("\033[31mtiming [adv_panphasia_cell_properties2]: %f s \033[0m", get_wtime() - tp);
      // tp = get_wtime();

      /////////////////////////////////////////////////////////////////////////
      // transform and convolve with Legendres
      g1.FourierTransformForward();
      g2.FourierTransformForward();
      g3.FourierTransformForward();
      g4.FourierTransformForward();

#pragma omp parallel for
      for (size_t i = 0; i < g0.size(0); i++)
      {
        for (size_t j = 0; j < g0.size(1); j++)
        {
          for (size_t k = 0; k < g0.size(2); k++)
          {
            if (!g0.is_nyquist_mode(i, j, k))
            {
              auto kvec = g0.get_k<real_t>(i, j, k);

              auto argx = 0.5 * M_PI * kvec[0] / g0.kny_[0];
              auto argy = 0.5 * M_PI * kvec[1] / g0.kny_[1];
              auto argz = 0.5 * M_PI * kvec[2] / g0.kny_[2];

              auto fx = real_t(sinc(argx));
              auto gx = ccomplex_t(0.0, dsinc(argx));
              auto fy = real_t(sinc(argy));
              auto gy = ccomplex_t(0.0, dsinc(argy));
              auto fz = real_t(sinc(argz));
              auto gz = ccomplex_t(0.0, dsinc(argz));

              auto y1(g1.kelem(i, j, k)), y2(g2.kelem(i, j, k)), y3(g3.kelem(i, j, k)), y4(g4.kelem(i, j, k));

              g0.kelem(i, j, k) += real_t(3.0) * (y1 * gx * gy * fz + y2 * fx * gy * gz + y3 * gx * fy * gz) + sqrt27 * y4 * gx * gy * gz;

              // do final phase shift to account for corner centered coordinates vs. cell centers
              auto phase_shift = -0.5 * M_PI * (kvec[0] / g0.kny_[0] + kvec[1] / g0.kny_[1] + kvec[2] / g0.kny_[2]);

              g0.kelem(i, j, k) *= std::exp(ccomplex_t(0, phase_shift));
            }
          }
        }
      }

      g1.reset();
      g2.reset();
      g3.reset();
      g4.reset();

      g.allocate();
      g0.FourierInterpolateCopyTo(g);

      music::ilog.Print("time for calculating PANPHASIA field : %f s, %f µs/cell", get_wtime() - t1,
                        1e6 * (get_wtime() - t1) / g.global_size(0) / g.global_size(1) / g.global_size(2));
      music::ilog.Print("PANPHASIA k-space statistices: mean Re = %f, std = %f", g.mean(), g.std());
    }
  };
};