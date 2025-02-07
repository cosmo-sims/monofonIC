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
#include <random_plugin.hh>
#include "random_music_wnoise_generator.hh"

typedef music_wnoise_generator<real_t> rng;

class RNG_music : public RNG_plugin
{
protected:
  std::vector<long> rngseeds_;
  std::vector<std::string> rngfnames_;
  std::vector<rng *> randc_;
  unsigned ran_cube_size_;

  int levelmin_, levelmax_, levelmin_seed_;

  bool disk_cached_;
  bool restart_;
  bool initialized_;

  std::vector<std::vector<real_t> *> mem_cache_;

  //! checks if the specified string is numeric
  bool is_number(const std::string &s);

  //! parses the random number parameters in the conf file
  void parse_random_parameters(void);

  //! computes the white noise fields and keeps them either in memory or on disk
  void compute_random_numbers(void);

  //! adjusts averages
  //void correct_avg(int icoarse, int ifine);

  //! store the white noise fields in memory or on disk
  //void store_rnd(int ilevel, rng *prng);

  bool is_power_of_two(size_t x) const
  {
      return (x & (x - 1)) == 0;
  }

  unsigned int_log2( size_t v ) const
  {
    unsigned r{0}; // r will be lg(v)
    while (v >>= 1) r++;
    return r;
  }

public:
  explicit RNG_music(config_file &cf) : RNG_plugin(cf), initialized_(false) 
  {
    // we need to make sure that the chosen resolution is a power of 2 resolution
    size_t res = pcf_->get_value<size_t>("setup", "GridRes");
    if( !is_power_of_two(res) ){
      std::string errmsg("MUSIC random number plugin requires [setup]/GridRes to be a power of 2!");
      music::flog << errmsg << std::endl;
      throw std::runtime_error(errmsg.c_str());
    } 
    levelmin_ = int_log2(res);
    levelmax_ = levelmin_;
    music::ilog << "MUSIC1 RNG plugin: setting levelmin = levelmax = " << levelmin_ << std::endl;

    this->initialize_for_grid_structure( );
    
  }

  ~RNG_music() 
  {
    for( auto rr : randc_ ){
      if( rr != nullptr ) delete rr;
    }
  }

  bool isMultiscale() const { return true; }

  void Fill_Grid( Grid_FFT<real_t>& g ) 
  {
    // determine extent of grid to be filled (can be a slab with MPI)
    const size_t i0 = g.local_0_start_, j0{0}, k0{0};
    const size_t Ni = g.rsize(0), Nj = g.rsize(1), Nk = g.rsize(2);

    // make sure we're in real space
    g.FourierTransformBackward();

    // copy over
    #pragma omp parallel for
    for( auto i = i0; i<Ni; ++i )
    {
      size_t ip  = i-i0; // index in g
      for( auto j = j0; j<Nj; ++j )
      {
        auto   jp = j-j0; // index in g
        for( auto k = k0; k<Nk; ++k )
        {
          auto   kp = k-k0; // index in g
          g.relem(ip,jp,kp) = (*randc_[levelmin_])(i,j,k);
        } 
      }  
    }
  } 

  void initialize_for_grid_structure()
  {
    ran_cube_size_ = pcf_->get_value_safe<unsigned>("random", "cubesize", DEF_RAN_CUBE_SIZE);
    disk_cached_ = pcf_->get_value_safe<bool>("random", "disk_cached", true);
    restart_ = pcf_->get_value_safe<bool>("random", "restart", false);

    mem_cache_.assign(levelmax_ - levelmin_ + 1, (std::vector<real_t> *)NULL);

    if (restart_ && !disk_cached_)
    {
      music::elog.Print("Cannot restart from mem cached random numbers.");
      throw std::runtime_error("Cannot restart from mem cached random numbers.");
    }

    //... determine seed/white noise file data to be applied
    parse_random_parameters();

    if (!restart_)
    {
      //... compute the actual random numbers
      compute_random_numbers();
    }

    initialized_ = true;
  }

  //void fill_grid(int level, DensityGrid<real_t> &R);
};

bool RNG_music::is_number(const std::string &s)
{
  for (size_t i = 0; i < s.length(); i++)
    if (!std::isdigit(s[i]) && s[i] != '-')
      return false;

  return true;
}

void RNG_music::parse_random_parameters(void)
{
  //... parse random number options
  for (int i = 0; i <= 100; ++i)
  {
    char seedstr[128];
    std::string tempstr;
    bool noseed = false;
    sprintf(seedstr, "seed[%d]", i);
    if (pcf_->contains_key("random", seedstr))
      tempstr = pcf_->get_value<std::string>("random", seedstr);
    else
    {
      // "-2" means that no seed entry was found for that level
      tempstr = std::string("-2");
      noseed = true;
    }

    if (is_number(tempstr))
    {
      long ltemp;
      pcf_->convert(tempstr, ltemp);
      rngfnames_.push_back("");
      if (noseed) // ltemp < 0 )
        //... generate some dummy seed which only depends on the level, negative so we know it's not
        //... an actual seed and thus should not be used as a constraint for coarse levels
        // rngseeds_.push_back( -abs((unsigned)(ltemp-i)%123+(unsigned)(ltemp+827342523521*i)%123456789) );
        rngseeds_.push_back(-std::abs((long)(ltemp - i) % 123 + (long)(ltemp + 7342523521 * i) % 123456789));
      else
      {
        if (ltemp <= 0)
        {
          music::elog.Print("Specified seed [random]/%s needs to be a number >0!", seedstr);
          throw std::runtime_error("Seed values need to be >0");
        }
        rngseeds_.push_back(ltemp);
      }
    }
    else
    {
      rngfnames_.push_back(tempstr);
      rngseeds_.push_back(-1);
      music::ilog.Print("Random numbers for level %3d will be read from file.", i);
    }
  }

  //.. determine for which levels random seeds/random number files are given
  levelmin_seed_ = -1;
  for (unsigned ilevel = 0; ilevel < rngseeds_.size(); ++ilevel)
  {
    if (levelmin_seed_ < 0 && (rngfnames_[ilevel].size() > 0 || rngseeds_[ilevel] >= 0))
      levelmin_seed_ = ilevel;
  }

  if( levelmin_seed_ < 0 ){
    music::elog.Print("No seed specified for MUSIC1 RNG plugin!");
    throw std::runtime_error("No seed specified for MUSIC1 RNG plugin!");
  }
}

void RNG_music::compute_random_numbers(void)
{
  bool rndsign = pcf_->get_value_safe<bool>("random", "grafic_sign", false);
  bool music2_rng = pcf_->get_value_safe<bool>("random", "music2_rng", false);
  if (music2_rng)
  {
    music::ilog.Print("MUSIC1 RNG plugin: using MUSIC2 compatible seeds");
  }

  //--- FILL ALL WHITE NOISE ARRAYS FOR WHICH WE NEED THE FULL FIELD ---//

  randc_.assign(std::max(levelmax_, levelmin_seed_) + 1, nullptr);


  //... seeds are given for a level coarser than levelmin
  if (levelmin_seed_ < levelmin_)
  {
    if (rngfnames_[levelmin_seed_].size() > 0)
      randc_[levelmin_seed_] = new rng(1 << levelmin_seed_, rngfnames_[levelmin_seed_], rndsign);
    else
      randc_[levelmin_seed_] = new rng(1 << levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true);

    for (int i = levelmin_seed_ + 1; i <= levelmin_; ++i)
    {
      //#warning add possibility to read noise from file also here!

      if (rngfnames_[i].size() > 0)
        music::wlog.Print("Warning: Cannot use filenames for higher levels currently! Ignoring!");

      randc_[i] = new rng(*randc_[i - 1], ran_cube_size_, rngseeds_[i], music2_rng, true);
      delete randc_[i - 1];
      randc_[i - 1] = NULL;
    }
  }

  //... seeds are given for a level finer than levelmin, obtain by averaging
  if (levelmin_seed_ > levelmin_)
  {
    if (rngfnames_[levelmin_seed_].size() > 0)
      randc_[levelmin_seed_] = new rng(1 << levelmin_seed_, rngfnames_[levelmin_seed_], rndsign);
    else
      randc_[levelmin_seed_] =
          new rng(1 << levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true); //, x0, lx );

    for (int ilevel = levelmin_seed_ - 1; ilevel >= (int)levelmin_; --ilevel)
    {
      if (rngseeds_[ilevel - levelmin_] > 0)
        music::ilog.Print("Warning: random seed for level %d will be ignored.\n"
                "            consistency requires that it is obtained by restriction from level %d",
                ilevel, levelmin_seed_);

      // if( brealspace_tf && ilevel < levelmax_ )
      //  randc_[ilevel] = new rng( *randc_[ilevel+1], false );
      // else // do k-space averaging
      randc_[ilevel] = new rng(*randc_[ilevel + 1], true);

      if (ilevel + 1 > levelmax_)
      {
        delete randc_[ilevel + 1];
        randc_[ilevel + 1] = NULL;
      }
    }
  }

  //--- GENERATE AND STORE ALL LEVELS, INCLUDING REFINEMENTS ---//

  //... levelmin
  if (randc_[levelmin_] == NULL)
  {
    if (rngfnames_[levelmin_].size() > 0)
      randc_[levelmin_] = new rng(1 << levelmin_, rngfnames_[levelmin_], rndsign);
    else
      randc_[levelmin_] = new rng(1 << levelmin_, ran_cube_size_, rngseeds_[levelmin_], true);
  }

// if( levelmax_ == levelmin_ )
#if 0
  {
    //... apply constraints to coarse grid
    //... if no constraints are specified, or not for this level
    //... constraints.apply will return without doing anything
    int x0[3] = { 0, 0, 0 };
    int lx[3] = { 1<<levelmin_, 1<<levelmin_, 1<<levelmin_ };
    constraints.apply( levelmin_, x0, lx, randc_[levelmin_] );
  }
#endif

  // store_rnd(levelmin_, randc_[levelmin_]);

  //... refinement levels
  // for (int ilevel = levelmin_ + 1; ilevel <= levelmax_; ++ilevel)
  // {
  //   int lx[3], x0[3];
  //   int shift[3], levelmin_poisson;
  //   shift[0] = pcf_->get_value<int>("setup", "shift_x");
  //   shift[1] = pcf_->get_value<int>("setup", "shift_y");
  //   shift[2] = pcf_->get_value<int>("setup", "shift_z");

  //   levelmin_poisson = pcf_->get_value<unsigned>("setup", "levelmin");

  //   int lfac = 1 << (ilevel - levelmin_poisson);

  //   lx[0] = 2 * prefh_->size(ilevel, 0);
  //   lx[1] = 2 * prefh_->size(ilevel, 1);
  //   lx[2] = 2 * prefh_->size(ilevel, 2);
  //   x0[0] = prefh_->offset_abs(ilevel, 0) - lfac * shift[0] - lx[0] / 4;
  //   x0[1] = prefh_->offset_abs(ilevel, 1) - lfac * shift[1] - lx[1] / 4;
  //   x0[2] = prefh_->offset_abs(ilevel, 2) - lfac * shift[2] - lx[2] / 4;

  //   if (randc_[ilevel] == NULL)
  //     randc_[ilevel] =
  //         new rng(*randc_[ilevel - 1], ran_cube_size_, rngseeds_[ilevel], kavg, ilevel == levelmin_ + 1, x0, lx);
  //   delete randc_[ilevel - 1];
  //   randc_[ilevel - 1] = NULL;

  //   //... apply constraints to this level, if any
  //   // if( ilevel == levelmax_ )
  //   // constraints.apply( ilevel, x0, lx, randc_[ilevel] );

  //   //... store numbers
  //   store_rnd(ilevel, randc_[ilevel]);
  // }

  // delete randc_[levelmax_];
  // randc_[levelmax_] = NULL;

  //... make sure that the coarse grid contains oct averages where it overlaps with a fine grid
  //... this also ensures that constraints enforced on fine grids are carried to the coarser grids
  // if (brealspace_tf)
  // {
  //   for (int ilevel = levelmax_; ilevel > levelmin_; --ilevel)
  //     correct_avg(ilevel - 1, ilevel);
  // }

  //.. we do not have random numbers for a coarse level, generate them
  /*if( levelmax_rand_ >= (int)levelmin_ )
    {
    std::cerr << "lmaxread >= (int)levelmin\n";
    randc_[levelmax_rand_] = new rng( (unsigned)pow(2,levelmax_rand_), rngfnames_[levelmax_rand_] );
    for( int ilevel = levelmax_rand_-1; ilevel >= (int)levelmin_; --ilevel )
    randc_[ilevel] = new rng( *randc_[ilevel+1] );
         }*/
}



namespace
{
RNG_plugin_creator_concrete<RNG_music> creator("MUSIC1");
}
