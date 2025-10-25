
// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn & Michael Buehlmann (this file)
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
#ifdef USE_HDF5
#include <unistd.h> // for unlink
#include <array>
#include <output_plugin.hh>
#include "HDF_IO.hh"

template <typename T>
std::vector<T> from_7array(const std::array<T,7>& a)
{
  return std::vector<T>{{a[0], a[1], a[2], a[3], a[4], a[5], a[6]}};
}

template <typename T>
std::vector<T> from_value(const T a)
{
  return std::vector<T>{{a}};
}

template <typename write_real_t>
class swift_output_plugin : public output_plugin
{

protected:
  int num_ranks_, this_rank_;
  int num_files_;

  real_t lunit_, vunit_, munit_;
  real_t boxsize_, hubble_param_, astart_, zstart_;
  bool blongids_, bdobaryons_;

  std::array<uint64_t,7> npart_;
  std::array<uint32_t,7> npartTotal_, npartTotalHighWord_;
  std::array<double,7> mass_;
  double time_;
  double ceint_, h_;

public:
  //! constructor
  explicit swift_output_plugin(config_file &cf, std::unique_ptr<cosmology::calculator> &pcc)
      : output_plugin(cf, pcc, "SWIFT")
  {

    // SWIFT uses a single file as IC */
    num_files_ = 1;
    this_rank_ = 0;
    num_ranks_ = 1;

    real_t astart = 1.0 / (1.0 + cf_.get_value<double>("setup", "zstart"));
    const double rhoc = 27.7536609198; // in h^2 1e10 M_sol / Mpc^3 assuming SWIFT's internal constants

    hubble_param_ = pcc->cosmo_param_["h"];
    zstart_ = cf_.get_value<double>("setup","zstart");
    astart_ = 1.0/(1.0 + zstart_);
    boxsize_ = cf_.get_value<double>("setup", "BoxLength");

    lunit_ = boxsize_ / hubble_param_; // final units will be in Mpc (without h)
    vunit_ = boxsize_; // final units will be in km/s
    munit_ = rhoc * std::pow(boxsize_, 3) / hubble_param_; // final units will be in 1e10 M_sol

    blongids_ = cf_.get_value_safe<bool>("output", "UseLongids", true);
    bdobaryons_ = cf_.get_value<bool>("setup","DoBaryons");

    for (int i = 0; i < 7; ++i)
    {
      npart_[i] = 0;
      npartTotal_[i] = 0;
      npartTotalHighWord_[i] = 0;
      mass_[i] = 0.0;
    }

    time_ = astart;

#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks_);
#endif

    if (bdobaryons_) {

      const double gamma  = cf_.get_value_safe<double>("cosmology", "gamma", 5.0 / 3.0);
      const double YHe    = pcc_->cosmo_param_["YHe"];
      const double omegab = pcc_->cosmo_param_["Omega_b"];
      const double Tcmb0  = pcc_->cosmo_param_["Tcmb"];

      // compute gas internal energy
      const double npol = (fabs(1.0 - gamma) > 1e-7) ? 1.0 / (gamma - 1.) : 1.0;
      const double unitv = 1e5;
      const double adec = 1.0 / (160. * std::pow(omegab * hubble_param_ * hubble_param_ / 0.022, 2.0 / 5.0));
      const double Tini = astart_ < adec ? Tcmb0 / astart_ : Tcmb0 / astart_ / astart_ * adec;
      const double mu = (Tini > 1.e4) ? 4.0 / (8. - 5. * YHe) : 4.0 / (1. + 3. * (1. - YHe));
      ceint_ = 1.3806e-16 / 1.6726e-24 * Tini * npol / mu / unitv / unitv;

      music::ilog.Print("Swift : Calculated initial gas temperature: %.2f K/mu", Tini / mu);
      music::ilog.Print("Swift : set initial internal energy to %.2e km^2/s^2", ceint_);

      h_ = boxsize_ / hubble_param_ / cf_.get_value<double>("setup","GridRes");
      music::ilog.Print("Swift : set initial smoothing length to mean inter-part separation: %.2f Mpc", h_);
    }

    // Only ranks 0 writes the header
    if (this_rank_ != 0) return;

    // delete output file if it exists
    unlink(fname_.c_str());

    // create output HDF5 file
    HDFCreateFile(fname_);

    // Write UNITS header using the physical constants assumed internally by SWIFT
    HDFCreateGroup(fname_, "Units");
    HDFWriteGroupAttribute(fname_, "Units", "Unit mass in cgs (U_M)", 1.98841e43);         // 10^10 Msun in grams
    HDFWriteGroupAttribute(fname_, "Units", "Unit length in cgs (U_L)", 3.08567758149e24); // 1 Mpc in cm
    HDFWriteGroupAttribute(fname_, "Units", "Unit time in cgs (U_t)", 3.08567758149e19);   // so that unit vel is 1 km/s
    HDFWriteGroupAttribute(fname_, "Units", "Unit current in cgs (U_I)", 1.0);             // 1 Ampere
    HDFWriteGroupAttribute(fname_, "Units", "Unit temperature in cgs (U_T)", 1.0);         // 1 Kelvin

    // Write MUSIC configuration header
    int order = cf_.get_value<int>("setup", "LPTorder");
    std::string load = cf_.get_value<std::string>("setup", "ParticleLoad");
    std::string tf = cf_.get_value<std::string>("cosmology", "transfer");
    std::string cosmo_set = cf_.get_value<std::string>("cosmology", "ParameterSet");
    std::string rng = cf_.get_value<std::string>("random", "generator");
    int do_fixing = cf_.get_value<bool>("setup", "DoFixing");
    int do_invert = cf_.get_value<bool>("setup", "DoInversion");
    int do_baryons = cf_.get_value<bool>("setup", "DoBaryons");
    int do_baryonsVrel = cf_.get_value<bool>("setup", "DoBaryonVrel");
    int L = cf_.get_value<int>("setup", "GridRes");

    HDFCreateGroup(fname_, "ICs_parameters");
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Code", std::string("MUSIC2 - monofonIC"));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Git Revision", std::string(GIT_REV));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Git Tag", std::string(GIT_TAG));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Git Branch", std::string(GIT_BRANCH));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Precision", std::string(CMAKE_PRECISION_STR));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Convolutions", std::string(CMAKE_CONVOLVER_STR));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "PLT", std::string(CMAKE_PLT_STR));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "LPT Order", order);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Particle Load", load);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Transfer Function", tf);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Cosmology Parameter Set", cosmo_set);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Random Generator", rng);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Mode Fixing", do_fixing);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Mode inversion", do_invert);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Baryons", do_baryons);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Baryons Relative Velocity", do_baryonsVrel);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Grid Resolution", L);

    if (tf == "CLASS") {
      double ztarget = cf_.get_value<double>("cosmology", "ztarget");
      HDFWriteGroupAttribute(fname_, "ICs_parameters", "Target Redshift", ztarget);
    }
    if (rng == "PANPHASIA") {
      std::string desc = cf_.get_value<std::string>("random", "descriptor");
      HDFWriteGroupAttribute(fname_, "ICs_parameters", "Descriptor", desc);
    }
  }

  // use destructor to write header post factum
  ~swift_output_plugin()
  {
    if (!std::uncaught_exceptions())
    {
      if (this_rank_  == 0) {

	// Write Standard Gadget / SWIFT hdf5 header
	HDFCreateGroup(fname_, "Header");
	HDFWriteGroupAttribute(fname_, "Header", "Dimension", 3);
	HDFWriteGroupAttribute(fname_, "Header", "BoxSize", boxsize_ / hubble_param_);  // in Mpc, not Mpc/h

	HDFWriteGroupAttribute(fname_, "Header", "NumPart_Total", from_7array<unsigned>(npartTotal_));
	HDFWriteGroupAttribute(fname_, "Header", "NumPart_Total_HighWord", from_7array<unsigned>(npartTotalHighWord_));
	HDFWriteGroupAttribute(fname_, "Header", "NumPart_ThisFile", from_7array<uint64_t>(npart_));
	HDFWriteGroupAttribute(fname_, "Header", "MassTable", from_7array<double>(mass_));

	HDFWriteGroupAttribute(fname_, "Header", "Time", from_value<double>(time_));
	HDFWriteGroupAttribute(fname_, "Header", "Redshift", from_value<double>(zstart_));
	HDFWriteGroupAttribute(fname_, "Header", "Flag_Entropy_ICs", from_value<int>(0));

	HDFWriteGroupAttribute(fname_, "Header", "NumFilesPerSnapshot", from_value<int>(num_files_));

	music::ilog << "Done writing SWIFT IC file to " << fname_ << std::endl;
	music::ilog << "Note that the IC file does not contain any h-factors nor any extra sqrt(a)-factors for the velocities" << std::endl;
	music::ilog << "The SWIFT parameters 'InitialConditions:cleanup_h_factors' and 'InitialConditions:cleanup_velocity_factors' should hence *NOT* be set!" << std::endl;
      }
    }
  }

  output_type write_species_as(const cosmo_species &) const { return output_type::particles; }

  real_t position_unit() const { return lunit_; }

  real_t velocity_unit() const { return vunit_; }

  real_t mass_unit() const { return munit_; }

  bool has_64bit_reals() const
  {
    if (typeid(write_real_t) == typeid(double))
      return true;
    return false;
  }

  bool has_64bit_ids() const
  {
    if (blongids_)
      return true;
    return false;
  }

  int get_species_idx(const cosmo_species &s) const
  {
    switch (s)
    {
    case cosmo_species::dm:
      return 1;
    case cosmo_species::baryon:
      return 0;
    case cosmo_species::neutrino:
      return 6;
    }
    return -1;
  }

  void write_particle_data(const particle::container &pc, const cosmo_species &s, double Omega_species)
  {
    int sid = get_species_idx(s);
    assert(sid != -1);

    npart_[sid] = (pc.get_global_num_particles());
    npartTotal_[sid] = (uint32_t)(pc.get_global_num_particles());
    npartTotalHighWord_[sid] = (uint32_t)((pc.get_global_num_particles()) >> 32);

    std::array<double,7> particle_masses;
    if( pc.bhas_individual_masses_ )
      particle_masses[sid] = 0.0;
    else
      particle_masses[sid] = Omega_species * munit_ / pc.get_global_num_particles();

    const size_t global_num_particles =  pc.get_global_num_particles();

    // start by creating the full empty datasets in the file
    if (this_rank_ == 0) {

      HDFCreateGroup(fname_, std::string("PartType") + std::to_string(sid));

      if (this->has_64bit_reals())
      {
	HDFCreateEmptyDatasetVector<double>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), global_num_particles, true);
	HDFCreateEmptyDatasetVector<double>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), global_num_particles, true);
      }
      else
      {
	HDFCreateEmptyDatasetVector<float>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), global_num_particles, true);
	HDFCreateEmptyDatasetVector<float>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), global_num_particles, true);
      }

      if (this->has_64bit_ids())
	HDFCreateEmptyDataset<uint64_t>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), global_num_particles, true);
      else
	HDFCreateEmptyDataset<uint32_t>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), global_num_particles, true);

      if (this->has_64bit_reals())
	HDFCreateEmptyDataset<double>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), global_num_particles, true);
      else
	HDFCreateEmptyDataset<float>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), global_num_particles, true);

      if( bdobaryons_ && s == cosmo_species::baryon) {

	// note: despite this being a constant array we still need to handle it in a distributed way
	HDFCreateEmptyDataset<write_real_t>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/InternalEnergy"), global_num_particles, true);
	HDFCreateEmptyDataset<write_real_t>(fname_, std::string("PartType") + std::to_string(sid) + std::string("/SmoothingLength"), global_num_particles, true);
      }

      music::ilog << "Created empty arrays for PartType" << std::to_string(sid) << " into file " << fname_ << "." << std::endl;
    }

    // compute each rank's offset in the global array
    const size_t n_local = pc.get_local_num_particles();
    size_t offset = 0;
#ifdef USE_MPI
    MPI_Exscan(&n_local, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

    // now each node writes its own chunk in a round-robin fashion, appending at the end of the currently existing data 
    for (int rank = 0; rank < num_ranks_; ++rank) {

#ifdef USE_MPI
      // wait until the initialisation or the previous rank in the loop is done
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      if (rank == this_rank_) {

	//... write positions and velocities.....
	if (this->has_64bit_reals())
	  {
	    HDFWriteDatasetVectorChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), pc.positions64_, offset);
	    HDFWriteDatasetVectorChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), pc.velocities64_, offset);
	  }
	else
	  {
	    HDFWriteDatasetVectorChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), pc.positions32_, offset);
	    HDFWriteDatasetVectorChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), pc.velocities32_, offset);
	  }

	//... write ids.....
	if (this->has_64bit_ids())
	  HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), pc.ids64_, offset);
	else
	  HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), pc.ids32_, offset);

	//... write masses.....
	if( pc.bhas_individual_masses_ )
	  {
	    if (this->has_64bit_reals())
	      HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), pc.mass64_, offset);
	    else
	      HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), pc.mass32_, offset);
	  }
	else
	  {
	    std::vector<write_real_t> data( pc.get_local_num_particles(), particle_masses[sid]);

	    if (this->has_64bit_reals())
	      HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), data, offset);
	    else
	      HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), data, offset);
	  }

	// write GAS internal energy and smoothing length if baryons are enabled
	if(bdobaryons_ && s == cosmo_species::baryon) {

	  std::vector<write_real_t> data( pc.get_local_num_particles(), ceint_ );
	  HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/InternalEnergy"), data, offset);

	  data.assign( pc.get_local_num_particles(), h_);
	  HDFWriteDatasetChunk(fname_, std::string("PartType") + std::to_string(sid) + std::string("/SmoothingLength"), data, offset);
	}
      }

      if (this_rank_ == 0) {
	music::ilog << "Rank " << rank << " wrote its PartType" << std::to_string(sid) << " data to the IC file." << std::endl;
      }
    }

#ifdef USE_MPI
    // end with a barrier to make sure everyone is done before the destructor does its job
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
};

namespace
{
output_plugin_creator_concrete<swift_output_plugin<double>> creator301("SWIFT");
} // namespace

#endif
