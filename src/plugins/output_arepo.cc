
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
#ifdef USE_HDF5
#include <unistd.h> // for unlink
#include <output_plugin.hh>
#include "HDF_IO.hh"

template <typename T>
std::vector<T> from_6array(const T *a)
{
  return std::vector<T>{{a[0], a[1], a[2], a[3], a[4], a[5]}};
}

template <typename T>
std::vector<T> from_value(const T a)
{
  return std::vector<T>{{a}};
}

template <typename write_real_t>
class arepo_output_plugin : public output_plugin
{
  struct header_t
  {
    unsigned npart[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    unsigned int npartTotal[6];
    int flag_cooling;
    int num_files;
    double BoxSize;
    double Omega0;
    double OmegaBaryon;
    double OmegaLambda;
    double HubbleParam;
    int flag_stellarage;
    int flag_metals;
    unsigned int npartTotalHighWord[6];
    int flag_entropy_instead_u;
    int flag_doubleprecision;
  };

protected:
  int num_files_, num_simultaneous_writers_;
  header_t header_;
  real_t lunit_, vunit_, munit_;
  bool blongids_;
  unsigned hdf5_compression_level_;
  std::string this_fname_;
  double Tini_;
  unsigned pmgrid_;
  unsigned gridboost_;
  int doublePrec_;
  int doBaryons_;
  double softening_;

public:
  //! constructor
  explicit arepo_output_plugin(config_file &cf, std::unique_ptr<cosmology::calculator>& pcc)
      : output_plugin(cf, pcc, "AREPO-HDF5")
  {
    num_files_ = 1;
#ifdef USE_MPI
    // use as many output files as we have MPI tasks
    MPI_Comm_size(MPI_COMM_WORLD, &num_files_);
#endif
    real_t astart = 1.0 / (1.0 + cf_.get_value<double>("setup", "zstart"));
    const double rhoc = 27.7519737; // in h^2 1e10 M_sol / Mpc^3

    lunit_ = cf_.get_value<double>("setup", "BoxLength");
    vunit_ = lunit_ / std::sqrt(astart);
    munit_ = rhoc * std::pow(cf_.get_value<double>("setup", "BoxLength"), 3); // in 1e10 h^-1 M_sol
    
    blongids_ = cf_.get_value_safe<bool>("output", "UseLongids", false);
    num_simultaneous_writers_ = cf_.get_value_safe<int>("output", "NumSimWriters", num_files_);

    hdf5_compression_level_ = cf_.get_value_safe<unsigned>("output", "hdf5CompressionLevel", 0);

    for (int i = 0; i < 6; ++i)
    {
      header_.npart[i] = 0;
      header_.npartTotal[i] = 0;
      header_.npartTotalHighWord[i] = 0;
      header_.mass[i] = 0.0;
    }

    header_.time = astart;
    header_.redshift = 1.0 / astart - 1.0;
    header_.flag_sfr = 0;
    header_.flag_feedback = 0;
    header_.flag_cooling = 0;
    header_.num_files = num_files_;
    header_.BoxSize = lunit_;
    header_.Omega0 = pcc->cosmo_param_["Omega_m"];
    header_.OmegaLambda = pcc->cosmo_param_["Omega_DE"];
    header_.OmegaBaryon = pcc->cosmo_param_["Omega_b"];
    header_.HubbleParam = pcc->cosmo_param_["h"];
    header_.flag_stellarage = 0;
    header_.flag_metals = 0;
    header_.flag_entropy_instead_u = 0;
    header_.flag_doubleprecision = (typeid(write_real_t) == typeid(double)) ? true : false;

    // initial gas temperature
    double Tcmb0 = pcc->cosmo_param_["Tcmb"];
    double Omegab = pcc->cosmo_param_["Omega_b"];
    double h = pcc->cosmo_param_["h"], h2 = h*h;
    double adec = 1.0 / (160.0 * pow(Omegab * h2 / 0.022, 2.0 / 5.0));
    Tini_ = astart < adec ? Tcmb0 / astart : Tcmb0 / astart / astart * adec;

    music::ilog << "Setting initial gas temperature to T = " << Tini_ << "K/µ" << std::endl;

    // suggested PM res
    pmgrid_ = 2*cf_.get_value<double>("setup", "GridRes");
    gridboost_ = 1;
    softening_ = cf_.get_value<double>("setup", "BoxLength")/pmgrid_/20;
    doBaryons_ = cf_.get_value<bool>("setup", "DoBaryons");
#if !defined(USE_SINGLEPRECISION)
    doublePrec_ = 1;
#else
    doublePrec_ = 0;
#endif

    this_fname_ = fname_;
#ifdef USE_MPI
    int thisrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
    if (num_files_ > 1)
      this_fname_ += "." + std::to_string(thisrank);
#endif

    unlink(this_fname_.c_str());
    HDFCreateFile(this_fname_);
  }

  // use destructor to write header post factum
  ~arepo_output_plugin()
  {
    HDFCreateGroup(this_fname_, "Header");
    HDFWriteGroupAttribute(this_fname_, "Header", "NumPart_ThisFile", from_6array<unsigned>(header_.npart));
    HDFWriteGroupAttribute(this_fname_, "Header", "MassTable", from_6array<double>(header_.mass));
    HDFWriteGroupAttribute(this_fname_, "Header", "Time", from_value<double>(header_.time));
    HDFWriteGroupAttribute(this_fname_, "Header", "Redshift", from_value<double>(header_.redshift));
    HDFWriteGroupAttribute(this_fname_, "Header", "NumPart_Total", from_6array<unsigned>(header_.npartTotal));
    HDFWriteGroupAttribute(this_fname_, "Header", "NumPart_Total_HighWord", from_6array<unsigned>(header_.npartTotalHighWord));
    HDFWriteGroupAttribute(this_fname_, "Header", "NumFilesPerSnapshot", from_value<int>(header_.num_files));
    HDFWriteGroupAttribute(this_fname_, "Header", "BoxSize", from_value<double>(header_.BoxSize));
    HDFWriteGroupAttribute(this_fname_, "Header", "Omega0", from_value<double>(header_.Omega0));
    HDFWriteGroupAttribute(this_fname_, "Header", "OmegaBaryon", from_value<double>(header_.OmegaBaryon));
    HDFWriteGroupAttribute(this_fname_, "Header", "OmegaLambda", from_value<double>(header_.OmegaLambda));
    HDFWriteGroupAttribute(this_fname_, "Header", "HubbleParam", from_value<double>(header_.HubbleParam));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Sfr", from_value<int>(0));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Cooling", from_value<int>(0));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_StellarAge", from_value<int>(0));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Metals", from_value<int>(0));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Feedback", from_value<int>(0));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_DoublePrecision", (int)doublePrec_);
    // HDFWriteGroupAttribute(this_fname_, "Header", "Music_levelmin", levelmin_);
    // HDFWriteGroupAttribute(this_fname_, "Header", "Music_levelmax", levelmax_);
    // HDFWriteGroupAttribute(this_fname_, "Header", "Music_levelcounts", levelcounts);
    HDFWriteGroupAttribute(this_fname_, "Header", "haveBaryons", from_value<int>((int)doBaryons_));
    HDFWriteGroupAttribute(this_fname_, "Header", "longIDs", from_value<int>((int)blongids_));
    HDFWriteGroupAttribute(this_fname_, "Header", "suggested_pmgrid", from_value<int>(pmgrid_));
    HDFWriteGroupAttribute(this_fname_, "Header", "suggested_gridboost", from_value<int>(gridboost_));
    HDFWriteGroupAttribute(this_fname_, "Header", "suggested_highressoft", from_value<double>(softening_));
    HDFWriteGroupAttribute(this_fname_, "Header", "suggested_gas_Tinit", from_value<double>(Tini_));

    music::ilog << this->interface_name_ << " output plugin wrote to \'" << this_fname_ << "\'" << std::endl;
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
      return 3;
    }
    return -1;
  }

  void write_particle_data(const particle::container &pc, const cosmo_species &s, double Omega_species)
  {
    int sid = get_species_idx(s);

    assert(sid != -1);

    header_.npart[sid] = (pc.get_local_num_particles());
    header_.npartTotal[sid] = (uint32_t)(pc.get_global_num_particles());
    header_.npartTotalHighWord[sid] = (uint32_t)((pc.get_global_num_particles()) >> 32);

    if( pc.bhas_individual_masses_ ){
      header_.mass[sid] = 0.0;
    }else{
      header_.mass[sid] = Omega_species * munit_ / pc.get_global_num_particles();
    }

    HDFCreateGroup(this_fname_, std::string("PartType") + std::to_string(sid));

    //... write positions and velocities.....
    if (this->has_64bit_reals())
    {
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), pc.positions64_);
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), pc.velocities64_);
    }
    else
    {
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), pc.positions32_);
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), pc.velocities32_);
    }

    //... write ids.....
    if (this->has_64bit_ids()){
      HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), pc.ids64_);
    }else{
      HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), pc.ids32_);
    }

    //... write masses.....
    if( pc.bhas_individual_masses_ ){
      if (this->has_64bit_reals()){
        HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), pc.mass64_);
      }else{
        HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), pc.mass32_);
      }
    }

    // std::cout << ">>>A> " << header_.npart[sid] << std::endl;
  }
};

namespace
{
#if !defined(USE_SINGLEPRECISION)
output_plugin_creator_concrete<arepo_output_plugin<double>> creator880("AREPO");
#else
output_plugin_creator_concrete<arepo_output_plugin<float>> creator881("AREPO");
#endif
} // namespace

#endif
