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

#include <cassert>
#include <cmath>
#include <string>
#include <typeinfo>

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

template <typename T>
std::vector<T> recenter_pkd_coordinates(const std::vector<T> &coords)
{
  std::vector<T> shifted(coords);
  for (size_t i = 0; i < shifted.size(); ++i)
  {
    double x = std::fmod(static_cast<double>(shifted[i]), 1.0);
    if (x < 0.0)
      x += 1.0;
    shifted[i] = static_cast<T>(x - 0.5);
  }
  return shifted;
}

template <typename write_real_t>
class pkdgrav3_hdf5_output_plugin : public output_plugin
{
  struct header_t
  {
    size_t npart[6];
    size_t npart_total[6];
    unsigned npart_total_highword[6];
    double mass[6];

    double time;
    double redshift;
    double box_size;
    double omega_m;
    double omega_lambda;
    double omega_b;
    double hubble_param;

    int num_files;
    int flag_sfr;
    int flag_feedback;
    int flag_cooling;
    int flag_stellarage;
    int flag_metals;
    int flag_entropy_ics;
    int flag_doubleprecision;
  };

protected:
  int num_files_;
  header_t header_;
  real_t lunit_, vunit_, munit_;
  bool blongids_;
  std::string this_fname_;

  double dMsolUnit_;
  double dKpcUnit_;
  double dGasConst_;
  double dErgPerGmUnit_;
  double dGmPerCcUnit_;
  double dSecUnit_;
  double dKmPerSecUnit_;
  double dComovingGmPerCcUnit_;
  double astart_;

public:
  explicit pkdgrav3_hdf5_output_plugin(config_file &cf, std::unique_ptr<cosmology::calculator> &pcc)
      : output_plugin(cf, pcc, (std::string("PKDGRAV3-HDF5 ") + typeid(write_real_t).name()).c_str())
  {
    num_files_ = 1;
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &num_files_);
#endif

    const double boxlen = cf_.get_value<double>("setup", "BoxLength"); // Mpc/h
    const double zstart = cf_.get_value<double>("setup", "zstart");
    astart_ = 1.0 / (1.0 + zstart);
    const double h = pcc_->cosmo_param_["h"];

    // PKDGRAV3 snapshot coordinates are normalized to box units in [-0.5, 0.5).
    lunit_ = 1.0;
    // In PKDGRAV3 code units, the total mass in the box is Omega_m.
    munit_ = 1.0;

    // Unit definitions follow PKDGRAV3's cosmological unit conventions.
    constexpr double KBOLTZ = 1.3806485e-16;
    constexpr double MHYDR = 1.6735575e-24;
    constexpr double MSOLG = 1.98847e33;
    constexpr double GCGS = 6.67408e-8;
    constexpr double KPCCM = 3.085678e21;

    dKpcUnit_ = boxlen * 1.0e3 / h;
    const double H0_cgs = 100.0 * h * 1.0e5 / (1.0e3 * KPCCM);
    constexpr double PI = 3.14159265358979323846;
    const double critDens0_cgs = 3.0 * H0_cgs * H0_cgs / (8.0 * PI * GCGS);
    const double critDens0 = critDens0_cgs * (KPCCM * KPCCM * KPCCM / MSOLG);
    dMsolUnit_ = critDens0 * dKpcUnit_ * dKpcUnit_ * dKpcUnit_;

    dGasConst_ = dKpcUnit_ * KPCCM * KBOLTZ / (MHYDR * GCGS * dMsolUnit_ * MSOLG);
    dErgPerGmUnit_ = GCGS * dMsolUnit_ * MSOLG / (dKpcUnit_ * KPCCM);
    dGmPerCcUnit_ = (dMsolUnit_ * MSOLG) / std::pow(dKpcUnit_ * KPCCM, 3.0);
    dSecUnit_ = std::sqrt(1.0 / (dGmPerCcUnit_ * GCGS));
    dKmPerSecUnit_ = std::sqrt(GCGS * dMsolUnit_ * MSOLG / (dKpcUnit_ * KPCCM)) / 1.0e5;
    dComovingGmPerCcUnit_ = dGmPerCcUnit_;
 
    vunit_ = boxlen / (astart_ * dKmPerSecUnit_);

    blongids_ = cf_.get_value_safe<bool>("output", "UseLongids", true);

    for (int i = 0; i < 6; ++i)
    {
      header_.npart[i] = 0;
      header_.npart_total[i] = 0;
      header_.npart_total_highword[i] = 0;
      header_.mass[i] = 0.0;
    }

    header_.time = astart_;
    header_.redshift = zstart;
    header_.box_size = boxlen;
    header_.omega_m = pcc_->cosmo_param_["Omega_m"];
    header_.omega_lambda = pcc_->cosmo_param_["Omega_DE"];
    header_.omega_b = pcc_->cosmo_param_["Omega_b"];
    header_.hubble_param = h;
    header_.num_files = num_files_;
    header_.flag_sfr = 0;
    header_.flag_feedback = 0;
    header_.flag_cooling = 0;
    header_.flag_stellarage = 0;
    header_.flag_metals = 0;
    header_.flag_entropy_ics = 0;
    header_.flag_doubleprecision = (typeid(write_real_t) == typeid(double)) ? 1 : 0;

    std::string::size_type pos = fname_.find_last_of(".");
    std::string fname_prefix = (pos == std::string::npos) ? fname_ : fname_.substr(0, pos);
    std::string fname_suffix = (pos == std::string::npos) ? "hdf5" : fname_.substr(pos + 1);

    this_fname_ = fname_prefix;
#ifdef USE_MPI
    int thisrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisrank);
    if (num_files_ > 1)
      this_fname_ += "." + std::to_string(thisrank);
#endif
    this_fname_ += "." + fname_suffix;

    unlink(this_fname_.c_str());
    HDFCreateFile(this_fname_);
  }

  ~pkdgrav3_hdf5_output_plugin()
  {
    if (std::uncaught_exceptions())
      return;

    HDFCreateGroup(this_fname_, "Header");
    HDFWriteGroupAttribute(this_fname_, "Header", "NumPart_ThisFile", from_6array<size_t>(header_.npart));
    HDFWriteGroupAttribute(this_fname_, "Header", "NumPart_Total", from_6array<size_t>(header_.npart_total));
    HDFWriteGroupAttribute(this_fname_, "Header", "NumPart_Total_HighWord", from_6array<unsigned>(header_.npart_total_highword));
    HDFWriteGroupAttribute(this_fname_, "Header", "MassTable", from_6array<double>(header_.mass));
    HDFWriteGroupAttribute(this_fname_, "Header", "Time", from_value<double>(header_.time));
    HDFWriteGroupAttribute(this_fname_, "Header", "Redshift", from_value<double>(header_.redshift));
    HDFWriteGroupAttribute(this_fname_, "Header", "NumFilesPerSnapshot", from_value<int>(header_.num_files));
    HDFWriteGroupAttribute(this_fname_, "Header", "BoxSize", from_value<double>(header_.box_size));
    HDFWriteGroupAttribute(this_fname_, "Header", "Omega0", from_value<double>(header_.omega_m));
    HDFWriteGroupAttribute(this_fname_, "Header", "OmegaLambda", from_value<double>(header_.omega_lambda));
    HDFWriteGroupAttribute(this_fname_, "Header", "OmegaBaryon", from_value<double>(header_.omega_b));
    HDFWriteGroupAttribute(this_fname_, "Header", "HubbleParam", from_value<double>(header_.hubble_param));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Sfr", from_value<int>(header_.flag_sfr));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Feedback", from_value<int>(header_.flag_feedback));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Cooling", from_value<int>(header_.flag_cooling));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_StellarAge", from_value<int>(header_.flag_stellarage));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Metals", from_value<int>(header_.flag_metals));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_Entropy_ICs", from_value<int>(header_.flag_entropy_ics));
    HDFWriteGroupAttribute(this_fname_, "Header", "Flag_DoublePrecision", from_value<int>(header_.flag_doubleprecision));

    HDFCreateGroup(this_fname_, "Cosmology");
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "Omega_m", header_.omega_m);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "Omega_lambda", header_.omega_lambda);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "Omega_b", header_.omega_b);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "h", header_.hubble_param);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "Omega0", header_.omega_m);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "OmegaLambda", header_.omega_lambda);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "HubbleParam", header_.hubble_param);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "dOmega0", header_.omega_m);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "dLambda", header_.omega_lambda);
    HDFWriteGroupAttribute(this_fname_, "Cosmology", "dHubble0", header_.hubble_param);

    HDFCreateGroup(this_fname_, "Units");
    HDFWriteGroupAttribute(this_fname_, "Units", "MsolUnit", dMsolUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "KpcUnit", dKpcUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "GasConst", dGasConst_);
    HDFWriteGroupAttribute(this_fname_, "Units", "ErgPerGmUnit", dErgPerGmUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "GmPerCcUnit", dGmPerCcUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "SecUnit", dSecUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "KmPerSecUnit", dKmPerSecUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "ComovingGmPerCcUnit", dComovingGmPerCcUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dMsolUnit", dMsolUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dKpcUnit", dKpcUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dGasConst", dGasConst_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dErgPerGmUnit", dErgPerGmUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dGmPerCcUnit", dGmPerCcUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dSecUnit", dSecUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dKmPerSecUnit", dKmPerSecUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "dComovingGmPerCcUnit", dComovingGmPerCcUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "Unit mass in solar masses", dMsolUnit_);
    HDFWriteGroupAttribute(this_fname_, "Units", "Unit velocity in km/sec", dKmPerSecUnit_);

    music::ilog << this->interface_name_ << " output plugin wrote to '" << this_fname_ << "'" << std::endl;
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

    header_.npart[sid] = pc.get_local_num_particles();
    header_.npart_total[sid] = pc.get_global_num_particles();
    header_.npart_total_highword[sid] = static_cast<unsigned>(pc.get_global_num_particles() >> 32);

    if (pc.bhas_individual_masses_)
      header_.mass[sid] = 0.0;
    else
      header_.mass[sid] = Omega_species * munit_ / pc.get_global_num_particles();

    HDFCreateGroup(this_fname_, std::string("PartType") + std::to_string(sid));

    if (this->has_64bit_reals())
    {
      const auto shifted_pos = recenter_pkd_coordinates(pc.positions64_);
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), shifted_pos);
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), pc.velocities64_);
    }
    else
    {
      const auto shifted_pos = recenter_pkd_coordinates(pc.positions32_);
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Coordinates"), shifted_pos);
      HDFWriteDatasetVector(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Velocities"), pc.velocities32_);
    }

    if (this->has_64bit_ids())
      HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), pc.ids64_);
    else
      HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/ParticleIDs"), pc.ids32_);

    if (pc.bhas_individual_masses_)
    {
      if (this->has_64bit_reals())
        HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), pc.mass64_);
      else
        HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), pc.mass32_);
    }
    else
    {
      std::vector<write_real_t> masses(pc.get_local_num_particles(), static_cast<write_real_t>(header_.mass[sid]));
      HDFWriteDataset(this_fname_, std::string("PartType") + std::to_string(sid) + std::string("/Masses"), masses);
    }
  }
};

namespace
{
output_plugin_creator_concrete<pkdgrav3_hdf5_output_plugin<float>> creator993("pkdgrav3_hdf5");
output_plugin_creator_concrete<pkdgrav3_hdf5_output_plugin<float>> creator994("pkdgrav3");
output_plugin_creator_concrete<pkdgrav3_hdf5_output_plugin<float>> creator995("PKDGRAV3");
#if !defined(USE_SINGLEPRECISION)
output_plugin_creator_concrete<pkdgrav3_hdf5_output_plugin<double>> creator996("pkdgrav3_hdf5_double");
output_plugin_creator_concrete<pkdgrav3_hdf5_output_plugin<double>> creator997("pkdgrav3_double");
#endif
} // namespace

#endif
