// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn & Michael Buehlmann (this file)
//		 2026 by Matthieu Schaller (this file)
//		 2026 by Will Roper (this file)
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
#include "HDF_IO.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <output_plugin.hh>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#if defined(USE_CLASS)
// CLASS version information, recorded in the ICs metadata
extern "C" {
extern const char *CLASS_GIT_REV;
extern const char *CLASS_GIT_TAG;
}
#endif

/**
 * @brief Signed 64-bit type used for cell counts and particle offsets.
 *
 * This is deliberately long long, rather than int64_t: HDF_IO's GetDataType<>
 * dispatches on typeid and recognises long long, while int64_t is a different
 * type on some supported platforms.
 */
using swift_count_t = long long;

/** @brief Struct to track a contiguous section of a particle dataset written by
 * one MPI rank. */
struct swift_write_range {
  //! First particle index in the global dataset.
  hsize_t file_start;

  //! Number of particles in the section.
  hsize_t length;
};

/**
 * @brief Convert a seven-element particle-type array to an HDF5 attribute.
 *
 * @param a Values indexed by Gadget/SWIFT particle type.
 *
 * @return Copy of @p a in the container expected by HDF_IO.
 */
template <typename T> std::vector<T> from_7array(const std::array<T, 7> &a) {
  return std::vector<T>{{a[0], a[1], a[2], a[3], a[4], a[5], a[6]}};
}

/**
 * @brief Wrap a scalar in the one-element array used by SWIFT attributes.
 *
 * @param a Scalar attribute value.
 *
 * @return One-element vector containing @p a.
 */
template <typename T> std::vector<T> from_value(const T a) {
  return std::vector<T>{{a}};
}

/**
 * @brief Write one SWIFT-compatible HDF5 initial-conditions file.
 *
 * Particles are stored in SWIFT top-level-cell order to both reduce the
 * movement of data during the initial redistribute of particles but also enable
 * the use of SWIFT tools such as swiftsimio on IC files generated with
 * Monofonic.
 *
 * With parallel HDF5, every rank opens the shared file collectively and writes
 * a union of hyperslabs. With serial HDF5, ranks use the same layout but take
 * turns opening and writing the file. Rank zero creates datasets and metadata
 * in both modes.
 *
 * @tparam write_real_t Floating-point type stored in particle datasets.
 */
template <typename write_real_t>
class swift_output_plugin : public output_plugin {

protected:
  //! Number of MPI ranks, or one without MPI.
  int num_ranks_;

  //! Rank in MPI_COMM_WORLD, or zero without MPI.
  int this_rank_;

  //! Number of output files; SWIFT output always uses one.
  int num_files_;

  //! Position conversion factor to Mpc.
  real_t lunit_;

  //! Velocity conversion factor to km/s.
  real_t vunit_;

  //! Mass conversion factor to 1e10 solar masses.
  real_t munit_;

  //! Comoving box length in Mpc/h from the configuration.
  real_t boxsize_;

  //! Dimensionless Hubble parameter.
  real_t hubble_param_;

  //! Initial scale factor.
  real_t astart_;

  //! Initial redshift.
  real_t zstart_;

  //! Whether particle IDs are stored as 64-bit integers.
  bool blongids_;

  //! Whether gas particle fields are required.
  bool bdobaryons_;

  //! Number of top-level cells along each dimension.
  size_t cdim_;

  //! Total number of top-level cells (cdim_^3).
  size_t ncells_;

  //! Box side length in output units (Mpc, no h); alias for lunit_.
  double dim_;

  //! Top-level cell side length in Mpc.
  double cell_width_;

  //! Inverse top-level cell width, matching SWIFT iwidth.
  double iwidth_;

  //! Particle counts in this output file.
  std::array<uint64_t, 7> npart_;

  //! Low words of global particle counts.
  std::array<uint32_t, 7> npartTotal_;

  //! High words of global counts.
  std::array<uint32_t, 7> npartTotalHighWord_;

  //! Initial gas internal energy per unit mass in (km/s)^2 (physical).
  double internal_energy_;

  //! Initial gas smoothing length in Mpc (comoving).
  double smoothing_length_;

  //! Initial gas temperature over the mean molecular weight, in K.
  double gas_temperature_over_mu_;

  //! Mean particle mass of each particle type, for the ICs metadata.
  std::array<double, 7> mean_particle_mass_;

public:
  /**
   * @brief Initialise units, cell geometry, and rank-independent metadata.
   *
   * Rank zero creates the file and writes /Units, /Code,
   * /ICs_parameters, and the species-independent part of /Cells.
   *
   * Configuration keys used from [output] are filename, UseLongids, and
   * optional top_level_cells (default 32), which must be in [1, 1290]; the
   * upper bound keeps cdim^3 within the int32 'nr_cells' attribute and the
   * int MPI counts of the cell count reductions.
   *
   * @param cf Complete monofonIC configuration.
   * @param pcc Cosmology calculator owned by the IC generator.
   *
   * @throws std::runtime_error If top_level_cells is out of range.
   */
  explicit swift_output_plugin(config_file &cf,
                               std::unique_ptr<cosmology::calculator> &pcc)
      : output_plugin(cf, pcc, "SWIFT") {

    // SWIFT uses a single IC file; rank info is filled in below under MPI
    num_files_ = 1;
    this_rank_ = 0;
    num_ranks_ = 1;

    // Critical density in h^2 1e10 M_sol / Mpc^3, using SWIFT's internal
    // physical constants
    const double rhoc = 27.7536609198;

    // Cosmology and box size from the configuration
    hubble_param_ = pcc->cosmo_param_["h"];
    zstart_ = cf_.get_value<double>("setup", "zstart");
    astart_ = 1.0 / (1.0 + zstart_);
    boxsize_ = cf_.get_value<double>("setup", "BoxLength");

    // Conversion factors to the output units: Mpc (no h), km/s and
    // 1e10 M_sol
    lunit_ = boxsize_ / hubble_param_;
    vunit_ = boxsize_;
    munit_ = rhoc * std::pow(boxsize_, 3) / hubble_param_;

    // Output options
    blongids_ = cf_.get_value_safe<bool>("output", "UseLongids", true);
    bdobaryons_ = cf_.get_value<bool>("setup", "DoBaryons");

    // SWIFT top-level cell grid
    cdim_ = cf_.get_value_safe<size_t>("output", "top_level_cells", 32);
    dim_ = lunit_;

    // Validate the grid. cdim^3 has to fit both the int32 'nr_cells'
    // attribute and the int MPI counts used for the per-cell count reductions,
    // and 1290^3 is the largest cube below INT_MAX
    if (cdim_ == 0) {
      throw std::runtime_error(
          "SWIFT output: 'top_level_cells' must be larger than zero.");
    }
    if (cdim_ > 1290) {
      throw std::runtime_error("SWIFT output: 'top_level_cells' must not "
                               "exceed 1290 (cell count overflow).");
    }

    // Derived cell geometry
    ncells_ = cdim_ * cdim_ * cdim_;
    cell_width_ = dim_ / double(cdim_);
    iwidth_ = double(cdim_) / dim_;

    // Particle counts and masses are filled in per species as they are
    // written
    for (int i = 0; i < 7; ++i) {
      npart_[i] = 0;
      npartTotal_[i] = 0;
      npartTotalHighWord_[i] = 0;
      mean_particle_mass_[i] = 0.0;
    }

#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks_);
#endif

    // Gas needs an initial internal energy and smoothing length
    if (bdobaryons_) {

      const double gamma =
          cf_.get_value_safe<double>("cosmology", "gamma", 5.0 / 3.0);
      const double YHe = pcc_->cosmo_param_["YHe"];
      const double omegab = pcc_->cosmo_param_["Omega_b"];
      const double Tcmb0 = pcc_->cosmo_param_["Tcmb"];

      // Gas temperature follows the CMB until thermal decoupling at adec and
      // then cools adiabatically as a^-2. The internal energy is
      // k_B T / (mu m_p (gamma - 1)), using SWIFT's cgs values of k_B and m_p
      // (physical_constants_cgs.h), converted from (cm/s)^2 to (km/s)^2. SWIFT
      // reads IC internal energies as physical
      const double npol = (fabs(1.0 - gamma) > 1e-7) ? 1.0 / (gamma - 1.) : 1.0;
      const double unitv = 1e5;
      const double adec =
          1.0 / (160. * std::pow(omegab * hubble_param_ * hubble_param_ / 0.022,
                                 2.0 / 5.0));
      const double Tini =
          astart_ < adec ? Tcmb0 / astart_ : Tcmb0 / astart_ / astart_ * adec;
      const double mu =
          (Tini > 1.e4) ? 4.0 / (8. - 5. * YHe) : 4.0 / (1. + 3. * (1. - YHe));
      const double k_B = 1.380649e-16;
      const double m_p = 1.67262192369e-24;
      internal_energy_ = k_B / m_p * Tini * npol / mu / unitv / unitv;
      gas_temperature_over_mu_ = Tini / mu;

      // SWIFT requires a SmoothingLength field in gas ICs but recomputes h in
      // its first density loop, so the grid spacing is only an initial guess.
      // bcc and fcc loads have a smaller mean particle separation
      smoothing_length_ =
          boxsize_ / hubble_param_ / cf_.get_value<double>("setup", "GridRes");
    }

    // Report the output setup, matching the aligned "key : value" style of
    // the rest of the log (only rank 0 logs at info level)
    this->log_setup();

    // Only rank 0 creates the file and writes the run metadata
    if (this_rank_ != 0) {
      return;
    }

    // Create the output file, replacing any existing one
    HDFCreateFile(fname_);

    // Units group, using the physical constants assumed internally by SWIFT:
    // 1e10 Msun in g, 1 Mpc in cm, a time unit giving 1 km/s velocities,
    // 1 Ampere and 1 Kelvin. Values are written as 1-element arrays, the
    // shape SWIFT uses and swiftsimio expects
    HDFCreateGroup(fname_, "Units");
    HDFWriteGroupAttribute(fname_, "Units", "Unit mass in cgs (U_M)",
                           from_value<double>(1.98841e43));
    HDFWriteGroupAttribute(fname_, "Units", "Unit length in cgs (U_L)",
                           from_value<double>(3.08567758149e24));
    HDFWriteGroupAttribute(fname_, "Units", "Unit time in cgs (U_t)",
                           from_value<double>(3.08567758149e19));
    HDFWriteGroupAttribute(fname_, "Units", "Unit current in cgs (U_I)",
                           from_value<double>(1.0));
    HDFWriteGroupAttribute(fname_, "Units", "Unit temperature in cgs (U_T)",
                           from_value<double>(1.0));

    // swiftsimio identifies SWIFT files through this group
    HDFCreateGroup(fname_, "Code");
    HDFWriteGroupAttribute(fname_, "Code", "Code", std::string("SWIFT"));

    this->write_common_cell_metadata();
  }

  /**
   * @brief Destructor: finalise the IC file by writing the /Header group.
   *
   * The header is written here because it needs the particle counts of every
   * species, which are only known once all write_particle_data() calls have
   * finished. Only rank zero writes it. The IC generator destroys output
   * plugins before finalising MPI, so this runs while MPI is still active.
   * Header creation is skipped during exception unwinding to avoid presenting
   * a partially written file as complete.
   */
  ~swift_output_plugin() {
    if (!std::uncaught_exceptions()) {
      if (this_rank_ == 0) {

        // Record everything needed to reproduce these ICs; SWIFT copies this
        // group into every snapshot
        this->write_ics_metadata();

        // Write Standard Gadget / SWIFT hdf5 header
        HDFCreateGroup(fname_, "Header");

        // BoxSize is in Mpc, not Mpc/h for SWIFT
        HDFWriteGroupAttribute(
            fname_, "Header", "BoxSize",
            std::vector<double>(3, boxsize_ / hubble_param_));

        // Now we can write the particle counts and masses, because we have them
        HDFWriteGroupAttribute(fname_, "Header", "NumPart_Total",
                               from_7array<unsigned>(npartTotal_));
        HDFWriteGroupAttribute(fname_, "Header", "NumPart_Total_HighWord",
                               from_7array<unsigned>(npartTotalHighWord_));
        HDFWriteGroupAttribute(fname_, "Header", "NumPart_ThisFile",
                               from_7array<uint64_t>(npart_));
        // Masses are always written per particle, so the mass table is zero
        HDFWriteGroupAttribute(fname_, "Header", "MassTable",
                               std::vector<double>(7, 0.0));

        // Cosmological parameters and initial time
        HDFWriteGroupAttribute(fname_, "Header", "Scale-factor",
                               from_value<double>(astart_));
        HDFWriteGroupAttribute(fname_, "Header", "Time",
                               from_value<double>(astart_));
        HDFWriteGroupAttribute(fname_, "Header", "Redshift",
                               from_value<double>(zstart_));

        HDFWriteGroupAttribute(fname_, "Header", "Flag_Entropy_ICs",
                               from_value<int>(0));
        HDFWriteGroupAttribute(fname_, "Header", "Dimension",
                               from_value<int>(3));

        HDFWriteGroupAttribute(fname_, "Header", "NumFilesPerSnapshot",
                               from_value<int>(num_files_));

        // Report the finished file and how SWIFT should read it
        music::ilog << "Done writing SWIFT IC file to " << fname_ << std::endl;
        music::ilog << "Note that the IC file does not contain any h-factors "
                       "nor any extra sqrt(a)-factors for the velocities"
                    << std::endl;
        music::ilog
            << "The SWIFT parameters 'InitialConditions:cleanup_h_factors' and "
               "'InitialConditions:cleanup_velocity_factors' should hence "
               "*NOT* be set!"
            << std::endl;
      }
    }
  }

  /**
   * @brief Tell the IC generator whether to produce particles or grid fields.
   *
   * Required by output_plugin. SWIFT ICs are always particle based, so every
   * species is requested as particles.
   *
   * @return output_type::particles for every species.
   */
  output_type write_species_as(const cosmo_species &) const {
    return output_type::particles;
  }

  /**
   * @brief Position unit of the output.
   *
   * @return Position conversion factor from box units to Mpc.
   */
  real_t position_unit() const { return lunit_; }

  /**
   * @brief Velocity unit of the output.
   *
   * @return Velocity conversion factor to km/s.
   */
  real_t velocity_unit() const { return vunit_; }

  /**
   * @brief Mass unit of the output.
   *
   * @return Mass conversion factor to 1e10 solar masses.
   */
  real_t mass_unit() const { return munit_; }

  /**
   * @brief Floating-point precision of the particle fields.
   *
   * @return Whether particle floating-point fields are stored as doubles.
   */
  bool has_64bit_reals() const {
    if (typeid(write_real_t) == typeid(double)) {
      return true;
    }
    return false;
  }

  /**
   * @brief Integer width of the particle IDs.
   *
   * @return Whether particle IDs are stored as unsigned 64-bit integers.
   */
  bool has_64bit_ids() const {
    if (blongids_) {
      return true;
    }
    return false;
  }

  /**
   * @brief Map a monofonIC species to its Gadget/SWIFT particle-type index.
   *
   * @param s Species to map.
   *
   * @return Particle-type index, or -1 for an unsupported species.
   */
  int get_species_idx(const cosmo_species &s) const {
    switch (s) {
    case cosmo_species::dm:
      return 1;
    case cosmo_species::baryon:
      return 0;
    case cosmo_species::neutrino:
      return 6;
    }
    return -1;
  }

protected:
  /**
   * @brief Synchronise all ranks; a no-op in non-MPI builds.
   */
  inline void barrier() const {
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  /**
   * @brief Periodically wrap one coordinate into [0, dim_) at write precision.
   *
   * The wrapped value is rounded to the type stored in the file before the
   * upper edge is checked, so that sorting and writing see exactly the same
   * coordinate and every particle is assigned to the cell its stored
   * position lies in.
   *
   * @tparam T Floating-point type of the stored coordinate.
   * @param x Coordinate in output length units.
   *
   * @return Wrapped coordinate; a rounded upper boundary maps to zero.
   */
  template <typename T> inline T wrap_coord(double x) const {
    double y = std::fmod(x, dim_);
    if (y < 0.0) {
      y += dim_;
    }
    // Rounding, in the wrap or to T, can land exactly on the upper edge
    const T yt = (T)y;
    if (!(yt < (T)dim_)) {
      return T(0);
    }
    return yt;
  }

  /**
   * @brief Compute SWIFT's flattened top-level-cell index.
   *
   * Coordinates are expected in [0, dim_). The upper clamp protects
   * against floating-point roundoff at cell boundaries.
   *
   * @param x Wrapped three-dimensional position in output length units.
   *
   * @return Cell index with z varying fastest and x slowest, matching SWIFT's
   * cell_getid(cdim, i, j, k) macro.
   */
  inline size_t cell_getid(const double x[3]) const {
    size_t ijk[3];
    for (int d = 0; d < 3; ++d) {
      const size_t i = (size_t)(x[d] * iwidth_);
      ijk[d] = (i >= cdim_) ? cdim_ - 1 : i;
    }
    return ijk[2] + cdim_ * (ijk[1] + cdim_ * ijk[0]);
  }

  /**
   * @brief Copy one particle field into top-level-cell order for writing.
   *
   * The particle container itself is never reordered. Instead,
   * sort_particles_into_cells() produces a permutation and each field is
   * copied through it here, so that the buffer handed to HDF5 matches the
   * cell-ordered file ranges. Using the same permutation for every field
   * keeps particle i at the same index in every dataset.
   *
   * @tparam T Field element type.
   * @param src Flat field in the container's original particle order.
   * @param perm Source particle index for each output slot, so that
   * out[i] = src[perm[i]].
   * @param width Values per particle: 3 for vectors, 1 for scalars.
   *
   * @return Cell-ordered copy of src.
   */
  template <typename T>
  static std::vector<T> apply_cell_order(const std::vector<T> &src,
                                         const std::vector<size_t> &perm,
                                         size_t width) {
    std::vector<T> out(perm.size() * width);
    for (size_t i = 0; i < perm.size(); ++i) {
      for (size_t d = 0; d < width; ++d) {
        out[i * width + d] = src[perm[i] * width + d];
      }
    }
    return out;
  }

  /**
   * @brief Sort this rank's particles into SWIFT's top-level cells.
   *
   * Builds a stable counting-sort permutation into cell order, together with
   * the local particle count of every cell. Positions are wrapped before
   * assigning cells.
   *
   * @tparam pos_t Floating-point type of the position array.
   * @param pos Flat local position array in (x, y, z) order.
   * @param n_local Number of local particles.
   * @param[out] perm Source indices ordered by cell, then original local order.
   * @param[out] counts Number of local particles in every cell.
   */
  template <typename pos_t>
  void sort_particles_into_cells(const std::vector<pos_t> &pos, size_t n_local,
                                 std::vector<size_t> &perm,
                                 std::vector<swift_count_t> &counts) const {

    // Initialize the per-cell counts
    counts.assign(ncells_, 0);

    // Reserve the permutation vector and a temporary cell ID array
    std::vector<uint32_t> cellid(n_local);

    // Loop over local particles finding their cell
    for (size_t i = 0; i < n_local; ++i) {

      // Get the wrapped position in output units
      double x[3];
      for (int d = 0; d < 3; ++d) {

        // Wrap to [0, dim_) if needs be
        const double xd = (double)pos[3 * i + d];
        x[d] = (double)this->wrap_coord<pos_t>(xd);
      }

      // Compute the cell index
      const size_t c = this->cell_getid(x);

      // Store the cell index
      cellid[i] = (uint32_t)c;

      // Update the count for this cell
      ++counts[c];
    }

    // Build the stable counting-sort permutation. next[c] is the next free slot
    // for cell c in the output array, and is incremented as each particle is
    // placed.
    std::vector<size_t> next(ncells_);
    size_t acc = 0;
    for (size_t c = 0; c < ncells_; ++c) {
      next[c] = acc;
      acc += (size_t)counts[c];
    }

    // Fill the permutation array in cell order
    perm.resize(n_local);
    for (size_t i = 0; i < n_local; ++i) {
      perm[next[cellid[i]]++] = i;
    }
  }

  /**
   * @brief Build the global top-level cell lookup table across all ranks.
   *
   * Combines the local cell counts across ranks and derives where
   * each cell, and this rank's share of it, starts in the particle species
   * datasets. rank_offset[c] equals the global start of cell c plus
   * contributions from lower-numbered ranks. This produces a deterministic
   * cell-major, rank-major layout without redistributing particles between
   * ranks.
   *
   * @param local_counts This rank's particle count in each cell.
   * @param[out] global_counts Global particle count in each cell.
   * @param[out] cell_offset Start of each cell in the species dataset.
   * @param[out] rank_offset Start of this rank's contribution to each cell.
   */
  void build_cell_lookup_table(const std::vector<swift_count_t> &local_counts,
                               std::vector<swift_count_t> &global_counts,
                               std::vector<swift_count_t> &cell_offset,
                               std::vector<swift_count_t> &rank_offset) const {

    // Compute the global counts and the prefix sum of local counts across ranks
    global_counts = local_counts;
    std::vector<swift_count_t> preceding(ncells_, 0);

#ifdef USE_MPI
    // Collect the counts from all ranks
    MPI_Allreduce(MPI_IN_PLACE, global_counts.data(), (int)ncells_,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Exscan(local_counts.data(), preceding.data(), (int)ncells_,
               MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    // MPI_Exscan leaves rank 0 undefined
    if (this_rank_ == 0) {
      std::fill(preceding.begin(), preceding.end(), 0);
    }
#endif

    // Build the cell offsets and this rank's contribution to each cell
    cell_offset.resize(ncells_);
    rank_offset.resize(ncells_);
    swift_count_t off = 0;
    for (size_t c = 0; c < ncells_; ++c) {
      cell_offset[c] = off;
      rank_offset[c] = off + preceding[c];
      off += global_counts[c];
    }
  }

  /**
   * @brief Turn this rank's cell offsets into contiguous file ranges.
   *
   * Cells this rank has no particles in are skipped, and ranges for adjacent
   * cells are merged so that each one becomes a single hyperslab.
   *
   * @param local_counts This rank's particle count in each cell.
   * @param rank_offset Start of this rank's contribution to each cell.
   *
   * @return File ranges in cell order, matching the cell-ordered particles.
   */
  std::vector<swift_write_range>
  build_write_ranges(const std::vector<swift_count_t> &local_counts,
                     const std::vector<swift_count_t> &rank_offset) const {

    // Merge adjacent cells into contiguous file ranges. Each range is a
    // hyperslab in the file, and the union of all hyperslabs is the selection
    // for this rank's write.
    std::vector<swift_write_range> ranges;
    for (size_t c = 0; c < ncells_; ++c) {

      // Nothing to do for empty cells
      if (local_counts[c] == 0) {
        continue;
      }

      // Get this rank's start in the file for this cell
      const hsize_t fs = (hsize_t)rank_offset[c];

      // If the last range ends where this one starts, merge them
      if (!ranges.empty() &&
          ranges.back().file_start + ranges.back().length == fs) {
        ranges.back().length += (hsize_t)local_counts[c];
      }

      // Otherwise, start a new range
      else {
        ranges.push_back({fs, (hsize_t)local_counts[c]});
      }
    }

    return ranges;
  }

  /**
   * @brief Attach SWIFT unit exponents to one particle dataset.
   *
   * Current and temperature exponents, plus the Hubble-factor exponent, are
   * zero for every field written by this plugin.
   *
   * @param dset Dataset path.
   * @param uM Mass-unit exponent.
   * @param uL Length-unit exponent.
   * @param ut Time-unit exponent.
   * @param a_exp Cosmological scale-factor exponent.
   */
  void write_unit_attributes(const std::string &dset, double uM, double uL,
                             double ut, double a_exp) const {
    HDFWriteDatasetAttribute(fname_, dset, "U_M exponent",
                             std::vector<double>{uM});
    HDFWriteDatasetAttribute(fname_, dset, "U_L exponent",
                             std::vector<double>{uL});
    HDFWriteDatasetAttribute(fname_, dset, "U_t exponent",
                             std::vector<double>{ut});
    HDFWriteDatasetAttribute(fname_, dset, "U_I exponent",
                             std::vector<double>{0.0});
    HDFWriteDatasetAttribute(fname_, dset, "U_T exponent",
                             std::vector<double>{0.0});
    HDFWriteDatasetAttribute(fname_, dset, "a-scale exponent",
                             std::vector<double>{a_exp});
    HDFWriteDatasetAttribute(fname_, dset, "h-scale exponent",
                             std::vector<double>{0.0});
  }

  /**
   * @brief Create empty particle datasets and their unit metadata.
   *
   * Called only by rank zero. Serial HDF5 output retains checksums and
   * compression. Parallel output disables filters because collective writes to
   * filtered datasets are not portable across supported HDF5 versions.
   *
   * @param grp Particle group name, for example PartType1.
   * @param global_num_particles Dataset length.
   * @param s Species used to decide whether gas-only fields are needed.
   */
  void create_species_datasets(const std::string &grp,
                               size_t global_num_particles,
                               const cosmo_species &s) const {
#ifdef USE_PARALLEL_HDF5
    const bool filter = false;
#else
    const bool filter = true;
#endif

    HDFCreateGroup(fname_, grp);

    // Call the appropriate template instantiation for the floating-point type
    if (this->has_64bit_reals()) {
      HDFCreateEmptyDatasetVector<double>(fname_, grp + "/Coordinates",
                                          global_num_particles, filter);
      HDFCreateEmptyDatasetVector<double>(fname_, grp + "/Velocities",
                                          global_num_particles, filter);
      HDFCreateEmptyDataset<double>(fname_, grp + "/Masses",
                                    global_num_particles, filter);
    } else {
      HDFCreateEmptyDatasetVector<float>(fname_, grp + "/Coordinates",
                                         global_num_particles, filter);
      HDFCreateEmptyDatasetVector<float>(fname_, grp + "/Velocities",
                                         global_num_particles, filter);
      HDFCreateEmptyDataset<float>(fname_, grp + "/Masses",
                                   global_num_particles, filter);
    }

    // Create the ParticleIDs dataset with the appropriate integer type
    if (this->has_64bit_ids()) {
      HDFCreateEmptyDataset<uint64_t>(fname_, grp + "/ParticleIDs",
                                      global_num_particles, filter);
    } else {
      HDFCreateEmptyDataset<uint32_t>(fname_, grp + "/ParticleIDs",
                                      global_num_particles, filter);
    }

    // Write the unit exponents for every dataset. SWIFT expects the units to be
    this->write_unit_attributes(grp + "/Coordinates", 0.0, 1.0, 0.0, 1.0);
    this->write_unit_attributes(grp + "/Velocities", 0.0, 1.0, -1.0, 0.0);
    this->write_unit_attributes(grp + "/Masses", 1.0, 0.0, 0.0, 0.0);
    this->write_unit_attributes(grp + "/ParticleIDs", 0.0, 0.0, 0.0, 0.0);

    // Handle baryons if appropriate. SWIFT expects an InternalEnergy and
    // SmoothingLength dataset for baryons
    if (bdobaryons_ && s == cosmo_species::baryon) {
      HDFCreateEmptyDataset<write_real_t>(fname_, grp + "/InternalEnergy",
                                          global_num_particles, filter);
      HDFCreateEmptyDataset<write_real_t>(fname_, grp + "/SmoothingLength",
                                          global_num_particles, filter);
      this->write_unit_attributes(grp + "/InternalEnergy", 0.0, 2.0, -2.0, 0.0);
      this->write_unit_attributes(grp + "/SmoothingLength", 0.0, 1.0, 0.0, 1.0);
    }
  }

  /**
   * @brief Write one cell-ordered field into this rank's HDF5 hyperslabs.
   *
   * All ranges are ORed into one file-space selection. Ranks with no local
   * particles use null/empty selections but still participate in collective
   * writes, as required by parallel HDF5.
   *
   * @tparam T Dataset element type.
   * @param fid Open HDF5 file handle.
   * @param name Dataset path.
   * @param ranges Disjoint particle ranges owned by this rank.
   * @param data Packed field data corresponding to @p ranges in order.
   * @param width Elements per particle, normally one or three.
   * @param dxpl HDF5 transfer property list (H5P_DEFAULT or collective).
   * @throws std::runtime_error If the dataset cannot be written.
   */
  template <typename T>
  static void write_ranges(hid_t fid, const std::string &name,
                           const std::vector<swift_write_range> &ranges,
                           const std::vector<T> &data, size_t width,
                           hid_t dxpl) {

    // Open the dataset created by rank 0
    hid_t dset = H5Dopen(fid, name.c_str());

    // Get the file space to select this rank's ranges in
    hid_t fspace = H5Dget_space(dset);
    hid_t mspace;

    // Ranks without local particles still take part in the (collective) write,
    // with empty memory and file selections
    if (ranges.empty()) {
      mspace = H5Screate(H5S_NULL);
      H5Sselect_none(fspace);
    }

    // Otherwise, select the union of this rank's ranges in the file
    else {

      // Memory space covering the packed, cell-ordered buffer
      hsize_t mdims[2] = {(hsize_t)(data.size() / width), (hsize_t)width};
      mspace = H5Screate_simple(width > 1 ? 2 : 1, mdims, NULL);

      // OR each range into the file selection, the first one replacing it
      bool first = true;
      for (const auto &r : ranges) {
        hsize_t start[2] = {r.file_start, 0};
        hsize_t count[2] = {r.length, (hsize_t)width};
        H5Sselect_hyperslab(fspace, first ? H5S_SELECT_SET : H5S_SELECT_OR,
                            start, NULL, count, NULL);
        first = false;
      }
    }

    // Write the buffer into the selected ranges
    if (H5Dwrite(dset, GetDataType<T>(), mspace, fspace, dxpl,
                 data.empty() ? NULL : &data[0]) < 0) {
      throw std::runtime_error("SWIFT output: failed to write dataset " + name);
    }

    // Clean up
    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(dset);
  }

  /**
   * @brief Pack and write every dataset for one species using a shared
   * ordering.
   *
   * Coordinates are periodically wrapped before writing so they lie inside
   * the top-level cells they were sorted into. Every other field is
   * permuted identically, preserving particle correspondence.
   *
   * @param fid Open HDF5 file handle.
   * @param dxpl Dataset transfer property list.
   * @param grp Particle group name.
   * @param pc Source particle container.
   * @param s Particle species.
   * @param perm Stable local cell-order permutation.
   * @param ranges Global file ranges owned by this rank.
   * @param particle_mass Constant mass used when no individual masses exist.
   */
  void write_all_fields(hid_t fid, hid_t dxpl, const std::string &grp,
                        const particle::container &pc, const cosmo_species &s,
                        const std::vector<size_t> &perm,
                        const std::vector<swift_write_range> &ranges,
                        double particle_mass) const {

    // Number of local particles, all of which are written
    const size_t n_local = perm.size();

    // Write the positions, wrapped exactly as when they were sorted into cells
    if (this->has_64bit_reals()) {
      auto buf = apply_cell_order(pc.positions64_, perm, 3);
      for (auto &x : buf) {
        x = this->wrap_coord<double>(x);
      }
      write_ranges(fid, grp + "/Coordinates", ranges, buf, 3, dxpl);
    } else {
      auto buf = apply_cell_order(pc.positions32_, perm, 3);
      for (auto &x : buf) {
        x = this->wrap_coord<float>((double)x);
      }
      write_ranges(fid, grp + "/Coordinates", ranges, buf, 3, dxpl);
    }

    // Write the velocities
    if (this->has_64bit_reals()) {
      write_ranges(fid, grp + "/Velocities", ranges,
                   apply_cell_order(pc.velocities64_, perm, 3), 3, dxpl);
    } else {
      write_ranges(fid, grp + "/Velocities", ranges,
                   apply_cell_order(pc.velocities32_, perm, 3), 3, dxpl);
    }

    // Write the particle IDs with the appropriate integer type
    if (this->has_64bit_ids()) {
      write_ranges(fid, grp + "/ParticleIDs", ranges,
                   apply_cell_order(pc.ids64_, perm, 1), 1, dxpl);
    } else {
      write_ranges(fid, grp + "/ParticleIDs", ranges,
                   apply_cell_order(pc.ids32_, perm, 1), 1, dxpl);
    }

    // Write the masses, either per particle from the container...
    if (pc.bhas_individual_masses_) {
      if (this->has_64bit_reals()) {
        write_ranges(fid, grp + "/Masses", ranges,
                     apply_cell_order(pc.mass64_, perm, 1), 1, dxpl);
      } else {
        write_ranges(fid, grp + "/Masses", ranges,
                     apply_cell_order(pc.mass32_, perm, 1), 1, dxpl);
      }
    }

    // ...or as a constant for every particle
    else {
      if (this->has_64bit_reals()) {
        write_ranges(fid, grp + "/Masses", ranges,
                     std::vector<double>(n_local, particle_mass), 1, dxpl);
      } else {
        write_ranges(fid, grp + "/Masses", ranges,
                     std::vector<float>(n_local, (float)particle_mass), 1,
                     dxpl);
      }
    }

    // Write the gas internal energy and smoothing length, which are constant
    // for every particle, if baryons are enabled
    if (bdobaryons_ && s == cosmo_species::baryon) {
      write_ranges(
          fid, grp + "/InternalEnergy", ranges,
          std::vector<write_real_t>(n_local, (write_real_t)internal_energy_), 1,
          dxpl);
      write_ranges(
          fid, grp + "/SmoothingLength", ranges,
          std::vector<write_real_t>(n_local, (write_real_t)smoothing_length_),
          1, dxpl);
    }
  }

  /**
   * @brief Open the shared file and write all fields through the selected
   * backend.
   *
   * Parallel-HDF5 builds perform one collective file session on all ranks.
   * Other builds serialise rank access and perform one file session per rank.
   *
   * @param grp Particle group name.
   * @param pc Source particle container.
   * @param s Particle species.
   * @param perm Stable local cell-order permutation.
   * @param ranges Global file ranges corresponding to the permuted particles.
   * @param particle_mass Constant particle mass, ignored for individual masses.
   */
  void write_species_fields(const std::string &grp,
                            const particle::container &pc,
                            const cosmo_species &s,
                            const std::vector<size_t> &perm,
                            const std::vector<swift_write_range> &ranges,
                            double particle_mass) const {
#ifdef USE_PARALLEL_HDF5
    // Open the shared file collectively on every rank
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t fid = H5Fopen(fname_.c_str(), H5F_ACC_RDWR, fapl);
    if (fid < 0) {
      throw std::runtime_error("SWIFT output: cannot open " + fname_ +
                               " for collective writing.");
    }

    // Use collective transfers for every dataset write
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

    // Write every field
    this->write_all_fields(fid, dxpl, grp, pc, s, perm, ranges, particle_mass);

    // Clean up
    H5Pclose(dxpl);
    H5Fclose(fid);
    H5Pclose(fapl);
#else
    // Without parallel HDF5 the ranks take turns on the shared file, each one
    // opening it once
    for (int rank = 0; rank < num_ranks_; ++rank) {

      // Wait for the previous rank to close the file
      this->barrier();
      if (rank != this_rank_) {
        continue;
      }

      // Open the file, write every field and close it again
      hid_t fid = H5Fopen(fname_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      if (fid < 0) {
        throw std::runtime_error("SWIFT output: cannot open " + fname_ +
                                 " for writing.");
      }
      this->write_all_fields(fid, H5P_DEFAULT, grp, pc, s, perm, ranges,
                             particle_mass);
      H5Fclose(fid);
    }
#endif
  }

  /**
   * @brief Position of the same point within every top-level cell.
   *
   * The point sits at fraction @p frac of the cell width along each axis, so
   * 0 gives the lower corner, 0.5 the centre and 1 the upper corner. Cells are
   * flattened in cell_getid() order, as in every /Cells dataset.
   *
   * @param frac Fractional position within each cell.
   *
   * @return Positions in Mpc, flattened as (cell, dimension).
   */
  std::vector<double> cell_positions(double frac) const {
    std::vector<double> pos(3 * ncells_);
    for (size_t ix = 0; ix < cdim_; ++ix) {
      for (size_t iy = 0; iy < cdim_; ++iy) {
        for (size_t iz = 0; iz < cdim_; ++iz) {
          const size_t c = (ix * cdim_ + iy) * cdim_ + iz;
          pos[3 * c + 0] = (double(ix) + frac) * cell_width_;
          pos[3 * c + 1] = (double(iy) + frac) * cell_width_;
          pos[3 * c + 2] = (double(iz) + frac) * cell_width_;
        }
      }
    }
    return pos;
  }

  /**
   * @brief Write species-independent /Cells groups, attributes, and centres.
   *
   * Called only by rank zero during construction. Cell centres use the same
   * flattened order as cell_getid() and all per-species cell datasets.
   */
  void write_common_cell_metadata() const {

    // Write the cell look up table metadata
    HDFCreateGroup(fname_, "Cells");
    HDFCreateGroup(fname_, "Cells/Meta-data");
    HDFWriteGroupAttribute(fname_, "Cells/Meta-data", "nr_cells",
                           from_value<int>((int)ncells_));
    HDFWriteGroupAttribute(fname_, "Cells/Meta-data", "size",
                           std::vector<double>(3, cell_width_));
    HDFWriteGroupAttribute(fname_, "Cells/Meta-data", "dimension",
                           std::vector<int>(3, (int)cdim_));
    HDFWriteGroupAttribute(fname_, "Cells/Meta-data", "Origin",
                           std::vector<double>(3, 0.0));

    // Write the cell centres in cell_getid() order
    HDFWriteDatasetVector(fname_, "Cells/Centres", this->cell_positions(0.5));

    // Create the groups the per-species cell datasets are written into
    HDFCreateGroup(fname_, "Cells/Files");
    HDFCreateGroup(fname_, "Cells/OffsetsInFile");
    HDFCreateGroup(fname_, "Cells/Counts");
    HDFCreateGroup(fname_, "Cells/MinPositions");
    HDFCreateGroup(fname_, "Cells/MaxPositions");
  }

  /**
   * @brief Write counts, offsets, file indices and bounds for one species.
   *
   * Called only by rank zero. Files is zero for every cell because this
   * plugin always writes one shared output file. The bounds are the geometric
   * cell edges, which contain every particle since positions are wrapped
   * into the cell they are sorted into; swiftsimio uses them to mask without
   * padding the requested region.
   *
   * @param grp Particle group name used as each metadata dataset name.
   * @param global_counts Global particle count in each cell.
   * @param cell_offset Start of each cell in the species datasets.
   */
  void write_species_cell_metadata(
      const std::string &grp, const std::vector<swift_count_t> &global_counts,
      const std::vector<swift_count_t> &cell_offset) const {

    // A single shared file, so every cell lives in file 0
    HDFWriteDataset(fname_, "Cells/Files/" + grp, std::vector<int>(ncells_, 0));

    // Where each cell starts in the species datasets, and how many particles
    // it holds
    HDFWriteDataset(fname_, "Cells/OffsetsInFile/" + grp, cell_offset);
    HDFWriteDataset(fname_, "Cells/Counts/" + grp, global_counts);

    // The geometric bounds of each cell
    HDFWriteDatasetVector(fname_, "Cells/MinPositions/" + grp,
                          this->cell_positions(0.0));
    HDFWriteDatasetVector(fname_, "Cells/MaxPositions/" + grp,
                          this->cell_positions(1.0));
  }

  /**
   * @brief Split "key = value" lines into trimmed key-value pairs.
   *
   * Lines are split at their first '='. Blank lines, comments starting with
   * '#' and lines without '=' are skipped.
   *
   * @param in Stream to read the lines from.
   *
   * @return Key-value pairs in the order they appear.
   */
  static std::vector<std::pair<std::string, std::string>>
  parse_key_values(std::istream &in) {

    // Strip surrounding whitespace
    auto trim = [](const std::string &str) {
      const size_t first = str.find_first_not_of(" \t\r");
      const size_t last = str.find_last_not_of(" \t\r");
      return (first == std::string::npos) ? std::string()
                                          : str.substr(first, last - first + 1);
    };

    // Keep each "key = value" line, split at its first '='
    std::vector<std::pair<std::string, std::string>> pairs;
    std::string line;
    while (std::getline(in, line)) {

      // Skip blank lines, comments and lines that are not key-value pairs
      const size_t eq = line.find('=');
      if (line.empty() || line[0] == '#' || eq == std::string::npos) {
        continue;
      }

      pairs.emplace_back(trim(line.substr(0, eq)), trim(line.substr(eq + 1)));
    }

    return pairs;
  }

  /**
   * @brief Copy one of monofonIC's input tables into the ICs metadata.
   *
   * The table is read from the text file monofonIC writes next to the
   * configuration, and each column becomes a dataset named after its header
   * (the k column is named "k", with its units as a group attribute).
   *
   * @param grp Group to write the columns into.
   * @param filename Table file name, relative to the configuration.
   */
  void write_input_table(const std::string &grp,
                         const std::string &filename) const {

    // The table is written by the IC generator, so it should always exist
    std::ifstream in(cf_.get_path_relative_to_config(filename));
    if (!in) {
      music::wlog << "SWIFT : could not read " << filename
                  << " to store in the ICs metadata" << std::endl;
      return;
    }

    // The column header is written with a field width of 20 characters per
    // column, the first starting after "# "; data lines are whitespace
    // separated numbers
    std::vector<std::string> names;
    std::vector<std::vector<double>> columns;
    std::string line;
    while (std::getline(in, line)) {

      // The header line naming the columns starts with the k column
      if (line.rfind("#", 0) == 0) {
        if (line.find("k [h/Mpc]") != std::string::npos) {
          names = {"k"};
          for (size_t pos = 20; pos < line.size(); pos += 20) {
            std::istringstream field(line.substr(pos, 20));
            std::string name;
            field >> name;
            names.push_back(name);
          }
        }
        continue;
      }

      // Append each number to its column
      std::istringstream values(line);
      double value;
      size_t column = 0;
      while (values >> value) {
        if (columns.size() <= column) {
          columns.emplace_back();
        }
        columns[column++].push_back(value);
      }
    }

    // Write one dataset per column. The header can name more columns than
    // are written, so the first names are used; generic names are the
    // fallback if it names too few
    HDFCreateGroup(fname_, grp);
    HDFWriteGroupAttribute(fname_, grp, "k units", std::string("h/Mpc"));
    HDFWriteGroupAttribute(fname_, grp, "Scale-factor", (double)astart_);
    for (size_t c = 0; c < columns.size(); ++c) {
      const std::string name = (names.size() >= columns.size())
                                   ? names[c]
                                   : "column_" + std::to_string(c);
      HDFWriteDataset(fname_, grp + "/" + name, columns[c]);
    }
  }

  /**
   * @brief Write everything needed to reproduce these ICs to /ICs_parameters.
   *
   * SWIFT copies this group unchanged into every snapshot, so the run can be
   * reproduced from any SWIFT output. Called by rank zero from the
   * destructor, once the input tables have been written and every species'
   * particle mass is known. Alongside the summary attributes of the group
   * itself, it holds:
   *  - Config: every configuration option, one sub-group per section;
   *  - Cosmology: the cosmological parameters actually used, after defaults;
   *  - Build: compiler, libraries, CLASS version and parallel layout;
   *  - Derived: growth factors, sigma_8, particle masses and gas state;
   *  - CLASS: the CLASS input parameters, when CLASS is the transfer function;
   *  - InputPowerSpectrum and InputTransferFunctions: the input tables.
   */
  void write_ics_metadata() const {

    // Read the main choices from the configuration
    int order = cf_.get_value<int>("setup", "LPTorder");
    std::string load = cf_.get_value<std::string>("setup", "ParticleLoad");
    std::string tf = cf_.get_value<std::string>("cosmology", "transfer");
    std::string cosmo_set =
        cf_.get_value<std::string>("cosmology", "ParameterSet");
    std::string rng = cf_.get_value<std::string>("random", "generator");
    int do_fixing = cf_.get_value<bool>("setup", "DoFixing");
    int do_invert = cf_.get_value<bool>("setup", "DoInversion");
    int do_baryons = cf_.get_value<bool>("setup", "DoBaryons");
    int do_baryonsVrel = cf_.get_value<bool>("setup", "DoBaryonVrel");
    int L = cf_.get_value<int>("setup", "GridRes");

    // Write them as a summary, as attributes of the group itself
    HDFCreateGroup(fname_, "ICs_parameters");
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Code",
                           std::string("MUSIC2 - monofonIC"));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Git Revision",
                           std::string(GIT_REV));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Git Tag",
                           std::string(GIT_TAG));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Git Branch",
                           std::string(GIT_BRANCH));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Precision",
                           std::string(CMAKE_PRECISION_STR));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Convolutions",
                           std::string(CMAKE_CONVOLVER_STR));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "PLT",
                           std::string(CMAKE_PLT_STR));
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "LPT Order", order);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Particle Load", load);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Transfer Function", tf);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Cosmology Parameter Set",
                           cosmo_set);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Random Generator", rng);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Mode Fixing", do_fixing);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Mode inversion",
                           do_invert);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Baryons", do_baryons);
    HDFWriteGroupAttribute(fname_, "ICs_parameters",
                           "Baryons Relative Velocity", do_baryonsVrel);
    HDFWriteGroupAttribute(fname_, "ICs_parameters", "Grid Resolution", L);

    // Only CLASS back-scales from a target redshift
    if (tf == "CLASS") {
      double ztarget = cf_.get_value<double>("cosmology", "ztarget");
      HDFWriteGroupAttribute(fname_, "ICs_parameters", "Target Redshift",
                             ztarget);
    }

    // PANPHASIA's phases are set by its descriptor
    if (rng == "PANPHASIA") {
      std::string desc = cf_.get_value<std::string>("random", "descriptor");
      HDFWriteGroupAttribute(fname_, "ICs_parameters", "Descriptor", desc);
    }

    // Dump every configuration option exactly as read, as
    // "section/key = value" lines
    std::ostringstream dump;
    cf_.dump(dump);
    std::istringstream config_lines(dump.str());

    // Write each option into its section's sub-group, creating the sub-group
    // the first time the section appears
    HDFCreateGroup(fname_, "ICs_parameters/Config");
    std::set<std::string> sections;
    for (const auto &kv : parse_key_values(config_lines)) {

      // Split "section/key"; options without a section go into "general"
      const size_t slash = kv.first.find('/');
      const std::string section =
          (slash == std::string::npos) ? "general" : kv.first.substr(0, slash);
      const std::string key =
          (slash == std::string::npos) ? kv.first : kv.first.substr(slash + 1);
      const std::string section_grp = "ICs_parameters/Config/" + section;

      // Create the section's sub-group on first use
      if (sections.insert(section).second) {
        HDFCreateGroup(fname_, section_grp);
      }

      HDFWriteGroupAttribute(fname_, section_grp, key, kv.second);
    }

    // The cosmological parameters actually used, including defaults and
    // derived values
    HDFCreateGroup(fname_, "ICs_parameters/Cosmology");
    for (const auto &param : pcc_->cosmo_param_.get_all()) {
      HDFWriteGroupAttribute(fname_, "ICs_parameters/Cosmology", param.first,
                             param.second);
    }

    // How the code was built and run
    const std::string build = "ICs_parameters/Build";
    HDFCreateGroup(fname_, build);
    HDFWriteGroupAttribute(fname_, build, "Build Type",
                           std::string(CMAKE_BUILDTYPE_STR));
    HDFWriteGroupAttribute(fname_, build, "Compiler", std::string(__VERSION__));
    HDFWriteGroupAttribute(fname_, build, "FFTW Version",
                           std::string(FFTW_API(version)));

    // HDF5 reports its version as three numbers
    unsigned hdf5_major, hdf5_minor, hdf5_release;
    H5get_libversion(&hdf5_major, &hdf5_minor, &hdf5_release);
    HDFWriteGroupAttribute(fname_, build, "HDF5 Version",
                           std::to_string(hdf5_major) + "." +
                               std::to_string(hdf5_minor) + "." +
                               std::to_string(hdf5_release));
#if defined(USE_CLASS)
    // Only builds with CLASS know its version
    HDFWriteGroupAttribute(fname_, build, "CLASS Revision",
                           std::string(CLASS_GIT_REV));
    HDFWriteGroupAttribute(fname_, build, "CLASS Tag",
                           std::string(CLASS_GIT_TAG));
#endif

    // The parallel layout the ICs were generated with
    HDFWriteGroupAttribute(fname_, build, "MPI Ranks", num_ranks_);
    HDFWriteGroupAttribute(fname_, build, "Threads", CONFIG::num_threads);

    // Derived quantities, to check a reproduction against
    const std::string derived = "ICs_parameters/Derived";
    HDFCreateGroup(fname_, derived);
    HDFWriteGroupAttribute(fname_, derived, "D+ start", pcc_->Dplus_start_);
    HDFWriteGroupAttribute(fname_, derived, "D+ target", pcc_->Dplus_target_);
    HDFWriteGroupAttribute(fname_, derived, "sigma_8",
                           (double)pcc_->compute_sigma8());

    // Only species that were written have a particle mass
    for (int i = 0; i < 7; ++i) {
      if (npart_[i] > 0) {
        HDFWriteGroupAttribute(fname_, derived,
                               "Mean particle mass PartType" +
                                   std::to_string(i) + " [1e10 Msun]",
                               mean_particle_mass_[i]);
      }
    }

    // The initial gas state, when baryons are enabled
    if (bdobaryons_) {
      HDFWriteGroupAttribute(fname_, derived, "Gas temperature over mu [K]",
                             gas_temperature_over_mu_);
      HDFWriteGroupAttribute(fname_, derived, "Gas internal energy [(km/s)^2]",
                             internal_energy_);
      HDFWriteGroupAttribute(fname_, derived, "Gas smoothing length [Mpc]",
                             smoothing_length_);
    }

    // The CLASS input parameters, as written by the CLASS plugin
    if (tf == "CLASS") {
      std::ifstream ini(
          cf_.get_path_relative_to_config("input_class_parameters.ini"));
      HDFCreateGroup(fname_, "ICs_parameters/CLASS");
      for (const auto &kv : parse_key_values(ini)) {
        HDFWriteGroupAttribute(fname_, "ICs_parameters/CLASS", kv.first,
                               kv.second);
      }
    }

    // The input power spectra and transfer functions at the starting time
    this->write_input_table("ICs_parameters/InputPowerSpectrum",
                            "input_powerspec.txt");
    this->write_input_table("ICs_parameters/InputTransferFunctions",
                            "input_transfer.txt");
  }

  /**
   * @brief Log the cell grid, HDF5 backend and initial gas state.
   *
   * Called from the constructor on every rank; only rank 0 logs at info
   * level. Follows the "SWIFT : ..." sentence style of the other plugins.
   */
  void log_setup() const {

    // The top-level cell grid the particles are sorted into
    music::ilog.Print("SWIFT : sorting particles into %zu^3 top-level cells of "
                      "width %.4g Mpc",
                      cdim_, cell_width_);

    // How the IC file is written
#ifdef USE_PARALLEL_HDF5
    music::ilog.Print("SWIFT : writing with parallel HDF5 from %d task(s)",
                      num_ranks_);
#else
    music::ilog.Print("SWIFT : writing with serial HDF5, %d task(s) in turn",
                      num_ranks_);
#endif

    // Initial gas state
    if (bdobaryons_) {
      music::ilog.Print("SWIFT : set initial gas temperature to %.2f K/mu",
                        gas_temperature_over_mu_);
      music::ilog.Print("SWIFT : set initial internal energy to %.2e km^2/s^2",
                        internal_energy_);
      music::ilog.Print("SWIFT : set initial smoothing length to the grid "
                        "spacing: %.4g Mpc",
                        smoothing_length_);
    }
  }

  /**
   * @brief Start a timed stage in the same style as the LPT steps.
   *
   * @param task Name of the stage, printed after the task symbol.
   */
  static void log_task_start(const std::string &task) {
    music::ilog << colors::SYM_CHECK << " " << colors::TASK_NAME << task
                << colors::RESET << std::setw(77 - (int)task.size())
                << std::setfill('.') << std::left << "" << std::endl;
  }

  /**
   * @brief Finish a timed stage with its right-aligned "took" line.
   *
   * @param t_start Wall-clock time at which the stage started.
   */
  static void log_task_took(double t_start) {
    music::ilog << std::setw(70) << std::setfill(' ') << std::right
                << "took : " << std::setw(8) << get_wtime() - t_start << "s"
                << std::endl;
  }

  /**
   * @brief Log how particles are spread over the top-level cells and how many
   * hyperslabs the busiest task will write.
   *
   * Called on every rank because of the hyperslab reduction; only rank 0
   * logs at info level.
   *
   * @param n_global Global particle count for the species.
   * @param global_counts Global particle count in each cell.
   * @param n_ranges Number of hyperslabs this rank will write.
   */
  void log_cell_occupancy(size_t n_global,
                          const std::vector<swift_count_t> &global_counts,
                          size_t n_ranges) const {

    // Largest hyperslab count across ranks
    unsigned long long ranges_max = n_ranges;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &ranges_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                  MPI_COMM_WORLD);
#endif

    // Occupancy of the global cell lookup table
    const auto minmax =
        std::minmax_element(global_counts.begin(), global_counts.end());
    const size_t n_filled =
        ncells_ -
        (size_t)std::count(global_counts.begin(), global_counts.end(), 0);

    music::ilog.Print("SWIFT : %zu of %zu cells filled, %lld to %lld particles "
                      "per cell (mean %.3g)",
                      n_filled, ncells_, *minmax.first, *minmax.second,
                      double(n_global) / double(ncells_));
    music::ilog.Print("SWIFT : at most %llu hyperslab(s) per task", ranges_max);
  }

  /**
   * @brief Log the amount of particle data written and the write rate.
   *
   * @param s Species, used to include gas fields in the byte count.
   * @param n_global Global particle count for the species.
   * @param t_write Time taken to write the particle data in seconds.
   */
  void log_write_rate(const cosmo_species &s, size_t n_global,
                      double t_write) const {

    // Bytes written: positions, velocities, ID and mass, plus internal
    // energy and smoothing length for gas
    const size_t nreal = this->has_64bit_reals() ? 8 : 4;
    const size_t nid = this->has_64bit_ids() ? 8 : 4;
    size_t per_particle = 7 * nreal + nid;
    if (bdobaryons_ && s == cosmo_species::baryon) {
      per_particle += 2 * sizeof(write_real_t);
    }
    const double gbytes =
        double(n_global) * double(per_particle) / (1024. * 1024. * 1024.);

    music::ilog.Print("SWIFT : wrote %.3g GB at %.3g GB/s", gbytes,
                      t_write > 0. ? gbytes / t_write : 0.);
  }

public:
  /**
   * @brief Sort and write one species, then emit its /Cells metadata.
   *
   * Entry point called by the IC generator once per species. Particles are
   * sorted into top-level cells (sort_particles_into_cells), the global cell
   * lookup table is built (build_cell_lookup_table) and turned into this rank's
   * file ranges (build_write_ranges). Rank zero then creates the datasets,
   * every rank writes its particles (write_species_fields), and rank zero
   * records the species' /Cells metadata.
   *
   * @param pc Local particle arrays plus local and global particle counts.
   * @param s Species represented by @p pc.
   * @param Omega_species Cosmological density parameter used for constant mass.
   */
  void write_particle_data(const particle::container &pc,
                           const cosmo_species &s, double Omega_species) {

    // Get the SWIFT particle type and its group name
    const int sid = get_species_idx(s);
    assert(sid != -1);
    const std::string grp = std::string("PartType") + std::to_string(sid);

    // Get the local and global particle counts
    const size_t n_local = pc.get_local_num_particles();
    const size_t global_num_particles = pc.get_global_num_particles();

    // Record the global count for the header, split into the 32-bit words the
    // Gadget-style header expects
    npart_[sid] = global_num_particles;
    npartTotal_[sid] = (uint32_t)(global_num_particles);
    npartTotalHighWord_[sid] = (uint32_t)((global_num_particles) >> 32);

    // Constant particle mass, unused when the container has individual masses
    const double particle_mass =
        pc.bhas_individual_masses_
            ? 0.0
            : Omega_species * munit_ / double(global_num_particles);
    mean_particle_mass_[sid] =
        Omega_species * munit_ / double(global_num_particles);

    // Sort the local particles into SWIFT's top-level cells
    log_task_start("Sorting " + grp + " into top-level cells");
    double t_start = get_wtime();
    std::vector<size_t> perm;
    std::vector<swift_count_t> local_counts;
    if (this->has_64bit_reals()) {
      this->sort_particles_into_cells(pc.positions64_, n_local, perm,
                                      local_counts);
    } else {
      this->sort_particles_into_cells(pc.positions32_, n_local, perm,
                                      local_counts);
    }
    log_task_took(t_start);

    // Combine every rank's cells into the global cell lookup table and get
    // the file ranges this rank contributes
    log_task_start("Building the " + grp + " cell lookup table");
    t_start = get_wtime();
    std::vector<swift_count_t> global_counts, cell_offset, rank_offset;
    this->build_cell_lookup_table(local_counts, global_counts, cell_offset,
                                  rank_offset);
    const std::vector<swift_write_range> ranges =
        this->build_write_ranges(local_counts, rank_offset);
    this->log_cell_occupancy(global_num_particles, global_counts,
                             ranges.size());
    log_task_took(t_start);

    // Rank 0 creates the full, empty datasets in the file, then every rank
    // writes its own particles into their cells' slots
    log_task_start("Writing " + grp + " particle data");
    t_start = get_wtime();
    if (this_rank_ == 0) {
      this->create_species_datasets(grp, global_num_particles, s);
    }
    this->barrier();
    this->write_species_fields(grp, pc, s, perm, ranges, particle_mass);
    this->barrier();
    this->log_write_rate(s, global_num_particles, get_wtime() - t_start);
    log_task_took(t_start);

    // Rank 0 records the cell metadata for this species
    log_task_start("Writing " + grp + " cell metadata");
    t_start = get_wtime();
    if (this_rank_ == 0) {
      this->write_species_cell_metadata(grp, global_counts, cell_offset);
    }
    this->barrier();
    log_task_took(t_start);
  }
};

namespace {
// SWIFT expects double-precision IC data, so particle data are always written
// as doubles, whatever precision monofonIC was built with (CODE_PRECISION)
output_plugin_creator_concrete<swift_output_plugin<double>> creator301("SWIFT");
} // namespace

#endif
