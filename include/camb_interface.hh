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

#ifdef USE_CAMB

#include <vector>
#include <string>
#include <memory>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/**
 * @brief C++ wrapper for Python CAMB interface
 *
 * This class manages the Python interpreter and provides methods to:
 * - Initialize CAMB with cosmological parameters
 * - Compute matter power spectra and transfer functions at arbitrary redshifts
 * - Extract transfer function data for interpolation
 */
class CAMBInterface {
public:
    /**
     * @brief Construct CAMB interface
     * @throws std::runtime_error if Python CAMB module cannot be imported
     */
    CAMBInterface();

    /**
     * @brief Destructor - ensures proper Python cleanup
     */
    ~CAMBInterface();

    /**
     * @brief Set cosmological parameters
     * @param H0 Hubble parameter at z=0 [km/s/Mpc]
     * @param ombh2 Omega_baryon * h^2
     * @param omch2 Omega_CDM * h^2
     * @param omnuh2 Omega_neutrino * h^2
     * @param omk Omega_curvature
     * @param As Primordial amplitude at k=0.05/Mpc
     * @param ns Primordial spectral index
     * @param tau Optical depth to reionization
     * @param w Dark energy equation of state (constant)
     * @param wa Dark energy equation of state evolution parameter
     * @param Tcmb CMB temperature today [K]
     * @param Neff Effective number of relativistic species
     */
    void set_cosmology(
        double H0, double ombh2, double omch2, double omnuh2, double omk,
        double As, double ns, double tau, double w, double wa,
        double Tcmb, double Neff);

    /**
     * @brief Set accuracy parameters for CAMB
     * @param accuracy_boost Boost factor for accuracy (default 1.0)
     * @param k_eta_max_scalar Max k*eta for scalar perturbations
     * @param lmax Maximum multipole for CMB
     */
    void set_accuracy(double accuracy_boost = 1.0, double k_eta_max_scalar = -1.0, int lmax = 2500);

    /**
     * @brief Set range of wavenumbers and redshifts for matter power
     * @param kmin Minimum k [h/Mpc]
     * @param kmax Maximum k [h/Mpc]
     * @param nk Number of k points (log-spaced)
     * @param zmax Maximum redshift
     */
    void set_matter_power(double kmin, double kmax, int nk, double zmax);

    /**
     * @brief Compute transfer functions and matter power
     * @throws std::runtime_error if CAMB computation fails
     */
    void compute();

    /**
     * @brief Get transfer function data at specified redshift
     * @param z Redshift
     * @param k Output: wavenumber array [h/Mpc]
     * @param T_cdm Output: CDM transfer function (dimensionless)
     * @param T_b Output: Baryon transfer function (dimensionless)
     * @param T_g Output: Photon transfer function (dimensionless)
     * @param T_nu Output: Massless neutrino transfer function (dimensionless)
     * @param T_mnu Output: Massive neutrino transfer function (dimensionless)
     * @param T_tot Output: Total matter transfer function (dimensionless)
     * @param v_cdm Output: CDM velocity transfer function [km/s]
     * @param v_b Output: Baryon velocity transfer function [km/s]
     *
     * Transfer functions are in synchronous gauge and follow CAMB conventions
     */
    void get_transfers(double z,
        std::vector<double>& k,
        std::vector<double>& T_cdm,
        std::vector<double>& T_b,
        std::vector<double>& T_g,
        std::vector<double>& T_nu,
        std::vector<double>& T_mnu,
        std::vector<double>& T_tot,
        std::vector<double>& v_cdm,
        std::vector<double>& v_b);

    /**
     * @brief Get CAMB version string
     * @return Version string
     */
    std::string get_version() const;

    /**
     * @brief Get primordial power spectrum normalization
     * @return Primordial amplitude As
     */
    double get_As() const { return As_; }

    /**
     * @brief Check if CAMB has been successfully computed
     * @return true if results are available
     */
    bool is_computed() const { return computed_; }

private:
    py::object camb_module_;      ///< Python CAMB module
    py::object camb_params_;      ///< CAMB parameters object
    py::object camb_results_;     ///< CAMB results object

    double As_;                   ///< Stored primordial amplitude
    bool computed_;               ///< Flag indicating if CAMB has been run

    /**
     * @brief Import CAMB module and check version
     */
    void import_camb();
};

#endif // USE_CAMB
