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

#ifdef USE_CAMB

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>

#include <general.hh>
#include <config_file.hh>
#include <transfer_function_plugin.hh>
#include <ic_generator.hh>
#include <math/interpolate.hh>
#include <logger.hh>

#include <pybind11/embed.h>
#include <camb_interface.hh>

namespace py = pybind11;

class transfer_CAMB_plugin : public TransferFunction_plugin
{
protected:
    // Python interpreter guard (ensure single initialization)
    static py::scoped_interpreter* guard_;
    static int instance_count_;

    std::unique_ptr<CAMBInterface> camb_;

private:
    using TransferFunction_plugin::cosmo_params_;

    // Interpolation objects for transfer functions at target and z=0
    interpolated_function_1d<true, true, false> delta_c_, delta_b_, delta_m_;
    interpolated_function_1d<true, true, false> theta_c_, theta_b_, theta_m_;
    interpolated_function_1d<true, true, false> delta_c0_, delta_b0_, delta_m0_;
    interpolated_function_1d<true, true, false> theta_c0_, theta_b0_, theta_m0_;

    double zstart_, ztarget_, kmax_, kmin_, h_, tnorm_, f_b_, f_c_;
    std::ofstream ofs_camb_input_;

    /**
     * @brief Compute weighted total matter velocity from CDM and baryon velocities
     *
     * @param k Wavenumber array [h/Mpc]
     * @param v_cdm CDM velocity from CAMB
     * @param v_b Baryon velocity from CAMB
     * @param v_tot Output total velocity
     */
    void compute_total_velocity(
        const std::vector<double>& k,
        const std::vector<double>& v_cdm,
        const std::vector<double>& v_b,
        std::vector<double>& v_tot)
    {
        v_tot.resize(k.size());
        for (size_t i = 0; i < k.size(); ++i) {
            // Weighted average by density fractions (excluding neutrinos)
            v_tot[i] = f_c_ * v_cdm[i] + f_b_ * v_b[i];
        }
    }

public:
    explicit transfer_CAMB_plugin(config_file &cf, const cosmology::parameters& cosmo_params)
        : TransferFunction_plugin(cf, cosmo_params)
    {
        // Initialize Python interpreter once for all instances
        if (instance_count_ == 0) {
            guard_ = new py::scoped_interpreter();
        }
        instance_count_++;

        this->tf_isnormalised_ = true;

        // Open output file for CAMB parameters (for debugging)
        ofs_camb_input_.open(cf.get_path_relative_to_config("input_camb_parameters.ini"),
                            std::ios::trunc);

        // Get configuration parameters
        ztarget_ = pcf_->get_value_safe<double>("cosmology", "ztarget", 0.0);
        zstart_ = pcf_->get_value<double>("setup", "zstart");

        f_b_ = cosmo_params_["Omega_b"] / (cosmo_params_["Omega_b"] + cosmo_params_["Omega_c"]);
        f_c_ = cosmo_params_["Omega_c"] / (cosmo_params_["Omega_b"] + cosmo_params_["Omega_c"]);
        h_ = cosmo_params_["h"];

        music::ilog << "CAMB: Initializing Python CAMB interface..." << std::endl;

        // Create CAMB interface
        camb_ = std::make_unique<CAMBInterface>();

        // Set cosmological parameters
        double H0 = 100.0 * h_;
        double ombh2 = cosmo_params_["Omega_b"] * h_ * h_;
        double omch2 = cosmo_params_["Omega_c"] * h_ * h_;
        double omnuh2 = cosmo_params_["Omega_nu_massive"] * h_ * h_;
        double omk = 1.0 - cosmo_params_["Omega_m"] - cosmo_params_["Omega_DE"];

        double A_s = cosmo_params_["A_s"];
        double n_s = cosmo_params_["n_s"];
        double tau = 0.054;  // CAMB default for tau_reio (not in monofonIC params)

        double w0 = cosmo_params_["w_0"];
        double wa = cosmo_params_["w_a"];

        double Tcmb = cosmo_params_["Tcmb"];
        double Neff = cosmo_params_["N_ur"];

        // Log parameters
        ofs_camb_input_ << "# CAMB parameters for monofonIC run" << std::endl;
        ofs_camb_input_ << "H0 = " << H0 << std::endl;
        ofs_camb_input_ << "ombh2 = " << ombh2 << std::endl;
        ofs_camb_input_ << "omch2 = " << omch2 << std::endl;
        ofs_camb_input_ << "omnuh2 = " << omnuh2 << std::endl;
        ofs_camb_input_ << "omk = " << omk << std::endl;
        ofs_camb_input_ << "A_s = " << A_s << std::endl;
        ofs_camb_input_ << "n_s = " << n_s << std::endl;
        ofs_camb_input_ << "tau = " << tau << std::endl;
        ofs_camb_input_ << "w0 = " << w0 << std::endl;
        ofs_camb_input_ << "wa = " << wa << std::endl;
        ofs_camb_input_ << "Tcmb = " << Tcmb << std::endl;
        ofs_camb_input_ << "Neff = " << Neff << std::endl;

        if (A_s > 0.0) {
            music::ilog << "CAMB: Using A_s=" << colors::CONFIG_VALUE << A_s
                       << colors::RESET << " to normalise the transfer function." << std::endl;
        } else {
            throw std::runtime_error("CAMB plugin requires A_s to be specified (sigma8 not yet supported)");
        }

        // Set cosmology in CAMB
        camb_->set_cosmology(H0, ombh2, omch2, omnuh2, omk, A_s, n_s, tau, w0, wa, Tcmb, Neff);

        // Set accuracy
        double accuracy_boost = pcf_->get_value_safe<double>("cosmology", "camb_accuracy_boost", 1.0);
        camb_->set_accuracy(accuracy_boost);

        // Determine k range
        double lbox = pcf_->get_value<double>("setup", "BoxLength");
        int nres = pcf_->get_value<double>("setup", "GridRes");
        kmax_ = std::max(20.0, 2.0 * M_PI / lbox * nres / 2 * std::sqrt(3.0) * 2.0);
        double kmin = 0.001;  // Similar to CLASS

        ofs_camb_input_ << "kmin = " << kmin << std::endl;
        ofs_camb_input_ << "kmax = " << kmax_ << std::endl;
        ofs_camb_input_ << "ztarget = " << ztarget_ << std::endl;
        ofs_camb_input_ << "zstart = " << zstart_ << std::endl;

        // Set matter power parameters
        double zmax = std::max(zstart_, ztarget_);
        int nk = 1000;  // Number of k points
        camb_->set_matter_power(kmin, kmax_, nk, zmax);

        // Compute CAMB results
        music::ilog << "CAMB: Computing transfer functions..." << std::endl;
        camb_->compute();

        // Get normalization from CAMB
        double A_s_camb = camb_->get_As();
        double k_p = cosmo_params["k_p"] / cosmo_params["h"];
        tnorm_ = std::sqrt(2.0 * M_PI * M_PI * A_s_camb *
                          std::pow(1.0 / k_p, cosmo_params["n_s"] - 1) /
                          std::pow(2.0 * M_PI, 3.0));

        music::ilog << "CAMB: Extracting transfer functions at z=0..." << std::endl;

        std::vector<double> k;
        // Extract transfer functions at z=0
        std::vector<double> dc0, db0, dg0, dn0, dmn0, dm0, tc0, tb0, tm0;
        camb_->get_transfers(0.0, k, dc0, db0, dg0, dn0, dmn0, dm0, tc0, tb0);

        // Compute weighted total matter velocity
        compute_total_velocity(k, tc0, tb0, tm0);

        // Store in interpolation objects
        delta_c0_.set_data(k, dc0);
        delta_b0_.set_data(k, db0);
        delta_m0_.set_data(k, dm0);
        theta_c0_.set_data(k, tc0);
        theta_b0_.set_data(k, tb0);
        theta_m0_.set_data(k, tm0);

        music::ilog << "CAMB: Extracting transfer functions at z=" << ztarget_ << "..." << std::endl;

        // Extract transfer functions at ztarget
        std::vector<double> dc, db, dg, dn, dmn, dm, tc, tb, tm;
        camb_->get_transfers(ztarget_, k, dc, db, dg, dn, dmn, dm, tc, tb);

        // Compute weighted total matter velocity
        compute_total_velocity(k, tc, tb, tm);

        // Store in interpolation objects
        delta_c_.set_data(k, dc);
        delta_b_.set_data(k, db);
        delta_m_.set_data(k, dm);
        theta_c_.set_data(k, tc);
        theta_b_.set_data(k, tb);
        theta_m_.set_data(k, tm);

        kmin_ = k[0];
        kmax_ = k.back();

        music::ilog << "CAMB table contains k = " << colors::CONFIG_VALUE << this->get_kmin()
                   << colors::RESET << " to " << colors::CONFIG_VALUE << this->get_kmax()
                   << colors::RESET << " h Mpc-1." << std::endl;

        // Set plugin capabilities
        tf_distinct_ = true;
        tf_withvel_ = true;
        tf_withtotal0_ = true;

        music::ilog << "CAMB: Transfer function plugin initialized successfully" << std::endl;

        // Clear CAMB results to free memory (similar to CLASS)
        camb_.reset();
    }

    ~transfer_CAMB_plugin()
    {
        instance_count_--;
        if (instance_count_ == 0 && guard_ != nullptr) {
            delete guard_;
            guard_ = nullptr;
        }
    }

    inline double compute(double k, tf_type type) const
    {
        k *= h_;

        if (k < kmin_ || k > kmax_) {
            return 0.0;
        }

        real_t val(0.0);
        switch (type) {
            // Values at ztarget:
            case delta_matter:
                val = delta_m_(k); break;
            case delta_cdm:
                val = delta_c_(k); break;
            case delta_baryon:
                val = delta_b_(k); break;
            case theta_matter:
                val = theta_m_(k); break;
            case theta_cdm:
                val = theta_c_(k); break;
            case theta_baryon:
                val = theta_b_(k); break;
            case delta_bc:
                val = delta_b_(k) - delta_c_(k); break;
            case theta_bc:
                val = theta_b_(k) - theta_c_(k); break;

            // Values at z=0:
            case delta_matter0:
                val = delta_m0_(k); break;
            case delta_cdm0:
                val = delta_c0_(k); break;
            case delta_baryon0:
                val = delta_b0_(k); break;
            case theta_matter0:
                val = theta_m0_(k); break;
            case theta_cdm0:
                val = theta_c0_(k); break;
            case theta_baryon0:
                val = theta_b0_(k); break;
            default:
                throw std::runtime_error("Invalid type requested in transfer function evaluation");
        }
        return val * tnorm_;
    }

    inline double get_kmin(void) const { return kmin_ / h_; }
    inline double get_kmax(void) const { return kmax_ / h_; }
};

// Static member initialization
py::scoped_interpreter* transfer_CAMB_plugin::guard_ = nullptr;
int transfer_CAMB_plugin::instance_count_ = 0;

namespace
{
    TransferFunction_plugin_creator_concrete<transfer_CAMB_plugin> creator("CAMB");
}

#endif // USE_CAMB
