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

#include <camb_interface.hh>
#include <stdexcept>
#include <sstream>
#include <logger.hh>

CAMBInterface::CAMBInterface()
    : As_(0.0), computed_(false)
{
    import_camb();
}

CAMBInterface::~CAMBInterface()
{
    // pybind11 handles cleanup automatically via RAII
}

void CAMBInterface::import_camb()
{
    try {
        // Import CAMB module
        camb_module_ = py::module_::import("camb");

        // Get version
        std::string version = get_version();
        music::ilog << "CAMB Python interface initialized (version " << version << ")" << std::endl;

    } catch (const py::error_already_set& e) {
        std::stringstream ss;
        ss << "Failed to import CAMB Python module: " << e.what();
        throw std::runtime_error(ss.str());
    }
}

void CAMBInterface::set_cosmology(
    double H0, double ombh2, double omch2, double omnuh2, double omk,
    double As, double ns, double tau, double w, double wa,
    double Tcmb, double Neff)
{
    try {
        // Create CAMB parameter object
        camb_params_ = camb_module_.attr("CAMBparams")();

        // Set cosmological parameters
        // CAMB uses H0, ombh2, omch2, etc.
        // Note: omnuh2 is not a parameter in set_cosmology, neutrino mass is set via mnu
        // Convert omnuh2 to mnu: Omega_nu * h^2 = sum(m_nu) / 93.14 eV
        double mnu = (omnuh2 > 0.0) ? omnuh2 * 93.14 : 0.0;  // Total neutrino mass in eV
        int num_massive_nu = (mnu > 1e-6) ? 1 : 0;

        camb_params_.attr("set_cosmology")(
            py::arg("H0") = H0,
            py::arg("ombh2") = ombh2,
            py::arg("omch2") = omch2,
            py::arg("omk") = omk,
            py::arg("neutrino_hierarchy") = "degenerate",
            py::arg("num_massive_neutrinos") = num_massive_nu,
            py::arg("mnu") = mnu,
            py::arg("nnu") = Neff,
            py::arg("YHe") = py::none(),  // Use CAMB default
            py::arg("meffsterile") = 0.0,
            py::arg("standard_neutrino_neff") = 3.046,
            py::arg("TCMB") = Tcmb,
            py::arg("tau") = tau
        );

        // Set dark energy parameters
        camb_params_.attr("set_dark_energy")(
            py::arg("w") = w,
            py::arg("wa") = wa,
            py::arg("dark_energy_model") = "ppf"
        );

        // Set initial power spectrum (primordial)
        auto initial_power = camb_params_.attr("InitPower");
        initial_power.attr("set_params")(
            py::arg("As") = As,
            py::arg("ns") = ns,
            py::arg("r") = 0.0
        );

        // Store As for normalization
        As_ = As;

        music::ilog << "CAMB cosmology parameters set: H0=" << H0
                   << ", ombh2=" << ombh2 << ", omch2=" << omch2 << std::endl;

    } catch (const py::error_already_set& e) {
        std::stringstream ss;
        ss << "Failed to set CAMB cosmology: " << e.what();
        throw std::runtime_error(ss.str());
    }
}

void CAMBInterface::set_accuracy(double accuracy_boost, double k_eta_max_scalar, int lmax)
{
    try {
        // Set accuracy parameters
        auto accuracy = camb_params_.attr("Accuracy");
        accuracy.attr("AccuracyBoost") = accuracy_boost;

        if (k_eta_max_scalar > 0.0) {
            accuracy.attr("lSampleBoost") = k_eta_max_scalar;
        }

        // Set max_l for CMB calculations (not critical for transfer functions)
        camb_params_.attr("max_l") = lmax;

        music::dlog << "CAMB accuracy boost set to " << accuracy_boost << std::endl;

    } catch (const py::error_already_set& e) {
        std::stringstream ss;
        ss << "Failed to set CAMB accuracy: " << e.what();
        throw std::runtime_error(ss.str());
    }
}

void CAMBInterface::set_matter_power(double kmin, double kmax, int nk, double zmax)
{
    try {
        // Set redshifts for matter power spectrum
        // CAMB needs a list of redshifts where we want output
        py::list redshifts;
        redshifts.append(0.0);  // Always include z=0

        // Add zmax if different from 0
        if (zmax > 0.001) {
            redshifts.append(zmax);
        }

        // Set matter power parameters
        camb_params_.attr("set_matter_power")(
            py::arg("redshifts") = redshifts,
            py::arg("kmax") = kmax,
            py::arg("k_per_logint") = 0,  // Use linear spacing option
            py::arg("nonlinear") = false
        );

        // Set k range - CAMB uses kmax in set_matter_power
        // and we'll interpolate the outputs
        camb_params_.attr("NonLinear") = camb_module_.attr("model").attr("NonLinear_none");

        // Want transfer functions
        camb_params_.attr("WantTransfer") = true;
        camb_params_.attr("WantCls") = false;  // Don't need CMB power spectra

        // Set transfer function k range
        auto transfer = camb_params_.attr("Transfer");
        transfer.attr("high_precision") = true;
        transfer.attr("accurate_massive_neutrinos") = true;
        transfer.attr("kmax") = kmax;
        transfer.attr("k_per_logint") = 0;

        music::dlog << "CAMB matter power set: kmax=" << kmax << ", zmax=" << zmax << std::endl;

    } catch (const py::error_already_set& e) {
        std::stringstream ss;
        ss << "Failed to set CAMB matter power: " << e.what();
        throw std::runtime_error(ss.str());
    }
}

void CAMBInterface::compute()
{
    try {
        // Validate parameters
        camb_params_.attr("validate")();

        music::ilog << "Computing CAMB transfer functions..." << std::endl;

        // Compute results
        camb_results_ = camb_module_.attr("get_results")(camb_params_);

        computed_ = true;

        music::ilog << "CAMB computation completed successfully" << std::endl;

    } catch (const py::error_already_set& e) {
        std::stringstream ss;
        ss << "CAMB computation failed: " << e.what();
        throw std::runtime_error(ss.str());
    }
}

void CAMBInterface::get_transfers(double z,
    std::vector<double>& k,
    std::vector<double>& T_cdm,
    std::vector<double>& T_b,
    std::vector<double>& T_g,
    std::vector<double>& T_nu,
    std::vector<double>& T_mnu,
    std::vector<double>& T_tot,
    std::vector<double>& v_cdm,
    std::vector<double>& v_b)
{
    if (!computed_) {
        throw std::runtime_error("CAMB must be computed before extracting transfer functions");
    }

    try {
        // Get transfer function data
        auto transfer_data = camb_results_.attr("get_matter_transfer_data")();

        // Get k values (in h/Mpc) - these are stored in transfer_data.q
        auto k_array = transfer_data.attr("q");
        py::array_t<double> k_np = k_array.cast<py::array_t<double>>();
        auto k_buf = k_np.request();
        double* k_ptr = static_cast<double*>(k_buf.ptr);
        size_t nk = k_buf.shape[0];

        k.resize(nk);
        T_cdm.resize(nk);
        T_b.resize(nk);
        T_g.resize(nk);
        T_nu.resize(nk);
        T_mnu.resize(nk);
        T_tot.resize(nk);
        v_cdm.resize(nk);
        v_b.resize(nk);

        // Copy k values
        for (size_t i = 0; i < nk; ++i) {
            k[i] = k_ptr[i];
        }

        // Get the full transfer_data array (shape: [n_types, nk, n_redshifts])
        auto transfer_array = transfer_data.attr("transfer_data");

        // Import transfer type indices from CAMB
        auto model = py::module_::import("camb.model");
        int idx_cdm = model.attr("Transfer_cdm").cast<int>()-1;
        int idx_b = model.attr("Transfer_b").cast<int>()-1;
        int idx_g = model.attr("Transfer_g").cast<int>()-1;
        int idx_r = model.attr("Transfer_r").cast<int>()-1;
        int idx_nu = model.attr("Transfer_nu").cast<int>()-1;
        int idx_tot = model.attr("Transfer_tot").cast<int>()-1;
        int idx_v_cdm = model.attr("Transfer_Newt_vel_cdm").cast<int>()-1;
        int idx_v_b = model.attr("Transfer_Newt_vel_baryon").cast<int>()-1;

        music::ilog << "CAMB transfer indices: cdm=" << idx_cdm << ", b=" << idx_b
                   << ", nu=" << idx_nu << ", tot=" << idx_tot << std::endl;

        // Redshift index (we only have one redshift)
        int iz = 0;

        // Extract transfer functions using Python's array indexing
        // This avoids any buffer/stride issues
        for (size_t i = 0; i < nk; ++i) {
            T_cdm[i] = transfer_array.attr("__getitem__")(py::make_tuple(idx_cdm, i, iz)).cast<double>();
            T_b[i] = transfer_array.attr("__getitem__")(py::make_tuple(idx_b, i, iz)).cast<double>();
            T_g[i] = transfer_array.attr("__getitem__")(py::make_tuple(idx_g, i, iz)).cast<double>();
            T_nu[i] = transfer_array.attr("__getitem__")(py::make_tuple(idx_nu, i, iz)).cast<double>();
            T_mnu[i] = T_nu[i];  // Use same for now (CAMB combines them)
            T_tot[i] = transfer_array.attr("__getitem__")(py::make_tuple(idx_tot, i, iz)).cast<double>();

            // Velocities
            v_cdm[i] = transfer_array.attr("__getitem__")(py::make_tuple(idx_v_cdm, i, iz)).cast<double>();
            v_b[i] = transfer_array.attr("__getitem__")(py::make_tuple(idx_v_b, i, iz)).cast<double>();
        }

        // Debug: print a few sample values at low k
        if (nk > 0) {
            music::ilog << "Sample values at k[0]=" << k[0] << ": T_cdm=" << T_cdm[0]
                       << ", T_b=" << T_b[0] << ", T_nu=" << T_nu[0] << std::endl;
        }

        music::dlog << "Extracted " << nk << " transfer function points at z=" << z << std::endl;

    } catch (const py::error_already_set& e) {
        std::stringstream ss;
        ss << "Failed to extract CAMB transfer functions: " << e.what();
        throw std::runtime_error(ss.str());
    }
}

std::string CAMBInterface::get_version() const
{
    try {
        auto version = camb_module_.attr("__version__");
        return version.cast<std::string>();
    } catch (...) {
        return "unknown";
    }
}

#endif // USE_CAMB
