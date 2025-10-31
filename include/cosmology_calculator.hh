// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
//   Thomas Montandon (toma): 2nd and 3rd author growth functions
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

#include <array>
#include <vec.hh>

#include <cosmology_parameters.hh>
#include <physical_constants.hh>
#include <transfer_function_plugin.hh>
#include <math/ode_integrate.hh>
#include <logger.hh>

#include <math/interpolate.hh>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace cosmology
{
//! structure for ODE parameters required to evaluate RHS of ODEs
struct ODEParams {
    double Omega_r;
    double Omega_m;
    double Omega_DE;
    double H0;
    double w0;
    double wa;
};

/*!
 * @brief ODE system for growth factors to third order
 *
 * This function defines the system of ODEs to be solved to obtain the growth factors
 * up to third order in perturbation theory for a single fluid in an expanding universe.
 * The cosmological model is specified through the parameters Omega_m, w0, wa, and H0.
 *
 * @param t The independent variable (conformal time)
 * @param y The dependent variable array containing [a, D, D', E, E', F_a, F_a', F_b, F_b', F_c]
 * @param dydt The output array for the derivatives
 * @param params Pointer to ODEParams structure containing cosmological parameters
 * @return GSL_SUCCESS on successful evaluation
 */
inline int ode_system(double t, const double y[], double dydt[], void *params) {
    ODEParams *p = static_cast<ODEParams *>(params);
    double Omega_r  = p->Omega_r;
    double Omega_m  = p->Omega_m;
    double Omega_DE = p->Omega_DE;
    
    double w0 = p->w0;
    double wa = p->wa;
    double H0 = p->H0;

    double a = y[0];
    double D = y[1];
    double Dprime = y[2];
    double E = y[3];
    double Eprime = y[4];
    double Fa = y[5];
    double Faprime = y[6];
    double Fb = y[7];
    double Fbprime = y[8];
    // double Fc = y[9];

    double H_of_a = H0 * std::sqrt(Omega_r / (a*a*a*a) + Omega_m / (a * a * a) + Omega_DE * std::pow(a, -3 * (1 + w0 + wa)) * std::exp(-3 * wa * (1 - a)));

    dydt[0] = a * a * H_of_a;
    dydt[1] = Dprime;
    dydt[2] = -a * H_of_a * Dprime + 3.0 / 2.0 * Omega_m * H0 * H0 * D / a;
    dydt[3] = Eprime;
    dydt[4] = -a * H_of_a * Eprime + 3.0 / 2.0 * Omega_m * H0 * H0 * (E - D * D) / a;
    dydt[5] = Faprime;
    dydt[6] = -a * H_of_a * Faprime + 3.0 / 2.0 * Omega_m * H0 * H0 * (Fa - 2.0 * D * D * D) / a;
    dydt[7] = Fbprime;
    dydt[8] = -a * H_of_a * Fbprime + 3.0 / 2.0 * Omega_m * H0 * H0 * (Fb + 2.0 * D * D * D - 2.0 * E * D) / a;
    dydt[9] = D * Eprime - E * Dprime;

    return GSL_SUCCESS;
}

/*!
 * @class cosmology::calculator
 * @brief provides functions to compute cosmological quantities
 *
 * This class provides member functions to compute cosmological quantities
 * related to the Friedmann equations and linear perturbation theory, it also
 * provides the functionality to work with back-scaled cosmological fields
 */
class calculator
{
public:
    //! data structure to store cosmological parameters
    cosmology::parameters cosmo_param_;

    //! pointer to an instance of a transfer function plugin
    std::unique_ptr<TransferFunction_plugin> transfer_function_;
double Dnow_, Dplus_start_, Dplus_target_, astart_, atarget_;

private:
    static constexpr double REL_PRECISION = 1e-10;
    
    double m_n_s_, m_sqrtpnorm_, tnorm_;

    //! interpolation functions of growth functions
    interpolated_function_1d<true, true, false> D_of_a_, dotD_of_a_, f_of_a_, a_of_D_, E_of_a_, dotE_of_a_, Fa_of_a_, Fb_of_a_, dotFa_of_a_, dotFb_of_a_, Fc_of_a_, dotFc_of_a_; // toma


    //! wrapper for GSL adaptive integration routine, do not use if many integrations need to be done as it allocates and deallocates memory
    //! set to 61-point Gauss-Kronrod and large workspace, used for sigma_8 normalisation
    real_t integrate(double (*func)(double x, void *params), double a, double b, void *params) const
    {
        constexpr size_t wspace_size{100000};

        double result{0.0};
        double error{0.0};

        gsl_function F;
        F.function = func;
        F.params = params;
        
        auto errh = gsl_set_error_handler_off();
        gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(wspace_size);
        gsl_integration_qag(&F, a, b, 0, REL_PRECISION, wspace_size, GSL_INTEG_GAUSS61, wspace, &result, &error);
        gsl_integration_workspace_free(wspace);
        gsl_set_error_handler(errh);

        if (error / result > REL_PRECISION)
            music::wlog << "no convergence in function 'integrate', rel. error=" << error / result << std::endl;

        return static_cast<real_t>(result);
    }

    /// @brief compute the perturbation theory growth factors up to third order D(a) E(a) F(a) by solving the single fluid ODE, returns tables
    /// @param tab_a (out) table of scale factors
    /// @param tab_D (out) table of linear growth factors D(a)
    /// @param tab_f (out) table of linear growth rates f(a)
    /// @param tab_E (out) table of second order growth factors E(a)
    /// @param tab_dotE (out) table of second order growth rates E'(a)
    /// @param tab_Fa (out) table of third order growth factors Fa(a)
    /// @param tab_dotFa (out) table of third order growth rates Fa'(a)
    /// @param tab_Fb (out) table of third order growth factors Fb(a)
    /// @param tab_dotFb (out) table of third order growth rates Fb'(a)
    /// @param tab_Fc (out) table of third order growth factors Fc(a)
    /// @param tab_dotFc (out) table of third order growth rates Fc'(a)
    void compute_growth(std::vector<double> &tab_a, std::vector<double> &tab_D, std::vector<double> &tab_dotD, std::vector<double> &tab_f,
                        std::vector<double> &tab_E, std::vector<double> &tab_dotE,
                        std::vector<double> &tab_Fa, std::vector<double> &tab_dotFa,
                        std::vector<double> &tab_Fb, std::vector<double> &tab_dotFb,
                        std::vector<double> &tab_Fc, std::vector<double> &tab_dotFc)
    {
        // using v_t = vec_t<10, double>;
        double Omega_m = cosmo_param_["Omega_m"];
        double Omega_r = cosmo_param_["Omega_r"];
        double Omega_DE = cosmo_param_["Omega_DE"];
        double H0 = cosmo_param_["H0"];
        double w0 = cosmo_param_["w_0"];
        double wa = cosmo_param_["w_a"];
        

        // set ICs, very deep in radiation domination
        const double a0 = 1e-6;
        
        // double H_a0 = H0 * std::sqrt(Omega_m / (a0 * a0 * a0) + (1 - Omega_m));
        double H_a0 = H_of_a(a0);
        const double t0 = 1.0 / (a0 * H_of_a(a0));

        // first order growth
        const double D0 = a0;
        const double Dprime0 = a0 * a0 * H_a0;
        
        // second order growth
        const double E0 = -3.0 / 7.0 * a0 * a0;
        const double Eprime0 = -6.0 / 7.0 * std::pow(a0,3) * H_a0;

        // third order growth
        const double Fa0 = -1.0 / 3.0 * a0 * a0 * a0;
        const double Faprime0 = -2.0 * std::pow(a0, 5.0 / 2.0);

        const double Fb0 = 10. / 21 * a0 * a0 * a0;
        const double Fbprime0 = 60. / 21 * std::pow(a0, 5.0 / 2.0);

        const double Fc0 = -1.0 / 7.0 * a0 * a0 * a0;

        double y[10] = {a0, D0, Dprime0, E0, Eprime0, Fa0, Faprime0, Fb0, Fbprime0, Fc0};

        const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
        gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 10);
        gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-8, 0.0);
        gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (10);
        // gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, T, 1e-6, 0.0, 0.0);

        ODEParams params = {Omega_r, Omega_m, Omega_DE, H0, w0, wa};
        gsl_odeiv2_system sys = {ode_system, nullptr, 10, &params};

        double h = 1e-6;
        double t = t0;
        double t1 = 0.1;
        double amax = 1.1;

        while (y[0] < amax) {
            
            int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);

            if (status != GSL_SUCCESS) {
                std::cerr << "Error: " << gsl_strerror(status) << std::endl;
                break;
            }

            tab_a.push_back(y[0]);
            tab_D.push_back(y[1]);
            tab_f.push_back(y[2]); // temporarily store D' in table
            tab_dotD.push_back(y[2]);
            tab_E.push_back(y[3]);
            tab_dotE.push_back(y[4]);
            tab_Fa.push_back(y[5]);
            tab_dotFa.push_back(y[6]);
            tab_Fb.push_back(y[7]);
            tab_dotFb.push_back(y[8]);
            tab_Fc.push_back(y[9]);
            tab_dotFc.push_back(y[9]); // temporarily store un-important random value
        }

        gsl_odeiv2_evolve_free (e);
        gsl_odeiv2_control_free (c);
        gsl_odeiv2_step_free (s);
        
        if (CONFIG::MPI_task_rank == 0)
        {
            // output growth factors to file for debugging
            std::ofstream output_file;
            output_file.open("GrowthFactors.txt");
            // compute f, before we stored here D'
            output_file << "#"
                        << "a"
                        << " "
                        << "D"
                        << " "
                        << "f"
                        << " "
                        << "E"
                        << " "
                        << "dotE"
                        << " "
                        << "Fa"
                        << " "
                        << "dotFa"
                        << " "
                        << "Fb"
                        << " "
                        << "dotFb"
                        << " "
                        << "Fc"
                        << " "
                        << "dotFc"
                        << "\n";
            for (size_t i = 0; i < tab_a.size(); ++i)
            {

                tab_dotFc[i] = (tab_D[i] * tab_dotE[i] - tab_E[i] * tab_f[i]); // toma
                tab_f[i] = tab_f[i] / (tab_a[i] * H_of_a(tab_a[i]) * tab_D[i]);

                // toma
                output_file << tab_a[i] << " " << tab_D[i] << " " << tab_dotD[i] << " " << tab_E[i] << " " << tab_dotE[i] << " " << tab_Fa[i] << " " << tab_dotFa[i] << " " << tab_Fb[i] << " " << tab_dotFb[i] << " " << tab_Fc[i] << " " << tab_dotFc[i] << "\n";
            }
            output_file.close();
        }
    }

public:
    
    //! default constructor [deleted]
    calculator() = delete;
    
    //! copy constructor [deleted]
    calculator(const calculator& c) = delete;

    //! constructor for a cosmology calculator object
    /*!
	 * @param acosmo a cosmological parameters structure
	 * @param pTransferFunction pointer to an instance of a transfer function object
	 */
    explicit calculator(config_file &cf)
        : cosmo_param_(cf), astart_( 1.0/(1.0+cf.get_value<double>("setup","zstart")) ),
            atarget_( 1.0/(1.0+cf.get_value_safe<double>("cosmology","ztarget",0.0)) )
    {
        // pre-compute growth factors and store for interpolation
        std::vector<double> tab_a, tab_D, tab_dotD, tab_f, tab_E, tab_dotE, tab_Fa, tab_dotFa, tab_Fb, tab_dotFb, tab_Fc, tab_dotFc;   // toma
        this->compute_growth(tab_a, tab_D, tab_dotD, tab_f, tab_E, tab_dotE, tab_Fa, tab_dotFa, tab_Fb, tab_dotFb, tab_Fc, tab_dotFc); // toma
        D_of_a_.set_data(tab_a, tab_D);
        dotD_of_a_.set_data(tab_a, tab_dotD);
        f_of_a_.set_data(tab_a, tab_f);
        a_of_D_.set_data(tab_D, tab_a);

        // toma
        E_of_a_.set_data(tab_a, tab_E);
        dotE_of_a_.set_data(tab_a, tab_dotE);

        Fa_of_a_.set_data(tab_a, tab_Fa);
        Fb_of_a_.set_data(tab_a, tab_Fb);
        Fc_of_a_.set_data(tab_a, tab_Fc);
        dotFa_of_a_.set_data(tab_a, tab_dotFa);
        dotFb_of_a_.set_data(tab_a, tab_dotFb);
        dotFc_of_a_.set_data(tab_a, tab_dotFc);

        Dnow_ = D_of_a_(1.0);

        Dplus_start_ = D_of_a_(astart_) / Dnow_;
        Dplus_target_ = D_of_a_(atarget_) / Dnow_;

        music::ilog << "Linear growth factors: D+_target = " << colors::CONFIG_VALUE << Dplus_target_ << colors::RESET << ", D+_start = " << colors::CONFIG_VALUE << Dplus_start_ << colors::RESET << std::endl;

        // set up transfer functions and compute normalisation
        transfer_function_ = std::move(select_TransferFunction_plugin(cf, cosmo_param_));
        transfer_function_->intialise();
        if( !transfer_function_->tf_isnormalised_ ){
            cosmo_param_.set("pnorm", this->compute_pnorm_from_sigma8() );
        }else{
            cosmo_param_.set("pnorm", 1.0/Dplus_target_/Dplus_target_);
            auto sigma8 = this->compute_sigma8();
            music::ilog << "Measured sigma_8 for given PS normalisation is " << colors::CONFIG_VALUE << sigma8 << colors::RESET << std::endl;
        }
        cosmo_param_.set("sqrtpnorm", std::sqrt(cosmo_param_["pnorm"]));

        music::ilog << std::setw(32) << std::left << "TF supports distinct CDM+baryons"
                    << " : " << colors::CONFIG_VALUE << (transfer_function_->tf_is_distinct() ? "yes" : "no") << colors::RESET << std::endl;
        music::ilog << std::setw(32) << std::left << "TF minimum wave number"
                    << " : " << colors::CONFIG_VALUE << transfer_function_->get_kmin() << colors::RESET << " h/Mpc" << std::endl;
        music::ilog << std::setw(32) << std::left << "TF maximum wave number"
                    << " : " << colors::CONFIG_VALUE << transfer_function_->get_kmax() << colors::RESET << " h/Mpc" << std::endl;
        if( std::sqrt(3.0)* 2.0*M_PI / cf.get_value<double>("setup","BoxLength") * cf.get_value<double>("setup","GridRes")/2  >= transfer_function_->get_kmax() ){
            music::elog << "Simulation nyquist mode kny = " << std::sqrt(3.0)* 2.0*M_PI / cf.get_value<double>("setup","BoxLength") * cf.get_value<double>("setup","GridRes")/2 << " h/Mpc is beyond valid range of transfer function!" << std::endl;
        }

        m_n_s_ = cosmo_param_["n_s"];
        m_sqrtpnorm_ = cosmo_param_["sqrtpnorm"];
        
        double k_p = cosmo_param_["k_p"] / cosmo_param_["h"];
        tnorm_ = std::sqrt(2.0 * M_PI * M_PI * cosmo_param_["A_s"] * std::pow(1.0 / k_p, cosmo_param_["n_s"] - 1) / std::pow(2.0 * M_PI, 3.0));
    }

    //! destructor
    ~calculator() { }

    /// @brief Write out a correctly scaled power spectrum at time a
    /// @param a scale factor
    /// @param fname file name
    void write_powerspectrum(real_t a, std::string fname) const
    {
        // const real_t Dplus0 = this->get_growth_factor(a);

        if (CONFIG::MPI_task_rank == 0)
        {
            double kmin = std::max(1e-4, transfer_function_->get_kmin());

            // write power spectrum to a file
            std::ofstream ofs(fname.c_str());
            std::stringstream ss;
            ss << " ,ap=" << a << "";
            ofs << "# " << std::setw(18) << "k [h/Mpc]"
                << std::setw(20) << ("P_dtot(k,a=ap)")
                << std::setw(20) << ("P_dcdm(k,a=ap)")
                << std::setw(20) << ("P_dbar(k,a=ap)")
                << std::setw(20) << ("P_tcdm(k,a=ap)")
                << std::setw(20) << ("P_tbar(k,a=ap)")
                << std::setw(20) << ("P_dtot(k,a=1)")
                << std::setw(20) << ("P_dcdm(k,a=1)")
                << std::setw(20) << ("P_dbar(k,a=1)")
                << std::setw(20) << ("P_tcdm(k,a=1)")
                << std::setw(20) << ("P_tbar(k,a=1)")
                << std::endl;
            for (double k = kmin; k < transfer_function_->get_kmax(); k *= 1.01)
            {
                ofs << std::setw(20) << std::setprecision(10) << k
                    << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_matter)*Dplus_start_, 2.0)
                    << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_cdm)*Dplus_start_, 2.0)
                    << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_baryon)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_matter)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_cdm)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_baryon)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_cdm)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_baryon)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_matter0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_cdm0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_baryon0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_cdm0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_baryon0)* Dplus_start_ / Dplus_target_, 2.0)
                    << std::endl;
            }
            if (ofs.fail()) {
                std::string message = "Could not write to file: " + fname;
                music::elog << message << std::endl;
                throw std::runtime_error(message);
            }
            ofs.close();
        }
        music::ilog << "Wrote power spectrum at a=" << a << " to file \'" << fname << "\'" << std::endl;
    }

    /// @brief Write out a correctly scaled transfer function at time a
    /// @param[in] fname filename to write to
    void write_transfer( std::string fname ) const
    {
        // const real_t Dplus0 = this->get_growth_factor(a);

        if (CONFIG::MPI_task_rank == 0)
        {
            double kmin = std::max(1e-4, transfer_function_->get_kmin());

            // write power spectrum to a file
            std::ofstream ofs(fname.c_str());
            std::stringstream ss;
            ss << " ,ap=" << astart_ << "";
            ofs << "# " << std::setw(18) << "k [h/Mpc]"
                << std::setw(20) << ("delta_c(k,a=ap)")
                << std::setw(20) << ("delta_b(k,a=ap)")
                << std::setw(20) << ("delta_m(k,a=ap)")
                << std::setw(20) << ("delta_bc(k,a=ap)")
                << std::setw(20) << ("theta_c(k,a=ap)")
                << std::setw(20) << ("theta_b(k,a=ap)")
                << std::setw(20) << ("theta_m(k,a=ap)")
                << std::setw(20) << ("theta_bc(k,a=ap)")
                << std::endl;
            double fb = cosmo_param_["f_b"], fc = cosmo_param_["f_c"];
            for (double k = kmin; k < transfer_function_->get_kmax(); k *= 1.01)
            {
                const double dm  = this->get_amplitude(k, delta_matter) * Dplus_start_ / Dplus_target_;
                const double dbc = this->get_amplitude(k, delta_bc);
                const double db  = dm + fc * dbc;
                const double dc  = dm - fb * dbc;
                const double tm  = this->get_amplitude(k, delta_matter) * Dplus_start_ / Dplus_target_;
                const double tbc = this->get_amplitude(k, theta_bc);
                const double tb  = dm + fc * dbc;
                const double tc  = dm - fb * dbc;
                
                ofs << std::setw(20) << std::setprecision(10) << k
                    << std::setw(20) << std::setprecision(10) << dc
                    << std::setw(20) << std::setprecision(10) << db
                    << std::setw(20) << std::setprecision(10) << dm
                    << std::setw(20) << std::setprecision(10) << dbc + 2 * tbc * (std::sqrt( Dplus_target_ / Dplus_start_ ) - 1.0)
                    << std::setw(20) << std::setprecision(10) << tc / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::setw(20) << std::setprecision(10) << tb / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::setw(20) << std::setprecision(10) << tm / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::setw(20) << std::setprecision(10) << tbc / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::endl;
            }
            if (ofs.fail()) {
                std::string message = "Could not write to file: " + fname;
                music::elog << message << std::endl;
                throw std::runtime_error(message);
            }
            ofs.close();
        }
        music::ilog << "Wrote input transfer functions at a=" << astart_ << " to file \'" << fname << "\'" << std::endl;
    }

    /// @brief return the cosmological parameter object
    /// @return cosmological parameter object
    const cosmology::parameters &get_parameters(void) const noexcept
    {
        return cosmo_param_;
    }

    /// @brief return the value of the Hubble function H(a) = dloga/dt 
    /// @param[in] a scale factor
    /// @return H(a)
    inline double H_of_a(double a) const noexcept
    {
        double HH2 = 0.0;
        HH2 += cosmo_param_["Omega_r"] / (a * a * a * a);
        HH2 += cosmo_param_["Omega_m"] / (a * a * a);
        HH2 += cosmo_param_["Omega_k"] / (a * a);
        HH2 += cosmo_param_["Omega_DE"] * std::pow(a, -3. * (1. + cosmo_param_["w_0"] + cosmo_param_["w_a"])) * exp(-3. * (1.0 - a) * cosmo_param_["w_a"]);
        return cosmo_param_["H0"] * std::sqrt(HH2);
    }

    /// @brief Computes the linear theory growth factor D+, normalised to D+(a=1)=1
    /// @param[in] a scale factor
    /// @return D+(a)
    real_t get_growth_factor(real_t a) const noexcept
    {
        return D_of_a_(a) / Dnow_;
    }

    /// @brief Computes the inverse of get_growth_factor, i.e. a(D+)
    /// @param[in] Dplus growth factor
    /// @return a(D+)
    real_t get_a( real_t Dplus ) const noexcept
    {
        return a_of_D_( Dplus * Dnow_ );
    }

    //! Computes the linear theory growth rate f
    /*! Function computes (by interpolating on precalculated table)
        *   f = dlog D+ / dlog a
        */
    real_t get_f(real_t a) const noexcept
    {
        return f_of_a_(a);
    }

    //! Computes the time derivative w.r.t. conformal time of the linear growth factor D+
    real_t get_dotD(real_t a) const noexcept
    {
        return dotD_of_a_(a) / Dnow_;
    }

    // toma
    //! Computes the second-order growth factor E,
    //! see definition in 2205.11347 (eq. 28) (or 1602.05933 eq. 2.21)
    real_t get_2growth_factor(real_t a) const noexcept
    {
        return E_of_a_(a) / Dnow_ / Dnow_;
    }

    // toma
    //! Computes the time derivative w.r.t. conformal time of the second-order growth factor E
    //! see definition in 2205.11347 (eq. 28) (or 1602.05933 eq. 2.21)
    real_t get_dotE(real_t a) const noexcept
    {
        return dotE_of_a_(a) / Dnow_ / Dnow_ ;
    }

    //! Computes the third-order growth factor Fa,
    //! see definition in 2205.11347
    real_t get_3growthA_factor(real_t a) const noexcept
    {
        return Fa_of_a_(a) / Dnow_ / Dnow_ / Dnow_;
    }

    //! Computes the third-order growth factor Fb,
    //! see definition in 2205.11347
    real_t get_3growthB_factor(real_t a) const noexcept
    {
        return Fb_of_a_(a) / Dnow_ / Dnow_ / Dnow_;
    }

    //! Computes the third-order growth factor Fc,
    //! see definition in 2205.11347
    real_t get_3growthC_factor(real_t a) const noexcept
    {
        return Fc_of_a_(a) / Dnow_ / Dnow_ / Dnow_;
    }

    //! Computes the conformal time derivative of the third-order growth factor Fa,
    //! see definition in 2205.11347
    real_t get_dotFa(real_t a) const noexcept
    {
        return dotFa_of_a_(a) / Dnow_ / Dnow_ / Dnow_;
    }

    //! computes the conformal time derivative of the third-order growth factor Fb,
    //! see definition in 2205.11347
    real_t get_dotFb(real_t a) const noexcept
    {
        return dotFb_of_a_(a) / Dnow_ / Dnow_ / Dnow_ ;
    }

    //! computes the conformal time derivative of the third-order growth factor Fc,
    //! see definition in 2205.11347
    real_t get_dotFc(real_t a) const noexcept
    {
        return dotFc_of_a_(a) / Dnow_ / Dnow_ / Dnow_;
    }

    //! Compute the factor relating particle displacement and velocity
    /*! Function computes
        *  vfac = a * (H(a)/h) * dlogD+ / dlog a
        */
    real_t get_vfacD(real_t a) const noexcept
    {
        return get_dotD(a) /  get_growth_factor(a) / cosmo_param_["h"];
    }

    //! Compute the factor relating second-order particle displacement and velocity
    real_t get_vfacE(real_t a) const noexcept
    {
        return get_dotE(a) /  get_2growth_factor(a) / cosmo_param_["h"];
    }

    //! Compute the factor relating third-order (a) particle displacement and velocity
    real_t get_vfacFa(real_t a) const noexcept
    {
        return get_dotFa(a) /  get_3growthA_factor(a) / cosmo_param_["h"];
    }

    //! Compute the factor relating third-order (b) particle displacement and velocity
    real_t get_vfacFb(real_t a) const noexcept
    {
        return get_dotFb(a) /  get_3growthB_factor(a) / cosmo_param_["h"];
    }

    //! Compute the factor relating third-order (c) particle displacement and velocity
    real_t get_vfacFc(real_t a) const noexcept
    {
        return get_dotFc(a) /  get_3growthC_factor(a) / cosmo_param_["h"];
    }

    //! Integrand for the sigma_8 normalization of the power spectrum
    /*! Returns the value of the primordial power spectrum multiplied with 
     the transfer function and the window function of 8 Mpc/h at wave number k */
    static double dSigma8(double k, void *pParams)
    {
        cosmology::calculator *pcc = reinterpret_cast<cosmology::calculator *>(pParams);

        const double x = k * 8.0;
        const double w = (x < 0.001)? 1.0-0.1*x*x : 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
            
        static double nspect = (double)pcc->cosmo_param_["n_s"];
        double tf = pcc->transfer_function_->compute(k, delta_matter);

        //... no growth factor since we compute at z=0 and normalize so that D+(z=0)=1
        return k * k * w * w * pow((double)k, (double)nspect) * tf * tf;
    }

    //! Integrand for the sigma_8 normalization of the power spectrum
    /*! Returns the value of the primordial power spectrum multiplied with 
	 the transfer function and the window function of 8 Mpc/h at wave number k */
    static double dSigma8_0(double k, void *pParams)
    {
        cosmology::calculator *pcc = reinterpret_cast<cosmology::calculator *>(pParams);

        const double x = k * 8.0;
        const double w = (x < 0.001)? 1.0-0.1*x*x : 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);

        static double nspect = static_cast<double>(pcc->cosmo_param_["n_s"]);
        double tf = pcc->transfer_function_->compute(k, delta_matter0);

        //... no growth factor since we compute at z=0 and normalize so that D+(z=0)=1
        return k * k * w * w * std::pow(k, nspect) * tf * tf;
    }

    //! Computes the amplitude of a mode from the power spectrum
    /*! Function evaluates the supplied transfer function transfer_function_
	 * and returns the amplitude of fluctuations (norm * T(k) * sqrt(k^n_s)) at wave 
     * number k (in h/Mpc) back-scaled to z=z_start
	 * @param k wave number at which to evaluate
     * @param type one of the species: {delta,theta}_{matter,cdm,baryon,neutrino}
	 */
    inline real_t get_amplitude( const real_t k, const tf_type type) const
    {
        return std::pow(k, 0.5 * m_n_s_) * transfer_function_->compute(k, type) * m_sqrtpnorm_;
    }

    //! Computes the normalised transfer function T(k)
    /*! Function evaluates the supplied transfer function transfer_function_
        * and returns the transfer function (-norm * T(k) * k^2) at wave number k 
        * (in h/Mpc) back-scaled to z=z_start
        * @param k wave number at which to evaluate
        * @param type one of the species: {delta,theta}_{matter,cdm,baryon,neutrino}
     */
    inline real_t get_transfer( const real_t k, const tf_type type) const
    {
        return -transfer_function_->compute(k, type)*k*k / tnorm_ * m_sqrtpnorm_;
    }

    //! Compute amplitude of the back-scaled delta_bc mode, with decaying velocity v_bc included or not (in which case delta_bc=const)
    inline real_t get_amplitude_delta_bc( const real_t k, bool withvbc ) const
    {
        const real_t Dratio = Dplus_target_ / Dplus_start_;
        const real_t dbc = transfer_function_->compute(k, delta_bc) + (withvbc? 2 * transfer_function_->compute(k, theta_bc) * (std::sqrt(Dratio) - 1.0) : 0.0);
        // need to multiply with Dplus_target since sqrtpnorm rescales like that
        return std::pow(k, 0.5 * m_n_s_) * dbc * (m_sqrtpnorm_ * Dplus_target_);
    }

    //! Compute amplitude of the back-scaled relative velocity theta_bc mode if withvbc==true, otherwise return zero
    inline real_t get_amplitude_theta_bc( const real_t k, bool withvbc ) const
    {
        const real_t Dratio = Dplus_target_ / Dplus_start_;
        const real_t tbc = transfer_function_->compute(k, theta_bc) * std::sqrt(Dratio);
        // need to multiply with Dplus_target since sqrtpnorm rescales like that
        return withvbc ? std::pow(k, 0.5 * m_n_s_) * tbc * (m_sqrtpnorm_ * Dplus_target_) : 0.0;
    }


    //! Computes the normalization for the power spectrum
    /*!
	 * integrates the power spectrum to fix the normalization to that given
	 * by the sigma_8 parameter
	 */
    real_t compute_sigma8(void)
    {
        real_t sigma0, kmin, kmax;
        kmax = transfer_function_->get_kmax();
        kmin = transfer_function_->get_kmin();

        // if (!transfer_function_->tf_has_total0())
        sigma0 = 4.0 * M_PI * integrate(&dSigma8, static_cast<double>(kmin), static_cast<double>(kmax), this);
        // else{
        //     sigma0 = 4.0 * M_PI * integrate(&dSigma8_0, static_cast<double>(kmin), static_cast<double>(kmax), this);
        // }

        return std::sqrt(sigma0);
    }

    //! Computes the normalization for the power spectrum
    /*!
	 * integrates the power spectrum to fix the normalization to that given
	 * by the sigma_8 parameter
	 */
    real_t compute_pnorm_from_sigma8(void)
    {
        auto measured_sigma8 = this->compute_sigma8();
        return cosmo_param_["sigma_8"] * cosmo_param_["sigma_8"] / (measured_sigma8  * measured_sigma8);
    }
};

//! compute the jeans sound speed
/*! given a density in g/cm^-3 and a mass in g it gives back the sound
 *  speed in cm/s for which the input mass is equal to the jeans mass
 *  @param rho density 
 *  @param mass mass scale
 *  @returns jeans sound speed
 */
// inline double jeans_sound_speed(double rho, double mass)
// {
//     const double G = 6.67e-8;
//     return pow(6.0 * mass / M_PI * std::sqrt(rho) * std::pow(G, 1.5), 1.0 / 3.0);
// }

} // namespace cosmology
