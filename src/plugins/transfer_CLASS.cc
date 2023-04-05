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

#ifdef USE_CLASS

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <sstream>

#include <ClassEngine.hh>

#include <general.hh>
#include <config_file.hh>
#include <transfer_function_plugin.hh>
#include <ic_generator.hh>

#include <math/interpolate.hh>

class transfer_CLASS_plugin : public TransferFunction_plugin
{
private:
  using TransferFunction_plugin::cosmo_params_;

  interpolated_function_1d<true, false, false> delta_c_, delta_b_, delta_n_, delta_t_, theta_c_, theta_b_, theta_n_, theta_t_;
  interpolated_function_1d<true, false, false> delta_c0_, delta_b0_, delta_n0_, delta_t0_, theta_c0_, theta_b0_, theta_n0_, theta_t0_;

  double omega_ncdm_total_, omega_nu_total_, omega_c_;
  std::vector<double> omega_ncdm_, m_ncdm_, T_ncdm_, omega_nu_;
  size_t N_ncdm_, N_nu_massive_;

  double zstart_, ztarget_, astart_, atarget_, kmax_, kmin_, h_, tnorm_, N_nu_massive;

  ClassParams pars_;
  std::unique_ptr<ClassEngine> the_ClassEngine_;
  std::ofstream ofs_class_input_;

  template <typename T>
  void add_class_parameter(std::string parameter_name, const T parameter_value)
  {
    pars_.add(parameter_name, parameter_value);
    ofs_class_input_ << parameter_name << " = " << parameter_value << std::endl;
  }

  //! Set up class parameters from MUSIC cosmological parameters
  void init_ClassEngine(void)
  {
    //--- general parameters ------------------------------------------
    add_class_parameter("z_max_pk", std::max(std::max(zstart_, ztarget_), 199.0)); // use 1.2 as safety
    add_class_parameter("P_k_max_h/Mpc", std::max(2.0, kmax_));
    add_class_parameter("output", "dTk,vTk");
    add_class_parameter("extra metric transfer functions", "yes");
    // add_class_parameter("lensing", "no");

    //--- choose gauge ------------------------------------------------
    // add_class_parameter("extra metric transfer functions", "yes");
    add_class_parameter("gauge", "synchronous");

    //--- cosmological parameters, densities --------------------------
    add_class_parameter("h", cosmo_params_.get("h"));

    add_class_parameter("Omega_b", cosmo_params_.get("Omega_b"));
    add_class_parameter("Omega_k", cosmo_params_.get("Omega_k"));
    
    //--- dark energy EOS ---------------------------------------------
    // CLASS docu:
    // 8a) Dark energy contributions. At least one out of three conditions must be satisfied:
    // i)   'Omega_Lambda' unspecified.
    // ii)  'Omega_fld' unspecified.
    // iii) 'Omega_scf' set to a negative value. [Will be refered to as
    //      unspecified in the following text.]
    // The code will then use the first unspecified component to satisfy the
    // closure equation (sum_i Omega_i) equals (1 + Omega_k)
    // (default: 'Omega_fld' and 'Omega_scf' set to 0 and 'Omega_Lambda' inferred by code)

    add_class_parameter("fluid_equation_of_state","CLP");
    add_class_parameter("Omega_Lambda", 0.0);
    add_class_parameter("Omega_scf", 0.0);
    add_class_parameter("w0_fld", cosmo_params_.get("w_0") );
    add_class_parameter("wa_fld", cosmo_params_.get("w_a") );
    add_class_parameter("cs2_fld", 1);

    //--- massive neutrinos -------------------------------------------
    // TBD: add explanations of massive neutrinos and nCDM species
    add_class_parameter("N_ur", cosmo_params_.get("N_ur"));

    // treating of massive trans-relativistic species
    std::stringstream sstr_m, sstr_omega, sstr_T;
    const double sum_m_nu = cosmo_params_.get("m_nu1") + cosmo_params_.get("m_nu2") + cosmo_params_.get("m_nu3");

    if (N_ncdm_ > 0)
    {
      for (size_t i = 0; i < omega_ncdm_.size(); ++i)
      {
        sstr_m << m_ncdm_[i];
        sstr_omega << omega_ncdm_[i];
        sstr_T << T_ncdm_[i];
        omega_ncdm_total_ += omega_ncdm_[i];

        if( i<omega_ncdm_.size()-1 )
        {
          sstr_m << ", ";
          sstr_omega << ", ";
          sstr_T << ", ";
        }

        music::ilog.Print(" * Found extra nCDM species: m=%g eV, T=%g (T_gamma), Omega=%g", m_ncdm_[i], T_ncdm_[i], omega_ncdm_[i]);
      }
      omega_c_ -= omega_ncdm_total_;

      music::ilog << "Energy density of all DM species:" << std::endl;
      music::ilog.Print(" Omega_ncdm = %8g,     Omega_c = %8g\n", omega_ncdm_total_, omega_c_);
      music::ilog << std::endl;
    }

    N_nu_massive_ = 0;
    if (cosmo_params_.get("N_nu_massive") > 0)
    {
      // if already nCDM species were written, need to add a comma
      if (N_ncdm_ > 0)
      {
        sstr_m << ", ";
        sstr_omega << ", ";
        sstr_T << ", ";
      }

      if (cosmo_params_.get("m_nu1") > 1e-9)
      {
        const double omeganu1 = cosmo_params_.get("m_nu1") / sum_m_nu * cosmo_params_.get("Omega_nu_massive");
        sstr_m << cosmo_params_.get("m_nu1");
        sstr_omega << omeganu1;
        sstr_T << 0.71611; // CLASS recommended value for active massive neutrinos, yields m/omega = 93.14 eV
        omega_nu_.push_back(omeganu1);
        omega_nu_total_ += omeganu1;
        N_nu_massive_ += 1;
        if (cosmo_params_.get("m_nu2") > 1e-9)
        {
          const double omeganu2 = cosmo_params_.get("m_nu2") / sum_m_nu * cosmo_params_.get("Omega_nu_massive");
          sstr_m << ", " << cosmo_params_.get("m_nu2");
          sstr_omega << ", " << omeganu2;
          sstr_T << ", " << 0.71611; // CLASS recommended value for active massive neutrinos, yields m/omega = 93.14 eV
          omega_nu_.push_back(omeganu2);
          omega_nu_total_ += omeganu2;
          N_nu_massive_ += 1;
          if (cosmo_params_.get("m_nu3") > 1e-9)
          {
            const double omeganu3 = cosmo_params_.get("m_nu3") / sum_m_nu * cosmo_params_.get("Omega_nu_massive");
            sstr_m << ", " << cosmo_params_.get("m_nu3");
            sstr_omega << ", " << omeganu3;
            sstr_T << ", " << 0.71611; // CLASS recommended value for active massive neutrinos, yields m/omega = 93.14 eV
            omega_nu_.push_back(omeganu3);
            omega_nu_total_ += omeganu3;
            N_nu_massive_ += 1;
          }
        }
      }
    }

    add_class_parameter("Omega_cdm", std::max(omega_c_, 1e-9));

    add_class_parameter("Omega_ncdm", sstr_omega.str().c_str());
    add_class_parameter("m_ncdm", sstr_m.str().c_str());
    add_class_parameter("T_ncdm", sstr_T.str().c_str());

    add_class_parameter("N_ncdm", cosmo_params_.get("N_nu_massive") + N_ncdm_);

    //--- cosmological parameters, primordial -------------------------
    add_class_parameter("P_k_ini type", "analytic_Pk");

    if (cosmo_params_.get("A_s") > 0.0)
    {
      add_class_parameter("A_s", cosmo_params_.get("A_s"));
    }
    else
    {
      add_class_parameter("sigma8", cosmo_params_.get("sigma_8"));
    }
    add_class_parameter("n_s", cosmo_params_.get("n_s"));
    add_class_parameter("alpha_s", 0.0);
    add_class_parameter("T_cmb", cosmo_params_.get("Tcmb"));
    add_class_parameter("YHe", cosmo_params_.get("YHe"));

    // additional parameters
    add_class_parameter("reio_parametrization", "reio_none");

    // precision parameters
    add_class_parameter("k_per_decade_for_pk", 100);
    add_class_parameter("k_per_decade_for_bao", 100);
    add_class_parameter("compute damping scale", "yes");
    add_class_parameter("tol_perturbations_integration", 1.e-8);
    add_class_parameter("tol_background_integration", 1e-9);

    // high precision options from cl_permille.pre:
    // precision file to be passed as input in order to achieve at least percent precision on scalar Cls
    add_class_parameter("hyper_flat_approximation_nu", 7000.);
    add_class_parameter("transfer_neglect_delta_k_S_t0", 0.17);
    add_class_parameter("transfer_neglect_delta_k_S_t1", 0.05);
    add_class_parameter("transfer_neglect_delta_k_S_t2", 0.17);
    add_class_parameter("transfer_neglect_delta_k_S_e", 0.13);
    add_class_parameter("delta_l_max", 1000);

    int class_verbosity = 0;

    add_class_parameter("background_verbose", class_verbosity);
    add_class_parameter("thermodynamics_verbose", class_verbosity);
    add_class_parameter("perturbations_verbose", class_verbosity);
    add_class_parameter("transfer_verbose", class_verbosity);
    add_class_parameter("primordial_verbose", class_verbosity);
    add_class_parameter("harmonic_verbose", class_verbosity);
    add_class_parameter("fourier_verbose", class_verbosity);
    add_class_parameter("lensing_verbose", class_verbosity);
    add_class_parameter("output_verbose", class_verbosity);

    // output parameters, only needed for the control CLASS .ini file that we output
    std::stringstream zlist;
    if (ztarget_ == zstart_)
      zlist << ztarget_ << ((ztarget_ != 0.0) ? ", 0.0" : "");
    else
      zlist << std::max(ztarget_, zstart_) << ", " << std::min(ztarget_, zstart_) << ", 0.0";
    add_class_parameter("z_pk", zlist.str());

    music::ilog << ">>> Computing transfer function via ClassEngine..." << std::endl;
    double wtime = get_wtime();

    the_ClassEngine_ = std::make_unique<ClassEngine>(pars_, false);

    music::ilog << std::setw(70) << std::setfill(' ') << std::right << "took : " << std::setw(8) << get_wtime() - wtime << "s" << std::endl;
  }

  //! run ClassEngine with parameters set up
  void run_ClassEngine(double z, int gauge, std::vector<double> &k, std::vector<double> &dc, std::vector<double> &tc, std::vector<double> &db, std::vector<double> &tb,
                       std::vector<double> &dn, std::vector<double> &tn, std::vector<double> &dm, std::vector<double> &tm, std::vector<double> &dt, std::vector<double> &tt,
                       std::vector<double> &phi_or_h_prime, std::vector<double> &psi_or_eta_prime)
  {
    k.clear();
    dc.clear();
    db.clear();
    dn.clear();
    dm.clear();
    dt.clear();
    tc.clear();
    tb.clear();
    tn.clear();
    tm.clear();
    tt.clear();

    std::vector<std::vector<double>> dncdm, tncdm;

    // if gauge < 0  : run class in synchronous gauge, and convert *only* velocities to conformal Newtonian gauge below, phi_or_h_prime is h_prime, psi_or_eta_prime is eta_prime
    // if gauge == 0 : run class in synchronous gauge, do not convert velocities, phi_or_h_prime is h_prime, psi_or_eta_prime is eta_prime
    // if gauge == 1 : run class in Newtonian gauge, do not convert velocities, phi_or_h is phi, psi_or_eta_prime is psi
    // ...
    const int class_gauge = (gauge < 0) ? 0 : gauge;
    double fHa{0.0};
    the_ClassEngine_->getTk(z, class_gauge, k, dc, db, dncdm, dm, dt, tc, tb, tncdm, tm, tt, phi_or_h_prime, psi_or_eta_prime, fHa);

    // allocate memory for neutrinos
    dn.assign(k.size(), 0.0);
    tn.assign(k.size(), 0.0);

    for (size_t index_k = 0; index_k < k.size(); ++index_k)
    {
      // add all neutrino species into one single massive neutrino TF
      for (size_t i = N_ncdm_, j = 0; i < N_ncdm_ + N_nu_massive_; ++i, ++j)
      {
        dn[index_k] += omega_nu_[j] / omega_nu_total_ * dncdm[i][index_k];
        tn[index_k] += omega_nu_[j] / omega_nu_total_ * tncdm[i][index_k];
      }
      // add all dark matter species into one single mixed-dark-matter TF (we do not support multiple DM species yet)
      dc[index_k] *= omega_c_ / (omega_c_ + omega_ncdm_total_);
      tc[index_k] *= omega_c_ / (omega_c_ + omega_ncdm_total_);
      for (size_t i = 0; i < N_ncdm_; ++i)
      {
        dc[index_k] += omega_ncdm_[i] / (omega_c_ + omega_ncdm_total_) * dncdm[i][index_k];
        tc[index_k] += omega_ncdm_[i] / (omega_c_ + omega_ncdm_total_) * tncdm[i][index_k];
      }
    }

    // convert velocities to conformal Newtonian gauge if gauge < 0
    if (gauge < 0)
    {
      for (size_t index_k = 0; index_k < k.size(); ++index_k)
      {
        // gauge transformation to conformal Newtonian velocity, Ma & Bertschinger eq.(27b)
        double alphak2 = (phi_or_h_prime[index_k] + 6 * psi_or_eta_prime[index_k]) / 2;
        tc[index_k] = (-alphak2) / fHa;
        tb[index_k] = (-alphak2 + tb[index_k]) / fHa;
        tn[index_k] = (-alphak2 + tn[index_k]) / fHa;
        tm[index_k] = (-alphak2 + tm[index_k]) / fHa;
        tt[index_k] = (-alphak2 + tt[index_k]) / fHa;
      }
    }
    else
    {
      _unused(fHa);
    }

    const double h = cosmo_params_.get("h");

    for (size_t i = 0; i < k.size(); ++i)
    {
      // convert to 'CAMB' format, since we interpolate loglog and
      // don't want negative numbers...
      auto ik2 = 1.0 / (k[i] * k[i]) * h * h;
      dc[i] = -dc[i] * ik2;
      db[i] = -db[i] * ik2;
      dn[i] = -dn[i] * ik2;
      dm[i] = -dm[i] * ik2;
      dt[i] = -dt[i] * ik2;
      tc[i] = -tc[i] * ik2;
      tb[i] = -tb[i] * ik2;
      tn[i] = -tn[i] * ik2;
      tm[i] = -tm[i] * ik2;
      tt[i] = -tt[i] * ik2;
      phi_or_h_prime[i] = -phi_or_h_prime[i] * ik2;
      psi_or_eta_prime[i] = -psi_or_eta_prime[i] * ik2;
    }
  }

public:
  explicit transfer_CLASS_plugin(config_file &cf, const cosmology::parameters &cosmo_params)
      : TransferFunction_plugin(cf, cosmo_params)
  {
    this->tf_isnormalised_ = true;

    ofs_class_input_.open("input_class_parameters.ini", std::ios::trunc);

    // all cosmological parameters need to be passed through the_cosmo_calc

    ztarget_ = pcf_->get_value_safe<double>("cosmology", "ztarget", 0.0);
    atarget_ = 1.0 / (1.0 + ztarget_);
    zstart_ = pcf_->get_value<double>("setup", "zstart");
    astart_ = 1.0 / (1.0 + zstart_);

    h_ = cosmo_params_["h"];

    if (cosmo_params_["A_s"] > 0.0)
    {
      music::ilog << "CLASS: Using A_s=" << cosmo_params_["A_s"] << " to normalise the transfer function." << std::endl;
    }
    else
    {
      double sigma8 = cosmo_params_["sigma_8"];
      if (sigma8 < 0)
      {
        throw std::runtime_error("Need to specify either A_s or sigma_8 for CLASS plugin...");
      }
      music::ilog << "CLASS: Using sigma8_ =" << sigma8 << " to normalise the transfer function." << std::endl;
    }

    // determine highest k we will need for the resolution selected
    double lbox = pcf_->get_value<double>("setup", "BoxLength");
    int nres = pcf_->get_value<double>("setup", "GridRes");
    kmax_ = std::max(20.0, 2.0 * M_PI / lbox * nres / 2 * sqrt(3) * 2.0); // 120% of spatial diagonal, or k=10h Mpc-1

    // initialize parameters for alternative dark matter models
    omega_c_ = cosmo_params.get("Omega_c");
    omega_ncdm_total_ = 0.0;
    omega_nu_total_ = 0.0;
    if (pcf_->contains_key("cosmology", "m_ncdm"))
    {
      m_ncdm_ = pcf_->get_value<std::vector<double>>("cosmology", "m_ncdm");
      if (m_ncdm_.size() > 1 && !pcf_->contains_key("cosmology", "Omega_ncdm"))
      {
        throw std::runtime_error("Need to specify \'cosmology/Omega_ncdm\' if multiple m_ncdm are specified");
      }
      else
      {
        if (m_ncdm_.size() == 1)
        {
          // if m_ncdm given for one non-CDM species, default to Omega_ncdm = Omega_c if not given explicitly
          omega_ncdm_.push_back(pcf_->get_value_safe<double>("cosmology", "Omega_ncdm", omega_c_));
          T_ncdm_.push_back(pcf_->get_value_safe<double>("cosmology", "T_ncdm", 0.7137658555));
        }
        else
        {

          omega_ncdm_ = pcf_->get_value<std::vector<double>>("cosmology", "Omega_ncdm");
          if (!pcf_->contains_key("cosmology", "T_ncdm"))
          {
            T_ncdm_.assign(m_ncdm_.size(), 0.7137658555); // (4/11)**(1/3)
          }
          else
          {
            T_ncdm_ = pcf_->get_value<std::vector<double>>("cosmology", "T_ncdm");
          }
        }
      }
      if (m_ncdm_.size() != omega_ncdm_.size())
      {
        throw std::runtime_error("Need to specify as many \'cosmology/Omega_ncdm\' as \'cosmology/m_ncdm\' (if using more than one).");
      }
      if (m_ncdm_.size() != T_ncdm_.size())
      {
        throw std::runtime_error("Need to specify as many \'cosmology/T_ncdm\' as \'cosmology/m_ncdm\' (if specifying at least one).");
      }
    }
    N_ncdm_ = m_ncdm_.size();

    // initialise CLASS and get the normalisation
    this->init_ClassEngine();
    double A_s_ = the_ClassEngine_->get_A_s(); // this either the input one, or the one computed from sigma8

    // compute the normalisation to interface with MUSIC
    double k_p = cosmo_params["k_p"] / cosmo_params["h"];
    tnorm_ = std::sqrt(2.0 * M_PI * M_PI * A_s_ * std::pow(1.0 / k_p, cosmo_params["n_s"] - 1) / std::pow(2.0 * M_PI, 3.0));

    // compute the transfer function at z=0 using CLASS engine
    constexpr int gauge{-1}; // always use synchronous gauge
    std::vector<double> k, dc, tc, db, tb, dn, tn, dm, tm, dt, tt, phi_or_h, psi_or_eta;
    this->run_ClassEngine(0.0, gauge, k, dc, tc, db, tb, dn, tn, dm, tm, dt, tt, phi_or_h, psi_or_eta);

    delta_c0_.set_data(k, dc);
    theta_c0_.set_data(k, tc);
    delta_b0_.set_data(k, db);
    theta_b0_.set_data(k, tb);
    delta_n0_.set_data(k, dn);
    theta_n0_.set_data(k, tn);
    delta_t0_.set_data(k, dm);
    theta_t0_.set_data(k, tm);

    // compute the transfer function at z=z_target using CLASS engine
    this->run_ClassEngine(ztarget_, gauge, k, dc, tc, db, tb, dn, tn, dm, tm, dt, tt, phi_or_h, psi_or_eta);

    delta_c_.set_data(k, dc);
    theta_c_.set_data(k, tc);
    delta_b_.set_data(k, db);
    theta_b_.set_data(k, tb);
    delta_n_.set_data(k, dn);
    theta_n_.set_data(k, tn);
    delta_t_.set_data(k, dm);
    theta_t_.set_data(k, tm);

    kmin_ = k[0];
    kmax_ = k.back();

    music::ilog << "CLASS table contains k = " << this->get_kmin() << " to " << this->get_kmax() << " h Mpc-1." << std::endl;

    tf_distinct_ = true;
    tf_withvel_ = true;
    tf_withtotal0_ = true;
  }

  ~transfer_CLASS_plugin()
  {
  }

  inline double compute(double k, tf_type type) const
  {
    k *= h_;

    if (k < kmin_ || k > kmax_)
    {
      return 0.0;
    }

    real_t val(0.0);
    switch (type)
    {
      // values at ztarget:
    case delta_matter:
      val = delta_t_(k);
      break;
    case delta_cdm:
      val = delta_c_(k);
      break;
    case delta_baryon:
      val = delta_b_(k);
      break;
    case theta_matter:
      val = theta_t_(k);
      break;
    case theta_cdm:
      val = theta_c_(k);
      break;
    case theta_baryon:
      val = theta_b_(k);
      break;
    case delta_bc:
      val = delta_b_(k) - delta_c_(k);
      break;
    case theta_bc:
      val = theta_b_(k) - theta_c_(k);
      break;

      // values at zstart:
    case delta_matter0:
      val = delta_t0_(k);
      break;
    case delta_cdm0:
      val = delta_c0_(k);
      break;
    case delta_baryon0:
      val = delta_b0_(k);
      break;
    case theta_matter0:
      val = theta_t0_(k);
      break;
    case theta_cdm0:
      val = theta_c0_(k);
      break;
    case theta_baryon0:
      val = theta_b0_(k);
      break;
    default:
      throw std::runtime_error("Invalid type requested in transfer function evaluation");
    }
    return val * tnorm_;
  }

  inline double get_kmin(void) const { return kmin_ / h_; }
  inline double get_kmax(void) const { return kmax_ / h_; }
};

namespace
{
  TransferFunction_plugin_creator_concrete<transfer_CLASS_plugin> creator("CLASS");
}

#endif // USE_CLASS
