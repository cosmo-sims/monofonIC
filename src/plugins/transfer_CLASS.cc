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

#include <general.hh>
#include <config_file.hh>
#include <transfer_function_plugin.hh>
#include <ic_generator.hh>

#include <math/interpolate.hh>

// disable warnings that pop up due to CLASS headers non-compliance.
#pragma GCC diagnostic ignored "-Wvariadic-macros"
#include "class.h"
#pragma GCC diagnostic pop

#include "transfer_CLASS.hh"


class transfer_CLASS_plugin : public TransferFunction_plugin
{
protected:
    // structs containing all the CLASS parameters
    struct file_content fc_;
    struct precision pr_;
    struct background ba_;
    struct thermodynamics th_;           /* for thermodynamics */
    struct perturbations pt_;         /* for source functions */
    struct transfer tr_;        /* for transfer functions */
    struct primordial pm_;       /* for primordial spectra */
    struct harmonic sp_;          /* for output spectra */
    struct fourier nl_;        /* for non-linear spectra */
    struct lensing le_;          /* for lensed spectra */
    struct distortions sd_;          /* for distortions */
    struct output op_;           /* for output files */

    ClassParams pars_;
    bool class_dofree_;

private:
  
  using TransferFunction_plugin::cosmo_params_;

  interpolated_function_1d<true, true, false> delta_c_, delta_b_, delta_n_, delta_m_, theta_c_, theta_b_, theta_n_, theta_m_;
  interpolated_function_1d<true, true, false> delta_c0_, delta_b0_, delta_n0_, delta_m0_, theta_c0_, theta_b0_, theta_n0_, theta_m0_;

  double zstart_, ztarget_, astart_, atarget_, kmax_, kmin_, h_, tnorm_;

  
  // std::unique_ptr<ClassEngine> the_ClassEngine_;
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
    add_class_parameter("z_max_pk", std::max(std::max(zstart_, ztarget_),199.0)); // use 1.2 as safety
    add_class_parameter("P_k_max_h/Mpc", std::max(2.0,kmax_));
    add_class_parameter("output", "dTk,vTk");
    add_class_parameter("extra metric transfer functions","yes");
    // add_class_parameter("lensing", "no");

    //--- choose gauge ------------------------------------------------
    // add_class_parameter("extra metric transfer functions", "yes");
    add_class_parameter("gauge", "synchronous");

    //--- cosmological parameters, densities --------------------------
    add_class_parameter("h", cosmo_params_.get("h"));

    add_class_parameter("Omega_b", cosmo_params_.get("Omega_b"));
    add_class_parameter("Omega_cdm", cosmo_params_.get("Omega_c"));
    add_class_parameter("Omega_k", cosmo_params_.get("Omega_k"));
    // add_class_parameter("Omega_fld", 0.0);
    // add_class_parameter("Omega_scf", 0.0);


    // add_class_parameter("fluid_equation_of_state","CLP");
    add_class_parameter("Omega_Lambda", 0.0);
    add_class_parameter("w0_fld", cosmo_params_.get("w_0") );
    add_class_parameter("wa_fld", cosmo_params_.get("w_a") );
    // add_class_parameter("cs2_fld", 1);

    //--- massive neutrinos -------------------------------------------
#if 0
    //default off
    // add_class_parameter("Omega_ur",0.0);
    add_class_parameter("N_ur", cosmo_params_.get("N_ur"));
    add_class_parameter("N_ncdm", 0);

#else
    
    add_class_parameter("N_ur", cosmo_params_.get("N_ur"));
    add_class_parameter("N_ncdm", cosmo_params_.get("N_nu_massive"));
    if( cosmo_params_.get("N_nu_massive") > 0 ){
      std::stringstream sstr;
      if( cosmo_params_.get("m_nu1") > 1e-9 ) sstr << cosmo_params_.get("m_nu1");
      if( cosmo_params_.get("m_nu2") > 1e-9 ) sstr << ", " << cosmo_params_.get("m_nu2");
      if( cosmo_params_.get("m_nu3") > 1e-9 ) sstr << ", " << cosmo_params_.get("m_nu3");
      add_class_parameter("m_ncdm", sstr.str().c_str());
    }
    
    // change above to enable
    //add_class_parameter("omega_ncdm", 0.0006451439);
    //add_class_parameter("m_ncdm", "0.4");
    //add_class_parameter("T_ncdm", 0.71611);
#endif

    //--- cosmological parameters, primordial -------------------------
    add_class_parameter("P_k_ini type", "analytic_Pk");

    if( cosmo_params_.get("A_s") > 0.0 ){
      add_class_parameter("A_s", cosmo_params_.get("A_s"));
    }else{
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
      zlist << ztarget_ << ((ztarget_!=0.0)? ", 0.0" : "");
    else
      zlist << std::max(ztarget_, zstart_) << ", " << std::min(ztarget_, zstart_) << ", 0.0";
    add_class_parameter("z_pk", zlist.str());

    music::ilog << "Computing transfer function via ClassEngine..." << std::endl;
    double wtime = get_wtime();

    //.........................................................................
    this->setup_class( pars_ );
    //.........................................................................

    wtime = get_wtime() - wtime;
    music::ilog << "CLASS took " << wtime << " s." << std::endl;
  }

  /*! \brief Set up CLASS engine with parameters
   *
   *  \param pars ClassParams object containing the parameters
   *  \return _SUCCESS_ or _FAILURE_
   */
  int setup_class( const ClassParams &pars )
  {
    
    ErrorMsg errmsg;

    //prepare fp structure
    std::vector<std::string> parNames;
    size_t n=pars.size();
    char dummystr[] = "dummy";
    parser_init(&fc_,n,dummystr, errmsg);
    
    //config
    for (size_t i=0;i<pars.size();i++){
      strcpy(fc_.name[i],pars.key(i).c_str());
      strcpy(fc_.value[i],pars.value(i).c_str());
      parNames.push_back(pars.key(i));
    }

    // if (input_init(&fc_,&pr_,&ba_,&th_,&pt_,&tr_,&pm_,&sp_,&nl_,&le_,&op_,errmsg) == _FAILURE_) {
    if (input_read_from_file(&fc_, &pr_, &ba_, &th_, &pt_, &tr_, &pm_, &sp_, &nl_, &le_, &sd_, &op_, errmsg) == _FAILURE_){
      music::elog.Print(" running input_init_from_arguments \n=>%s\n",errmsg);
      class_dofree_=false;
      return _FAILURE_;
    }

    // proetction parametres mal defini
    for (size_t i = 0; i < pars.size(); i++)
    {
      if (fc_.read[i] != _TRUE_)
        throw std::runtime_error(std::string("invalid CLASS parameter: ") + fc_.name[i]);
    }

    if (background_init(&pr_,&ba_) == _FAILURE_) {
      music::elog.Print(" running background_init \n=>%s\n",ba_.error_message);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (thermodynamics_init(&pr_,&ba_,&th_) == _FAILURE_) {
      music::elog.Print(" in thermodynamics_init \n=>%s\n",th_.error_message);
      background_free(&ba_);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (perturbations_init(&pr_,&ba_,&th_,&pt_) == _FAILURE_) {
      music::elog.Print(" in perturb_init \n=>%s\n",pt_.error_message);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (primordial_init(&pr_,&pt_,&pm_) == _FAILURE_) {
      music::elog.Print(" in primordial_init \n=>%s\n",pm_.error_message);
      perturbations_free(&pt_);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (fourier_init(&pr_,&ba_,&th_,&pt_,&pm_,&nl_) == _FAILURE_)  {
      music::elog.Print(" in fourier_init \n=>%s\n",nl_.error_message);
      primordial_free(&pm_);
      perturbations_free(&pt_);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (transfer_init(&pr_,&ba_,&th_,&pt_,&nl_,&tr_) == _FAILURE_) {
      music::elog.Print(" in transfer_init \n=>%s\n",tr_.error_message);
      fourier_free(&nl_);
      primordial_free(&pm_);
      perturbations_free(&pt_);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (harmonic_init(&pr_,&ba_,&pt_,&pm_,&nl_,&tr_,&sp_) == _FAILURE_) {
      music::elog.Print(" in harmonic_init \n=>%s\n",sp_.error_message);
      transfer_free(&tr_);
      fourier_free(&nl_);
      primordial_free(&pm_);
      perturbations_free(&pt_);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (lensing_init(&pr_,&pt_,&sp_,&nl_,&le_) == _FAILURE_) {
      music::elog.Print(" in lensing_init \n=>%s\n",le_.error_message);
      harmonic_free(&sp_);
      transfer_free(&tr_);
      fourier_free(&nl_);
      primordial_free(&pm_);
      perturbations_free(&pt_);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_=false;
      return _FAILURE_;
    }

    if (distortions_init(&pr_, &ba_, &th_, &pt_, &pm_, &sd_) == _FAILURE_)
    {
      music::elog.Print(" in distortions_init \n=>%s\n", sd_.error_message);
      lensing_free(&le_);
      harmonic_free(&sp_);
      transfer_free(&tr_);
      fourier_free(&nl_);
      primordial_free(&pm_);
      perturbations_free(&pt_);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_ = false;
      return _FAILURE_;
    }

    class_dofree_=true;
    return _SUCCESS_;
  }

  int free_class()
  {
    if (class_dofree_) {
      distortions_free(&sd_);
      lensing_free(&le_);
      harmonic_free(&sp_);
      transfer_free(&tr_);
      fourier_free(&nl_);
      primordial_free(&pm_);
      perturbations_free(&pt_);
      thermodynamics_free(&th_);
      background_free(&ba_);
      class_dofree_=false;
    }
    return _SUCCESS_;
  }

  void call_perturb_sources_at_tau(
                            int index_md,
                            int index_ic,
                            int index_tp,
                            double tau,
                            double * psource
                            ) {
    if (perturbations_sources_at_tau(&pt_, index_md, index_ic, index_tp, tau, psource) == _FAILURE_)
    {
      std::cerr << ">>>fail getting Tk type=" << (int)index_tp <<std::endl; 
      throw std::runtime_error(pt_.error_message);
    }
  }

  //! run ClassEngine with parameters set up
  int run_ClassEngine(double z, std::vector<double> &k, std::vector<double> &d_cdm, std::vector<double> &t_cdm, std::vector<double> &d_b, std::vector<double> &t_b,
                       std::vector<double> &d_ncdm, std::vector<double> &t_ncdm, std::vector<double> &d_tot, std::vector<double> &t_tot)
  {
    k.clear(); 
    d_cdm.clear(); d_b.clear(); d_ncdm.clear(); d_tot.clear();
    t_cdm.clear(); t_b.clear(); t_ncdm.clear(); t_tot.clear();
    
    z = std::max(z,1e-10);
    // the_ClassEngine_->getTk(z, k, dc, db, dn, dm, tc, tb, tn, tm);

    //transform redshift in conformal time
    double tau;
    int index;
    background_tau_of_z(&ba_,z,&tau);
    if(std::log(tau) < pt_.ln_tau[0]){
      music::elog << "Asking sources at a z bigger than z_max_pk, something probably went wrong\n";
      throw std::runtime_error(pt_.error_message);
    }

    double *pvecback=new double[ba_.bg_size];
    background_at_tau(&ba_,tau,long_info,inter_normal, &index, pvecback);
    double fHa = pvecback[ba_.index_bg_f] * (pvecback[ba_.index_bg_a]*pvecback[ba_.index_bg_H]);
    delete[] pvecback;

    //...
    // copy transfer func data to temporary
    const size_t index_md = pt_.index_md_scalars;
    d_cdm.assign( pt_.k_size[index_md], 0.0 );
    d_b.assign( pt_.k_size[index_md], 0.0 );
    d_ncdm.assign( pt_.k_size[index_md], 0.0 );
    d_tot.assign( pt_.k_size[index_md], 0.0 );
    t_cdm.assign( pt_.k_size[index_md], 0.0 );
    t_b.assign( pt_.k_size[index_md], 0.0 );
    t_ncdm.assign( pt_.k_size[index_md], 0.0 );
    t_tot.assign( pt_.k_size[index_md], 0.0 );

    if( pt_.ic_size[index_md] > 1 ){
      std::cerr << ">>>have more than 1 ICs, will use first and ignore others" << std::endl;
    }

    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_delta_cdm, tau, &d_cdm[0]);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_delta_b, tau, &d_b[0]);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_delta_ncdm1, tau, &d_ncdm[0]);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_delta_m, tau, &d_tot[0]);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_theta_b, tau, &t_b[0]);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_theta_ncdm1, tau, &t_ncdm[0]);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_theta_m, tau, &t_tot[0]);

    //
    std::vector<double> h_prime(pt_.k_size[index_md],0.0), eta_prime(pt_.k_size[index_md],0.0);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_eta_prime, tau, &eta_prime[0]);
    call_perturb_sources_at_tau(index_md, 0, pt_.index_tp_h_prime, tau, &h_prime[0]);

    // gauge trafo velocities, store k-vector
    for (int index_k=0; index_k<pt_.k_size[index_md]; index_k++) 
    {
      auto ak = pt_.k[index_md][index_k];

      // write data to vectors
      k.push_back( ak );

      // use the conformal Newtonian gauge for velocities
      // not correct, but N-body gauge currently not implemented
      double alphak2 = (h_prime[index_k]+6*eta_prime[index_k])/2; 

      t_cdm[index_k]  = (-alphak2) / fHa;
      t_b[index_k]    = (-alphak2 + t_b[index_k]) / fHa;
      t_ncdm[index_k] = (-alphak2 + t_ncdm[index_k]) / fHa;
      t_tot[index_k]  = (-alphak2 + t_tot[index_k]) / fHa;
    }

    //...

    const double h  = cosmo_params_.get("h");

    for (size_t i = 0; i < k.size(); ++i)
    {
      // convert to 'CAMB' format, since we interpolate loglog and
      // don't want negative numbers...
      auto ik2  = 1.0 / (k[i] * k[i]) * h * h;
      d_cdm[i]  = -d_cdm[i] * ik2;
      d_b[i]    = -d_b[i] * ik2;
      d_ncdm[i] = -d_ncdm[i] * ik2;
      d_tot[i]  = -d_tot[i] * ik2;
      t_cdm[i]  = -t_cdm[i] * ik2;
      t_b[i]    = -t_b[i] * ik2;
      t_ncdm[i] = -t_ncdm[i] * ik2;
      t_tot[i]  = -t_tot[i] * ik2;
    }

    return _SUCCESS_;
  }

public:
  explicit transfer_CLASS_plugin(config_file &cf, const cosmology::parameters& cosmo_params)
      : TransferFunction_plugin(cf,cosmo_params)
  {
    this->tf_isnormalised_ = true;

    ofs_class_input_.open(cf.get_path_relative_to_config("input_class_parameters.ini"), std::ios::trunc);

    // all cosmological parameters need to be passed through the_cosmo_calc
    
    ztarget_ = pcf_->get_value_safe<double>("cosmology", "ztarget", 0.0);
    atarget_ = 1.0 / (1.0 + ztarget_);
    zstart_ = pcf_->get_value<double>("setup", "zstart");
    astart_ = 1.0 / (1.0 + zstart_);

    h_ = cosmo_params_["h"];
    
    if (cosmo_params_["A_s"] > 0.0) {
      music::ilog << "CLASS: Using A_s=" << cosmo_params_["A_s"] << " to normalise the transfer function." << std::endl;
    }else{
      double sigma8 = cosmo_params_["sigma_8"];
      if( sigma8 < 0 ){
        throw std::runtime_error("Need to specify either A_s or sigma_8 for CLASS plugin...");
      }
      music::ilog << "CLASS: Using sigma8_ =" << sigma8<< " to normalise the transfer function." << std::endl;
    }

    // determine highest k we will need for the resolution selected
    double lbox = pcf_->get_value<double>("setup", "BoxLength");
    int nres = pcf_->get_value<double>("setup", "GridRes");
    kmax_ = std::max(20.0, 2.0 * M_PI / lbox * nres / 2 * sqrt(3) * 2.0); // 120% of spatial diagonal, or k=10h Mpc-1

    // initialise CLASS and get the normalisation
    this->init_ClassEngine();
    double A_s_ = pm_.A_s;//the_ClassEngine_->get_A_s(); // this either the input one, or the one computed from sigma8
    
    // compute the normalisation to interface with MUSIC
    double k_p = cosmo_params["k_p"] / cosmo_params["h"];
    tnorm_ = std::sqrt(2.0 * M_PI * M_PI * A_s_ * std::pow(1.0 / k_p, cosmo_params["n_s"] - 1) / std::pow(2.0 * M_PI, 3.0));

    // compute the transfer function at z=0 using CLASS engine
    std::vector<double> k, dc, tc, db, tb, dn, tn, dm, tm;
    this->run_ClassEngine(0.0, k, dc, tc, db, tb, dn, tn, dm, tm);

    delta_c0_.set_data(k, dc);
    theta_c0_.set_data(k, tc);
    delta_b0_.set_data(k, db);
    theta_b0_.set_data(k, tb);
    delta_n0_.set_data(k, dn);
    theta_n0_.set_data(k, tn);
    delta_m0_.set_data(k, dm);
    theta_m0_.set_data(k, tm);

     // compute the transfer function at z=z_target using CLASS engine
    this->run_ClassEngine(ztarget_, k, dc, tc, db, tb, dn, tn, dm, tm);
    delta_c_.set_data(k, dc);
    theta_c_.set_data(k, tc);
    delta_b_.set_data(k, db);
    theta_b_.set_data(k, tb);
    delta_n_.set_data(k, dn);
    theta_n_.set_data(k, tn);
    delta_m_.set_data(k, dm);
    theta_m_.set_data(k, tm);

    kmin_ = k[0];
    kmax_ = k.back();

    music::ilog << "CLASS table contains k = " << this->get_kmin() << " to " << this->get_kmax() << " h Mpc-1." << std::endl;

    tf_distinct_ = true;
    tf_withvel_ = true;
    tf_withtotal0_ = true;
  }

  ~transfer_CLASS_plugin()
  {
    this->free_class();
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
      val = delta_b_(k)-delta_c_(k); break;
    case theta_bc:
      val = theta_b_(k)-theta_c_(k); break;

      // values at zstart:
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

namespace
{
TransferFunction_plugin_creator_concrete<transfer_CLASS_plugin> creator("CLASS");
}

#endif // USE_CLASS
