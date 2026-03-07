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

#include <vector>
#include <transfer_function_plugin.hh>
#include <math/interpolate.hh>

class transfer_CAMB_file_plugin : public TransferFunction_plugin
{

private:

  using TransferFunction_plugin::cosmo_params_;

  interpolated_function_1d<true, true, false> delta_c_, delta_b_, delta_n_, delta_m_, theta_c_, theta_b_, theta_n_, theta_m_;

  double m_kmin, m_kmax;

  // bool m_linbaryoninterp;

  
  bool read_table( const std::string& filename )
  {

    size_t nlines{0};
    // m_linbaryoninterp = false;

    std::vector<double> k, dc, tc, db, tb, dn, tn, dm, tm;
    
    // Detect file format (9 or 13 columns)
    size_t num_columns = 0;
    bool has_velocities = false;

    if( CONFIG::MPI_task_rank == 0 )
    {
      music::ilog << "Reading tabulated transfer function data from file:" << std::endl
                  << "  \'" << filename << "\'" << std::endl;

      std::string line;
      std::ifstream ifs(filename.c_str());

      if (!ifs.good())
        throw std::runtime_error("Could not find transfer function file \'" + filename + "\'");
      
      // Detect number of columns by reading first non-comment line
      while (!ifs.eof() && num_columns == 0)
      {
        getline(ifs, line);
        if (ifs.eof())
          break;
        if (line.find("#") != std::string::npos || line.empty())
          continue;
        
        std::stringstream ss(line);
        double val;
        while (ss >> val)
          num_columns++;
      }
      
      ifs.clear();
      ifs.seekg(0, std::ios::beg);
      
      if (num_columns == 13)
      {
        has_velocities = true;
        music::ilog << "Detected 13-column CAMB format (with velocity fields)" << std::endl;
      }
      else if (num_columns == 9)
      {
        has_velocities = false;
        music::ilog << "Detected 9-column CAMB format (old format without velocity fields)" << std::endl;
        music::wlog << "Velocity transfer functions will be approximated from density transfer functions" << std::endl;
      }
      else
      {
        music::elog << "CAMB file has " << num_columns << " columns, expected 9 or 13" << std::endl;
        throw std::runtime_error("CAMB transfer file has unexpected number of columns");
      }
      
      while (!ifs.eof())
      {
        getline(ifs, line);
        if (ifs.eof())
          break;

        // OH: ignore line if it has a comment:
        if (line.find("#") != std::string::npos)
          continue;

        std::stringstream ss(line);

        double Tk, Tdc, Tdb, Tdn, Tdm, Tvb, Tvc, dummy;

        if (has_velocities)
        {
          // 13-column format: k, CDM, baryon, photon, nu, mass_nu, total, no_nu, total_de, Weyl, v_cdm, v_b, v_b-v_cdm
          ss >> Tk;    // k
          ss >> Tdc;   // cdm
          ss >> Tdb;   // baryon
          ss >> dummy; // photon
          ss >> dummy; // nu
          ss >> Tdn;   // mass_nu
          ss >> Tdm;   // total
          ss >> dummy; // no_nu
          ss >> dummy; // total_de
          ss >> dummy; // Weyl
          ss >> Tvc;   // v_cdm
          ss >> Tvb;   // v_b
          ss >> dummy; // v_b-v_cdm
        }
        else
        {
          // 9-column format: k, CDM, baryon, photon, nu, mass_nu, total, no_nu, total_de
          ss >> Tk;    // k
          ss >> Tdc;   // cdm
          ss >> Tdb;   // baryon
          ss >> dummy; // photon
          ss >> dummy; // nu
          ss >> Tdn;   // mass_nu
          ss >> Tdm;   // total
          ss >> dummy; // no_nu
          ss >> dummy; // total_de
          
          // Approximate velocities from density (simple approximation: theta ≈ -delta)
          Tvc = Tdc;
          Tvb = Tdb;
        }

        if (ss.bad() || ss.fail())
        {
          music::elog.Print("error reading the transfer function file (corrupt or not in expected format)!");
          throw std::runtime_error("error reading transfer function file \'" + filename + "\'");
        }

        // if (cosmo_params_["Omega_b"] < 1e-6)
        //   Tkvtot = Tktot;
        // else
        //   Tkvtot = cosmo_params_["f_c"] * Tkvc + cosmo_params_["f_b"]* Tkvb; 

        // m_linbaryoninterp |= Tkb < 0.0 || Tkvb < 0.0;

        k.push_back(Tk);
        dc.push_back(Tdc);
        db.push_back(Tdb);
        dn.push_back(Tdn);
        dm.push_back(Tdm);
        tc.push_back(Tvc);
        tb.push_back(Tvb);
        tn.push_back(Tdm);
        tm.push_back(Tdm);
        ++nlines;  
      }

      ifs.close();
      music::ilog.Print("Read CAMB transfer function table with %zu rows", nlines);

      // if (m_linbaryoninterp)
      //   music::ilog.Print("Using log-lin interpolation for baryons\n    (TF is not positive definite)");
    }

#if defined(USE_MPI)
    unsigned n = k.size();
    MPI_Bcast(&n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (CONFIG::MPI_task_rank > 0)
    {
      k.assign(n, 0);
      dc.assign(n, 0);
      tc.assign(n, 0);
      db.assign(n, 0);
      tb.assign(n, 0);
      dn.assign(n, 0);
      tn.assign(n, 0);
      dm.assign(n, 0);
      tm.assign(n, 0);
    }

    MPI_Bcast(&k[0],  n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dc[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tc[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&db[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tb[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dn[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tn[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dm[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tm[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    delta_c_.set_data(k, dc);
    theta_c_.set_data(k, tc);
    delta_b_.set_data(k, db);
    theta_b_.set_data(k, tb);
    delta_n_.set_data(k, dn);
    theta_n_.set_data(k, tn);
    delta_m_.set_data(k, dm);
    theta_m_.set_data(k, tm);

    // do not use first and last value since interpolation becomes lower order
    m_kmin = k[1];
    m_kmax = k[k.size()-2];

    return has_velocities;
  }

public:
  transfer_CAMB_file_plugin(config_file &cf, const cosmology::parameters& cosmo_params)
      : TransferFunction_plugin(cf, cosmo_params)
  {
    bool has_velocities = false;

    music::wlog << "The CAMB file plugin is not well tested! Proceed with checks of correctness of output before running a simulation!" << std::endl;

    std::string filename = pcf_->get_value<std::string>("cosmology", "transfer_file");

    has_velocities = this->read_table( filename );

    // Check if transfer function k-range is sufficient for simulation
    const double gridres = pcf_->get_value<double>("setup", "GridRes");
    const double boxlength = pcf_->get_value<double>("setup", "BoxLength");
    const double k_nyquist = M_PI * gridres / boxlength * std::sqrt(3.0) * 1.2; // 1.2 safety * space diagonal Nyquist
    
    if (m_kmax < k_nyquist)
    {
      music::elog << "Transfer function file k-range is insufficient!" << std::endl;
      music::elog << "  Required k_max (Nyquist): " << k_nyquist << " h/Mpc" << std::endl;
      music::elog << "  Available k_max in file:  " << m_kmax << " h/Mpc" << std::endl;
      music::elog << "  Grid resolution: " << gridres << ", Box length: " << boxlength << " Mpc/h" << std::endl;
      music::elog << "Please provide a transfer function file with larger k-range." << std::endl;
      throw std::runtime_error("Transfer function k-range insufficient for simulation parameters");
    }
    
    music::ilog << "Transfer function k-range check: OK" << std::endl;
    music::ilog << "  k_min = " << m_kmin << " h/Mpc, k_max = " << m_kmax << " h/Mpc" << std::endl;
    music::ilog << "  Required k_Nyquist = " << k_nyquist << " h/Mpc" << std::endl;

    // set properties of this transfer function plugin:
    tf_distinct_ = true; // different density between CDM v.s. Baryon
    tf_withvel_ = has_velocities;  // using velocity transfer function
    tf_withtotal0_ = false; // only have 1 file for the moment
  }

  ~transfer_CAMB_file_plugin()
  {
  }

  //!< return log-log interpolated values for transfer funtion 'type'
  inline double compute(double k, tf_type type) const
  {
    // use constant interpolation on the left side of the tabulated values, i.e.
    // set values k<k_min to value at k_min! (since transfer functions asymptote to constant)
    // also clamp to k_max to prevent interpolation errors
    k = std::max(k,m_kmin);
    k = std::min(k,m_kmax);

    switch (type)
    {
    case delta_matter0:
    case delta_matter:  
      return delta_m_(k);

    case delta_cdm0:
    case delta_cdm:
      return delta_c_(k);

    case delta_baryon0:
    case delta_baryon:
      return delta_b_(k);

    case theta_matter0:
    case theta_matter:
      return theta_m_(k);

    case theta_cdm0:
    case theta_cdm:
      return theta_c_(k);

    case theta_baryon0:
    case theta_baryon:
      return theta_b_(k);

    case delta_bc:
      return delta_b_(k)-delta_c_(k);
    
    case theta_bc:
      return theta_b_(k)-theta_c_(k);

    default:
      throw std::runtime_error("Invalid type requested in transfer function evaluation");
    }
  }

  //!< Return minimum k for which we can interpolate
  inline double get_kmin(void) const { return m_kmin; }

  //!< Return maximum k for which we can interpolate
  inline double get_kmax(void) const { return m_kmax; }
};

namespace
{
  TransferFunction_plugin_creator_concrete<transfer_CAMB_file_plugin> creator("CAMB_file");
}
