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

#include "HDF_IO.hh"
#include <logger.hh>
#include <output_plugin.hh>

class generic_output_plugin : public output_plugin
{
private:
	std::string get_field_name( const cosmo_species &s, const fluid_component &c );
protected:
	bool out_eulerian_;
    unsigned hdf5_compression_level_;
public:
	//! constructor
	explicit generic_output_plugin(config_file &cf, std::unique_ptr<cosmology::calculator> &pcc )
	: output_plugin(cf, pcc, "Generic HDF5")
	{
		real_t astart   = 1.0/(1.0+cf_.get_value<double>("setup", "zstart"));
		real_t boxsize  = cf_.get_value<double>("setup", "BoxLength");
		real_t omegab   = pcc->cosmo_param_["Omega_b"];
		real_t omegam   = pcc->cosmo_param_["Omega_m"];
		real_t omegal   = pcc->cosmo_param_["Omega_DE"];
		

		out_eulerian_   = cf_.get_value_safe<bool>("output", "generic_out_eulerian",false);

        hdf5_compression_level_ = cf_.get_value_safe<unsigned>("output", "hdf5CompressionLevel", 0);

		if( CONFIG::MPI_task_rank == 0 )
		{
			unlink(fname_.c_str());
			HDFCreateFile( fname_ );
			HDFCreateGroup( fname_, "Header" );
			HDFWriteGroupAttribute<double>( fname_, "Header", "Boxsize", boxsize );
			HDFWriteGroupAttribute<double>( fname_, "Header", "astart", astart );
			HDFWriteGroupAttribute<double>( fname_, "Header", "Omega_b", omegab );
			HDFWriteGroupAttribute<double>( fname_, "Header", "Omega_m", omegam );
			HDFWriteGroupAttribute<double>( fname_, "Header", "Omega_L", omegal );
			
		}

#if defined(USE_MPI)
		MPI_Barrier( MPI_COMM_WORLD );
#endif
	}

    output_type write_species_as( const cosmo_species &s ) const
	{ 
		if( out_eulerian_ )
			return output_type::field_eulerian;
		return output_type::field_lagrangian;
	}

	bool has_64bit_reals() const{ return true; }

	bool has_64bit_ids() const{ return true; }

	real_t position_unit() const { return 1.0; }
	
	real_t velocity_unit() const { return 1.0; }

	real_t mass_unit() const { return 1.0; }
	
	void write_grid_data(const Grid_FFT<real_t> &g, const cosmo_species &s, const fluid_component &c );
};


std::string generic_output_plugin::get_field_name( const cosmo_species &s, const fluid_component &c )
{
	std::string field_name;
	switch( s ){
		case cosmo_species::dm: 
			field_name += "DM"; break;
		case cosmo_species::baryon: 
			field_name += "BA"; break;
		case cosmo_species::neutrino: 
			field_name += "NU"; break;
		default: break;
	}
	field_name += "_";
	switch( c ){
		case fluid_component::density:
			field_name += "delta"; break;
		case fluid_component::vx:
			field_name += "vx"; break;
		case fluid_component::vy:
			field_name += "vy"; break;
		case fluid_component::vz:
			field_name += "vz"; break;
		case fluid_component::dx:
			field_name += "dx"; break;
		case fluid_component::dy:
			field_name += "dy"; break;
		case fluid_component::dz:
			field_name += "dz"; break;
		default: break;
	}
	return field_name;
}

void generic_output_plugin::write_grid_data(const Grid_FFT<real_t> &g, const cosmo_species &s, const fluid_component &c ) 
{
	std::string field_name = this->get_field_name( s, c );
	g.Write_to_HDF5(fname_, field_name, hdf5_compression_level_);
	music::ilog << interface_name_ << " : Wrote field \'" << field_name << "\' to file \'" << fname_ << "\'" << std::endl;
}

namespace
{
   output_plugin_creator_concrete<generic_output_plugin> creator001("generic"); 
} // namespace

#endif
