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

#include "output_plugin.hh"

/**
 * @brief Get the output plugin map object
 * 
 * @return std::map< std::string, output_plugin_creator *>& 
 */
std::map< std::string, output_plugin_creator *>& get_output_plugin_map()
{
	static std::map< std::string, output_plugin_creator* > output_plugin_map;
	return output_plugin_map;
}

/**
 * @brief Print out the names of all output plugins compiled in
 * 
 */
void print_output_plugins()
{
	std::map< std::string, output_plugin_creator *>& m = get_output_plugin_map();
	
	std::map< std::string, output_plugin_creator *>::iterator it;
	it = m.begin();
	music::ilog << "Available output plug-ins:\n";
	while( it!=m.end() )
	{
		if( it->second )
			music::ilog << "\t\'" << it->first << "\'\n";
		++it;
	}
	music::ilog << std::endl;
}

/**
 * @brief Return a pointer to the desired output plugin as given in the config file
 * 
 * Implements the abstract factory pattern (https://en.wikipedia.org/wiki/Abstract_factory_pattern)
 * 
 * @param cf reference to config_file object
 * @param pcc reference to cosmology::calculator object
 * @return std::unique_ptr<output_plugin> 
 */
std::unique_ptr<output_plugin> select_output_plugin( config_file& cf, std::unique_ptr<cosmology::calculator>& pcc )
{
	std::string formatname = cf.get_value<std::string>( "output", "format" );
	
	output_plugin_creator *the_output_plugin_creator = get_output_plugin_map()[ formatname ];
	
	if( !the_output_plugin_creator )
	{	
		music::elog << "Output plug-in \'" << formatname << "\' not found." << std::endl;
		print_output_plugins();
		throw std::runtime_error("Unknown output plug-in");
		
	}else{
		music::ilog << music::HRULE << std::endl;
		music::ilog << std::setw(32) << std::left << "Output plugin" << " : " << formatname << std::endl;
	}
	
	return the_output_plugin_creator->create( cf, pcc );
}



