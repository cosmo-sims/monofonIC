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

#include <transfer_function_plugin.hh>

/**
 * @brief Get the TransferFunction plugin map object
 * 
 * @return std::map<std::string, TransferFunction_plugin_creator *>& 
 */
std::map<std::string, TransferFunction_plugin_creator *> &
get_TransferFunction_plugin_map()
{
    static std::map<std::string, TransferFunction_plugin_creator *> TransferFunction_plugin_map;
    return TransferFunction_plugin_map;
}

/**
 * @brief Print out the names of all transfer function plugins compiled in
 * 
 */
void print_TransferFunction_plugins()
{
    std::map<std::string, TransferFunction_plugin_creator *> &m = get_TransferFunction_plugin_map();
    std::map<std::string, TransferFunction_plugin_creator *>::iterator it;
    it = m.begin();
    music::ilog << "Available transfer function plug-ins:" << std::endl;
    while (it != m.end())
    {
        if ((*it).second)
            music::ilog << "\t\'" << (*it).first << "\'" << std::endl;
        ++it;
    }
    music::ilog << std::endl;
}

/**
 * @brief Return a pointer to the desired transfer function plugin as given in the config file
 * 
 * Implements the abstract factory pattern (https://en.wikipedia.org/wiki/Abstract_factory_pattern)
 * 
 * @param cf  reference to config_file object
 * @param cosmo_param reference to cosmology::parameters object holding cosmological parameter values
 * @return std::unique_ptr<TransferFunction_plugin> 
 */
std::unique_ptr<TransferFunction_plugin> select_TransferFunction_plugin(config_file &cf, const cosmology::parameters& cosmo_param)
{
    std::string tfname = cf.get_value<std::string>("cosmology", "transfer");

    TransferFunction_plugin_creator *the_TransferFunction_plugin_creator = get_TransferFunction_plugin_map()[tfname];

    if (!the_TransferFunction_plugin_creator)
    {
        music::elog << "Invalid/Unregistered transfer function plug-in encountered : " << tfname << std::endl;
        print_TransferFunction_plugins();
        throw std::runtime_error("Unknown transfer function plug-in");
    }
    else
    {
        music::ilog << music::HRULE << std::endl;
        music::ilog << std::setw(32) << std::left << "Transfer function plugin" << " : " << tfname << std::endl;
    }

    return the_TransferFunction_plugin_creator->create(cf, cosmo_param);
}
