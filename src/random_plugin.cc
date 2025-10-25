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

#include <general.hh>
#include <random_plugin.hh>

/**
 * @brief Get the RNG plugin map object
 * 
 * @return std::map<std::string, RNG_plugin_creator *>& 
 */
std::map<std::string, RNG_plugin_creator *> &
get_RNG_plugin_map()
{
    static std::map<std::string, RNG_plugin_creator *> RNG_plugin_map;
    return RNG_plugin_map;
}

/**
 * @brief Print out the names of all RNG plugins compiled in
 * 
 */
void print_RNG_plugins()
{
    std::map<std::string, RNG_plugin_creator *> &m = get_RNG_plugin_map();
    std::map<std::string, RNG_plugin_creator *>::iterator it;
    it = m.begin();
    music::ilog << "Available random number generator plug-ins:" << std::endl;
    while (it != m.end())
    {
        if ((*it).second){
            music::ilog.Print("\t\'%s\'\n", (*it).first.c_str());
        }
        ++it;
    }
    music::ilog << std::endl;
}


/**
 * @brief Return a pointer to the desired random number generator plugin as given in the config file
 * 
 * Implements the abstract factory pattern (https://en.wikipedia.org/wiki/Abstract_factory_pattern)
 * 
 * @param cf reference to config_file object
 * @return std::unique_ptr<RNG_plugin> unique pointer to plugin
 */
std::unique_ptr<RNG_plugin> select_RNG_plugin(config_file &cf)
{
    std::string rngname = cf.get_value<std::string>("random", "generator");

    RNG_plugin_creator *the_RNG_plugin_creator = get_RNG_plugin_map()[rngname];

    if (!the_RNG_plugin_creator)
    {
        music::ilog.Print("Invalid/Unregistered random number generator plug-in encountered : %s", rngname.c_str());
        print_RNG_plugins();
        throw std::runtime_error("Unknown random number generator plug-in");
    }
    else
    {
        music::ilog << music::HRULE << std::endl;
        music::ilog << std::setw(32) << std::left << "Random number generator plugin" << " : " << colors::CONFIG_VALUE << rngname << colors::RESET << std::endl;
    }

    return the_RNG_plugin_creator->Create(cf);
}
