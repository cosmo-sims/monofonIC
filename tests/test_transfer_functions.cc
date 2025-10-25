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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <general.hh>
#include <config_file.hh>
#include <cosmology_parameters.hh>
#include <transfer_function_plugin.hh>

// Define CONFIG namespace globals (normally defined in main.cc)
namespace CONFIG {
    int MPI_thread_support = -1;
    int MPI_task_rank = 0;
    int MPI_task_size = 1;
    bool MPI_ok = false;
    bool MPI_threads_ok = false;
}

// Define memory tracking globals (normally defined in main.cc)
size_t global_mem_high_mark = 0;
size_t local_mem_high_mark = 0;

/**
 * @brief Regression test program for transfer function plugins
 *
 * This program evaluates transfer functions from different plugins at a range
 * of k values and either generates reference files or compares against them.
 *
 * Usage:
 *   test_transfer_functions --generate <plugin_name> <config_file> <output_file>
 *   test_transfer_functions --test <plugin_name> <config_file> <reference_file>
 */

// Transfer function types to test
const std::vector<tf_type> TF_TYPES = {
    tf_type::delta_matter,
    tf_type::delta_cdm,
    tf_type::delta_baryon,
    tf_type::theta_matter,
    tf_type::theta_cdm,
    tf_type::theta_baryon
};

const std::vector<std::string> TF_NAMES = {
    "delta_matter",
    "delta_cdm",
    "delta_baryon",
    "theta_matter",
    "theta_cdm",
    "theta_baryon"
};

// k range for sampling: 1e-4 to 10 h/Mpc, logarithmically spaced
constexpr int N_K_SAMPLES = 100;
constexpr double K_MIN = 1.0e-4;
constexpr double K_MAX = 10.0;

/**
 * @brief Generate logarithmically-spaced k values
 */
std::vector<double> generate_k_values()
{
    std::vector<double> k_values;
    k_values.reserve(N_K_SAMPLES);

    const double log_k_min = std::log10(K_MIN);
    const double log_k_max = std::log10(K_MAX);
    const double delta_log_k = (log_k_max - log_k_min) / (N_K_SAMPLES - 1);

    for (int i = 0; i < N_K_SAMPLES; ++i) {
        const double log_k = log_k_min + i * delta_log_k;
        k_values.push_back(std::pow(10.0, log_k));
    }

    return k_values;
}

/**
 * @brief Evaluate transfer functions at all k values
 */
void evaluate_transfer_functions(
    const TransferFunction_plugin& tf_plugin,
    const std::vector<double>& k_values,
    std::vector<std::vector<double>>& tf_values)
{
    const size_t n_types = TF_TYPES.size();
    tf_values.resize(n_types);

    for (size_t i = 0; i < n_types; ++i) {
        tf_values[i].reserve(k_values.size());
    }

    for (const double k : k_values) {
        // Check if k is within valid range
        if (k < tf_plugin.get_kmin() || k > tf_plugin.get_kmax()) {
            // Outside valid range, use NaN
            for (size_t i = 0; i < n_types; ++i) {
                tf_values[i].push_back(std::nan(""));
            }
            continue;
        }

        for (size_t i = 0; i < n_types; ++i) {
            try {
                const double value = tf_plugin.compute(k, TF_TYPES[i]);
                tf_values[i].push_back(value);
            } catch (...) {
                // If computation fails, use NaN
                tf_values[i].push_back(std::nan(""));
            }
        }
    }
}

/**
 * @brief Write transfer function values to a text file
 */
bool write_transfer_file(
    const std::string& filename,
    const std::vector<double>& k_values,
    const std::vector<std::vector<double>>& tf_values)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open output file: " << filename << std::endl;
        return false;
    }

    // Write header
    file << "# Transfer function regression test reference data\n";
    file << "# Columns: k[h/Mpc]";
    for (const auto& name : TF_NAMES) {
        file << " " << name;
    }
    file << "\n";

    // Write data with high precision
    file << std::scientific << std::setprecision(15);

    for (size_t j = 0; j < k_values.size(); ++j) {
        file << k_values[j];
        for (size_t i = 0; i < tf_values.size(); ++i) {
            file << " " << tf_values[i][j];
        }
        file << "\n";
    }

    file.close();
    return true;
}

/**
 * @brief Read transfer function values from a text file
 */
bool read_transfer_file(
    const std::string& filename,
    std::vector<double>& k_values,
    std::vector<std::vector<double>>& tf_values)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open reference file: " << filename << std::endl;
        return false;
    }

    k_values.clear();
    tf_values.clear();
    tf_values.resize(TF_TYPES.size());

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        double k;
        iss >> k;
        k_values.push_back(k);

        for (size_t i = 0; i < TF_TYPES.size(); ++i) {
            double value;
            iss >> value;
            tf_values[i].push_back(value);
        }
    }

    file.close();
    return true;
}

/**
 * @brief Compare two sets of transfer function values
 */
bool compare_transfer_functions(
    const std::vector<double>& k_values,
    const std::vector<std::vector<double>>& tf_ref,
    const std::vector<std::vector<double>>& tf_test,
    double rtol)
{
    bool all_passed = true;
    int total_failures = 0;

    std::cout << "\nComparing transfer functions with rtol = " << rtol << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (size_t i = 0; i < TF_TYPES.size(); ++i) {
        int n_failures = 0;
        double max_rel_diff = 0.0;
        double max_rel_diff_k = 0.0;

        for (size_t j = 0; j < k_values.size(); ++j) {
            const double ref = tf_ref[i][j];
            const double test = tf_test[i][j];

            // Skip NaN values (both must be NaN or both must be valid)
            if (std::isnan(ref) && std::isnan(test)) {
                continue;
            }
            if (std::isnan(ref) || std::isnan(test)) {
                std::cout << "ERROR: NaN mismatch at k = " << k_values[j]
                          << " for " << TF_NAMES[i] << "\n";
                n_failures++;
                continue;
            }

            // Compute relative difference
            const double abs_ref = std::abs(ref);
            const double abs_diff = std::abs(test - ref);
            const double rel_diff = (abs_ref > 1e-30) ? abs_diff / abs_ref : abs_diff;

            if (rel_diff > max_rel_diff) {
                max_rel_diff = rel_diff;
                max_rel_diff_k = k_values[j];
            }

            if (rel_diff > rtol) {
                n_failures++;
                if (n_failures <= 5) {  // Only print first 5 failures per type
                    std::cout << "FAIL: " << TF_NAMES[i]
                              << " at k = " << std::scientific << k_values[j]
                              << ", ref = " << ref
                              << ", test = " << test
                              << ", rel_diff = " << rel_diff << "\n";
                }
            }
        }

        if (n_failures > 0) {
            std::cout << "FAILED: " << TF_NAMES[i]
                      << " (" << n_failures << "/" << k_values.size() << " points failed)\n";
            std::cout << "        Max rel_diff = " << std::scientific << max_rel_diff
                      << " at k = " << max_rel_diff_k << "\n";
            all_passed = false;
            total_failures += n_failures;
        } else {
            std::cout << "PASSED: " << TF_NAMES[i]
                      << " (max rel_diff = " << std::scientific << max_rel_diff << ")\n";
        }
    }

    std::cout << std::string(80, '-') << "\n";

    if (all_passed) {
        std::cout << "SUCCESS: All transfer function comparisons passed!\n";
    } else {
        std::cout << "FAILURE: " << total_failures << " total comparison failures\n";
    }

    return all_passed;
}

void print_usage(const char* prog_name)
{
    std::cout << "Transfer Function Regression Test\n\n";
    std::cout << "Usage:\n";
    std::cout << "  Generate reference:\n";
    std::cout << "    " << prog_name << " --generate <plugin_name> <config_file> <output_file>\n\n";
    std::cout << "  Run test:\n";
    std::cout << "    " << prog_name << " --test <plugin_name> <config_file> <reference_file> [rtol]\n\n";
    std::cout << "Arguments:\n";
    std::cout << "  plugin_name    : Transfer function plugin (CLASS, eisenstein, etc.)\n";
    std::cout << "  config_file    : Configuration file with cosmology parameters\n";
    std::cout << "  output_file    : Output file for reference data (--generate mode)\n";
    std::cout << "  reference_file : Reference file to compare against (--test mode)\n";
    std::cout << "  rtol           : Relative tolerance (default: 1e-6)\n";
}

int main(int argc, char** argv)
{
    if (argc < 5) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string mode(argv[1]);
    const std::string plugin_name(argv[2]);
    const std::string config_filename(argv[3]);
    const std::string data_file(argv[4]);

    const double rtol = (argc > 5 && mode == "--test") ? std::atof(argv[5]) : 1.0e-6;

    if (mode != "--generate" && mode != "--test") {
        std::cerr << "ERROR: Mode must be --generate or --test\n";
        print_usage(argv[0]);
        return 1;
    }

    std::cout << "Transfer Function Regression Test\n";
    std::cout << "==================================\n";
    std::cout << "Mode          : " << mode << "\n";
    std::cout << "Plugin        : " << plugin_name << "\n";
    std::cout << "Config file   : " << config_filename << "\n";
    std::cout << "Data file     : " << data_file << "\n";
    if (mode == "--test") {
        std::cout << "Tolerance     : " << rtol << "\n";
    }
    std::cout << "\n";

    try {
        // Load configuration file
        config_file cf(config_filename);

        // Override the transfer function plugin name in config
        cf.insert_value("cosmology", "transfer", plugin_name);

        // Initialize cosmology parameters
        cosmology::parameters cosmo_params(cf);

        // Create and initialize transfer function plugin
        auto tf_plugin = select_TransferFunction_plugin(cf, cosmo_params);
        if (!tf_plugin) {
            std::cerr << "ERROR: Failed to create transfer function plugin: "
                      << plugin_name << "\n";
            return 1;
        }

        std::cout << "Initializing transfer function plugin...\n";
        tf_plugin->intialise();

        std::cout << "Plugin initialized successfully\n";
        std::cout << "  k range: [" << tf_plugin->get_kmin() << ", "
                  << tf_plugin->get_kmax() << "] h/Mpc\n";
        std::cout << "  Distinct CDM/baryon: " << (tf_plugin->tf_is_distinct() ? "yes" : "no") << "\n";
        std::cout << "  Has velocities: " << (tf_plugin->tf_has_velocities() ? "yes" : "no") << "\n";

        // Generate k values
        const auto k_values = generate_k_values();
        std::cout << "\nEvaluating transfer functions at " << k_values.size()
                  << " k values...\n";

        // Evaluate transfer functions
        std::vector<std::vector<double>> tf_values;
        evaluate_transfer_functions(*tf_plugin, k_values, tf_values);

        if (mode == "--generate") {
            // Write reference file
            std::cout << "Writing reference file: " << data_file << "\n";
            if (!write_transfer_file(data_file, k_values, tf_values)) {
                return 1;
            }
            std::cout << "Reference file written successfully\n";
            return 0;
        }
        else {  // mode == "--test"
            // Read reference file
            std::cout << "Reading reference file: " << data_file << "\n";
            std::vector<double> k_ref;
            std::vector<std::vector<double>> tf_ref;
            if (!read_transfer_file(data_file, k_ref, tf_ref)) {
                return 1;
            }

            // Check that k values match
            if (k_ref.size() != k_values.size()) {
                std::cerr << "ERROR: k value count mismatch (reference: "
                          << k_ref.size() << ", test: " << k_values.size() << ")\n";
                return 1;
            }

            // Compare transfer functions
            const bool passed = compare_transfer_functions(k_values, tf_ref, tf_values, rtol);
            return passed ? 0 : 1;
        }

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
