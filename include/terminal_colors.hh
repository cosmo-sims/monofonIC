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
#pragma once

/**
 * @file terminal_colors.hh
 * @brief ANSI color codes and symbols for terminal output
 *
 * Uses the Catppuccin Mocha color palette for modern, pleasing aesthetics
 * with excellent readability on dark terminal backgrounds.
 *
 * Color scheme: https://github.com/catppuccin/catppuccin
 */

 #include <string>

namespace colors {

// ============================================================================
// Catppuccin Mocha Color Palette (256-color mode)
// ============================================================================

// Task status colors
constexpr const char* TASK_EXECUTING = "\033[38;5;114m";  // Green (#a6e3a1) - Active tasks
constexpr const char* TASK_SKIPPED   = "\033[38;5;243m";  // Gray (#6c7086) - Skipped tasks
constexpr const char* SUCCESS        = "\033[38;5;115m";  // Mint (#94e2d5) - Success messages

// Information type colors
constexpr const char* HEADER         = "\033[38;5;147m";  // Lavender (#b4befe) - Section headers
constexpr const char* TASK_NAME      = "\033[1;38;5;117m"; // Bold Sky (#89dceb) - Task names
constexpr const char* TIMING         = "\033[38;5;186m";  // Peach (#fab387) - Timing information
constexpr const char* CONFIG_VALUE   = "\033[38;5;150m";  // Teal (#94e2d5) - Configuration values
constexpr const char* SPECIES        = "\033[38;5;183m";  // Mauve (#cba6f7) - Species/physics info
constexpr const char* HIGHLIGHT      = "\033[38;5;116m";  // Sky (#89dceb) - General highlights
constexpr const char* LOGO           = "\033[38;5;147m";  // Lavender (#b4befe) - ASCII art logo

// Existing logger colors (preserved for compatibility)
constexpr const char* ERROR          = "\033[31m";        // Red - Errors
constexpr const char* WARNING        = "\033[33m";        // Yellow - Warnings

// Text formatting
constexpr const char* BOLD           = "\033[1m";         // Bold text
constexpr const char* DIM            = "\033[2m";         // Dimmed text
constexpr const char* RESET          = "\033[0m";         // Reset to default

// ============================================================================
// Decorative symbols (UTF-8)
// ============================================================================

constexpr const char* SYM_CHECK      = "▸";               // Task executing/completed
constexpr const char* SYM_SKIP       = "○";               // Task skipped
constexpr const char* SYM_DIAMOND    = "✦";               // Section headers
constexpr const char* SYM_ATOM       = "⚛";               // Physics/species
constexpr const char* SYM_ARROW      = "→";               // Progress/direction
constexpr const char* SYM_DOT        = "•";               // List items

// ============================================================================
// Convenience functions
// ============================================================================

/**
 * @brief Create a colored string with automatic reset
 */
inline std::string colored(const std::string& text, const char* color) {
    return std::string(color) + text + RESET;
}

/**
 * @brief Create a colored and bold string
 */
inline std::string colored_bold(const std::string& text, const char* color) {
    return std::string(BOLD) + color + text + RESET;
}

} // namespace colors
