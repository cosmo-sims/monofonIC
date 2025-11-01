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

inline bool use_color = (std::getenv("NO_COLOR") == nullptr);

inline const char* maybe(const char* code) {
    return use_color ? code : "";
}

// Task status colors
inline const char* TASK_EXECUTING = maybe("\033[38;5;114m");  // Green (#a6e3a1) - Active tasks
inline const char* TASK_SKIPPED   = maybe("\033[38;5;243m");  // Gray (#6c7086) - Skipped tasks
inline const char* SUCCESS        = maybe("\033[38;5;115m");  // Mint (#94e2d5) - Success messages

// Information type colors
inline const char* HEADER         = maybe("\033[38;5;147m");  // Lavender (#b4befe) - Section headers
inline const char* TASK_NAME      = maybe("\033[1;38;5;117m"); // Bold Sky (#89dceb) - Task names
inline const char* TIMING         = maybe("\033[38;5;186m");  // Peach (#fab387) - Timing information
inline const char* CONFIG_VALUE   = maybe("\033[38;5;150m");  // Teal (#94e2d5) - Configuration values
inline const char* SPECIES        = maybe("\033[38;5;183m");  // Mauve (#cba6f7) - Species/physics info
inline const char* HIGHLIGHT      = maybe("\033[38;5;116m");  // Sky (#89dceb) - General highlights
inline const char* LOGO           = maybe("\033[38;5;147m");  // Lavender (#b4befe) - ASCII art logo

// Existing logger colors (preserved for compatibility)
inline const char* ERROR          = maybe("\033[31m");        // Red - Errors
inline const char* WARNING        = maybe("\033[33m");        // Yellow - Warnings

// Text formatting
inline const char* BOLD           = maybe("\033[1m");         // Bold text
inline const char* DIM            = maybe("\033[2m");         // Dimmed text
inline const char* RESET          = maybe("\033[0m");         // Reset to default

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
