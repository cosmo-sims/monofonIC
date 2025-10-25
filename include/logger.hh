// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn & Michael Michaux (this file)
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

#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <sstream>
#include <regex>
#include <type_traits>
#include <terminal_colors.hh>

namespace music {

// Helper function to strip ANSI color codes from strings
inline std::string strip_ansi_codes(const std::string& str) {
  static const std::regex ansi_regex("\033\\[[0-9;]*m");
  return std::regex_replace(str, ansi_regex, "");
}

enum log_level : int {
  off     = 0,
  fatal   = 1,
  error   = 2,
  warning = 3,
  info    = 4,
  debug   = 5
};

inline const std::string HRULE = std::string(colors::HEADER) + "───────────────────────────────────────────────────────────────────────────────" + colors::RESET;

class logger {
private:
  static log_level log_level_;
  static std::ofstream output_file_;

public:
  logger()  = default;
  ~logger() = default;

  static void set_level(const log_level &level);
  static log_level get_level();

  static void set_output(const std::string filename);
  static void unset_output();

  static std::ofstream &get_output();

  template <typename T> logger &operator<<(const T &item) {
    std::cout << item;
    if (output_file_.is_open()) {
      // Only strip ANSI codes from string types, pass everything else through
      // to preserve stream formatting (e.g., std::setw)
      using bare_type = typename std::decay<T>::type;
      if constexpr (std::is_same_v<bare_type, std::string> ||
                    std::is_same_v<bare_type, const char*> ||
                    std::is_same_v<bare_type, char*>) {
        std::string str_item = item;
        output_file_ << strip_ansi_codes(str_item);
      } else {
        output_file_ << item;
      }
    }
    return *this;
  }

  logger &operator<<(std::ostream &(*fp)(std::ostream &)) {
    std::cout << fp;
    if (output_file_.is_open()) {
      output_file_ << fp;
    }
    return *this;
  }
};

class log_stream {
private:
  logger &logger_;
  log_level stream_level_;
  std::string line_prefix_, line_postfix_;

  bool newline;

public:
  log_stream(logger &logger, const log_level &level)
    : logger_(logger), stream_level_(level), newline(true) {
    switch (stream_level_) {
      case log_level::fatal:
        line_prefix_ = colors::ERROR + std::string("Fatal : ");
        break;
      case log_level::error:
        line_prefix_ = colors::ERROR + std::string("Error : ");
        break;
      case log_level::warning:
        line_prefix_ = colors::WARNING + std::string("Warning : ");
        break;
      case log_level::info:
        line_prefix_ = "";
        break;
      case log_level::debug:
        line_prefix_ = colors::HIGHLIGHT + std::string("Debug : ") + colors::RESET;
        break;
      default:
        line_prefix_ = colors::RESET + std::string("");
        break;
    }
    line_postfix_ = colors::RESET;
  }
  ~log_stream() = default;

  inline std::string GetPrefix() const {
    return line_prefix_;
  }

  template <typename T> log_stream &operator<<(const T &item) {
    if (logger::get_level() >= stream_level_) {
      if (newline) {
        logger_ << line_prefix_;
        newline = false;
      }
      logger_ << item;
    }
    return *this;
  }

  log_stream &operator<<(std::ostream &(*fp)(std::ostream &)) {
    if (logger::get_level() >= stream_level_) {
      logger_ << fp;
      logger_ << line_postfix_;
      newline = true;
    }
    return *this;
  }

  inline void Print(const char *str, ...) __attribute__ ((format (printf, 2, 3))) {
    char out[1024];
    va_list argptr;
    va_start(argptr, str);
    vsprintf(out, str, argptr);
    va_end(argptr);
    std::string out_string = std::string(out);
    // out_string.erase(std::remove(out_string.begin(), out_string.end(), '\n'),
    //                  out_string.end());
    (*this) << out_string << std::endl;
  }
};

// global instantiations for different levels
extern logger glogger;
extern log_stream flog;
extern log_stream elog;
extern log_stream wlog;
extern log_stream ilog;
extern log_stream dlog;

} // namespace music
