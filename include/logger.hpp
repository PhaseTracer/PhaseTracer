// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef PHASETRACER_LOGGER_HPP_
#define PHASETRACER_LOGGER_HPP_
/**
   Convenience macros for logging and setting log level
*/

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>
#define LOG BOOST_LOG_TRIVIAL
#include <boost/log/expressions.hpp>

#define LOGGER(level)                  \
  boost::log::core::get()->set_filter( \
      boost::log::trivial::severity >= \
      boost::log::trivial::level);

#endif // PHASETRACER_LOGGER_HPP_
