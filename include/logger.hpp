#ifndef PHASETRACER_LOGGER_HPP_INCLUDED
#define PHASETRACER_LOGGER_HPP_INCLUDED
/**
   Convenience macros for logging and setting log level
*/

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>
#define LOG BOOST_LOG_TRIVIAL
#include <boost/log/expressions.hpp>

#define LOGGER(level) \
  boost::log::core::get()->set_filter( \
    boost::log::trivial::severity >= \
    boost::log::trivial::level);

#endif
