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

#ifndef EFFECTIVEPOTENTIAL_PROPERTY_HPP_
#define EFFECTIVEPOTENTIAL_PROPERTY_HPP_
/**
   Convenience macro for adding getters and setters to our objects
*/

#define PROPERTY(type, name, default) \
 public: \
  void set_##name(type _##name) {name = _##name;} \
  type get_##name() const {return name;} \
 private: \
  type name = default;

#define PROTECTED_PROPERTY(type, name, default) \
 public: \
  void set_##name(type _##name) {name = _##name;} \
  type get_##name() const {return name;} \
 protected: \
  type name = default;

#define PROPERTY_CUSTOM_SETTER(type, name, default) \
 public: \
  type get_##name() const {return name;} \
 private: \
  type name = default; \

#define PROTECTED_PROPERTY_CUSTOM_SETTER(type, name, default) \
 public: \
  type get_##name() const {return name;} \
 protected: \
  type name = default;

#endif  // EFFECTIVEPOTENTIAL_PROPERTY_HPP_
