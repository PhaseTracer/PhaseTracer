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

#ifndef PHASETRACER_SCALE_HPP_
#define PHASETRACER_SCALE_HPP_

namespace PhaseTracer {

class Scale 
{
public:
  static const Scale MEV;
  static const Scale GEV;
  static const Scale TEV;

  double operator()() const { return factor_; }
  const char* name() const { return name_; }

  bool operator==(const Scale& other) const { return factor_ == other.factor_; }
  bool operator!=(const Scale& other) const { return !(*this == other); }

private:
  constexpr Scale(double factor, const char* name) : factor_(factor), name_(name) {}
  double factor_;
  const char* name_;
};

inline constexpr Scale Scale::MEV{1e3, "MeV"};
inline constexpr Scale Scale::GEV{1.0, "GeV"};
inline constexpr Scale Scale::TEV{1e-3, "TeV"};

inline Scale scale = Scale::GEV;

} // namespace PhaseTracer

#endif // PHASETRACER_SCALE_HPP_
