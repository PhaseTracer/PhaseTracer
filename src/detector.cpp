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

#include <cmath>
#include <stdexcept>
#include <utility>
#include <iomanip>

#include "detector.hpp"

namespace PhaseTracer {

std::ostream &operator<<(std::ostream &o, const Detector &d)
{
	o << d.name() << " (threshold = " << d.threshold() << ", T_obs = " << d.T_obs()
		<< " yr, band = [" << d.f_min() << ", " << d.f_max() << "] Hz)";
	return o;
}

void check_detector(const Detector &d)
{
	auto require = [&d](bool condition, const std::string &what) {
		if (!condition) throw std::invalid_argument("Detector " + d.name() + ": " + what);
	};
	require(!d.name().empty(), "name must not be empty");
	require(d.threshold() > 0., "threshold must be > 0");
	require(d.T_obs() > 0., "T_obs must be > 0");
	require(d.f_min() > 0. && d.f_min() < d.f_max(), "f_min must be in (0, f_max)");
	require(d.omega_noise(std::sqrt(d.f_min() * d.f_max())) > 0., "noise must be > 0 inside its band");
}

std::ostream &operator<<(std::ostream &o, const SNRResult &r)
{
	o<< r.detector << " = " << std::setprecision(6) << r.snr << " (threshold " << r.threshold << ")";
	return o;
}

} // namespace PhaseTracer
