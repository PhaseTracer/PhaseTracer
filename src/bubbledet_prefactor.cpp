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

#include <algorithm>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

#include <Eigen/Core>
#include <interpolation.h>
#include <nlohmann/json.hpp>

#include "bubbledet_prefactor.hpp"
#include "logger.hpp"

namespace PhaseTracer {

using json = nlohmann::json;
using EffectivePotential::ParticleSpec;
using EffectivePotential::ZeroModeType;

namespace {

/** Below this peak m^2 (GeV^2) along the bounce a field is treated as massless
 *  and dropped: it does not change across the wall and gives a degenerate
 *  determinant. */
constexpr double kMasslessThreshold = 1e-6;

/** Standard analytic thermal prefactor T^4 (S/2pi)^{3/2}, used as a fallback. */
double analytic_prefactor(double temperature, double action_on_T) {
  const double t4 = temperature * temperature * temperature * temperature;
  return t4 * std::pow(action_on_T / (2. * M_PI), 1.5);
}

const char *zero_mode_string(ZeroModeType z) {
  switch (z) {
  case ZeroModeType::Higgs:
    return "Higgs";
  case ZeroModeType::Goldstone:
    return "Goldstone";
  default:
    return "None";
  }
}

/**
 * Maps the 1D bounce coordinate s (Profile1D.Phi -- a physical field for a
 * single-field bounce, or the path arc-length for a multifield one) to the
 * field-space point phi(s) and the path tangent dphi/ds. For multifield bounces
 * this is built by splining phi_for_profile1D against the path coordinate; for
 * single-field bounces it is the identity s -> (s).
 */
class PathMap {
public:
  PathMap(const PhaseTracer::Profile1D &profile,
          const std::vector<Eigen::VectorXd> &phi_for_profile, int n_scalars)
      : n_(n_scalars), identity_(phi_for_profile.empty()) {
    if (identity_) {
      if (n_ != 1) {
        throw std::runtime_error("multifield bounce is missing its field-space path (phi_for_profile)");
      }
      return;
    }
    // Sort (coordinate, phi) pairs by coordinate and drop duplicates so ALGLIB
    // gets a strictly increasing abscissa.
    const int m = static_cast<int>(profile.Phi.size());
    std::vector<int> order(m);
    for (int i = 0; i < m; ++i) {
      order[i] = i;
    }
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return profile.Phi[a] < profile.Phi[b]; });

    std::vector<double> coord;
    std::vector<std::vector<double>> comp(n_);
    coord.reserve(m);
    for (int c = 0; c < n_; ++c) {
      comp[c].reserve(m);
    }
    double last = std::numeric_limits<double>::quiet_NaN();
    for (int k = 0; k < m; ++k) {
      const int i = order[k];
      const double s = profile.Phi[i];
      if (k > 0 && std::abs(s - last) < 1e-12) {
        continue; // skip duplicate coordinate
      }
      last = s;
      coord.push_back(s);
      for (int c = 0; c < n_; ++c) {
        comp[c].push_back(phi_for_profile[i](c));
      }
    }
    if (coord.size() < 2) {
      throw std::runtime_error("degenerate tunneling path");
    }
    splines_.resize(n_);
    alglib::real_1d_array x;
    x.setcontent(coord.size(), coord.data());
    for (int c = 0; c < n_; ++c) {
      alglib::real_1d_array y;
      y.setcontent(comp[c].size(), comp[c].data());
      alglib::spline1dbuildcubic(x, y, splines_[c]);
    }
  }

  Eigen::VectorXd at(double s) const {
    Eigen::VectorXd v(n_);
    if (identity_) {
      v(0) = s;
      return v;
    }
    for (int c = 0; c < n_; ++c) {
      v(c) = alglib::spline1dcalc(splines_[c], s);
    }
    return v;
  }

  Eigen::VectorXd tangent(double s) const {
    Eigen::VectorXd t(n_);
    if (identity_) {
      t(0) = 1.0;
      return t;
    }
    for (int c = 0; c < n_; ++c) {
      double val, d, dd;
      alglib::spline1ddiff(splines_[c], s, val, d, dd);
      t(c) = d;
    }
    return t;
  }

private:
  int n_;
  bool identity_;
  std::vector<alglib::spline1dinterpolant> splines_;
};

/**
 * Minimal persistent child process with line-oriented bidirectional pipes.
 * The parent writes one request line to the child's stdin and reads one
 * response line from its stdout. Closing stdin tells the worker to exit.
 */
class PipeProcess {
public:
  PipeProcess(const std::string &exe, const std::string &script) {
    // Writing to a worker that has died (e.g. exec failed, or it crashed) must
    // not kill us with SIGPIPE -- we want fwrite/fflush to fail so we can fall
    // back to the analytic prefactor. Ignore it process-wide (idempotent).
    std::signal(SIGPIPE, SIG_IGN);

    int in_pipe[2];  // parent -> child stdin
    int out_pipe[2]; // child stdout -> parent
    if (pipe(in_pipe) != 0 || pipe(out_pipe) != 0) {
      return;
    }
    const pid_t pid = fork();
    if (pid < 0) {
      return;
    }
    if (pid == 0) {
      // Child: wire pipes to stdio and exec the worker.
      dup2(in_pipe[0], STDIN_FILENO);
      dup2(out_pipe[1], STDOUT_FILENO);
      close(in_pipe[0]);
      close(in_pipe[1]);
      close(out_pipe[0]);
      close(out_pipe[1]);
      execlp(exe.c_str(), exe.c_str(), script.c_str(), static_cast<char *>(nullptr));
      _exit(127); // exec failed
    }
    // Parent.
    close(in_pipe[0]);
    close(out_pipe[1]);
    pid_ = pid;
    to_child_ = fdopen(in_pipe[1], "w");
    from_child_ = fdopen(out_pipe[0], "r");
  }

  ~PipeProcess() { shutdown(); }

  PipeProcess(const PipeProcess &) = delete;
  PipeProcess &operator=(const PipeProcess &) = delete;

  bool ok() const { return pid_ > 0 && to_child_ != nullptr && from_child_ != nullptr; }

  bool write_line(const std::string &s) {
    if (to_child_ == nullptr) {
      return false;
    }
    if (fwrite(s.data(), 1, s.size(), to_child_) != s.size()) {
      return false;
    }
    return fputc('\n', to_child_) != EOF && fflush(to_child_) == 0;
  }

  bool read_line(std::string &line) {
    if (from_child_ == nullptr) {
      return false;
    }
    line.clear();
    int c;
    while ((c = fgetc(from_child_)) != EOF) {
      if (c == '\n') {
        return true;
      }
      line.push_back(static_cast<char>(c));
    }
    return !line.empty(); // last line may lack a trailing newline
  }

  void shutdown() {
    if (to_child_ != nullptr) {
      fclose(to_child_); // EOF on child's stdin -> worker exits its read loop
      to_child_ = nullptr;
    }
    if (from_child_ != nullptr) {
      fclose(from_child_);
      from_child_ = nullptr;
    }
    if (pid_ > 0) {
      int status = 0;
      waitpid(pid_, &status, 0);
      pid_ = -1;
    }
  }

private:
  pid_t pid_ = -1;
  FILE *to_child_ = nullptr;
  FILE *from_child_ = nullptr;
};

} // namespace

// ----------------------------------------------------------------------------
//  Implementation (held behind a shared_ptr so BubbleDetPrefactor is copyable
//  and storable in a std::function).
// ----------------------------------------------------------------------------
class BubbleDetPrefactor::Impl {
public:
  Impl(const EffectivePotential::Potential &potential, BubbleDetConfig config)
      : potential_(potential), config_(std::move(config)) {
    if (config_.worker_script.empty()) {
      LOG(warning) << "BubbleDetPrefactor: no worker_script configured; will always "
                      "fall back to the analytic prefactor.";
      return;
    }
    proc_ = std::make_unique<PipeProcess>(config_.python_executable, config_.worker_script);
    if (!proc_->ok()) {
      LOG(warning) << "BubbleDetPrefactor: failed to start worker '"
                   << config_.python_executable << " " << config_.worker_script
                   << "'; will fall back to the analytic prefactor.";
      proc_.reset();
    }
  }

  double compute(double T, double action_on_T, const ActionResult &bounce) {
    try {
      const json request = build_request(T, action_on_T, bounce);

      std::string response_line;
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!proc_ || !proc_->ok()) {
          throw std::runtime_error("BubbleDet worker is not running");
        }
        if (!proc_->write_line(request.dump())) {
          throw std::runtime_error("failed to send request to BubbleDet worker");
        }
        if (!proc_->read_line(response_line)) {
          throw std::runtime_error("no response from BubbleDet worker (it may have crashed)");
        }
      }

      const json response = json::parse(response_line);
      if (!response.value("ok", false)) {
        throw std::runtime_error("BubbleDet worker: " + response.value("error", std::string("unknown error")));
      }
      const double S1 = response.at("S1").get<double>();
      const double A_bubbledet = std::exp(-S1); // A = exp(-S1)
      const double A_analytic = analytic_prefactor(T, action_on_T);
      LOG(debug) << "BubbleDetPrefactor at T=" << T << ", S/T=" << action_on_T
                 << ": A_BubbleDet=" << A_bubbledet << " (S1=" << S1 << "), "
                 << "A_analytic=" << A_analytic
                 << ", ratio=" << A_bubbledet / A_analytic;
      return A_bubbledet;
    } catch (const std::exception &e) {
      LOG(warning) << "BubbleDetPrefactor falling back to analytic prefactor at T=" << T
                   << ": " << e.what();
      return analytic_prefactor(T, action_on_T);
    }
  }

private:
  json build_request(double T, double /*action_on_T*/, const ActionResult &bounce) const {
    const int n = static_cast<int>(potential_.get_n_scalars());

    const Profile1D &prof = bounce.bubble_profile;
    const int m = static_cast<int>(prof.R.size());
    if (m < 2) {
      throw std::runtime_error("bounce profile has too few points");
    }

    // BubbleDet evaluates the fluctuation operator at r=0, but PhaseTracer's
    // grid starts at a small finite rmin. Prepend the origin using spherical
    // regularity: dPhi(0)=0 and Phi(0)=Phi(rmin) to leading order. Here Phi is
    // the 1D bounce coordinate (a field for n=1, the path arc-length for n>1).
    std::vector<double> R, Phi, dPhi;
    R.reserve(m + 1);
    Phi.reserve(m + 1);
    dPhi.reserve(m + 1);
    if (prof.R[0] > 0.0) {
      R.push_back(0.0);
      Phi.push_back(prof.Phi[0]);
      dPhi.push_back(0.0);
    }
    for (int i = 0; i < m; ++i) {
      R.push_back(prof.R[i]);
      Phi.push_back(prof.Phi[i]);
      dPhi.push_back(prof.dPhi[i]);
    }

    const double phi_metaMin = Phi.back(); // false vacuum (large r)

    // Map the bounce coordinate to the field-space point (and tangent).
    PathMap path(prof, bounce.phi_for_profile, n);

    // Sampling grid over the bounce coordinate, spanning the profile range plus
    // a margin (all of BubbleDet's evaluations stay within it).
    double lo = *std::min_element(Phi.begin(), Phi.end());
    double hi = *std::max_element(Phi.begin(), Phi.end());
    const double span = hi - lo;
    const double margin = config_.phi_grid_margin * (span > 0. ? span : 1.0);
    lo -= margin;
    hi += margin;
    const int G = std::max(config_.phi_grid_points, 4);

    std::vector<double> phi_grid(G), Vv(G);
    std::vector<Eigen::VectorXd> field(G);
    for (int i = 0; i < G; ++i) {
      const double s = lo + (hi - lo) * i / (G - 1);
      phi_grid[i] = s;
      field[i] = path.at(s);
      Vv[i] = potential_.V(field[i], T);
    }

    json request;
    request["T"] = T;
    request["dim"] = config_.dim;
    request["thermal"] = config_.thermal;
    request["phi_metaMin"] = phi_metaMin;
    request["R"] = R;
    request["Phi"] = Phi;
    request["dPhi"] = dPhi;
    request["phi_grid"] = phi_grid;
    request["V"] = Vv;

    if (n == 1) {
      // Single field: the bounce coordinate is the physical field, so exact
      // field-space derivatives are the coordinate derivatives. For n>1 we omit
      // them and let BubbleDet finite-difference V(coord), which gives the
      // correct longitudinal (path) derivatives.
      std::vector<double> dVv(G), d2Vv(G);
      for (int i = 0; i < G; ++i) {
        dVv[i] = potential_.dV_dx(field[i], T)(0);
        d2Vv[i] = potential_.d2V_dx2(field[i], T)(0, 0);
      }
      request["dV"] = dVv;
      request["d2V"] = d2Vv;
    }

    request["particles"] = build_particles(T, phi_grid, field, path);
    return request;
  }

  json build_particles(double T, const std::vector<double> &phi_grid,
                       const std::vector<Eigen::VectorXd> &field, const PathMap &path) const {
    const int n = static_cast<int>(potential_.get_n_scalars());
    const int G = static_cast<int>(phi_grid.size());
    json particles = json::array();

    // 1. Longitudinal (Higgs) mode: W = d2V along the path. The worker reads
    //    this off BubbleConfig.d2V, so it stays consistent with the V used for
    //    the bounce (exact for n=1, finite-differenced from V(coord) for n>1).
    particles.push_back({
        {"name", "Higgs_longitudinal"},
        {"spin", 0},
        {"dof", config_.higgs_dof},
        {"zero_modes", "Higgs"},
        {"W_source", "d2V"},
    });

    // 2. Transverse CP-even scalar modes: W = n_hat^T (d2V/dphi^2) n_hat, the
    //    effective-potential Hessian projected onto the path-normal direction.
    const int n_trans = config_.transverse_scalar_modes < 0 ? (n - 1) : config_.transverse_scalar_modes;
    if (n_trans > 0) {
      if (n != 2) {
        throw std::runtime_error(
            "transverse CP-even projection is only implemented for 2 fields; set transverse_scalar_modes=0");
      }
      std::vector<double> Wt(G);
      double wmax = 0.;
      for (int i = 0; i < G; ++i) {
        Eigen::VectorXd t = path.tangent(phi_grid[i]);
        const double tn = t.norm();
        if (tn > 0.) {
          t /= tn;
        }
        Eigen::VectorXd nhat(2);
        nhat << -t(1), t(0); // unit normal to the path in field space
        const Eigen::MatrixXd H = potential_.d2V_dx2(field[i], T);
        Wt[i] = nhat.dot(H * nhat);
        wmax = std::max(wmax, std::abs(Wt[i]));
      }
      if (wmax >= kMasslessThreshold) {
        particles.push_back({
            {"name", "transverse_cp_even"},
            {"spin", 0},
            {"dof", config_.transverse_dof},
            {"zero_modes", "None"},
            {"W", Wt},
        });
      }
    }

    // 3. Remaining fields (gauge bosons, Goldstones, ...) from the spectrum,
    //    evaluated as m^2(phi(s)) along the path. Fermions (spin 1/2) and any
    //    Higgs-tagged entry are skipped -- the functor owns the longitudinal
    //    Higgs, and BubbleDet's ParticleConfig accepts only spin 0 or 1.
    const std::vector<ParticleSpec> spectrum =
        config_.spectrum_provider ? config_.spectrum_provider(T) : potential_.get_fluctuation_spectrum(T);
    for (const auto &s : spectrum) {
      if (!s.mass_sq) {
        LOG(warning) << "BubbleDetPrefactor: particle '" << s.name << "' has no mass function; skipping.";
        continue;
      }
      if (s.spin != 0.0 && s.spin != 1.0) {
        LOG(warning) << "BubbleDetPrefactor: skipping spin-" << s.spin << " particle '" << s.name
                     << "' (BubbleDet supports spin 0 and 1 only).";
        continue;
      }
      if (s.zero_mode == ZeroModeType::Higgs) {
        LOG(warning) << "BubbleDetPrefactor: ignoring spectrum particle '" << s.name
                     << "' tagged Higgs; the longitudinal Higgs is supplied internally.";
        continue;
      }
      std::vector<double> W(G);
      double wmax = 0.;
      for (int i = 0; i < G; ++i) {
        W[i] = s.mass_sq(field[i]);
        wmax = std::max(wmax, std::abs(W[i]));
      }
      if (wmax < kMasslessThreshold) {
        // A field that stays massless along the bounce (e.g. the photon) does
        // not change across the wall and gives a degenerate determinant that
        // BubbleDet cannot extrapolate -- drop it.
        LOG(debug) << "BubbleDetPrefactor: skipping massless particle '" << s.name << "'.";
        continue;
      }
      particles.push_back({
          {"name", s.name},
          {"spin", static_cast<int>(s.spin)},
          {"dof", s.dof},
          {"zero_modes", zero_mode_string(s.zero_mode)},
          {"W", W},
      });
    }
    return particles;
  }

  const EffectivePotential::Potential &potential_;
  BubbleDetConfig config_;
  std::unique_ptr<PipeProcess> proc_;
  std::mutex mutex_;
};

// ----------------------------------------------------------------------------
//  Public surface
// ----------------------------------------------------------------------------
BubbleDetPrefactor::BubbleDetPrefactor(const EffectivePotential::Potential &potential, BubbleDetConfig config)
    : impl_(std::make_shared<Impl>(potential, std::move(config))) {}

double BubbleDetPrefactor::operator()(double temperature, double action_on_T, const ActionResult &bounce) const {
  return impl_->compute(temperature, action_on_T, bounce);
}

} // namespace PhaseTracer
