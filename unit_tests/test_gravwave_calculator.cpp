#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>

#include "catch/catch.hpp"
#include "models/1D_test_model.hpp"
#include "models/2D_test_model.hpp"
#include "models/Z2_scalar_singlet_model.hpp"
#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "action_calculator.hpp"
#include "thermo_finder.hpp"
#include "detector.hpp"
#include "gravwave_calculator.hpp"
#include "scale.hpp"
#include "logger.hpp"

TEST_CASE("Compute gravitational wave spectrum for a one-dimensional model", "[1DTestModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::OneDimModel model;

  const double rel_tol = 1.e-3;
  const double abs_tol = 1.e-3;

  // Make PhaseFinder object and find the transitions
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.set_find_min_x_tol_rel(1.e-8);
  pf.set_find_min_x_tol_abs(1.e-8);
  pf.find_phases();

  PhaseTracer::ActionCalculator ac(pf);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.set_TC_tol_rel(1e-16);
  tf.find_transitions();

  // PhaseTracer::GravWaveCalculator gc(tf);
  // gc.calc_spectrums();
  // auto spec = gc.get_spectrums();

  // SECTION("Check number of spectra")

  // REQUIRE(spec.size() == 1);

  // SECTION("Check value of spectra")

  // REQUIRE(spec[0].Tref == Approx(57.3603638993).epsilon(rel_tol));
  // REQUIRE(spec[0].alpha == Approx(0.0013946917).epsilon(rel_tol));
  // REQUIRE(spec[0].beta_H == Approx(6873.221784433).epsilon(rel_tol));
  // REQUIRE(spec[0].peak_frequency == Approx(0.1191659843).epsilon(rel_tol));

  // Catch::StringMaker<double>::precision = 30;
  // REQUIRE(spec[0].peak_amplitude == Approx(3.8981023e-23).epsilon(rel_tol));
}

TEST_CASE("Compute gravitational wave spectrum for a two-dimensional model", "[2DTestModel]") {

  LOGGER(fatal);

  // Construct our model
  EffectivePotential::TwoDimModel model;

  // Make PhaseFinder object and find the transitions
  PhaseTracer::PhaseFinder pf(model);
  pf.set_seed(1);
  pf.find_phases();

  PhaseTracer::ActionCalculator ac(pf);
  PhaseTracer::TransitionFinder tf(pf, ac);
  tf.find_transitions();
  auto transitions = tf.get_transitions();

  // PhaseTracer::GravWaveCalculator gc(tf);
  // gc.calc_spectrums();
  // auto spec = gc.get_spectrums();

  // const double rel_tol = 1.e-3;

  // SECTION("Check number of spectra")

  // REQUIRE(spec.size() == 1);

  // SECTION("Check value of spectra")

  // REQUIRE(spec[0].Tref == Approx(84.1901920318).epsilon(rel_tol));
  // REQUIRE(spec[0].alpha == Approx(0.0970871041).epsilon(rel_tol));
  // REQUIRE(spec[0].beta_H == Approx(1781.3022934312).epsilon(rel_tol));
  // REQUIRE(spec[0].peak_frequency == Approx(0.0452179242).epsilon(rel_tol));

  // Catch::StringMaker<double>::precision = 30;
  // REQUIRE(spec[0].peak_amplitude == Approx(1.21499917010577e-16).epsilon(rel_tol));
}
TEST_CASE("GW fit formulas do not depend on the energy unit", "[GravWaveCalculator]") {

  LOGGER(fatal);

  struct ScaleGuard {
    ~ScaleGuard() { PhaseTracer::scale = PhaseTracer::Scale::GEV; }
  } guard;

  EffectivePotential::OneDimModel model;
  PhaseTracer::PhaseFinder pf(model);
  PhaseTracer::ActionCalculator ac(pf);
  PhaseTracer::TransitionFinder tf(pf);
  PhaseTracer::ThermoFinder tm(tf, ac);
  PhaseTracer::GravWaveCalculator gw(tm);

  const double alpha = 0.1;
  const double beta_H = 100.;
  const double vw = 0.7;
  const double cs = 1. / std::sqrt(3.);
  const double g_eff = 100.;
  const double h_eff = 100.;

  for (const bool legacy : {true, false}) {
    gw.set_use_legacy_gw_methods(legacy);
    // 0.05 GeV is below the 10 GeV bubble-collision threshold and 50 GeV is above it
    for (const double T_GeV : {0.05, 50.}) {
      for (const double f : {1e-5, 1e-3, 1e-1}) {
        PhaseTracer::scale = PhaseTracer::Scale::GEV;
        const double omega_GeV = gw.fit_omega(f, alpha, beta_H, T_GeV, T_GeV, vw, cs, g_eff, h_eff);
        PhaseTracer::scale = PhaseTracer::Scale::MEV;
        const double T_MeV = 1e3 * T_GeV;
        const double omega_MeV = gw.fit_omega(f, alpha, beta_H, T_MeV, T_MeV, vw, cs, g_eff, h_eff);
        CHECK(omega_GeV > 0.);
        CHECK(omega_MeV == Approx(omega_GeV).epsilon(1e-12));
      }
    }
  }
}

TEST_CASE("SNR is computed for every detector the calculator owns", "[Detector]") {

  LOGGER(fatal);

  EffectivePotential::OneDimModel model;
  PhaseTracer::PhaseFinder pf(model);
  PhaseTracer::ActionCalculator ac(pf);
  PhaseTracer::TransitionFinder tf(pf);
  PhaseTracer::ThermoFinder tm(tf, ac);
  PhaseTracer::GravWaveCalculator gw(tm);

  const auto &dets = gw.get_detectors();
  const double year = 365.25 * 86400;

  /* A broad bump in ln f, peaked at 1 mHz but still sizeable in the PTA band */
  auto omega = [](double f) {
    const double x = std::log(f / 1e-3) / 8.;
    return 1e-10 * std::exp(-0.5 * x * x);
  };
  auto tabulate = [&omega](double f_lo, double f_hi, int n, std::vector<double> &fs, std::vector<double> &os) {
    fs.resize(n);
    os.resize(n);
    for (int ii = 0; ii < n; ii++) {
      fs[ii] = f_lo * std::pow(f_hi / f_lo, static_cast<double>(ii) / (n - 1));
      os[ii] = omega(fs[ii]);
    }
  };

  // SECTION("Detector list") {
  //   REQUIRE(dets.size() >= 2);
  //   CHECK(dets[0]->name() == "LISA");
  //   CHECK(dets[1]->name() == "Taiji");
  //   std::set<std::string> names;
  //   for (const auto &d : dets) {
  //     names.insert(d->name());
  //   }
  //   CHECK(names.size() == dets.size());
  //   for (const double f : {1e-4, 1e-3, 1e-2}) {
  //     CHECK(dets[0]->omega_noise(f) == PhaseTracer::noise_omega_LISA(f));
  //     CHECK(dets[1]->omega_noise(f) == PhaseTracer::noise_omega_Taiji(f));
  //   }
  //   for (const auto &d : dets) {
  //     CHECK_NOTHROW(PhaseTracer::check_detector(*d));
  //   }
  // }

  SECTION("PTA band") {
    const PhaseTracer::PTA pta(34, 18, 10, 1e-6, 10, "EPTA");
    CHECK(pta.f_min() == Approx(1. / (18 * year)));
    CHECK(pta.f_max() == Approx(1. / (10 * 86400.)));
    CHECK(pta.threshold() == 10.);
    /* Omega_n grows as f^5 */
    CHECK(pta.omega_noise(2e-8) / pta.omega_noise(1e-8) == Approx(32.));
  }

  SECTION("Noise curves match TransitionListener") {
    /* h^2 Omega_noise from TransitionListener's tab_data/noise-dataset.csv (1811.11175) */
    const std::vector<std::tuple<std::shared_ptr<const PhaseTracer::Detector>, double, double>> refs = {
      {std::make_shared<PhaseTracer::B_DECIGO>(), 0.027542287033382, 7.27843840314172e-10},
      {std::make_shared<PhaseTracer::B_DECIGO>(), 0.275422870333816, 9.92252705038004e-11},
      {std::make_shared<PhaseTracer::B_DECIGO>(), 10.0, 1.46597444431235e-06},
      {std::make_shared<PhaseTracer::DECIGO>(), 0.003019951720402, 7.40313602614322e-12},
      {std::make_shared<PhaseTracer::DECIGO>(), 0.131825673855641, 1.32654062057266e-13},
      {std::make_shared<PhaseTracer::DECIGO>(), 10.0, 4.72503340258623e-08},
      {std::make_shared<PhaseTracer::BBO>(), 0.003019951720402, 3.80891660270254e-12},
      {std::make_shared<PhaseTracer::BBO>(), 0.173780082874938, 2.14401960604109e-14},
      {std::make_shared<PhaseTracer::BBO>(), 10.0, 5.72927408786253e-08},
      {std::make_shared<PhaseTracer::PTA>(1000, 20, 14, 1e-7, 4, "SKA_20"), 1e-8, 1.7596306780268e-12},
      {std::make_shared<PhaseTracer::PTA>(1000, 20, 14, 1e-7, 4, "SKA_20"), 1e-7, 1.7596306780268e-07},
    };
    for (const auto &[d, f, expected] : refs) {
      INFO(d->name() << " at f = " << f);
      CHECK(d->omega_noise(f) == Approx(expected).epsilon(1e-4));
    }
  }

  SECTION("SNR integral matches a direct calculation") {
    std::vector<double> fs, os;
    tabulate(1e-10, 1e3, 4001, fs, os);
    const auto results = gw.get_SNR_tabulated(fs, os);
    REQUIRE(results.size() == dets.size());

    for (size_t jj = 0; jj < dets.size(); jj++) {
      const auto &d = *dets[jj];
      const auto &r = results[jj];
      /* Trapezoid in ln f on a much finer grid over the part of the band the grid covers. Steep
         integrands (a PTA's falls like f^-9) are dominated by the lowest in-band grid point */
      const double f_lo = *std::find_if(fs.begin(), fs.end(), [&d](double f) { return f >= d.f_min() * (1 - 1e-9); });
      const double f_hi = *std::find_if(fs.rbegin(), fs.rend(), [&d](double f) { return f <= d.f_max() * (1 + 1e-9); });
      const int n = 200001;
      const double du = std::log(f_hi / f_lo) / (n - 1);
      double integral = 0.;
      for (int ii = 0; ii < n; ii++) {
        const double f = f_lo * std::exp(ii * du);
        const double g = f * std::pow(omega(f) / d.omega_noise(f), 2);
        integral += (ii == 0 || ii == n - 1 ? 0.5 : 1.) * g * du;
      }
      CHECK(r.detector == d.name());
      CHECK(r.snr > 0.);
      CHECK(r.snr == Approx(std::sqrt(d.ndet() * integral * d.T_obs() * year)).epsilon(1e-4));
      CHECK(r.threshold == d.threshold());
      CHECK(r.T_obs == d.T_obs());
      CHECK(r.f_min == d.f_min());
      CHECK(r.f_max == d.f_max());
      CHECK(r.band_covered);
      CHECK(r.detected == (r.snr >= r.threshold));
      CHECK(r.noise.empty());
    }
  }

  SECTION("Band coverage") {
    std::vector<double> fs, os;
    tabulate(1e-4, 1e-2, 201, fs, os);
    for (const auto &r : gw.get_SNR_tabulated(fs, os)) {
      CHECK_FALSE(r.band_covered);
      if (r.f_max < fs.front() || r.f_min >= fs.back()) {
        CHECK(r.snr == 0.);  // fewer than two grid points lie in the band, e.g. a PTA or B_DECIGO
      } else {
        CHECK(r.snr > 0.);
      }
    }

    /* Entirely above every band: no in-band points */
    tabulate(1e3, 1e4, 11, fs, os);
    for (const auto &r : gw.get_SNR_tabulated(fs, os)) {
      CHECK(r.snr == 0.);
      CHECK_FALSE(r.detected);
    }
  }

  // SECTION("check_detector rejects unphysical detectors") {
  //   struct BadBand : PhaseTracer::LISA {
  //     double f_min() const override { return 1e-1; }
  //     double f_max() const override { return 1e-5; }
  //   };
  //   struct BadTime : PhaseTracer::LISA {
  //     double T_obs() const override { return 0.; }
  //   };
  //   CHECK_THROWS_AS(PhaseTracer::check_detector(BadBand()), std::invalid_argument);
  //   CHECK_THROWS_AS(PhaseTracer::check_detector(BadTime()), std::invalid_argument);
  //   /* A negative observation time gives a negative band */
  //   CHECK_THROWS_AS(PhaseTracer::check_detector(PhaseTracer::PTA(34, -18, 10, 1e-6, 10, "bad")), std::invalid_argument);
  // }

  SECTION("Spectra carry one result and noise curve per detector") {
    gw.set_vw(0.7);
    PhaseTracer::TransitionMilestone m{};
    m.status = PhaseTracer::MilestoneStatus::YES;
    m.temperature = 50.;
    m.reheating_temperature = 52.;
    m.alpha = 0.2;
    m.betaH_eff = 50.;
    m.cs_plus = 1. / std::sqrt(3.);
    m.vw = 0.7;
    m.g_eff = 106.75;
    m.h_eff = 106.75;

    const auto sp = gw.calc_spectrum(m);
    REQUIRE(sp.snr.size() == dets.size());
    for (size_t jj = 0; jj < dets.size(); jj++) {
      REQUIRE(sp.snr[jj].noise.size() == sp.frequency.size());
      CHECK(sp.snr[jj].noise[7] == dets[jj]->omega_noise(sp.frequency[7]));
    }
    CHECK(&sp.snr_for("Taiji") == &sp.snr[1]);
    CHECK_THROWS_AS(sp.snr_for("no such detector"), std::out_of_range);

    /* The spectrum's SNR uses the fit formulae over each detector's full band, so it matches get_SNR
       and is non-zero even in bands the frequency grid misses (e.g. the PTAs) */
    const auto s = gw.get_SNR(m);
    REQUIRE(s.size() == dets.size());
    for (size_t jj = 0; jj < dets.size(); jj++) {
      CHECK(sp.snr[jj].snr == Approx(s[jj].snr).epsilon(1e-12));
      CHECK(sp.snr[jj].snr > 0.);
      CHECK(sp.snr[jj].band_covered);
    }
    CHECK(s[0].snr == Approx(136.21548820870757).epsilon(1e-9));
    CHECK(s[1].snr == Approx(145.33402812584717).epsilon(1e-9));

    /* Summing two copies doubles Omega, and so doubles the SNR */
    const auto summed = gw.sum_spectrums({sp, sp});
    REQUIRE(summed.snr.size() == dets.size());
    for (size_t jj = 0; jj < dets.size(); jj++) {
      CHECK(summed.snr[jj].snr == Approx(2. * sp.snr[jj].snr).epsilon(1e-12));
    }
  }
}
