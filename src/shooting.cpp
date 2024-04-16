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

#include "shooting.hpp"

namespace PhaseTracer {


/* Take one 5th-order Cash-Karp Runge-Kutta step */
std::pair<Eigen::Vector2d, Eigen::Vector2d> _rkck(Eigen::Vector2d y, Eigen::Vector2d dydt, double t, std::function<Eigen::Vector2d(Eigen::Vector2d, double)> f, double dt) {
  double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875, b21=0.2;
  double b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42=-0.9, b43=1.2;
  double b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0;
  double b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0, c1=37.0/378.0;
  double c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0, dc5=-277.00/14336.0;
  double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, dc6=c6-0.25;
    
  Eigen::Vector2d ytemp = y + b21 * dt * dydt;
  Eigen::Vector2d ak2 = f(ytemp, t + a2 * dt);
    
  ytemp = y + dt * (b31 * dydt + b32 * ak2);
  Eigen::Vector2d ak3 = f(ytemp, t + a3 * dt);
    
  ytemp = y + dt * (b41 * dydt + b42 * ak2 + b43 * ak3);
  Eigen::Vector2d ak4 = f(ytemp, t + a4 * dt);
    
  ytemp = y + dt * (b51 * dydt + b52 * ak2 + b53 * ak3 + b54 * ak4);
  Eigen::Vector2d ak5 = f(ytemp, t + a5 * dt);
    
  ytemp = y + dt * (b61 * dydt + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);
  Eigen::Vector2d ak6 = f(ytemp, t + a6 * dt);
    
  Eigen::Vector2d dyout = dt * (c1 * dydt + c3 * ak3 + c4 * ak4 + c6 * ak6);
  Eigen::Vector2d yerr = dt * (dc1 * dydt + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
    
  return {dyout, yerr};
}

/* Take a single 5th order Runge-Kutta step with error monitoring */
struct rkqs_rval {
  Eigen::Vector2d Delta_y;
  double Delta_t;
  double dtnxt;
  rkqs_rval(Eigen::Vector2d dy, double dt, double dtn) : Delta_y(dy), Delta_t(dt), dtnxt(dtn) {}
};
rkqs_rval rkqs(Eigen::Vector2d y, Eigen::Vector2d dydt, double t, std::function<Eigen::Vector2d(Eigen::Vector2d, double)> f, double dt_try, Eigen::Vector2d epsfrac, Eigen::Vector2d epsabs) {
  double dt = dt_try;
  double errmax;
  Eigen::Vector2d dy, yerr;
  while (true) {
    std::tie(dy, yerr) = _rkck(y, dydt, t, f, dt);
    errmax = (yerr.array() / epsabs.array()).abs().cwiseMin((yerr.array().abs() / (y.array().abs() + 1E-50) / epsfrac.array())).maxCoeff();
    if (errmax < 1.0) break;
    double dttemp = 0.9 * dt * pow(errmax, -0.25);
    dt = (dt > 0) ? std::max(dttemp, dt * 0.1) : std::min(dttemp, dt * 0.1);
    if (t + dt == t) {
        throw std::runtime_error("Stepsize rounds down to zero.");
    }
  }
  double dtnext = (errmax > 1.89e-4) ? 0.9 * dt * std::pow(errmax, -0.2) : 5 * dt;
  
  return rkqs_rval(dy, dt, dtnext);
}



double Shooting::dV_from_absMin(double delta_phi){
  const double phi = phi_absMin + delta_phi;
  double dV_ = ps.dV(phi);
  if (phi_eps_rel > 0){ // This will not be used if phi_eps_rel == 0
    const double phi_eps = phi_eps_rel * fabs(phi_absMin - phi_metaMin);
    const double dV_eps = ps.d2V(phi) * delta_phi;
    const double blend_factor = exp(-pow((delta_phi/phi_eps),2));
    dV_ = dV_eps*blend_factor + dV_*(1-blend_factor);
  }
  return dV_;
}

void Shooting::findBarrierLocation(){
  const double phi_tol = fabs(phi_metaMin - phi_absMin) * 1e-12;
  const double V_phimeta = ps.V(phi_metaMin);
  double phi1 = phi_metaMin;
  double phi2 = phi_absMin;
  double phi0 = 0.5 * (phi1 + phi2);
  while (std::abs(phi1 - phi2) > phi_tol) {
    const double V0 = ps.V(phi0);
    if (V0 > V_phimeta) {
      phi1 = phi0;
    } else {
      phi2 = phi0;
    }
    phi0 = 0.5 * (phi1 + phi2);
  }
  phi_bar = phi0;
  LOG(trace) << "phi_bar = "<< phi_bar;
  return;
}

void Shooting::findRScale() {
    double phi_tol = fabs(phi_bar - phi_metaMin) * 1e-6;
    double x1 = std::min(phi_bar, phi_metaMin);
    double x2 = std::max(phi_bar, phi_metaMin);
    
    int bits = std::ceil(std::log2(1.0 / phi_tol));
    std::pair<double, double>  result = boost::math::tools::brent_find_minima(
        [this](double x) { return -ps.V(x); }, x1, x2, bits);
    double phi_bar_top = result.first;
  
    if (phi_bar_top + phi_tol > x2 || phi_bar_top - phi_tol < x1) {
        throw std::runtime_error(
            "Minimization is placing the top of the potential barrier outside of the interval defined by phi_bar and phi_metaMin. Assume that the barrier does not exist.");
    }
    
    LOG(trace) << "phi_bar_top = "<< phi_bar_top;
  
    double Vtop = ps.V(phi_bar_top) - ps.V(phi_metaMin);
    double xtop = phi_bar_top - phi_metaMin;
    
    if (Vtop <= 0) {
        throw std::runtime_error("Barrier height is not positive, does not exist.");
    }
    
    rscale = std::abs(xtop) / std::sqrt(std::abs(6 * Vtop));
    LOG(trace) << "rscale = "<< rscale;
    return;
}

void Shooting::exactSolution(double r, double phi0, double dV_, double d2V_,
                   double* phi_r, double* dphi_r){
  const double beta = sqrt(fabs(d2V_));
  const double beta_r = beta*r;
  const double nu = 0.5*(alpha-1.);
  double phi = 0., dphi = 0.;
  if (beta_r < 1e-2){
    const int s = d2V_>0 ? +1 : -1;
    for (int k=1; k<4; k++){
      double temp = pow(0.5*beta_r, 2.*k-2.) * pow(s,k) / (std::tgamma(k+1.)*std::tgamma(k+1.+nu));
      phi += temp;
      dphi += temp * (2.*k);
    }
    phi *= 0.25 * std::tgamma(nu+1) * pow(r,2) * dV_ * s;
    dphi *= 0.25 * std::tgamma(nu+1) * r * dV_ * s;
    phi += phi0;
  }else if (d2V_>0){
    try {
      phi = (std::tgamma(nu+1)*pow(0.5*beta_r,-nu) * boost::math::cyl_bessel_i(nu, beta_r)-1) * dV_/d2V_;
      dphi = -nu*(pow(0.5*beta_r,-nu) / r) * boost::math::cyl_bessel_i(nu, beta_r);
      dphi += pow(0.5*beta_r,-nu) * 0.5*beta
              * (boost::math::cyl_bessel_i(nu-1, beta_r)+boost::math::cyl_bessel_i(nu+1, beta_r));
      dphi *= std::tgamma(nu+1) * dV_/d2V_;
      phi += phi0;
    } catch (const boost::wrapexcept<std::overflow_error>& e){
      phi = std::numeric_limits<double>::infinity();
      dphi = NAN;
    }
  }else{
    phi = (std::tgamma(nu+1)*pow(0.5*beta_r,-nu) * boost::math::cyl_bessel_j(nu, beta_r)-1) * dV_/d2V_;
    dphi = -nu*(pow(0.5*beta_r,-nu) / r) * boost::math::cyl_bessel_j(nu, beta_r);
    dphi += pow(0.5*beta_r,-nu) * 0.5*beta
            * (boost::math::cyl_bessel_j(nu-1, beta_r)-boost::math::cyl_bessel_j(nu+1, beta_r));
    dphi *= std::tgamma(nu+1) * dV_/d2V_;
    phi += phi0;
  }
//  LOG(debug) << "For phi(r=0) = " << phi0 << " : ";
//  LOG(debug) << "  beta_r = "<< beta_r << ", d2V[phi(r=0)] = " << d2V_;
//  LOG(debug) << "  phi(r=" << r <<") = " << phi << ", dphi(r=" << r <<") = " << dphi;
  
  *phi_r = phi;
  *dphi_r = dphi;
}

void Shooting::initialConditions(double delta_phi0, double rmin, double delta_phi_cutoff,
                       double* r0, double* phi_r0, double* dphi_r0){

  const double phi0 = phi_absMin + delta_phi0;
  const double dV_ = dV_from_absMin(delta_phi0);
  const double d2V_ = ps.d2V(phi0);
  
  *r0 = rmin;
  
  const auto deltaPhiDiff = [this, phi0, dV_, d2V_, phi_r0, dphi_r0, delta_phi_cutoff](double rtry) {
    this->exactSolution(rtry, phi0, dV_, d2V_, phi_r0, dphi_r0);
    return fabs(*phi_r0 - phi_absMin) - fabs(delta_phi_cutoff);
  };
  if ( deltaPhiDiff(rmin) > 0) { // phi_rmin - phi_absMin > delta_phi_cutoff
    LOG(trace) << "Initial point for phi(r=0) = " << phi0 << " : ";
    LOG(trace) << "  r_0 =  r_min = " << *r0 << ", phi = " << *phi_r0 << ", dphi = " << *dphi_r0;
    return;
  }
  if ( (*dphi_r0)*delta_phi0 < 0 ) { // wrong direction
    LOG(trace) << "Initial point for phi(r=0) = " << phi0 << " : ";
    LOG(trace) << "  r_0 = " << *r0 << ", phi = " << *phi_r0 << ", dphi = " << *dphi_r0;
    LOG(trace) << "  Wrong direction!";
    return;
  }

  // Find the smallest r0 such that delta_phi_r0 > delta_phi_cutoff
  double r = rmin;
  double rlast = r;
  while(std::isfinite(r)){
    rlast = r;
    r *= 10.;
    if ( deltaPhiDiff(r) > 0) break;
  }
  
  // Find phi - phi_absMin == delta_phi_cutoff exactly
  const auto result = boost::math::tools::bisect(deltaPhiDiff, rlast, r, boost::math::tools::eps_tolerance<double>(), max_iter);
  *r0 = (result.first + result.second) * 0.5;
  exactSolution(*r0, phi0, dV_, d2V_, phi_r0, dphi_r0);
  LOG(trace) << "Initial point for phi(r=0) = " << phi0 << " : ";
  LOG(trace) << "  r_0 = " << *r0 << ", phi = " << *phi_r0 << ", dphi = " << *dphi_r0;
  return;
}

int Shooting::integrateProfile(double r0, std::vector<double> y0_,
          double* rf, std::vector<double>* yf,
          double dr0, std::vector<double> epsabs_, std::vector<double> epsfrac_,
          double drmin, double rmax){
  
  Eigen::Map<Eigen::Vector2d> y0(y0_.data());
  Eigen::Map<Eigen::Vector2d> epsabs(epsabs_.data());
  Eigen::Map<Eigen::Vector2d> epsfrac(epsfrac_.data());
  
  double dr = dr0;
  std::function<Eigen::Vector2d(Eigen::Vector2d, double)> dY = [this](Eigen::Vector2d y, double t){
    return equationOfMotion(y, t);
  };
  Eigen::Vector2d dydr0 = dY(y0, r0);
  int ysign = y0[0] > phi_metaMin ? 1 : -1;
  rmax += r0;
  
  int convergence_type;
  while (true) {
    auto rval = rkqs(y0, dydr0, r0, dY, dr, epsfrac, epsabs);
    
    Eigen::Vector2d dy = rval.Delta_y;
    dr = rval.Delta_t;
    double drnext = rval.dtnxt;
    
    double r1 = r0 + dr;
    Eigen::Vector2d y1 = y0 + dy;
    Eigen::Vector2d dydr1 = dY(y1,r1);
    
//    LOG(trace) << "r = " << r1 << ", phi = " << y1[0] << ", dphi/dr = " << y1[1];
    if (r1 > rmax){
      throw std::runtime_error("r > rmax");
    } else if (dr < drmin){
      throw std::runtime_error("dr < drmin");
    } else if (fabs(y1[0]-phi_metaMin)< 3.*epsabs[0] and fabs(y1[1])< 3.*epsabs[1] ){
      convergence_type = 0; // converged
      *rf = r1;
      std::vector<double> y1_(y1.data(), y1.data() + y1.size()); // TODO
      *yf = y1_;
      break;
    } else if ( (y1[0]-phi_metaMin)*ysign < 0 or y1[1]*ysign > 0 ){
      CubicInterpFunction CF0(y0[0], dr*dydr0[0], y1[0], dr*dydr1[0], phi_metaMin);
      CubicInterpFunction CF1(y0[1], dr*dydr0[1], y1[1], dr*dydr1[1]);
      double x;
      if (y1[1]*ysign > 0){
        // Extrapolate to where dphi(r) = 0
        const auto result = boost::math::tools::bisect(CF1, 0., 1., boost::math::tools::eps_tolerance<double>(), max_iter);
        x = (result.first + result.second) * 0.5;
        convergence_type = -1; // "undershoot"
      } else {
        // Extrapolate to where phi(r) = phi_metaMin
        const auto result = boost::math::tools::bisect(CF0, 0., 1., boost::math::tools::eps_tolerance<double>(), max_iter);
        x = (result.first + result.second) * 0.5;
        convergence_type = 1; // "overshoot"
      }
      *rf = r0+dr*x;
      *yf = {phi_metaMin + CF0(x), CF1(x)} ;
      break;
    }
    r0 = r1;
    y0 = y1;
    dydr0 = dydr1;
    dr = drnext;
  }

  // Check convergence for a second time
  if (fabs((*yf)[0]-phi_metaMin)< 3*epsabs[0] and fabs((*yf)[1]) < 3*epsabs[1] ) {
    convergence_type = 0;
  }

  LOG(trace) << "Integrated to ";
  LOG(trace) << "  r = " << *rf << ", phi = " << (*yf)[0] << ", dphi/dr = " << (*yf)[1];
  LOG(trace) << "Convergence_type = " << ( (convergence_type ==0) ? "converged" : ( (convergence_type == 1) ? "overshoot":"undershoot") )  << std::endl;
//  std::exit(0);
  return convergence_type;
}

Profile1D Shooting::integrateAndSaveProfile(Eigen::VectorXd R, std::vector<double> y0_,
                      double dr, std::vector<double> epsabs_, std::vector<double> epsfrac_, double drmin){
  
  Eigen::Map<Eigen::Vector2d> y0(y0_.data());
  Eigen::Map<Eigen::Vector2d> epsabs(epsabs_.data());
  Eigen::Map<Eigen::Vector2d> epsfrac(epsfrac_.data());
  
  Profile1D profile;
  int N = R.size();
  Eigen::VectorXd Phi(N);
  Eigen::VectorXd dPhi(N);
  Phi(0) = y0[0];
  dPhi(0) = y0[1];
  
  double r0 = R[0];
  std::function<Eigen::Vector2d(Eigen::Vector2d, double)> dY = [this](Eigen::Vector2d y, double t){
    return equationOfMotion(y, t);
  };
  Eigen::Vector2d dydr0 = dY(y0, r0);
  double Rerr = NAN;
  
  int ii = 1;
  while (ii<N){
    auto rval = rkqs(y0, dydr0, r0, dY, dr, epsfrac, epsabs);
    
    Eigen::Vector2d dy = rval.Delta_y;
    dr = rval.Delta_t;
    double drnext = rval.dtnxt;

    double r1;
    Eigen::Vector2d y1;
    if (dr >= drmin) {
      r1 = r0 + dr;
      y1 = y0 + dy;
    } else {
      // TODO not checked
      y1 = y0 + dy*drmin/dr;
      drnext = drmin;
      dr = drnext;
      r1 = r0 + drmin;
      if ( not std::isnan(Rerr) ) Rerr = r1; // TODO
    }
    Eigen::Vector2d dydr1 = dY(y1,r1);
//    LOG(debug) << "r = " << r1 << ", phi = " << y1[0] << ", dphi/dr = " << y1[1];
    if (r0 < R[ii] and R[ii] <= r1) {
      CubicInterpFunction CF0(y0[0], dr*dydr0[0], y1[0], dr*dydr1[0]);
      CubicInterpFunction CF1(y0[1], dr*dydr0[1], y1[1], dr*dydr1[1]);
      while ( ii < N and r0 < R[ii] and R[ii] <= r1 ) {
        double x = (R[ii]-r0)/dr;
        Phi[ii] = CF0(x);
        dPhi[ii] = CF1(x);
        ii ++;
      }
    }

    r0 = r1;
    y0 = y1;
    dydr0 = dydr1;
    dr = drnext;
  }

  profile.R = R;
  profile.Phi = Phi;
  profile.dPhi = dPhi;
  profile.Rerr = Rerr;
  
  return profile;
}

Profile1D Shooting::findProfile(double metaMin, double absMin, double xguess, int max_interior_pts){
  phi_absMin = absMin;
  phi_metaMin = metaMin;
  LOG(trace) << "phi_absMin = " << phi_absMin;
  LOG(trace) << "phi_metaMin = " << phi_metaMin;
  findBarrierLocation();
  findRScale();
  
  // Set r parameters
  rmin *= rscale;
  double dr0 = rmin;
  double drmin = 0.01*rmin;
  rmax *= rscale;
  // Set the phi parameters
  const double delta_phi = phi_metaMin - phi_absMin;
  const std::vector<double> epsabs = {fabs(delta_phi*phitol), fabs(delta_phi*phitol/rscale)};
  const std::vector<double> epsfrac = {phitol, phitol};
  double delta_phi_cutoff = thinCutoff * delta_phi;
  // Set x parameters
  double xmin = xtol*10;
  double xmax = std::numeric_limits<double>::infinity();
  double x;

  if (std::isnan(xguess)){ // TODO check wheather this works
    x =  -log(fabs((phi_bar-phi_absMin)/delta_phi));
  } else {
    x = xguess;
  }
  
  double r0, rf=NAN;
  std::vector<double> y0(2);
  double delta_phi0;
  while (true){
    delta_phi0 = exp(-x)*delta_phi;
    double r0_, phi0, dphi0;
    initialConditions(delta_phi0, rmin, delta_phi_cutoff, &r0_, &phi0, &dphi0);
    
    if (not std::isfinite(r0_) or not std::isfinite(x)){
      if (std::isnan(rf)) {
        LOG(fatal) << "Failed to retrieve initial conditions on the first try.";
        std::terminate();
      }
      break;
    }
    
    r0 = r0_;
    y0[0] = phi0;
    y0[1] = dphi0;
    std::vector<double> yf(2);
    int shooting_type = integrateProfile(r0, y0, &rf, &yf, dr0, epsabs, epsfrac, drmin, rmax);

    if (shooting_type == 0){
      break;
    } else if (shooting_type == -1){ // undershoot
      xmin = x;
      x = std::isinf(xmax) ? x*5. : 0.5*(xmin+xmax);
    } else if (shooting_type == 1){ // overshoot
      xmax = x;
      x = .5*(xmin+xmax);
    }
    
    if ((xmax-xmin) < xtol){
      LOG(debug) << "Reached xtol";
      break;
    }
  
  }
  
  // Getting the points along the path
  Eigen::VectorXd R = Eigen::VectorXd::LinSpaced(npoints, r0, rf);
  Profile1D profile_final = integrateAndSaveProfile(R, y0, dr0, epsabs, epsfrac, drmin);
  
  // Make points interior to the bubble
  Profile1D profile_int;
  
  if (max_interior_pts <= 0) max_interior_pts = int(R.size()/2);
  double dx0 = R[1]-R[0];
  int n = std::ceil(R[0] / dx0);
  if (max_interior_pts>0 and n > 1){
    if (R[0] / dx0 <= max_interior_pts){
      profile_int.R = Eigen::VectorXd::LinSpaced(n, 0.0, R[0]).head(n-1);
    } else {
      n = max_interior_pts;
      double a = (R[0]/dx0 - n) * 2./(n*(n+1.));
      profile_int.R.resize(n);
      for (int i = 0; i < n; i++) {
          int N = n - i;
        profile_int.R[i] = R[0] - dx0 * (N + 0.5 * a * N * (N + 1.));
      }
      profile_int.R[0] = 0.0;
    }
    int int_size = profile_int.R.size();
    profile_int.Phi.resize(int_size);
    profile_int.dPhi.resize(int_size);
    profile_int.Phi[0] = phi_absMin + delta_phi0;
    profile_int.dPhi[0] = 0.0;
    
    double dV_ = dV_from_absMin(delta_phi0);
    double d2V_ = ps.d2V(profile_int.Phi[0]);
    
    for (int i=1; i < int_size; i++ ){
      exactSolution(profile_int.R[i], profile_int.Phi[0], dV_, d2V_, profile_int.Phi.data() + i, profile_int.dPhi.data() + i);
    }
  }
  Profile1D profile;
  int all_size = profile_int.R.size() + profile_final.R.size();
  profile.R.resize(all_size);
  profile.Phi.resize(all_size);
  profile.dPhi.resize(all_size);
  profile.R << profile_int.R, profile_final.R;
  profile.Phi << profile_int.Phi, profile_final.Phi;
  profile.dPhi << profile_int.dPhi, profile_final.dPhi;
  profile.Rerr = profile_final.Rerr;
  
  return profile;

}

double Shooting::calAction(Profile1D profile){
  
  const auto r = profile.R;
  const auto phi = profile.Phi;
  const auto dphi = profile.dPhi;
  size_t n = r.size();
  double d = alpha + 1;  // Number of dimensions in the integration
  // And integrate the profile
  alglib::real_1d_array x;
  alglib::real_1d_array integrand;
  x.setlength(n);
  integrand.setlength(n);
//  std::ofstream file("integrand.txt");
  for (size_t i = 0; i < n; ++i) {
    x[i] = r[i];
    // Find the area of an n-sphere (alpha=n):
    double area = std::pow(r[i], alpha) * 2 * pow(M_PI, d * 0.5) / tgamma(d * 0.5);
    integrand[i] = 0.5 * std::pow(dphi[i], 2) + ps.V(phi[i]) - ps.V(phi_metaMin);
    integrand[i] *= area;
//      file << r[i] << ", " << integrand[i] << std::endl;
//        std::cout << "integrand = " << integrand[i] << std::endl;
  }
//  file.close();

  alglib::spline1dinterpolant spline;
  alglib::spline1dbuildcubic(x, integrand, spline);
  
  double S = alglib::spline1dintegrate(spline, x[n-1]);
  
  // Find the bulk term in the bubble interior
  double volume = std::pow(r[0], d) * pow(M_PI, d * 0.5) / tgamma(d * 0.5 + 1);
  
  S += volume * (ps.V(phi[0]) - ps.V(phi_metaMin));
  return S;
}


void Shooting::evenlySpacedPhi(Profile1D pf, std::vector<double>* p,std::vector<double>* dp,
                     size_t npoints, bool fixAbs){
  
  Eigen::VectorXd phi = pf.Phi;
  Eigen::VectorXd dphi = pf.dPhi;
  
  // Excessively small intervals may cause errors in alglib::spline1dbuildcubic
  Eigen::VectorXd filtered_phi;
  Eigen::VectorXd filtered_dphi;
  filtered_phi.resize(phi.size());
  filtered_dphi.resize(dphi.size());
  int filtered_index = 0;
  for (int i = 1; i < phi.size(); i++){
    if (std::abs(phi(i) - phi(i - 1)) > 1E-12){
      filtered_phi(filtered_index) = phi(i);
      filtered_dphi(filtered_index) = dphi(i);
      filtered_index++;
    }
  }
  filtered_phi.conservativeResize(filtered_index);
  filtered_dphi.conservativeResize(filtered_index);

  alglib::real_1d_array phi_arr;
  alglib::real_1d_array dphi_arr;
  size_t arr_npoints = filtered_phi.size();
  if (fixAbs and phi_absMin<filtered_phi[0] ){
    arr_npoints++;
  }
  if (phi_metaMin>filtered_phi[filtered_phi.size()-1]){
    arr_npoints++;
  }
  phi_arr.setlength(arr_npoints);
  dphi_arr.setlength(arr_npoints);
  size_t ii=0;
  if (fixAbs and phi_absMin<filtered_phi[0] ){
    phi_arr[ii] = phi_absMin;
    dphi_arr[ii] = 0.;
    ii++;
  }
  for (size_t jj=0; jj<filtered_phi.size(); jj++){
    phi_arr[ii] = filtered_phi[jj];
    dphi_arr[ii] = filtered_dphi[jj];
    ii++;
  }
  if (phi_metaMin>filtered_phi[filtered_phi.size()-1]){
    phi_arr[ii] = phi_metaMin;
    dphi_arr[ii] = 0.;
  }

  alglib::spline1dinterpolant spl;
  alglib::spline1dbuildcubic(phi_arr, dphi_arr, spl);
  std::vector<double> evenly_phi(npoints);
  std::vector<double> evenly_dphi(npoints);
  double max = phi_metaMin;
  double min = fixAbs ? phi_absMin : phi_arr[0];
  double step = (max-min)/(npoints-1);
  for (size_t jj=0; jj<npoints; jj++){
    evenly_phi[jj] = min + jj * step;
    evenly_dphi[jj] = alglib::spline1dcalc(spl, evenly_phi[jj]);
  }
  
  *p = evenly_phi;
  *dp = evenly_dphi;
}

}  // namespace PhaseTracer
