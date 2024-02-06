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

double Shooting::dV_from_absMin(double delta_phi){
  const double phi = phi_absMin + delta_phi;
  double dV_ = dV(phi);
  if (phi_eps_rel > 0){ // This will not be used if phi_eps_rel == 0
    const double phi_eps = phi_eps_rel * fabs(phi_absMin - phi_metaMin);
    const double dV_eps = d2V(phi) * delta_phi;
    const double blend_factor = exp(-pow((delta_phi/phi_eps),2));
    dV_ = dV_eps*blend_factor + dV_*(1-blend_factor);
  }
  return dV_;
}

void Shooting::findBarrierLocation(){
  const double phi_tol = fabs(phi_metaMin - phi_absMin) * 1e-12;
  const double V_phimeta = V(phi_metaMin);
  double phi1 = phi_metaMin;
  double phi2 = phi_absMin;
  double phi0 = 0.5 * (phi1 + phi2);
  while (std::abs(phi1 - phi2) > phi_tol) {
    const double V0 = V(phi0);
    if (V0 > V_phimeta) {
      phi1 = phi0;
    } else {
      phi2 = phi0;
    }
    phi0 = 0.5 * (phi1 + phi2);
  }
  phi_bar = phi0;
  LOG(debug) << "phi_bar = "<< phi_bar;
  return;
}

void Shooting::findRScale() {
    double phi_tol = fabs(phi_bar - phi_metaMin) * 1e-6;
    double x1 = std::min(phi_bar, phi_metaMin);
    double x2 = std::max(phi_bar, phi_metaMin);
    
    int bits = std::ceil(std::log2(1.0 / phi_tol));
    std::pair<double, double>  result = boost::math::tools::brent_find_minima(
        [this](double x) { return -V(x); }, x1, x2, bits);
    double phi_bar_top = result.first;
  
    if (phi_bar_top + phi_tol > x2 || phi_bar_top - phi_tol < x1) {
        throw std::runtime_error(
            "Minimization is placing the top of the potential barrier outside of the interval defined by phi_bar and phi_metaMin. Assume that the barrier does not exist.");
    }
    
    LOG(debug) << "phi_bar_top = "<< phi_bar_top;
  
    double Vtop = V(phi_bar_top) - V(phi_metaMin);
    double xtop = phi_bar_top - phi_metaMin;
    
    if (Vtop <= 0) {
        throw std::runtime_error("Barrier height is not positive, does not exist.");
    }
    
    rscale = std::abs(xtop) / std::sqrt(std::abs(6 * Vtop));
    LOG(debug) << "rscale = "<< rscale;
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
    // TODO: ignore the warnings if there is
    phi = (std::tgamma(nu+1)*pow(0.5*beta_r,-nu) * boost::math::cyl_bessel_i(nu, beta_r)-1) * dV_/d2V_;
    dphi = -nu*(pow(0.5*beta_r,-nu) / r) * boost::math::cyl_bessel_i(nu, beta_r);
    dphi += pow(0.5*beta_r,-nu) * 0.5*beta
            * (boost::math::cyl_bessel_i(nu-1, beta_r)+boost::math::cyl_bessel_i(nu+1, beta_r));
    dphi *= std::tgamma(nu+1) * dV_/d2V_;
    phi += phi0;
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
  const double d2V_ = d2V(phi0);
  
  *r0 = rmin;
  
  const auto deltaPhiDiff = [this, phi0, dV_, d2V_, phi_r0, dphi_r0, delta_phi_cutoff](double rtry) {
    this->exactSolution(rtry, phi0, dV_, d2V_, phi_r0, dphi_r0);
    return fabs(*phi_r0 - phi_absMin) - fabs(delta_phi_cutoff);
  };
  if ( deltaPhiDiff(rmin) > 0) { // phi_rmin - phi_absMin > delta_phi_cutoff
    LOG(debug) << "Initial point for phi(r=0) = " << phi0 << " : ";
    LOG(debug) << "  r_0 =  r_min = " << *r0 << ", phi = " << *phi_r0 << ", dphi = " << *dphi_r0;
    return;
  }
  if ( (*dphi_r0)*delta_phi0 < 0 ) { // wrong direction
    LOG(debug) << "Initial point for phi(r=0) = " << phi0 << " : ";
    LOG(debug) << "  r_0 = " << *r0 << ", phi = " << *phi_r0 << ", dphi = " << *dphi_r0;
    LOG(debug) << "  Wrong direction!";
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
  LOG(debug) << "Initial point for phi(r=0) = " << phi0 << " : ";
  LOG(debug) << "  r_0 = " << *r0 << ", phi = " << *phi_r0 << ", dphi = " << *dphi_r0;
  return;
}

int Shooting::integrateProfile(double r0, std::vector<double> y0,
          double* rf, std::vector<double>* yf,
          double dr0, std::vector<double> epsabs, std::vector<double> epsfrac,
          double drmin, double rmax){
  
  int convergence_type;
  int ysign = y0[0] > phi_metaMin ? 1 : -1;

  using state_type = std::vector<double>;
  using error_stepper_type =
     boost::numeric::odeint::runge_kutta_dopri5<state_type>;
  using controlled_stepper_type =
     boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;
  state_type y1 = y0;
  double r1 = r0;
  double dr = dr0;
  state_type dydr1(2), dydr0(2);
  double eps_abs = epsabs[1];
  double eps_rel = epsfrac[1];
  controlled_stepper_type stepper
     = make_controlled(eps_rel, eps_abs,
                       error_stepper_type());
  while (true) {
    y0 = y1;
    r0 = r1;
    while ( stepper.try_step(
          [this](const state_type& y, state_type& dydr, double r) {equationOfMotion(y, dydr, r);},
                             y1, r1, dr) ) {};
    
//      LOG(debug) << "r = " << r1 << ", phi = " << y1[0] << ", dphi/dr = " << y1[1];
    if (r1>r0+rmax){
      throw std::runtime_error("r > rmax");
    } else if (dr < drmin){
      throw std::runtime_error("dr < drmin");
    } else if (fabs(y1[0]-phi_metaMin)< 3.*epsabs[0] and fabs(y1[1])< 3.*epsabs[1] ){
      convergence_type = 0; // converged
      *rf = r1;
      *yf = y1;
      break;
    } else if ( (y1[0]-phi_metaMin)*ysign < 0 or y1[1]*ysign > 0 ){
      dydr0 = dY(y0, r0);
      dydr1 = dY(y1, r1);
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
  }
  // Check convergence for a second time
  if (fabs((*yf)[0]-phi_metaMin)< 3*epsabs[0] and fabs((*yf)[1]) < 3*epsabs[1] ) {
    convergence_type = 0;
  }
  
  LOG(debug) << "Integrated to ";
  LOG(debug) << "  r = " << *rf << ", phi = " << (*yf)[0] << ", dphi/dr = " << (*yf)[1];
  LOG(debug) << "Convergence_type = " << ( (convergence_type ==0) ? "converged" : ( (convergence_type == 1) ? "overshoot":"undershoot") )  << std::endl;
  
  return convergence_type;
}

Profile1D Shooting::integrateAndSaveProfile(Eigen::VectorXd R, std::vector<double> y0,
                      double dr0, std::vector<double> epsabs, std::vector<double> epsfrac, double drmin){
  Profile1D profile;
  int N = R.size();
  double r0 = R[0];
  Eigen::VectorXd Phi(N);
  Eigen::VectorXd dPhi(N);
  Phi(0) = y0[0];
  dPhi(0) = y0[1];

  std::vector<double> y1 = y0;
  double r1 = r0;
  double dr_ = dr0;
  double Rerr = NAN;
  
  using state_type = std::vector<double>;
  using error_stepper_type =
     boost::numeric::odeint::runge_kutta_dopri5<state_type>;
  using controlled_stepper_type =
     boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;
  double eps_abs = epsabs[1];
  double eps_rel = epsfrac[1];
  controlled_stepper_type stepper
     = make_controlled(eps_rel, eps_abs,
                       error_stepper_type());
  int ii = 1;
  while (ii<N){
    y0 = y1;
    r0 = r1;
    while ( stepper.try_step(
          [this](const state_type& y, state_type& dydr, double r) {equationOfMotion(y, dydr, r);},
                             y1, r1, dr_) ) {};
    double dr = r1 - r0; // This may be different to dr_
    if (dr < drmin) {
      r1 = r0 + drmin;
      y1[0] = y0[0] + (y1[0]-y0[0])*drmin/dr;
      y1[1] = y0[1] + (y1[1]-y0[1])*drmin/dr;
      dr = drmin;
      if ( not std::isnan(Rerr) ) Rerr = r1; // TODO
    }
    if (r0 < R[ii] and R[ii] <= r1) {
      std::vector<double> dydr0 = dY(y0, r0);
      std::vector<double> dydr1 = dY(y1, r1);
      CubicInterpFunction CF0(y0[0], dr*dydr0[0], y1[0], dr*dydr1[0]);
      CubicInterpFunction CF1(y0[1], dr*dydr0[1], y1[1], dr*dydr1[1]);
      while ( ii < N and r0 < R[ii] and R[ii] <= r1 ) {
        double x = (R[ii]-r0)/dr;
        Phi[ii] = CF0(x);
        dPhi[ii] = CF1(x);
        ii ++;
      }
    }
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
      LOG(debug) << "Reached xtol" << std::endl;
      break;
    }
  
  }
  
  // Getting the points along the path
  Eigen::VectorXd R = Eigen::VectorXd::LinSpaced(npoints, r0, rf);
  Profile1D profile_final = integrateAndSaveProfile(R, y0, dr0, epsfrac, epsabs, drmin);
  
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
    double d2V_ = d2V(profile_int.Phi[0]);
    
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
  boost::multi_array<double, 1> integrand(boost::extents[n]);
  for (size_t i = 0; i < n; ++i) {
      // Find the area of an n-sphere (alpha=n):
      double area = std::pow(r[i], alpha) * 2 * pow(M_PI, d * 0.5) / tgamma(d * 0.5);
      integrand[i] = 0.5 * std::pow(dphi[i], 2) + V(phi[i]) - V(phi_metaMin);
      integrand[i] *= area;
//        std::cout << "integrand = " << integrand[i] << std::endl;
  }
  double lower = r[0], upper = r[n - 1];
  double S = boost::math::quadrature::gauss_kronrod<double, 15>::integrate([&](double x) {
          size_t index = static_cast<size_t>((x - lower) / (upper - lower) * (integrand.size() - 1));
          return integrand[index];
      }, lower, upper);
  // Find the bulk term in the bubble interior
  double volume = std::pow(r[0], d) * pow(M_PI, d * 0.5) / tgamma(d * 0.5 + 1);
  
  S += volume * (V(phi[0]) - V(phi_metaMin));
  return S;
}


}  // namespace PhaseTracer
