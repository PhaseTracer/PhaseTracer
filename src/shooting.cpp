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
  LOG(debug) << "phi_bar = "<< phi_bar;
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
    
    LOG(debug) << "phi_bar_top = "<< phi_bar_top;
  
    double Vtop = ps.V(phi_bar_top) - ps.V(phi_metaMin);
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
  LOG(debug) << "phi_absMin = " << phi_absMin;
  LOG(debug) << "phi_metaMin = " << phi_metaMin;
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
  boost::multi_array<double, 1> integrand(boost::extents[n]);
  for (size_t i = 0; i < n; ++i) {
      // Find the area of an n-sphere (alpha=n):
      double area = std::pow(r[i], alpha) * 2 * pow(M_PI, d * 0.5) / tgamma(d * 0.5);
      integrand[i] = 0.5 * std::pow(dphi[i], 2) + ps.V(phi[i]) - ps.V(phi_metaMin);
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
    dphi_arr[ii] = filtered_phi[jj];
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
  
  std::vector<double> _dphi = {0.        ,0.01140122,0.02280328,0.03420565,0.04560822,0.05701094,
    0.06841378,0.07981672,0.09104573,0.10210272,0.11314672,0.1241725 ,
    0.13518779,0.14632975,0.15759505,0.16892094,0.18026054,0.19158224,
    0.2028653 ,0.2140954 ,0.22526342,0.23636449,0.24739712,0.25836242,
    0.26925843,0.28008571,0.2908499 ,0.30154835,0.31218535,0.32276683,
    0.33328654,0.34375171,0.35416582,0.36453023,0.37484636,0.38511566,
    0.39533758,0.40551645,0.41565407,0.42575196,0.43581162,0.44583284,
    0.45581683,0.4657661 ,0.47568213,0.48556636,0.49542019,0.50524502,
    0.51504076,0.52480229,0.53453804,0.5442493 ,0.55393733,0.56359379,
    0.57322634,0.58283859,0.59242869,0.60198942,0.61153263,0.62105851,
    0.63055394,0.64003448,0.64949824,0.65893457,0.66835847,0.67776214,
    0.68714555,0.69651936,0.7058737 ,0.71521991,0.72455155,0.73386506,
    0.74316769,0.7524407 ,0.76170715,0.77093463,0.78015423,0.78933636,
    0.79850487,0.80763696,0.8167512 ,0.82582887,0.83488664,0.84390583,
    0.85290578,0.86186271,0.87080425,0.8796954 ,0.88857864,0.89740074,
    0.90622284,0.91497642,0.92372643,0.93242086,0.94109596,0.94973323,
    0.95833079,0.96691328,0.97543085,0.98394842,0.99239666,1.00083198,
    1.00922924,1.0175802 ,1.02593005,1.0341947 ,1.04245935,1.05067749,
    1.05885403,1.06703056,1.07511771,1.08320445,1.09125322,1.09924862,
    1.10724401,1.11516631,1.12306893,1.13096102,1.13876953,1.14657804,
    1.15434964,1.16206283,1.16977602,1.17743002,1.18504676,1.19266351,
    1.20020314,1.20772241,1.21524168,1.22267051,1.23009136,1.2375122 ,
    1.24483409,1.25215566,1.25947478,1.26669627,1.27391776,1.2811391 ,
    1.2882598 ,1.29538049,1.30250119,1.30952779,1.31654704,1.32356629,
    1.33050367,1.3374209 ,1.34433812,1.35119116,1.35800583,1.3648205 ,
    1.37159424,1.37830589,1.38501753,1.39171712,1.39832532,1.40493353,
    1.41154173,1.41806863,1.42457302,1.43107741,1.43754047,1.44394073,
    1.450341  ,1.45674126,1.46304155,1.46933741,1.47563327,1.48188052,
    1.48807175,1.49426297,1.50045419,1.50654922,1.5126356 ,1.51872199,
    1.52477516,1.53075656,1.53673796,1.54271936,1.54863131,1.5545076 ,
    1.56038389,1.56626018,1.57203652,1.57780761,1.5835787 ,1.5893304 ,
    1.59499624,1.60066208,1.60632791,1.61195555,1.6175161 ,1.62307666,
    1.62863721,1.63414666,1.63960193,1.6450572 ,1.65051247,1.65590966,
    1.66125967,1.66660969,1.67195971,1.67723705,1.68245768,1.68767831,
    1.69289894,1.6980371 ,1.70309735,1.70815761,1.71321787,1.71821382,
    1.72311706,1.7280203 ,1.73292354,1.73778845,1.74253792,1.7472874 ,
    1.75203687,1.75678108,1.76137975,1.76597842,1.77057709,1.77517576,
    1.77966063,1.78411118,1.78856173,1.79301227,1.79739716,1.80170201,
    1.80600686,1.81031171,1.81460451,1.81876586,1.82292721,1.82708855,
    1.8312499 ,1.83531647,1.83933628,1.84335608,1.84737589,1.85136645,
    1.85524646,1.85912648,1.8630065 ,1.86688652,1.87066917,1.87441096,
    1.87815274,1.88189453,1.88561487,1.88921979,1.89282472,1.89642964,
    1.90003456,1.90356246,1.90703171,1.91050096,1.91397021,1.91743946,
    1.92078224,1.92411684,1.92745144,1.93078604,1.93408417,1.93728498,
    1.94048578,1.94368659,1.9468874 ,1.95001241,1.95308013,1.95614785,
    1.95921557,1.96228328,1.9652408 ,1.96817599,1.97111117,1.97404636,
    1.97697359,1.97977665,1.98257971,1.98538278,1.98818584,1.99095498,
    1.9936262 ,1.99629741,1.99896863,2.00163984,2.00425518,2.00679469,
    2.00933419,2.0118737 ,2.01441321,2.01687861,2.01928641,2.02169421,
    2.02410201,2.02650982,2.02882874,2.03110472,2.0333807 ,2.03565668,
    2.03793266,2.0401082 ,2.04225211,2.04439602,2.04653993,2.04868383,
    2.05071875,2.05273021,2.05474168,2.05675315,2.05876461,2.06066131,
    2.06253984,2.06441837,2.0662969 ,2.06817544,2.06993601,2.071681  ,
    2.07342598,2.07517097,2.07691595,2.0785422 ,2.0801529 ,2.08176361,
    2.08337431,2.08498501,2.08647844,2.087954  ,2.08942957,2.09090513,
    2.0923807 ,2.09374252,2.09508198,2.09642143,2.09776088,2.09910033,
    2.10033149,2.10153373,2.10273598,2.10393823,2.10514047,2.10624161,
    2.10730543,2.10836926,2.10943308,2.1104969 ,2.1114684 ,2.11239246,
    2.11331652,2.11424058,2.11516465,2.1160066 ,2.11678944,2.11757228,
    2.11835512,2.11913795,2.11985018,2.12049021,2.12113023,2.12177026,
    2.12241029,2.12299232,2.12348781,2.12398331,2.12447881,2.12497431,
    2.12542538,2.1257745 ,2.12612362,2.12647274,2.12682186,2.1271409 ,
    2.12734166,2.12754243,2.12774319,2.12794395,2.12812958,2.12817986,
    2.12823015,2.12828044,2.12833072,2.12838101,2.12827877,2.12817632,
    2.12807387,2.12797142,2.12786897,2.12762713,2.12736953,2.12711194,
    2.12685435,2.12659675,2.12621274,2.12579745,2.12538215,2.12496686,
    2.12455157,2.12402241,2.12344671,2.12287101,2.1222953 ,2.1217196 ,
    2.12104191,2.12030293,2.11956394,2.11882495,2.11808596,2.11725591,
    2.1163506 ,2.11544529,2.11453997,2.11363466,2.11264794,2.11157308,
    2.11049822,2.10942337,2.10834851,2.1072003 ,2.1059525 ,2.1047047 ,
    2.1034569 ,2.1022091 ,2.10089403,2.0994697 ,2.09804537,2.09662104,
    2.09519671,2.0937088 ,2.09210415,2.0904995 ,2.08889485,2.08729019,
    2.08562284,2.08383387,2.0820449 ,2.08025592,2.07846695,2.0766155 ,
    2.07464571,2.07267592,2.07070613,2.06873633,2.06670806,2.06458147,
    2.06245488,2.06032828,2.05820169,2.05600877,2.05372893,2.05144909,
    2.04916925,2.0468894 ,2.04452974,2.04209826,2.03966678,2.0372353 ,
    2.03480382,2.03227417,2.02969264,2.02711111,2.02452958,2.02194804,
    2.0192452 ,2.01651519,2.01378518,2.01105517,2.00832293,2.00544599,
    2.00256906,1.99969212,1.99681519,1.99390199,1.99087967,1.98785735,
    1.98483502,1.9818127 ,1.97871554,1.97554935,1.97238316,1.96921697,
    1.96605078,1.9627667 ,1.95945815,1.95614959,1.95284104,1.94950803,
    1.9460586 ,1.94260917,1.93915974,1.93571032,1.93218318,1.92859435,
    1.92500553,1.9214167 ,1.91782787,1.91410386,1.91037709,1.90665033,
    1.90292356,1.89913646,1.8952732 ,1.89140994,1.88754668,1.88368343,
    1.87969268,1.87569437,1.87169605,1.86769773,1.86363442,1.85950247,
    1.85537052,1.85123857,1.84709739,1.84283321,1.83856904,1.83430486,
    1.83004068,1.82568559,1.82129059,1.81689559,1.81250059,1.80805871,
    1.80353428,1.79900984,1.79448541,1.78995172,1.78529923,1.78064674,
    1.77599425,1.77134176,1.76658468,1.76180551,1.75702633,1.75224716,
    1.74738994,1.74248544,1.73758095,1.73267646,1.72771438,1.72268593,
    1.71765747,1.71262901,1.7075575 ,1.70240643,1.69725536,1.69210429,
    1.68691886,1.68164651,1.67637417,1.67110183,1.66579811,1.66040583,
    1.65501356,1.64962129,1.64419501,1.63868414,1.63317328,1.62766241,
    1.62210945,1.61648131,1.61085318,1.60522505,1.59954138,1.59379732,
    1.58805326,1.5823092 ,1.57649093,1.57063226,1.5647736 ,1.55891494,
    1.5529583 ,1.54698637,1.54101443,1.53502774,1.52894387,1.52285999,
    1.51677612,1.51064261,1.50444814,1.49825366,1.49205919,1.48577551,
    1.47947177,1.47316804,1.46683892,1.46042726,1.45401561,1.44760396,
    1.4411171 ,1.43459889,1.42808067,1.42153841,1.41491499,1.40829158,
    1.40166817,1.39496124,1.38823401,1.38150677,1.37473512,1.36790545,
    1.36107578,1.35423432,1.34730361,1.34037291,1.3334422 ,1.32642639,
    1.31939608,1.31236577,1.30527191,1.29814344,1.29101496,1.28383853,
    1.27661336,1.26938819,1.26212487,1.2548045 ,1.24748413,1.24012982,
    1.23271578,1.22530174,1.21785258,1.21034643,1.20284028,1.19529266,
    1.187696  ,1.18009935,1.17244995,1.16476445,1.15707894,1.14932939,
    1.14156263,1.13379324,1.12596105,1.11812886,1.11027425,1.10237694,
    1.09447963,1.08653645,1.07857356,1.07060373,1.06257482,1.05454592,
    1.04647975,1.03838441,1.03028467,1.0221225 ,1.01396033,1.00575653,
    0.99752716,0.98928285,0.98098595,0.97268905,0.96433308,0.95596835,
    0.947565  ,0.93913215,0.93067816,0.92217698,0.91366909,0.90509939,
    0.89652968,0.88789604,0.87925766,0.87056367,0.86185654,0.85309915,
    0.84432323,0.83549947,0.82665479,0.8177618 ,0.80884845,0.7998835 ,
    0.79090163,0.78186213,0.772812  ,0.76369553,0.75456742,0.74538181,
    0.73617205,0.72691943,0.71762622,0.7083072 ,0.69892904,0.6895296 ,
    0.68008012,0.67059285,0.66107957,0.65150404,0.64190167,0.63226434,
    0.62257291,0.61285177,0.60309623,0.59328699,0.58344636,0.57357295,
    0.56365439,0.55369734,0.54370815,0.5336864 ,0.52362735,0.5135324 ,
    0.50340878,0.49325756,0.4830803 ,0.47287911,0.46265653,0.45241507,
    0.44216125,0.43189963,0.4216356 ,0.41137528,0.40110972,0.39082618,
    0.38052127,0.37019377,0.35984126,0.34946447,0.33906258,0.32863476,
    0.31817944,0.30769185,0.29717534,0.28662677,0.27604267,0.2654247 ,
    0.25476773,0.24407172,0.23333566,0.22255611,0.21173088,0.20085995,
    0.18994192,0.17897563,0.1679644 ,0.15691076,0.14582273,0.13471317,
    0.12360484,0.11251389,0.10140551,0.09027214,0.07910874,0.0679086 ,
    0.05666416,0.04536823,0.03403101,0.02268333,0.01126042,0.};
  *p = evenly_phi;
  *dp = evenly_dphi;
}

}  // namespace PhaseTracer
