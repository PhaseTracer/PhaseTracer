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
#include <boost/math/tools/roots.hpp>

#include "transition_finder.hpp"
#include "logger.hpp"

namespace PhaseTracer {

std::vector<Transition> TransitionFinder::find_transition(Phase phase1, Phase phase2, double T1, double T2) const {
  if (T1 > T2) {
    LOG(debug) << "Phases do not overlap in temperature - no critical temperature";
    return {{NON_OVERLAPPING_T}};
  }

  try {
    const auto f = [this, phase1, phase2](double T) {return this->pf.delta_potential_at_T(phase1, phase2, T);};

    const double root_bits = 1. - std::log2(TC_tol_rel);
    boost::math::tools::eps_tolerance<double> stop(root_bits);
    boost::uintmax_t non_const_max_iter = max_iter;
    const auto result = boost::math::tools::toms748_solve(f, T1, T2, stop, non_const_max_iter);
    const double TC = (result.first + result.second) * 0.5;

    LOG(debug) << "Found critical temperature = " << TC;

    const bool ordered = f(T1) < 0.;
    if (!ordered) {
      std::swap(phase1, phase2);
    }
    const auto phase1_at_critical = pf.phase_at_T(phase1, TC);
    const auto phase2_at_critical = pf.phase_at_T(phase2, TC);
    const auto delta_potential = phase1_at_critical.potential - phase2_at_critical.potential;
    const auto true_vacua = pf.symmetric_partners(phase1_at_critical.x);
    const auto false_vacua =  pf.symmetric_partners(phase2_at_critical.x);
    std::vector<Transition> unique_transitions;
    // Fix false_vacuum, and loop all possible symmetric partners of true_vacuum
    const auto false_vacuum = false_vacua[0];
    auto true_vacuum = true_vacua[0];
    for (size_t i_unique=0; i_unique < true_vacua.size(); i_unique++) {
      true_vacuum = true_vacua[i_unique];
      bool duplicate = false;
      const auto fv = pf.symmetric_partners(false_vacuum);
      const auto tv = pf.symmetric_partners(true_vacuum);
      LOG(trace) << "False vacuum: " << false_vacuum;
      LOG(trace) << "and true vacuum: " << true_vacuum;
      for (const auto tran : unique_transitions) {
        for (size_t i=0; i< fv.size(); ++i) {
          if ((fv[i]-tran.false_vacuum).norm() < change_abs_tol + change_rel_tol * fv[i].norm()
            && (tv[i]-tran.true_vacuum).norm() < change_abs_tol + change_rel_tol * tv[i].norm()) {
            LOG(trace) << "are duplicate";
            duplicate = true;
            break;
           }
        }
      }
      if (!duplicate) {
        LOG(trace) << "are not duplicate";
        const auto gamma_ = gamma(true_vacuum, false_vacuum, TC);
        const auto changed_ = changed(true_vacuum, false_vacuum);

        double TN = std::numeric_limits<double>::quiet_NaN();
        auto true_vacuum_TN = true_vacuum;
        auto false_vacuum_TN = false_vacuum;
        double action_TN = std::numeric_limits<double>::quiet_NaN();
        if (false){ // Need BubbleProfiler
          TN = get_Tnuc(phase1, phase2, i_unique, TC, T1);
          auto vacua = get_vacua_at_T(phase1, phase2, TN, i_unique);
          true_vacuum_TN = vacua[1];
          false_vacuum_TN = vacua[0];
//          false_vacuum_TN(0) = -1E-4;
          action_TN = get_action(true_vacuum_TN,false_vacuum_TN,TN);
        }
        unique_transitions.push_back({SUCCESS, TC, phase1, phase2,
          true_vacuum, false_vacuum, gamma_, changed_, delta_potential,
          TN, true_vacuum_TN, false_vacuum_TN, action_TN,
          unique_transitions.size()});


      }
    }


    return unique_transitions;
  } catch (const std::exception& e) {
    // TODO
    LOG(debug) << e.what() << " - probably no sign change between T = " << T1 << " and " << T2;
    return {{ERROR}};
  }
}

double TransitionFinder::get_Tnuc(Phase phase1, Phase phase2, size_t i_unique, double T_begin, double T_end) const{
  
//    print("Tunneling from phase %s to phase %s at T=%0.4g"
  
  if (T_begin < T_end) {
    LOG(fatal) << "T_begin < T_end, so swith the values. ";
    T_begin = T_begin + T_end;
    T_end   = T_begin - T_end;
    T_begin = T_begin - T_end;
  }
  
  LOG(debug) << "Find Tnuc between " << phase1.key << " and " << phase2.key << " in [" << T_begin << ", "<< T_end <<"]. ";
  
  
  const auto fun_nucleation = [this, phase1, phase2, i_unique](double Ttry) {
    return this->get_action(phase1, phase2, Ttry, i_unique)/Ttry - 140.;
  };
  
  LOG(debug) << "fun_nucleation(T_begin)= " << fun_nucleation(T_begin) ;
  LOG(debug) << "fun_nucleation(T_end)= " << fun_nucleation(T_end) ;
  
  if ( fun_nucleation(T_begin) < 0 ) {
    LOG(debug) << "The tunneling possibility at T_begin satisfys the nucleation condition." ;
    return T_begin;
  }
  double Tnuc = std::numeric_limits<double>::quiet_NaN();
  if ( fun_nucleation(T_end) > 0) {
    double Tnuc_try = T_begin;
    while (Tnuc_try>T_end) {
      Tnuc_try -= Tnuc_step;
      while (fun_nucleation(Tnuc_try) < 0) break;
    }
  }
  
  try{
    const double root_bits = 1. - std::log2(TC_tol_rel);
    boost::math::tools::eps_tolerance<double> stop(root_bits);
    boost::uintmax_t non_const_max_iter = max_iter;
    const auto result = boost::math::tools::bisect(fun_nucleation, T_end, T_begin, stop, non_const_max_iter);
    Tnuc = (result.first + result.second) * 0.5;
    LOG(debug) << "Found nucleation temperature = " << Tnuc;
    return Tnuc;
  } catch(char *str){
    std::cout << str << std::endl; // TODO
    // find the first
  }
  

//
//
//  for (double Ttry = T_begin ; Ttry > T_end; Ttry--){
//    LOG(debug) << "Cal action at " << Ttry << ". ";
//
//  }
}

std::vector<Transition> TransitionFinder::divide_and_find_transition(const Phase& phase1, const Phase& phase2, double T1, double T2) const {

#ifdef CAL_TNUC
  LOG(warning) << "When assume_only_one_transition is set to false, the calculation of Tnuc will not be performed. " << std::endl;
#endif
  
  std::vector<Transition> roots;
  for (double T = T1; T <= T2; T += separation) {
    const auto bs = find_transition(phase1, phase2, T, std::min(T + separation, T2));
    for (const auto &b : bs) {
      roots.push_back(b);
    }
  }
  return roots;
}

std::vector<Eigen::VectorXd> TransitionFinder::get_vacua_at_T(Phase phase1, Phase phase2, double T, size_t i_unique)const{
  const auto phase1_at_T = pf.phase_at_T(phase1, T);
  const auto phase2_at_T = pf.phase_at_T(phase2, T);
  const auto true_vacua_at_T = pf.symmetric_partners(phase1_at_T.x);
  const auto false_vacua_at_T =  pf.symmetric_partners(phase2_at_T.x);
  
  return {false_vacua_at_T[0], true_vacua_at_T[i_unique]};
}

double TransitionFinder::get_action(Phase phase1, Phase phase2, double T, size_t i_unique) const{
  const auto vacua = get_vacua_at_T(phase1, phase2, T, i_unique);
  return get_action(vacua[0], vacua[1], T);
}

double TransitionFinder::get_action(const Eigen::VectorXd& vacuum_1, const Eigen::VectorXd& vacuum_2, double T) const{
    double action=std::numeric_limits<double>::quiet_NaN();
    V_BubbleProfiler V_BP_=V_BP; // perturbative_profiler only accept non-const potential
    V_BP_.set_T(T); // This is necessary!!
    
    Eigen::VectorXd true_vacuum = vacuum_1;
    Eigen::VectorXd false_vacuum = vacuum_2;
    if ( pf.P.V(true_vacuum,T) > pf.P.V(false_vacuum,T) ) true_vacuum.swap(false_vacuum);
    
    size_t n_dims = 4;
    
    bool use_perturbative = false;
    if (n_fields == 1  && !use_perturbative) {

      double false_min = false_vacuum[0];
      double true_min = true_vacuum[0];
      auto barrier_ = V_BP_.find_one_dimensional_barrier( true_vacuum, false_vacuum, T);
      double barrier = barrier_[0];

      LOG(debug) << "Calculate action at T=" << T << ", with false, true, barrier = " << false_min << ", " << true_min << ", " << barrier;

      try{
        BubbleProfiler::Shooting one_dim;
        one_dim.solve(V_BP_, false_min,
                      true_min,barrier, n_dims, BubbleProfiler::Shooting::Solver_options::Compute_action);
        action = one_dim.get_euclidean_action();
      }catch (const std::exception& e) {
        LOG(warning) << "At T=" << T << ", between[" << false_min << "] and [" << true_min << "]: "   << e.what() << std::endl;
      }
    } else {
      LOG(debug) << "Calculate action at T=" << T << ", between [" << false_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "] and [" << true_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "]";
      BubbleProfiler::RK4_perturbative_profiler profiler;
      
//      profiler.set_domain_start(input.domain_start);
//      profiler.set_domain_end(input.domain_end);
      profiler.set_initial_step_size(1.e-2);
      profiler.set_interpolation_points_fraction(1.0);
      
      
      profiler.set_false_vacuum_loc(false_vacuum);
      profiler.set_true_vacuum_loc(true_vacuum);
      n_dims = 3;
      profiler.set_number_of_dimensions(n_dims);
      auto root_finder = std::make_shared<BubbleProfiler::GSL_root_finder<Eigen::Dynamic> >();
      profiler.set_root_finder(root_finder);
      std::shared_ptr<BubbleProfiler::Profile_guesser> guesser;
      guesser = std::make_shared<BubbleProfiler::Kink_profile_guesser>();
      profiler.set_initial_guesser(guesser);
      auto convergence_tester = std::make_shared<BubbleProfiler::Relative_convergence_tester>(
                                1.e-3, 1.e-3);
      profiler.set_convergence_tester(convergence_tester);
      
      try{
        profiler.calculate_bubble_profile(V_BP_);
        action = profiler.get_euclidean_action();
      }catch (const std::exception& e) {
        LOG(warning) << "At T=" << T <<  ", between [" << false_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "] and [" << true_vacuum.transpose().format(Eigen::IOFormat(4, Eigen::DontAlignCols, " ", " ")) << "]: "   << e.what() << std::endl;
      }
    }
    
    LOG(debug) << " S = " << action << std::endl;
    
    return action;
    
  }


double TransitionFinder::gamma(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd& false_vacuum, const double TC) const {
  const int b = true_vacuum.size() + 1;
  const int items = (b + (n_ew_scalars % b)) % b;
  return (true_vacuum-false_vacuum).head(items).norm()/TC;
}

std::vector<bool> TransitionFinder::changed(const Eigen::VectorXd& true_vacuum, const Eigen::VectorXd& false_vacuum) const {
  auto delta = true_vacuum-false_vacuum;
  std::vector<bool> changed_;
  for (unsigned int i = 0; i < delta.size(); i += 1) {
    changed_.push_back(std::abs(delta[i])> change_abs_tol + change_rel_tol * std::max(true_vacuum.norm(), false_vacuum.norm()));
  }
  return changed_;
}

bool TransitionFinder::phases_overlaped(const Phase& phase1, const Phase& phase2, double T) const {
  return pf.identical_within_tol(pf.phase_at_T(phase1, T).x, pf.phase_at_T(phase2, T).x);
}

std::vector<double> TransitionFinder::get_un_overlapped_T_range(const Phase& phase1, const Phase& phase2, double T1, double T2) const {
  std::vector<double> T_range;

  if (!(phases_overlaped(phase1, phase2, T1) || phases_overlaped(phase1, phase2, T2))) {
    T_range.push_back(T1);
    T_range.push_back(T2);
    return T_range;
  }

  if ( phases_overlaped(phase1, phase2, T1) && phases_overlaped(phase1, phase2, T2) ) {
    LOG(fatal) << "Phases " << phase1.key << " and " << phase2.key << " are redundant";
    return T_range;
  }

  double T_same = phases_overlaped(phase1, phase2, T1) ? T1 : T2;
  double T_diff = phases_overlaped(phase1, phase2, T2) ? T1 : T2;
  T_range.push_back(T_diff);

  while (std::abs(T_diff-T_same) > TC_tol_rel) {
    double T_mid = 0.5*(T_same+T_diff);
    if (phases_overlaped(phase1, phase2, T_mid)) {
      T_same = T_mid;
    } else {
      T_diff = T_mid;
    }
  }
  T_range.push_back(T_diff);
  T_range.push_back(T_same);

  return T_range;
}

void TransitionFinder::find_transitions() {
  if (calculated_transitions) {
    return;
  }

  transitions.empty();
  calculated_transitions = true;
  const auto phases = pf.get_phases();

  for (const auto &phase1 : phases) {
    for (const auto &phase2 : phases) {
      // Don't investigate n1 -> n1 or n1 -> n2 as well as n2 -> n1
      if (phase1.key >= phase2.key) {
        continue;
      }

      LOG(debug) << "Finding critical temperatures between phases "
                 << phase1.key << " and " << phase2.key;


      double tmax = std::min(phase1.T.back(), phase2.T.back());
      double tmin = std::max(phase1.T.front(), phase2.T.front());

      if (tmin > tmax) {
        LOG(debug) << "Phases do not overlap in temperature - no critical temperature";
        continue;
      }

      auto T_range = get_un_overlapped_T_range(phase1, phase2, tmin, tmax);

      if (T_range.size() == 0) continue;
      if (T_range.size() == 3) {
        const auto TC = T_range[2];
        const auto TN = T_range[2];
        LOG(debug) << "Phases " << phase1.key << " and " << phase2.key << " cross at T=" << TC;
        const auto phase1_at_critical = pf.phase_at_T(phase1, TC);
        const auto phase2_at_critical = pf.phase_at_T(phase2, TC);
        const double gamma_ = 0.;
        auto changed_ = changed(phase1_at_critical.x, phase2_at_critical.x);
        Transition f = {SUCCESS, TC, phase1, phase2,
                        phase1_at_critical.x, phase2_at_critical.x,
                        gamma_, changed_,
                        phase1_at_critical.potential - phase2_at_critical.potential,
                        TC, phase1_at_critical.x, phase2_at_critical.x, 0
        };
        transitions.push_back(f);
      }

      const auto found = assume_only_one_transition ?
        find_transition(phase1, phase2, T_range[0], T_range[1]) :
        divide_and_find_transition(phase1, phase2, T_range[0], T_range[1]);

      for (const auto &f : found) {
        if (f.message == SUCCESS) {
            transitions.push_back(f);
        }
      }
    }
  }

  // Finally sort transitions by critical temperature
  std::sort(transitions.begin(), transitions.end(),
    [](const Transition &a, const Transition &b) {return a.TC < b.TC;});
}

std::ostream& operator << (std::ostream& o, const TransitionFinder& a) {
  if (!a.calculated_transitions) {
    o << "no transitions calculated yet" << std::endl << std::endl;
  } else if (a.transitions.empty()) {
    o << "found no transitions" << std::endl << std::endl;
  } else {
    o << "found " << a.transitions.size() << " transition";
    if (a.transitions.size() > 1) {
      o << "s";
    }
    o << std::endl << std::endl;
    for (const auto &t : a.transitions) {
      o << t << std::endl;
    }
  }
  return o;
}
}  // namespace PhaseTracer
