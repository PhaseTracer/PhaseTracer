// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Tue 17 Nov 2020 15:32:58

#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_two_scale_matching.hpp"
#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching.hpp"
#include "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_two_scale_model.hpp"
#include "standard_model_two_scale_model.hpp"
#include "error.hpp"
#include "single_scale_constraint.hpp"

namespace flexiblesusy {

#define CLASSNAME ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching_up<Two_scale>

CLASSNAME::ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching_up(
   standard_model::StandardModel<Two_scale>* low_,
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>* high_,
   const Scale_getter& scale_getter_,
   int loop_order_,
   int higgs_idx_)
   : model(high_)
   , eft(low_)
   , scale_getter(scale_getter_)
   , scale(0.)
   , higgs_idx(higgs_idx_)
{
   set_loop_order(loop_order_);
}

double CLASSNAME::get_scale() const
{
   if (scale != 0.)
      return scale;

   if (!scale_getter)
      throw SetupError("Scale getter is not set!");

   return scale_getter();
}

void CLASSNAME::set_models(Model* low, Model* high)
{
   eft = cast_model<standard_model::StandardModel<Two_scale>*>(low);
   model = cast_model<ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>*>(high);
}

void CLASSNAME::match()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   if (model->get_thresholds() && loop_order)
      ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching::match_low_to_high_scale_model(*model, *eft, loop_order, higgs_idx);
   else
      ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching::match_low_to_high_scale_model_tree_level(*model, *eft);
}

void CLASSNAME::match_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching::match_low_to_high_scale_model_tree_level(*model, *eft);
}

int CLASSNAME::get_higgs_index() const
{
   return higgs_idx;
}

int CLASSNAME::get_loop_order() const
{
   return loop_order;
}

void CLASSNAME::set_higgs_index(int idx)
{
   higgs_idx = idx;
}

void CLASSNAME::set_loop_order(int loop_order_)
{
   loop_order = loop_order_;
}

void CLASSNAME::set_scale(const Scale_getter& scale_getter_)
{
   scale_getter = scale_getter_;
}

void CLASSNAME::set_scale(Scale_getter&& scale_getter_)
{
   scale_getter = std::move(scale_getter_);
}

void CLASSNAME::set_scale(double scale_)
{
   scale = scale_;
}

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define CLASSNAME ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching_down<Two_scale>

CLASSNAME::ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching_down(
   standard_model::StandardModel<Two_scale>* low_,
   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>* high_,
   const Scale_getter& scale_getter_,
   int loop_order_,
   int higgs_idx_)
   : model(high_)
   , eft(low_)
   , scale_getter(scale_getter_)
   , scale(0.)
   , higgs_idx(higgs_idx_)
{
   set_loop_order(loop_order_);
}

double CLASSNAME::get_scale() const
{
   if (scale != 0.)
      return scale;

   if (!scale_getter)
      throw SetupError("Scale getter is not set!");

   return scale_getter();
}

void CLASSNAME::set_models(Model* high, Model* low)
{
   eft = cast_model<standard_model::StandardModel<Two_scale>*>(low);
   model = cast_model<ScalarSingletZ2DMEWSBoutputlamHEFTHiggs<Two_scale>*>(high);
}

void CLASSNAME::match()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   if (model->get_thresholds() && loop_order)
      ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching::match_high_to_low_scale_model(*eft, *model, loop_order, higgs_idx);
   else
      ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching::match_high_to_low_scale_model_tree_level(*eft, *model, higgs_idx);
}

void CLASSNAME::match_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_standard_model_matching::match_high_to_low_scale_model_tree_level(*eft, *model, higgs_idx);
}

int CLASSNAME::get_higgs_index() const
{
   return higgs_idx;
}

int CLASSNAME::get_loop_order() const
{
   return loop_order;
}

void CLASSNAME::set_higgs_index(int idx)
{
   higgs_idx = idx;
}

void CLASSNAME::set_loop_order(int loop_order_)
{
   loop_order = loop_order_;
}

void CLASSNAME::set_scale(const Scale_getter& scale_getter_)
{
   scale_getter = scale_getter_;
}

void CLASSNAME::set_scale(Scale_getter&& scale_getter_)
{
   scale_getter = std::move(scale_getter_);
}

void CLASSNAME::set_scale(double scale_)
{
   scale = scale_;
}

} // namespace flexiblesusy
