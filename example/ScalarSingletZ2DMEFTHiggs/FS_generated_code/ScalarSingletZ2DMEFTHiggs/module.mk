DIR          := models/ScalarSingletZ2DMEFTHiggs
MODNAME      := ScalarSingletZ2DMEFTHiggs
SARAH_MODEL  := ScalarSingletZ2DM
WITH_$(MODNAME) := yes
MODScalarSingletZ2DMEFTHiggs_MOD := SM
MODScalarSingletZ2DMEFTHiggs_DEP := $(patsubst %,model_specific/%,$(MODScalarSingletZ2DMEFTHiggs_MOD))
MODScalarSingletZ2DMEFTHiggs_INC := $(patsubst %,-Imodel_specific/%,$(MODScalarSingletZ2DMEFTHiggs_MOD))
MODScalarSingletZ2DMEFTHiggs_LIB := $(foreach M,$(MODScalarSingletZ2DMEFTHiggs_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODScalarSingletZ2DMEFTHiggs_SUBMOD  := $(DIR)/cxx_qft
MODScalarSingletZ2DMEFTHiggs_SUBMOD_INC := $(patsubst %,-I%,$(MODScalarSingletZ2DMEFTHiggs_SUBMOD))

ScalarSingletZ2DMEFTHiggs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
ScalarSingletZ2DMEFTHiggs_INSTALL_CXXQFT_DIR := \
		$(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)/cxx_qft

ScalarSingletZ2DMEFTHiggs_MK     := \
		$(DIR)/module.mk

ScalarSingletZ2DMEFTHiggs_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

ScalarSingletZ2DMEFTHiggs_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

ScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(ScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_MK)
LIBScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_SRC ?= ''

ScalarSingletZ2DMEFTHiggs_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

ScalarSingletZ2DMEFTHiggs_INCLUDE_MK := \
		$(ScalarSingletZ2DMEFTHiggs_SUSY_BETAS_MK) \
		$(ScalarSingletZ2DMEFTHiggs_SOFT_BETAS_MK)

ScalarSingletZ2DMEFTHiggs_SLHA_INPUT := \
		$(DIR)/LesHouches.in.ScalarSingletZ2DMEFTHiggs_generated \
		$(DIR)/LesHouches.in.ScalarSingletZ2DMEFTHiggs

ScalarSingletZ2DMEFTHiggs_REFERENCES := \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_references.tex

ScalarSingletZ2DMEFTHiggs_GNUPLOT := \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_plot_rgflow.gnuplot \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_plot_spectrum.gnuplot

ScalarSingletZ2DMEFTHiggs_TARBALL := \
		$(MODNAME).tar.gz

LIBScalarSingletZ2DMEFTHiggs_SRC := \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_a_muon.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_edm.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_FFV_form_factors.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_f_to_f_conversion.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_l_to_lgamma.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_b_to_s_gamma.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_effective_couplings.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_info.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_input_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_mass_eigenstates.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_model_slha.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_observables.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_physical.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_slha_io.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_soft_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_susy_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_utilities.cpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_weinberg_angle.cpp

LIBScalarSingletZ2DMEFTHiggs_SRC += $(LIBScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_SRC)

EXEScalarSingletZ2DMEFTHiggs_SRC := \
		$(DIR)/run_ScalarSingletZ2DMEFTHiggs.cpp \
		$(DIR)/run_cmd_line_ScalarSingletZ2DMEFTHiggs.cpp \
		$(DIR)/scan_ScalarSingletZ2DMEFTHiggs.cpp
LLScalarSingletZ2DMEFTHiggs_LIB  :=
LLScalarSingletZ2DMEFTHiggs_OBJ  :=
LLScalarSingletZ2DMEFTHiggs_SRC  := \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_librarylink.cpp

LLScalarSingletZ2DMEFTHiggs_MMA  := \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_librarylink.m \
		$(DIR)/run_ScalarSingletZ2DMEFTHiggs.m

LIBScalarSingletZ2DMEFTHiggs_HDR := \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_a_muon.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_convergence_tester.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_edm.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_FFV_form_factors.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_f_to_f_conversion.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_l_to_lgamma.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_b_to_s_gamma.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_effective_couplings.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_ewsb_solver.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_ewsb_solver_interface.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_high_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_info.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_initial_guesser.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_input_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_low_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_mass_eigenstates.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_mass_eigenstates_interface.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_model.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_model_slha.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_observables.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_physical.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_slha_io.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_spectrum_generator.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_spectrum_generator_interface.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_soft_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_susy_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_susy_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_utilities.hpp \
		$(DIR)/ScalarSingletZ2DMEFTHiggs_weinberg_angle.hpp

LIBScalarSingletZ2DMEFTHiggs_CXXQFT_HDR := \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEFTHiggs_qft.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEFTHiggs_fields.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEFTHiggs_vertices.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEFTHiggs_context_base.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEFTHiggs_npointfunctions_wilsoncoeffs.hpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
-include $(DIR)/two_scale.mk
endif
ifneq ($(findstring lattice,$(SOLVERS)),)
-include $(DIR)/lattice.mk
endif
ifneq ($(findstring semi_analytic,$(SOLVERS)),)
-include $(DIR)/semi_analytic.mk
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(ScalarSingletZ2DMEFTHiggs_SUSY_BETAS_MK)
-include $(ScalarSingletZ2DMEFTHiggs_SOFT_BETAS_MK)
-include $(ScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_MK)
-include $(ScalarSingletZ2DMEFTHiggs_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(ScalarSingletZ2DMEFTHiggs_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMEFTHiggs_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMEFTHiggs_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

# remove duplicates in case all solvers are used
LIBScalarSingletZ2DMEFTHiggs_SRC := $(sort $(LIBScalarSingletZ2DMEFTHiggs_SRC))
EXEScalarSingletZ2DMEFTHiggs_SRC := $(sort $(EXEScalarSingletZ2DMEFTHiggs_SRC))

LIBScalarSingletZ2DMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBScalarSingletZ2DMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBScalarSingletZ2DMEFTHiggs_SRC)))

EXEScalarSingletZ2DMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEScalarSingletZ2DMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEScalarSingletZ2DMEFTHiggs_SRC)))

EXEScalarSingletZ2DMEFTHiggs_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEScalarSingletZ2DMEFTHiggs_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEScalarSingletZ2DMEFTHiggs_SRC)))

LIBScalarSingletZ2DMEFTHiggs_DEP := \
		$(LIBScalarSingletZ2DMEFTHiggs_OBJ:.o=.d)

EXEScalarSingletZ2DMEFTHiggs_DEP := \
		$(EXEScalarSingletZ2DMEFTHiggs_OBJ:.o=.d)

LLScalarSingletZ2DMEFTHiggs_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLScalarSingletZ2DMEFTHiggs_SRC)))

LLScalarSingletZ2DMEFTHiggs_OBJ  := $(LLScalarSingletZ2DMEFTHiggs_SRC:.cpp=.o)
LLScalarSingletZ2DMEFTHiggs_LIB  := $(LLScalarSingletZ2DMEFTHiggs_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBScalarSingletZ2DMEFTHiggs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_ScalarSingletZ2DMEFTHiggs := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_ScalarSingletZ2DMEFTHiggs := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBScalarSingletZ2DMEFTHiggs) $(EXEScalarSingletZ2DMEFTHiggs_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)install -d $(ScalarSingletZ2DMEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMEFTHiggs_SRC) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_SRC) $(ScalarSingletZ2DMEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMEFTHiggs_HDR) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMEFTHiggs_CXXQFT_HDR) $(ScalarSingletZ2DMEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEScalarSingletZ2DMEFTHiggs_SRC) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMEFTHiggs_SRC) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMEFTHiggs_MMA) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(ScalarSingletZ2DMEFTHiggs_MK) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEFTHiggs_INCLUDE_MK) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_MK) $(ScalarSingletZ2DMEFTHiggs_INSTALL_CXXQFT_DIR)

ifneq ($(ScalarSingletZ2DMEFTHiggs_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEFTHiggs_SLHA_INPUT) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEFTHiggs_REFERENCES) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEFTHiggs_GNUPLOT) $(ScalarSingletZ2DMEFTHiggs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEFTHiggs_DEP)
		$(Q)-rm -f $(EXEScalarSingletZ2DMEFTHiggs_DEP)
		$(Q)-rm -f $(LLScalarSingletZ2DMEFTHiggs_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEFTHiggs)
		$(Q)-rm -f $(LLScalarSingletZ2DMEFTHiggs_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEFTHiggs_OBJ)
		$(Q)-rm -f $(EXEScalarSingletZ2DMEFTHiggs_OBJ)
		$(Q)-rm -f $(LLScalarSingletZ2DMEFTHiggs_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEFTHiggs_SRC)
		$(Q)-rm -f $(LIBScalarSingletZ2DMEFTHiggs_HDR)
		$(Q)-rm -f $(LIBScalarSingletZ2DMEFTHiggs_CXXQFT_HDR)
		$(Q)-rm -f $(EXEScalarSingletZ2DMEFTHiggs_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMEFTHiggs_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMEFTHiggs_MMA)
		$(Q)-rm -f $(METACODE_STAMP_ScalarSingletZ2DMEFTHiggs)
		$(Q)-rm -f $(ScalarSingletZ2DMEFTHiggs_INCLUDE_MK)
		$(Q)-rm -f $(ScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(ScalarSingletZ2DMEFTHiggs_SLHA_INPUT)
		$(Q)-rm -f $(ScalarSingletZ2DMEFTHiggs_REFERENCES)
		$(Q)-rm -f $(ScalarSingletZ2DMEFTHiggs_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEScalarSingletZ2DMEFTHiggs_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(ScalarSingletZ2DMEFTHiggs_TARBALL) \
		$(LIBScalarSingletZ2DMEFTHiggs_SRC) $(LIBScalarSingletZ2DMEFTHiggs_HDR) $(LIBScalarSingletZ2DMEFTHiggs_CXXQFT_HDR) \
		$(EXEScalarSingletZ2DMEFTHiggs_SRC) \
		$(LLScalarSingletZ2DMEFTHiggs_SRC) $(LLScalarSingletZ2DMEFTHiggs_MMA) \
		$(ScalarSingletZ2DMEFTHiggs_MK) $(ScalarSingletZ2DMEFTHiggs_INCLUDE_MK) $(ScalarSingletZ2DMEFTHiggs_CXXQFT_VERTICES_MK) \
		$(ScalarSingletZ2DMEFTHiggs_SLHA_INPUT) $(ScalarSingletZ2DMEFTHiggs_REFERENCES) \
		$(ScalarSingletZ2DMEFTHiggs_GNUPLOT)

$(LIBScalarSingletZ2DMEFTHiggs_SRC) $(LIBScalarSingletZ2DMEFTHiggs_HDR) $(LIBScalarSingletZ2DMEFTHiggs_CXXQFT_HDR) $(EXEScalarSingletZ2DMEFTHiggs_SRC) $(LLScalarSingletZ2DMEFTHiggs_SRC) $(LLScalarSingletZ2DMEFTHiggs_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_ScalarSingletZ2DMEFTHiggs)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_ScalarSingletZ2DMEFTHiggs): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_ScalarSingletZ2DMEFTHiggs)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_ScalarSingletZ2DMEFTHiggs)"
		@echo "Note: to regenerate ScalarSingletZ2DMEFTHiggs source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_ScalarSingletZ2DMEFTHiggs)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_ScalarSingletZ2DMEFTHiggs):
		@true
endif

$(LIBScalarSingletZ2DMEFTHiggs_DEP) $(EXEScalarSingletZ2DMEFTHiggs_DEP) $(LLScalarSingletZ2DMEFTHiggs_DEP) $(LIBScalarSingletZ2DMEFTHiggs_OBJ) $(EXEScalarSingletZ2DMEFTHiggs_OBJ) $(LLScalarSingletZ2DMEFTHiggs_OBJ) $(LLScalarSingletZ2DMEFTHiggs_LIB): \
	CPPFLAGS += $(MODScalarSingletZ2DMEFTHiggs_SUBMOD_INC) $(MODScalarSingletZ2DMEFTHiggs_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)  $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBScalarSingletZ2DMEFTHiggs_DEP) $(EXEScalarSingletZ2DMEFTHiggs_DEP) $(LLScalarSingletZ2DMEFTHiggs_DEP) $(LIBScalarSingletZ2DMEFTHiggs_OBJ) $(EXEScalarSingletZ2DMEFTHiggs_OBJ) $(LLScalarSingletZ2DMEFTHiggs_OBJ) $(LLScalarSingletZ2DMEFTHiggs_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLScalarSingletZ2DMEFTHiggs_OBJ) $(LLScalarSingletZ2DMEFTHiggs_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBScalarSingletZ2DMEFTHiggs): $(LIBScalarSingletZ2DMEFTHiggs_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBScalarSingletZ2DMEFTHiggs) $(MODScalarSingletZ2DMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLScalarSingletZ2DMEFTHiggs_LIB): $(LLScalarSingletZ2DMEFTHiggs_OBJ) $(LIBScalarSingletZ2DMEFTHiggs) $(MODScalarSingletZ2DMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS)

ALLDEP += $(LIBScalarSingletZ2DMEFTHiggs_DEP) $(EXEScalarSingletZ2DMEFTHiggs_DEP)
ALLSRC += $(LIBScalarSingletZ2DMEFTHiggs_SRC) $(EXEScalarSingletZ2DMEFTHiggs_SRC)
ALLLIB += $(LIBScalarSingletZ2DMEFTHiggs)
ALLEXE += $(EXEScalarSingletZ2DMEFTHiggs_EXE)
ALLMODDEP += $(MODScalarSingletZ2DMEFTHiggs_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLScalarSingletZ2DMEFTHiggs_DEP)
ALLSRC += $(LLScalarSingletZ2DMEFTHiggs_SRC)
ALLLL  += $(LLScalarSingletZ2DMEFTHiggs_LIB)
endif
