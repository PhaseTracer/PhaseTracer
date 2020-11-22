DIR          := models/ScalarSingletZ2DMMhInput
MODNAME      := ScalarSingletZ2DMMhInput
SARAH_MODEL  := ScalarSingletZ2DM
WITH_$(MODNAME) := yes
MODScalarSingletZ2DMMhInput_MOD := SM
MODScalarSingletZ2DMMhInput_DEP := $(patsubst %,model_specific/%,$(MODScalarSingletZ2DMMhInput_MOD))
MODScalarSingletZ2DMMhInput_INC := $(patsubst %,-Imodel_specific/%,$(MODScalarSingletZ2DMMhInput_MOD))
MODScalarSingletZ2DMMhInput_LIB := $(foreach M,$(MODScalarSingletZ2DMMhInput_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODScalarSingletZ2DMMhInput_SUBMOD  := $(DIR)/cxx_qft
MODScalarSingletZ2DMMhInput_SUBMOD_INC := $(patsubst %,-I%,$(MODScalarSingletZ2DMMhInput_SUBMOD))

ScalarSingletZ2DMMhInput_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
ScalarSingletZ2DMMhInput_INSTALL_CXXQFT_DIR := \
		$(ScalarSingletZ2DMMhInput_INSTALL_DIR)/cxx_qft

ScalarSingletZ2DMMhInput_MK     := \
		$(DIR)/module.mk

ScalarSingletZ2DMMhInput_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

ScalarSingletZ2DMMhInput_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

ScalarSingletZ2DMMhInput_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(ScalarSingletZ2DMMhInput_CXXQFT_VERTICES_MK)
LIBScalarSingletZ2DMMhInput_CXXQFT_VERTICES_SRC ?= ''

ScalarSingletZ2DMMhInput_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

ScalarSingletZ2DMMhInput_INCLUDE_MK := \
		$(ScalarSingletZ2DMMhInput_SUSY_BETAS_MK) \
		$(ScalarSingletZ2DMMhInput_SOFT_BETAS_MK)

ScalarSingletZ2DMMhInput_SLHA_INPUT := \
		$(DIR)/LesHouches.in.ScalarSingletZ2DMMhInput_generated \
		$(DIR)/LesHouches.in.ScalarSingletZ2DM \
		$(DIR)/LesHouches.in.ScalarSingletZ2DMMhInput

ScalarSingletZ2DMMhInput_REFERENCES := \
		$(DIR)/ScalarSingletZ2DMMhInput_references.tex

ScalarSingletZ2DMMhInput_GNUPLOT := \
		$(DIR)/ScalarSingletZ2DMMhInput_plot_rgflow.gnuplot \
		$(DIR)/ScalarSingletZ2DMMhInput_plot_spectrum.gnuplot

ScalarSingletZ2DMMhInput_TARBALL := \
		$(MODNAME).tar.gz

LIBScalarSingletZ2DMMhInput_SRC := \
		$(DIR)/ScalarSingletZ2DMMhInput_a_muon.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_edm.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_FFV_form_factors.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_f_to_f_conversion.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_l_to_lgamma.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_b_to_s_gamma.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_effective_couplings.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_info.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_input_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_mass_eigenstates.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_model_slha.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_observables.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_physical.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_slha_io.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_soft_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_susy_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_utilities.cpp \
		$(DIR)/ScalarSingletZ2DMMhInput_weinberg_angle.cpp

LIBScalarSingletZ2DMMhInput_SRC += $(LIBScalarSingletZ2DMMhInput_CXXQFT_VERTICES_SRC)

EXEScalarSingletZ2DMMhInput_SRC := \
		$(DIR)/run_ScalarSingletZ2DMMhInput.cpp \
		$(DIR)/run_cmd_line_ScalarSingletZ2DMMhInput.cpp \
		$(DIR)/scan_ScalarSingletZ2DMMhInput.cpp
LLScalarSingletZ2DMMhInput_LIB  :=
LLScalarSingletZ2DMMhInput_OBJ  :=
LLScalarSingletZ2DMMhInput_SRC  := \
		$(DIR)/ScalarSingletZ2DMMhInput_librarylink.cpp

LLScalarSingletZ2DMMhInput_MMA  := \
		$(DIR)/ScalarSingletZ2DMMhInput_librarylink.m \
		$(DIR)/run_ScalarSingletZ2DMMhInput.m

LIBScalarSingletZ2DMMhInput_HDR := \
		$(DIR)/ScalarSingletZ2DMMhInput_a_muon.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_convergence_tester.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_edm.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_FFV_form_factors.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_f_to_f_conversion.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_l_to_lgamma.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_b_to_s_gamma.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_effective_couplings.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_ewsb_solver.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_ewsb_solver_interface.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_high_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_info.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_initial_guesser.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_input_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_low_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_mass_eigenstates.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_mass_eigenstates_interface.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_model.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_model_slha.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_observables.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_physical.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_slha_io.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_spectrum_generator.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_spectrum_generator_interface.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_soft_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_susy_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_susy_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_utilities.hpp \
		$(DIR)/ScalarSingletZ2DMMhInput_weinberg_angle.hpp

LIBScalarSingletZ2DMMhInput_CXXQFT_HDR := \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInput_qft.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInput_fields.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInput_vertices.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInput_context_base.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInput_npointfunctions_wilsoncoeffs.hpp

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
-include $(ScalarSingletZ2DMMhInput_SUSY_BETAS_MK)
-include $(ScalarSingletZ2DMMhInput_SOFT_BETAS_MK)
-include $(ScalarSingletZ2DMMhInput_CXXQFT_VERTICES_MK)
-include $(ScalarSingletZ2DMMhInput_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(ScalarSingletZ2DMMhInput_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMMhInput_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMMhInput_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMMhInput_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBScalarSingletZ2DMMhInput_SRC := $(sort $(LIBScalarSingletZ2DMMhInput_SRC))
EXEScalarSingletZ2DMMhInput_SRC := $(sort $(EXEScalarSingletZ2DMMhInput_SRC))

LIBScalarSingletZ2DMMhInput_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBScalarSingletZ2DMMhInput_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBScalarSingletZ2DMMhInput_SRC)))

EXEScalarSingletZ2DMMhInput_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEScalarSingletZ2DMMhInput_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEScalarSingletZ2DMMhInput_SRC)))

EXEScalarSingletZ2DMMhInput_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEScalarSingletZ2DMMhInput_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEScalarSingletZ2DMMhInput_SRC)))

LIBScalarSingletZ2DMMhInput_DEP := \
		$(LIBScalarSingletZ2DMMhInput_OBJ:.o=.d)

EXEScalarSingletZ2DMMhInput_DEP := \
		$(EXEScalarSingletZ2DMMhInput_OBJ:.o=.d)

LLScalarSingletZ2DMMhInput_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLScalarSingletZ2DMMhInput_SRC)))

LLScalarSingletZ2DMMhInput_OBJ  := $(LLScalarSingletZ2DMMhInput_SRC:.cpp=.o)
LLScalarSingletZ2DMMhInput_LIB  := $(LLScalarSingletZ2DMMhInput_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBScalarSingletZ2DMMhInput     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_ScalarSingletZ2DMMhInput := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_ScalarSingletZ2DMMhInput := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBScalarSingletZ2DMMhInput) $(EXEScalarSingletZ2DMMhInput_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)install -d $(ScalarSingletZ2DMMhInput_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInput_SRC) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInput_CXXQFT_VERTICES_SRC) $(ScalarSingletZ2DMMhInput_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInput_HDR) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInput_CXXQFT_HDR) $(ScalarSingletZ2DMMhInput_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEScalarSingletZ2DMMhInput_SRC) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMMhInput_SRC) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMMhInput_MMA) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(ScalarSingletZ2DMMhInput_MK) $(ScalarSingletZ2DMMhInput_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInput_INCLUDE_MK) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInput_CXXQFT_VERTICES_MK) $(ScalarSingletZ2DMMhInput_INSTALL_CXXQFT_DIR)

ifneq ($(ScalarSingletZ2DMMhInput_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInput_SLHA_INPUT) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInput_REFERENCES) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInput_GNUPLOT) $(ScalarSingletZ2DMMhInput_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInput_DEP)
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInput_DEP)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInput_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInput)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInput_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInput_OBJ)
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInput_OBJ)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInput_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInput_SRC)
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInput_HDR)
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInput_CXXQFT_HDR)
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInput_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInput_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInput_MMA)
		$(Q)-rm -f $(METACODE_STAMP_ScalarSingletZ2DMMhInput)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInput_INCLUDE_MK)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInput_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInput_SLHA_INPUT)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInput_REFERENCES)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInput_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInput_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(ScalarSingletZ2DMMhInput_TARBALL) \
		$(LIBScalarSingletZ2DMMhInput_SRC) $(LIBScalarSingletZ2DMMhInput_HDR) $(LIBScalarSingletZ2DMMhInput_CXXQFT_HDR) \
		$(EXEScalarSingletZ2DMMhInput_SRC) \
		$(LLScalarSingletZ2DMMhInput_SRC) $(LLScalarSingletZ2DMMhInput_MMA) \
		$(ScalarSingletZ2DMMhInput_MK) $(ScalarSingletZ2DMMhInput_INCLUDE_MK) $(ScalarSingletZ2DMMhInput_CXXQFT_VERTICES_MK) \
		$(ScalarSingletZ2DMMhInput_SLHA_INPUT) $(ScalarSingletZ2DMMhInput_REFERENCES) \
		$(ScalarSingletZ2DMMhInput_GNUPLOT)

$(LIBScalarSingletZ2DMMhInput_SRC) $(LIBScalarSingletZ2DMMhInput_HDR) $(LIBScalarSingletZ2DMMhInput_CXXQFT_HDR) $(EXEScalarSingletZ2DMMhInput_SRC) $(LLScalarSingletZ2DMMhInput_SRC) $(LLScalarSingletZ2DMMhInput_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_ScalarSingletZ2DMMhInput)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_ScalarSingletZ2DMMhInput): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_ScalarSingletZ2DMMhInput)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_ScalarSingletZ2DMMhInput)"
		@echo "Note: to regenerate ScalarSingletZ2DMMhInput source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_ScalarSingletZ2DMMhInput)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_ScalarSingletZ2DMMhInput):
		@true
endif

$(LIBScalarSingletZ2DMMhInput_DEP) $(EXEScalarSingletZ2DMMhInput_DEP) $(LLScalarSingletZ2DMMhInput_DEP) $(LIBScalarSingletZ2DMMhInput_OBJ) $(EXEScalarSingletZ2DMMhInput_OBJ) $(LLScalarSingletZ2DMMhInput_OBJ) $(LLScalarSingletZ2DMMhInput_LIB): \
	CPPFLAGS += $(MODScalarSingletZ2DMMhInput_SUBMOD_INC) $(MODScalarSingletZ2DMMhInput_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)  $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBScalarSingletZ2DMMhInput_DEP) $(EXEScalarSingletZ2DMMhInput_DEP) $(LLScalarSingletZ2DMMhInput_DEP) $(LIBScalarSingletZ2DMMhInput_OBJ) $(EXEScalarSingletZ2DMMhInput_OBJ) $(LLScalarSingletZ2DMMhInput_OBJ) $(LLScalarSingletZ2DMMhInput_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLScalarSingletZ2DMMhInput_OBJ) $(LLScalarSingletZ2DMMhInput_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBScalarSingletZ2DMMhInput): $(LIBScalarSingletZ2DMMhInput_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBScalarSingletZ2DMMhInput) $(MODScalarSingletZ2DMMhInput_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLScalarSingletZ2DMMhInput_LIB): $(LLScalarSingletZ2DMMhInput_OBJ) $(LIBScalarSingletZ2DMMhInput) $(MODScalarSingletZ2DMMhInput_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS)

ALLDEP += $(LIBScalarSingletZ2DMMhInput_DEP) $(EXEScalarSingletZ2DMMhInput_DEP)
ALLSRC += $(LIBScalarSingletZ2DMMhInput_SRC) $(EXEScalarSingletZ2DMMhInput_SRC)
ALLLIB += $(LIBScalarSingletZ2DMMhInput)
ALLEXE += $(EXEScalarSingletZ2DMMhInput_EXE)
ALLMODDEP += $(MODScalarSingletZ2DMMhInput_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLScalarSingletZ2DMMhInput_DEP)
ALLSRC += $(LLScalarSingletZ2DMMhInput_SRC)
ALLLL  += $(LLScalarSingletZ2DMMhInput_LIB)
endif
