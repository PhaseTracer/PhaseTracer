DIR          := models/ScalarSingletZ2DMMhInputMsInput
MODNAME      := ScalarSingletZ2DMMhInputMsInput
SARAH_MODEL  := ScalarSingletZ2DM
WITH_$(MODNAME) := yes
MODScalarSingletZ2DMMhInputMsInput_MOD := SM
MODScalarSingletZ2DMMhInputMsInput_DEP := $(patsubst %,model_specific/%,$(MODScalarSingletZ2DMMhInputMsInput_MOD))
MODScalarSingletZ2DMMhInputMsInput_INC := $(patsubst %,-Imodel_specific/%,$(MODScalarSingletZ2DMMhInputMsInput_MOD))
MODScalarSingletZ2DMMhInputMsInput_LIB := $(foreach M,$(MODScalarSingletZ2DMMhInputMsInput_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODScalarSingletZ2DMMhInputMsInput_SUBMOD  := $(DIR)/cxx_qft
MODScalarSingletZ2DMMhInputMsInput_SUBMOD_INC := $(patsubst %,-I%,$(MODScalarSingletZ2DMMhInputMsInput_SUBMOD))

ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
ScalarSingletZ2DMMhInputMsInput_INSTALL_CXXQFT_DIR := \
		$(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)/cxx_qft

ScalarSingletZ2DMMhInputMsInput_MK     := \
		$(DIR)/module.mk

ScalarSingletZ2DMMhInputMsInput_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

ScalarSingletZ2DMMhInputMsInput_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

ScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(ScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_MK)
LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_SRC ?= ''

ScalarSingletZ2DMMhInputMsInput_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

ScalarSingletZ2DMMhInputMsInput_INCLUDE_MK := \
		$(ScalarSingletZ2DMMhInputMsInput_SUSY_BETAS_MK) \
		$(ScalarSingletZ2DMMhInputMsInput_SOFT_BETAS_MK)

ScalarSingletZ2DMMhInputMsInput_SLHA_INPUT := \
		$(DIR)/LesHouches.in.ScalarSingletZ2DMMhInputMsInput_generated \


ScalarSingletZ2DMMhInputMsInput_REFERENCES := \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_references.tex

ScalarSingletZ2DMMhInputMsInput_GNUPLOT := \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_plot_rgflow.gnuplot \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_plot_spectrum.gnuplot

ScalarSingletZ2DMMhInputMsInput_TARBALL := \
		$(MODNAME).tar.gz

LIBScalarSingletZ2DMMhInputMsInput_SRC := \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_a_muon.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_edm.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_FFV_form_factors.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_f_to_f_conversion.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_l_to_lgamma.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_b_to_s_gamma.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_effective_couplings.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_info.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_input_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_model_slha.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_observables.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_physical.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_slha_io.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_soft_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_susy_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_utilities.cpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_weinberg_angle.cpp

LIBScalarSingletZ2DMMhInputMsInput_SRC += $(LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_SRC)

EXEScalarSingletZ2DMMhInputMsInput_SRC := \
		$(DIR)/run_ScalarSingletZ2DMMhInputMsInput.cpp \
		$(DIR)/run_cmd_line_ScalarSingletZ2DMMhInputMsInput.cpp \
		$(DIR)/scan_ScalarSingletZ2DMMhInputMsInput.cpp
LLScalarSingletZ2DMMhInputMsInput_LIB  :=
LLScalarSingletZ2DMMhInputMsInput_OBJ  :=
LLScalarSingletZ2DMMhInputMsInput_SRC  := \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_librarylink.cpp

LLScalarSingletZ2DMMhInputMsInput_MMA  := \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_librarylink.m \
		$(DIR)/run_ScalarSingletZ2DMMhInputMsInput.m

LIBScalarSingletZ2DMMhInputMsInput_HDR := \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_a_muon.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_convergence_tester.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_edm.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_FFV_form_factors.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_f_to_f_conversion.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_l_to_lgamma.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_b_to_s_gamma.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_effective_couplings.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_ewsb_solver.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_ewsb_solver_interface.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_high_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_info.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_initial_guesser.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_input_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_low_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_mass_eigenstates.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_interface.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_model.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_model_slha.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_observables.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_physical.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_slha_io.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_spectrum_generator.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_spectrum_generator_interface.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_soft_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_susy_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_susy_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_utilities.hpp \
		$(DIR)/ScalarSingletZ2DMMhInputMsInput_weinberg_angle.hpp

LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_HDR := \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInputMsInput_qft.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInputMsInput_fields.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInputMsInput_vertices.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInputMsInput_context_base.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMMhInputMsInput_npointfunctions_wilsoncoeffs.hpp

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
-include $(ScalarSingletZ2DMMhInputMsInput_SUSY_BETAS_MK)
-include $(ScalarSingletZ2DMMhInputMsInput_SOFT_BETAS_MK)
-include $(ScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_MK)
-include $(ScalarSingletZ2DMMhInputMsInput_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(ScalarSingletZ2DMMhInputMsInput_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMMhInputMsInput_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMMhInputMsInput_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBScalarSingletZ2DMMhInputMsInput_SRC := $(sort $(LIBScalarSingletZ2DMMhInputMsInput_SRC))
EXEScalarSingletZ2DMMhInputMsInput_SRC := $(sort $(EXEScalarSingletZ2DMMhInputMsInput_SRC))

LIBScalarSingletZ2DMMhInputMsInput_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBScalarSingletZ2DMMhInputMsInput_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBScalarSingletZ2DMMhInputMsInput_SRC)))

EXEScalarSingletZ2DMMhInputMsInput_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEScalarSingletZ2DMMhInputMsInput_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEScalarSingletZ2DMMhInputMsInput_SRC)))

EXEScalarSingletZ2DMMhInputMsInput_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEScalarSingletZ2DMMhInputMsInput_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEScalarSingletZ2DMMhInputMsInput_SRC)))

LIBScalarSingletZ2DMMhInputMsInput_DEP := \
		$(LIBScalarSingletZ2DMMhInputMsInput_OBJ:.o=.d)

EXEScalarSingletZ2DMMhInputMsInput_DEP := \
		$(EXEScalarSingletZ2DMMhInputMsInput_OBJ:.o=.d)

LLScalarSingletZ2DMMhInputMsInput_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLScalarSingletZ2DMMhInputMsInput_SRC)))

LLScalarSingletZ2DMMhInputMsInput_OBJ  := $(LLScalarSingletZ2DMMhInputMsInput_SRC:.cpp=.o)
LLScalarSingletZ2DMMhInputMsInput_LIB  := $(LLScalarSingletZ2DMMhInputMsInput_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBScalarSingletZ2DMMhInputMsInput     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_ScalarSingletZ2DMMhInputMsInput := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_ScalarSingletZ2DMMhInputMsInput := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBScalarSingletZ2DMMhInputMsInput) $(EXEScalarSingletZ2DMMhInputMsInput_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)install -d $(ScalarSingletZ2DMMhInputMsInput_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInputMsInput_SRC) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_SRC) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInputMsInput_HDR) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_HDR) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEScalarSingletZ2DMMhInputMsInput_SRC) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMMhInputMsInput_SRC) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMMhInputMsInput_MMA) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(ScalarSingletZ2DMMhInputMsInput_MK) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInputMsInput_INCLUDE_MK) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_MK) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_CXXQFT_DIR)

ifneq ($(ScalarSingletZ2DMMhInputMsInput_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInputMsInput_SLHA_INPUT) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInputMsInput_REFERENCES) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMMhInputMsInput_GNUPLOT) $(ScalarSingletZ2DMMhInputMsInput_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInputMsInput_DEP)
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInputMsInput_DEP)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInputMsInput_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInputMsInput)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInputMsInput_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInputMsInput_OBJ)
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInputMsInput_OBJ)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInputMsInput_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInputMsInput_SRC)
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInputMsInput_HDR)
		$(Q)-rm -f $(LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_HDR)
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInputMsInput_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInputMsInput_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMMhInputMsInput_MMA)
		$(Q)-rm -f $(METACODE_STAMP_ScalarSingletZ2DMMhInputMsInput)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInputMsInput_INCLUDE_MK)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInputMsInput_SLHA_INPUT)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInputMsInput_REFERENCES)
		$(Q)-rm -f $(ScalarSingletZ2DMMhInputMsInput_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEScalarSingletZ2DMMhInputMsInput_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(ScalarSingletZ2DMMhInputMsInput_TARBALL) \
		$(LIBScalarSingletZ2DMMhInputMsInput_SRC) $(LIBScalarSingletZ2DMMhInputMsInput_HDR) $(LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_HDR) \
		$(EXEScalarSingletZ2DMMhInputMsInput_SRC) \
		$(LLScalarSingletZ2DMMhInputMsInput_SRC) $(LLScalarSingletZ2DMMhInputMsInput_MMA) \
		$(ScalarSingletZ2DMMhInputMsInput_MK) $(ScalarSingletZ2DMMhInputMsInput_INCLUDE_MK) $(ScalarSingletZ2DMMhInputMsInput_CXXQFT_VERTICES_MK) \
		$(ScalarSingletZ2DMMhInputMsInput_SLHA_INPUT) $(ScalarSingletZ2DMMhInputMsInput_REFERENCES) \
		$(ScalarSingletZ2DMMhInputMsInput_GNUPLOT)

$(LIBScalarSingletZ2DMMhInputMsInput_SRC) $(LIBScalarSingletZ2DMMhInputMsInput_HDR) $(LIBScalarSingletZ2DMMhInputMsInput_CXXQFT_HDR) $(EXEScalarSingletZ2DMMhInputMsInput_SRC) $(LLScalarSingletZ2DMMhInputMsInput_SRC) $(LLScalarSingletZ2DMMhInputMsInput_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_ScalarSingletZ2DMMhInputMsInput)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_ScalarSingletZ2DMMhInputMsInput): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_ScalarSingletZ2DMMhInputMsInput)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_ScalarSingletZ2DMMhInputMsInput)"
		@echo "Note: to regenerate ScalarSingletZ2DMMhInputMsInput source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_ScalarSingletZ2DMMhInputMsInput)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_ScalarSingletZ2DMMhInputMsInput):
		@true
endif

$(LIBScalarSingletZ2DMMhInputMsInput_DEP) $(EXEScalarSingletZ2DMMhInputMsInput_DEP) $(LLScalarSingletZ2DMMhInputMsInput_DEP) $(LIBScalarSingletZ2DMMhInputMsInput_OBJ) $(EXEScalarSingletZ2DMMhInputMsInput_OBJ) $(LLScalarSingletZ2DMMhInputMsInput_OBJ) $(LLScalarSingletZ2DMMhInputMsInput_LIB): \
	CPPFLAGS += $(MODScalarSingletZ2DMMhInputMsInput_SUBMOD_INC) $(MODScalarSingletZ2DMMhInputMsInput_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBScalarSingletZ2DMMhInputMsInput_DEP) $(EXEScalarSingletZ2DMMhInputMsInput_DEP) $(LLScalarSingletZ2DMMhInputMsInput_DEP) $(LIBScalarSingletZ2DMMhInputMsInput_OBJ) $(EXEScalarSingletZ2DMMhInputMsInput_OBJ) $(LLScalarSingletZ2DMMhInputMsInput_OBJ) $(LLScalarSingletZ2DMMhInputMsInput_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLScalarSingletZ2DMMhInputMsInput_OBJ) $(LLScalarSingletZ2DMMhInputMsInput_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBScalarSingletZ2DMMhInputMsInput): $(LIBScalarSingletZ2DMMhInputMsInput_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBScalarSingletZ2DMMhInputMsInput) $(MODScalarSingletZ2DMMhInputMsInput_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLScalarSingletZ2DMMhInputMsInput_LIB): $(LLScalarSingletZ2DMMhInputMsInput_OBJ) $(LIBScalarSingletZ2DMMhInputMsInput) $(MODScalarSingletZ2DMMhInputMsInput_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBScalarSingletZ2DMMhInputMsInput_DEP) $(EXEScalarSingletZ2DMMhInputMsInput_DEP)
ALLSRC += $(LIBScalarSingletZ2DMMhInputMsInput_SRC) $(EXEScalarSingletZ2DMMhInputMsInput_SRC)
ALLLIB += $(LIBScalarSingletZ2DMMhInputMsInput)
ALLEXE += $(EXEScalarSingletZ2DMMhInputMsInput_EXE)
ALLMODDEP += $(MODScalarSingletZ2DMMhInputMsInput_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLScalarSingletZ2DMMhInputMsInput_DEP)
ALLSRC += $(LLScalarSingletZ2DMMhInputMsInput_SRC)
ALLLL  += $(LLScalarSingletZ2DMMhInputMsInput_LIB)
endif
