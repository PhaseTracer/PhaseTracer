DIR          := models/ScalarSingletZ2DM
MODNAME      := ScalarSingletZ2DM
SARAH_MODEL  := ScalarSingletZ2DM
WITH_$(MODNAME) := yes
MODScalarSingletZ2DM_MOD := SM
MODScalarSingletZ2DM_DEP := $(patsubst %,model_specific/%,$(MODScalarSingletZ2DM_MOD))
MODScalarSingletZ2DM_INC := $(patsubst %,-Imodel_specific/%,$(MODScalarSingletZ2DM_MOD))
MODScalarSingletZ2DM_LIB := $(foreach M,$(MODScalarSingletZ2DM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODScalarSingletZ2DM_SUBMOD  := $(DIR)/cxx_qft
MODScalarSingletZ2DM_SUBMOD_INC := $(patsubst %,-I%,$(MODScalarSingletZ2DM_SUBMOD))

ScalarSingletZ2DM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
ScalarSingletZ2DM_INSTALL_CXXQFT_DIR := \
		$(ScalarSingletZ2DM_INSTALL_DIR)/cxx_qft

ScalarSingletZ2DM_MK     := \
		$(DIR)/module.mk

ScalarSingletZ2DM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

ScalarSingletZ2DM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

ScalarSingletZ2DM_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

ScalarSingletZ2DM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

ScalarSingletZ2DM_INCLUDE_MK := \
		$(ScalarSingletZ2DM_SUSY_BETAS_MK) \
		$(ScalarSingletZ2DM_SOFT_BETAS_MK) \
		$(ScalarSingletZ2DM_CXX_QFT_VERTICES_MK)

ScalarSingletZ2DM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.ScalarSingletZ2DM_generated \
		$(DIR)/LesHouches.in.ScalarSingletZ2DM

ScalarSingletZ2DM_REFERENCES := \
		$(DIR)/ScalarSingletZ2DM_references.tex

ScalarSingletZ2DM_GNUPLOT := \
		$(DIR)/ScalarSingletZ2DM_plot_rgflow.gnuplot \
		$(DIR)/ScalarSingletZ2DM_plot_spectrum.gnuplot

ScalarSingletZ2DM_TARBALL := \
		$(MODNAME).tar.gz

LIBScalarSingletZ2DM_SRC := \
		$(DIR)/ScalarSingletZ2DM_a_muon.cpp \
		$(DIR)/ScalarSingletZ2DM_edm.cpp \
		$(DIR)/ScalarSingletZ2DM_FFV_form_factors.cpp \
		$(DIR)/ScalarSingletZ2DM_l_to_lgamma.cpp \
		$(DIR)/ScalarSingletZ2DM_effective_couplings.cpp \
		$(DIR)/ScalarSingletZ2DM_info.cpp \
		$(DIR)/ScalarSingletZ2DM_input_parameters.cpp \
		$(DIR)/ScalarSingletZ2DM_mass_eigenstates.cpp \
		$(DIR)/ScalarSingletZ2DM_observables.cpp \
		$(DIR)/ScalarSingletZ2DM_physical.cpp \
		$(DIR)/ScalarSingletZ2DM_slha_io.cpp \
		$(DIR)/ScalarSingletZ2DM_soft_parameters.cpp \
		$(DIR)/ScalarSingletZ2DM_susy_parameters.cpp \
		$(DIR)/ScalarSingletZ2DM_utilities.cpp \
		$(DIR)/ScalarSingletZ2DM_weinberg_angle.cpp

EXEScalarSingletZ2DM_SRC := \
		$(DIR)/run_ScalarSingletZ2DM.cpp \
		$(DIR)/run_cmd_line_ScalarSingletZ2DM.cpp \
		$(DIR)/scan_ScalarSingletZ2DM.cpp
LLScalarSingletZ2DM_LIB  :=
LLScalarSingletZ2DM_OBJ  :=
LLScalarSingletZ2DM_SRC  := \
		$(DIR)/ScalarSingletZ2DM_librarylink.cpp

LLScalarSingletZ2DM_MMA  := \
		$(DIR)/ScalarSingletZ2DM_librarylink.m \
		$(DIR)/run_ScalarSingletZ2DM.m

LIBScalarSingletZ2DM_HDR := \
		$(DIR)/ScalarSingletZ2DM_a_muon.hpp \
		$(DIR)/ScalarSingletZ2DM_convergence_tester.hpp \
		$(DIR)/ScalarSingletZ2DM_edm.hpp \
		$(DIR)/ScalarSingletZ2DM_FFV_form_factors.hpp \
		$(DIR)/ScalarSingletZ2DM_l_to_lgamma.hpp \
		$(DIR)/ScalarSingletZ2DM_effective_couplings.hpp \
		$(DIR)/ScalarSingletZ2DM_ewsb_solver.hpp \
		$(DIR)/ScalarSingletZ2DM_ewsb_solver_interface.hpp \
		$(DIR)/ScalarSingletZ2DM_high_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DM_info.hpp \
		$(DIR)/ScalarSingletZ2DM_initial_guesser.hpp \
		$(DIR)/ScalarSingletZ2DM_input_parameters.hpp \
		$(DIR)/ScalarSingletZ2DM_low_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DM_mass_eigenstates.hpp \
		$(DIR)/ScalarSingletZ2DM_model.hpp \
		$(DIR)/ScalarSingletZ2DM_model_slha.hpp \
		$(DIR)/ScalarSingletZ2DM_observables.hpp \
		$(DIR)/ScalarSingletZ2DM_physical.hpp \
		$(DIR)/ScalarSingletZ2DM_slha_io.hpp \
		$(DIR)/ScalarSingletZ2DM_spectrum_generator.hpp \
		$(DIR)/ScalarSingletZ2DM_spectrum_generator_interface.hpp \
		$(DIR)/ScalarSingletZ2DM_soft_parameters.hpp \
		$(DIR)/ScalarSingletZ2DM_susy_parameters.hpp \
		$(DIR)/ScalarSingletZ2DM_susy_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DM_utilities.hpp \
		$(DIR)/ScalarSingletZ2DM_weinberg_angle.hpp

LIBScalarSingletZ2DM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/ScalarSingletZ2DM_qft.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DM_fields.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DM_vertices.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DM_context_base.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DM_npointfunctions.hpp

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
-include $(ScalarSingletZ2DM_SUSY_BETAS_MK)
-include $(ScalarSingletZ2DM_SOFT_BETAS_MK)
-include $(ScalarSingletZ2DM_CXX_QFT_VERTICES_MK)
-include $(ScalarSingletZ2DM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(ScalarSingletZ2DM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DM_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBScalarSingletZ2DM_SRC := $(sort $(LIBScalarSingletZ2DM_SRC))
EXEScalarSingletZ2DM_SRC := $(sort $(EXEScalarSingletZ2DM_SRC))

LIBScalarSingletZ2DM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBScalarSingletZ2DM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBScalarSingletZ2DM_SRC)))

EXEScalarSingletZ2DM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEScalarSingletZ2DM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEScalarSingletZ2DM_SRC)))

EXEScalarSingletZ2DM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEScalarSingletZ2DM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEScalarSingletZ2DM_SRC)))

LIBScalarSingletZ2DM_DEP := \
		$(LIBScalarSingletZ2DM_OBJ:.o=.d)

EXEScalarSingletZ2DM_DEP := \
		$(EXEScalarSingletZ2DM_OBJ:.o=.d)

LLScalarSingletZ2DM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLScalarSingletZ2DM_SRC)))

LLScalarSingletZ2DM_OBJ  := $(LLScalarSingletZ2DM_SRC:.cpp=.o)
LLScalarSingletZ2DM_LIB  := $(LLScalarSingletZ2DM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBScalarSingletZ2DM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_ScalarSingletZ2DM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_ScalarSingletZ2DM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBScalarSingletZ2DM) $(EXEScalarSingletZ2DM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(ScalarSingletZ2DM_INSTALL_DIR)
		$(Q)install -d $(ScalarSingletZ2DM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DM_SRC) $(ScalarSingletZ2DM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DM_HDR) $(ScalarSingletZ2DM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DM_CXXQFT_HDR) $(ScalarSingletZ2DM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEScalarSingletZ2DM_SRC) $(ScalarSingletZ2DM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DM_SRC) $(ScalarSingletZ2DM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DM_MMA) $(ScalarSingletZ2DM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(ScalarSingletZ2DM_MK) $(ScalarSingletZ2DM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DM_INCLUDE_MK) $(ScalarSingletZ2DM_INSTALL_DIR)
ifneq ($(ScalarSingletZ2DM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DM_SLHA_INPUT) $(ScalarSingletZ2DM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DM_REFERENCES) $(ScalarSingletZ2DM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DM_GNUPLOT) $(ScalarSingletZ2DM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBScalarSingletZ2DM_DEP)
		$(Q)-rm -f $(EXEScalarSingletZ2DM_DEP)
		$(Q)-rm -f $(LLScalarSingletZ2DM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBScalarSingletZ2DM)
		$(Q)-rm -f $(LLScalarSingletZ2DM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBScalarSingletZ2DM_OBJ)
		$(Q)-rm -f $(EXEScalarSingletZ2DM_OBJ)
		$(Q)-rm -f $(LLScalarSingletZ2DM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBScalarSingletZ2DM_SRC)
		$(Q)-rm -f $(LIBScalarSingletZ2DM_HDR)
		$(Q)-rm -f $(LIBScalarSingletZ2DM_CXXQFT_HDR)
		$(Q)-rm -f $(EXEScalarSingletZ2DM_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DM_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_ScalarSingletZ2DM)
		$(Q)-rm -f $(ScalarSingletZ2DM_INCLUDE_MK)
		$(Q)-rm -f $(ScalarSingletZ2DM_SLHA_INPUT)
		$(Q)-rm -f $(ScalarSingletZ2DM_REFERENCES)
		$(Q)-rm -f $(ScalarSingletZ2DM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEScalarSingletZ2DM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(ScalarSingletZ2DM_TARBALL) \
		$(LIBScalarSingletZ2DM_SRC) $(LIBScalarSingletZ2DM_HDR) $(LIBScalarSingletZ2DM_CXXQFT_HDR) \
		$(EXEScalarSingletZ2DM_SRC) \
		$(LLScalarSingletZ2DM_SRC) $(LLScalarSingletZ2DM_MMA) \
		$(ScalarSingletZ2DM_MK) $(ScalarSingletZ2DM_INCLUDE_MK) \
		$(ScalarSingletZ2DM_SLHA_INPUT) $(ScalarSingletZ2DM_REFERENCES) \
		$(ScalarSingletZ2DM_GNUPLOT)

$(LIBScalarSingletZ2DM_SRC) $(LIBScalarSingletZ2DM_HDR) $(LIBScalarSingletZ2DM_CXXQFT_HDR) $(EXEScalarSingletZ2DM_SRC) $(LLScalarSingletZ2DM_SRC) $(LLScalarSingletZ2DM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_ScalarSingletZ2DM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_ScalarSingletZ2DM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_ScalarSingletZ2DM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_ScalarSingletZ2DM)"
		@echo "Note: to regenerate ScalarSingletZ2DM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_ScalarSingletZ2DM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_ScalarSingletZ2DM):
		@true
endif

$(LIBScalarSingletZ2DM_DEP) $(EXEScalarSingletZ2DM_DEP) $(LLScalarSingletZ2DM_DEP) $(LIBScalarSingletZ2DM_OBJ) $(EXEScalarSingletZ2DM_OBJ) $(LLScalarSingletZ2DM_OBJ) $(LLScalarSingletZ2DM_LIB): \
	CPPFLAGS += $(MODScalarSingletZ2DM_SUBMOD_INC) $(MODScalarSingletZ2DM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBScalarSingletZ2DM_DEP) $(EXEScalarSingletZ2DM_DEP) $(LLScalarSingletZ2DM_DEP) $(LIBScalarSingletZ2DM_OBJ) $(EXEScalarSingletZ2DM_OBJ) $(LLScalarSingletZ2DM_OBJ) $(LLScalarSingletZ2DM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLScalarSingletZ2DM_OBJ) $(LLScalarSingletZ2DM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBScalarSingletZ2DM): $(LIBScalarSingletZ2DM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBScalarSingletZ2DM) $(MODScalarSingletZ2DM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLScalarSingletZ2DM_LIB): $(LLScalarSingletZ2DM_OBJ) $(LIBScalarSingletZ2DM) $(MODScalarSingletZ2DM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBScalarSingletZ2DM_DEP) $(EXEScalarSingletZ2DM_DEP)
ALLSRC += $(LIBScalarSingletZ2DM_SRC) $(EXEScalarSingletZ2DM_SRC)
ALLLIB += $(LIBScalarSingletZ2DM)
ALLEXE += $(EXEScalarSingletZ2DM_EXE)
ALLMODDEP += $(MODScalarSingletZ2DM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLScalarSingletZ2DM_DEP)
ALLSRC += $(LLScalarSingletZ2DM_SRC)
ALLLL  += $(LLScalarSingletZ2DM_LIB)
endif
