DIR          := models/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs
MODNAME      := ScalarSingletZ2DMEWSBoutputlamHEFTHiggs
SARAH_MODEL  := ScalarSingletZ2DM
WITH_$(MODNAME) := yes
MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MOD := SM
MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP := $(patsubst %,model_specific/%,$(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MOD))
MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INC := $(patsubst %,-Imodel_specific/%,$(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MOD))
MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB := $(foreach M,$(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUBMOD  := $(DIR)/cxx_qft
MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUBMOD_INC := $(patsubst %,-I%,$(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUBMOD))

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_CXXQFT_DIR := \
		$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)/cxx_qft

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MK     := \
		$(DIR)/module.mk

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INCLUDE_MK := \
		$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUSY_BETAS_MK) \
		$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SOFT_BETAS_MK) \
		$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXX_QFT_VERTICES_MK)

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SLHA_INPUT := \
		$(DIR)/LesHouches.in.ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_generated \
		$(DIR)/LesHouches.in.ScalarSingletZ2DMEWSBoutputlamHEFTHiggs

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_REFERENCES := \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_references.tex

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_GNUPLOT := \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_plot_rgflow.gnuplot \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_plot_spectrum.gnuplot

ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_TARBALL := \
		$(MODNAME).tar.gz

LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC := \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_a_muon.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_edm.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_FFV_form_factors.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_l_to_lgamma.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_effective_couplings.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_info.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_mass_eigenstates.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_observables.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_physical.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_slha_io.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_utilities.cpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_weinberg_angle.cpp

EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC := \
		$(DIR)/run_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs.cpp \
		$(DIR)/run_cmd_line_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs.cpp \
		$(DIR)/scan_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs.cpp
LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB  :=
LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ  :=
LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC  := \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_librarylink.cpp

LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MMA  := \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_librarylink.m \
		$(DIR)/run_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs.m

LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_HDR := \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_a_muon.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_convergence_tester.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_edm.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_FFV_form_factors.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_l_to_lgamma.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_effective_couplings.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_ewsb_solver.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_ewsb_solver_interface.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_high_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_info.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_initial_guesser.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_input_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_low_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_mass_eigenstates.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_model.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_model_slha.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_observables.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_physical.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_slha_io.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_spectrum_generator.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_spectrum_generator_interface.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_soft_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_parameters.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_susy_scale_constraint.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_utilities.hpp \
		$(DIR)/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_weinberg_angle.hpp

LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXXQFT_HDR := \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_qft.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_fields.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_vertices.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_context_base.hpp \
		$(DIR)/cxx_qft/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_npointfunctions.hpp

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
-include $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUSY_BETAS_MK)
-include $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SOFT_BETAS_MK)
-include $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXX_QFT_VERTICES_MK)
-include $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC := $(sort $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC))
EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC := $(sort $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC))

LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)))

EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)))

EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)))

LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP := \
		$(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ:.o=.d)

EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP := \
		$(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ:.o=.d)

LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)))

LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ  := $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC:.cpp=.o)
LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB  := $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
		$(Q)install -d $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_HDR) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXXQFT_HDR) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MMA) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MK) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INCLUDE_MK) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
ifneq ($(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SLHA_INPUT) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_REFERENCES) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_GNUPLOT) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP)
		$(Q)-rm -f $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP)
		$(Q)-rm -f $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs)
		$(Q)-rm -f $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ)
		$(Q)-rm -f $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ)
		$(Q)-rm -f $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)
		$(Q)-rm -f $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_HDR)
		$(Q)-rm -f $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXXQFT_HDR)
		$(Q)-rm -f $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)
		$(Q)-rm -f $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MMA)
		$(Q)-rm -f $(METACODE_STAMP_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs)
		$(Q)-rm -f $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INCLUDE_MK)
		$(Q)-rm -f $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SLHA_INPUT)
		$(Q)-rm -f $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_REFERENCES)
		$(Q)-rm -f $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_TARBALL) \
		$(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_HDR) $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXXQFT_HDR) \
		$(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) \
		$(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MMA) \
		$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MK) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INCLUDE_MK) \
		$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SLHA_INPUT) $(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_REFERENCES) \
		$(ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_GNUPLOT)

$(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_HDR) $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_CXXQFT_HDR) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs)"
		@echo "Note: to regenerate ScalarSingletZ2DMEWSBoutputlamHEFTHiggs source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs):
		@true
endif

$(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP) $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB): \
	CPPFLAGS += $(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SUBMOD_INC) $(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP) $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs): $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs) $(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB): $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_OBJ) $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs) $(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP)
ALLSRC += $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC) $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)
ALLLIB += $(LIBScalarSingletZ2DMEWSBoutputlamHEFTHiggs)
ALLEXE += $(EXEScalarSingletZ2DMEWSBoutputlamHEFTHiggs_EXE)
ALLMODDEP += $(MODScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_DEP)
ALLSRC += $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_SRC)
ALLLL  += $(LLScalarSingletZ2DMEWSBoutputlamHEFTHiggs_LIB)
endif
