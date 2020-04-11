DIR          := models/THDMIISNMSSMBCsimple
MODNAME      := THDMIISNMSSMBCsimple
SARAH_MODEL  := THDMS
WITH_$(MODNAME) := yes
MODTHDMIISNMSSMBCsimple_MOD := SM
MODTHDMIISNMSSMBCsimple_DEP := $(patsubst %,model_specific/%,$(MODTHDMIISNMSSMBCsimple_MOD))
MODTHDMIISNMSSMBCsimple_INC := $(patsubst %,-Imodel_specific/%,$(MODTHDMIISNMSSMBCsimple_MOD))
MODTHDMIISNMSSMBCsimple_LIB := $(foreach M,$(MODTHDMIISNMSSMBCsimple_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODTHDMIISNMSSMBCsimple_SUBMOD  := $(DIR)/cxx_qft
MODTHDMIISNMSSMBCsimple_SUBMOD_INC := $(patsubst %,-I%,$(MODTHDMIISNMSSMBCsimple_SUBMOD))

THDMIISNMSSMBCsimple_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
THDMIISNMSSMBCsimple_INSTALL_CXXQFT_DIR := \
		$(THDMIISNMSSMBCsimple_INSTALL_DIR)/cxx_qft

THDMIISNMSSMBCsimple_MK     := \
		$(DIR)/module.mk

THDMIISNMSSMBCsimple_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

THDMIISNMSSMBCsimple_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

THDMIISNMSSMBCsimple_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

THDMIISNMSSMBCsimple_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

THDMIISNMSSMBCsimple_INCLUDE_MK := \
		$(THDMIISNMSSMBCsimple_SUSY_BETAS_MK) \
		$(THDMIISNMSSMBCsimple_SOFT_BETAS_MK) \
		$(THDMIISNMSSMBCsimple_CXX_QFT_VERTICES_MK)

THDMIISNMSSMBCsimple_SLHA_INPUT := \
		$(DIR)/LesHouches.in.THDMIISNMSSMBCsimple_generated \
		$(DIR)/LesHouches.in.THDMII \
		$(DIR)/LesHouches.in.THDMIISNMSSMBCsimple.BM1

THDMIISNMSSMBCsimple_REFERENCES := \
		$(DIR)/THDMIISNMSSMBCsimple_references.tex

THDMIISNMSSMBCsimple_GNUPLOT := \
		$(DIR)/THDMIISNMSSMBCsimple_plot_rgflow.gnuplot \
		$(DIR)/THDMIISNMSSMBCsimple_plot_spectrum.gnuplot

THDMIISNMSSMBCsimple_TARBALL := \
		$(MODNAME).tar.gz

LIBTHDMIISNMSSMBCsimple_SRC := \
		$(DIR)/THDMIISNMSSMBCsimple_a_muon.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_edm.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_FFV_form_factors.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_l_to_lgamma.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_effective_couplings.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_info.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_input_parameters.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_mass_eigenstates.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_observables.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_physical.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_slha_io.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_soft_parameters.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_susy_parameters.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_utilities.cpp \
		$(DIR)/THDMIISNMSSMBCsimple_weinberg_angle.cpp

EXETHDMIISNMSSMBCsimple_SRC := \
		$(DIR)/run_THDMIISNMSSMBCsimple.cpp \
		$(DIR)/run_cmd_line_THDMIISNMSSMBCsimple.cpp \
		$(DIR)/scan_THDMIISNMSSMBCsimple.cpp
LLTHDMIISNMSSMBCsimple_LIB  :=
LLTHDMIISNMSSMBCsimple_OBJ  :=
LLTHDMIISNMSSMBCsimple_SRC  := \
		$(DIR)/THDMIISNMSSMBCsimple_librarylink.cpp

LLTHDMIISNMSSMBCsimple_MMA  := \
		$(DIR)/THDMIISNMSSMBCsimple_librarylink.m \
		$(DIR)/run_THDMIISNMSSMBCsimple.m

LIBTHDMIISNMSSMBCsimple_HDR := \
		$(DIR)/THDMIISNMSSMBCsimple_a_muon.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_convergence_tester.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_edm.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_FFV_form_factors.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_l_to_lgamma.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_effective_couplings.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_ewsb_solver.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_ewsb_solver_interface.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_high_scale_constraint.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_info.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_initial_guesser.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_input_parameters.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_low_scale_constraint.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_mass_eigenstates.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_model.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_model_slha.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_observables.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_physical.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_slha_io.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_spectrum_generator.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_spectrum_generator_interface.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_soft_parameters.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_susy_parameters.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_susy_scale_constraint.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_utilities.hpp \
		$(DIR)/THDMIISNMSSMBCsimple_weinberg_angle.hpp

LIBTHDMIISNMSSMBCsimple_CXXQFT_HDR := \
		$(DIR)/cxx_qft/THDMIISNMSSMBCsimple_qft.hpp \
		$(DIR)/cxx_qft/THDMIISNMSSMBCsimple_fields.hpp \
		$(DIR)/cxx_qft/THDMIISNMSSMBCsimple_vertices.hpp \
		$(DIR)/cxx_qft/THDMIISNMSSMBCsimple_context_base.hpp \
		$(DIR)/cxx_qft/THDMIISNMSSMBCsimple_npointfunctions.hpp

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
-include $(THDMIISNMSSMBCsimple_SUSY_BETAS_MK)
-include $(THDMIISNMSSMBCsimple_SOFT_BETAS_MK)
-include $(THDMIISNMSSMBCsimple_CXX_QFT_VERTICES_MK)
-include $(THDMIISNMSSMBCsimple_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(THDMIISNMSSMBCsimple_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMIISNMSSMBCsimple_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMIISNMSSMBCsimple_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMIISNMSSMBCsimple_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBTHDMIISNMSSMBCsimple_SRC := $(sort $(LIBTHDMIISNMSSMBCsimple_SRC))
EXETHDMIISNMSSMBCsimple_SRC := $(sort $(EXETHDMIISNMSSMBCsimple_SRC))

LIBTHDMIISNMSSMBCsimple_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTHDMIISNMSSMBCsimple_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTHDMIISNMSSMBCsimple_SRC)))

EXETHDMIISNMSSMBCsimple_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXETHDMIISNMSSMBCsimple_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXETHDMIISNMSSMBCsimple_SRC)))

EXETHDMIISNMSSMBCsimple_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXETHDMIISNMSSMBCsimple_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXETHDMIISNMSSMBCsimple_SRC)))

LIBTHDMIISNMSSMBCsimple_DEP := \
		$(LIBTHDMIISNMSSMBCsimple_OBJ:.o=.d)

EXETHDMIISNMSSMBCsimple_DEP := \
		$(EXETHDMIISNMSSMBCsimple_OBJ:.o=.d)

LLTHDMIISNMSSMBCsimple_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLTHDMIISNMSSMBCsimple_SRC)))

LLTHDMIISNMSSMBCsimple_OBJ  := $(LLTHDMIISNMSSMBCsimple_SRC:.cpp=.o)
LLTHDMIISNMSSMBCsimple_LIB  := $(LLTHDMIISNMSSMBCsimple_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBTHDMIISNMSSMBCsimple     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_THDMIISNMSSMBCsimple := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_THDMIISNMSSMBCsimple := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBTHDMIISNMSSMBCsimple) $(EXETHDMIISNMSSMBCsimple_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(THDMIISNMSSMBCsimple_INSTALL_DIR)
		$(Q)install -d $(THDMIISNMSSMBCsimple_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBTHDMIISNMSSMBCsimple_SRC) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBTHDMIISNMSSMBCsimple_HDR) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBTHDMIISNMSSMBCsimple_CXXQFT_HDR) $(THDMIISNMSSMBCsimple_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXETHDMIISNMSSMBCsimple_SRC) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLTHDMIISNMSSMBCsimple_SRC) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLTHDMIISNMSSMBCsimple_MMA) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(THDMIISNMSSMBCsimple_MK) $(THDMIISNMSSMBCsimple_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(THDMIISNMSSMBCsimple_INCLUDE_MK) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
ifneq ($(THDMIISNMSSMBCsimple_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(THDMIISNMSSMBCsimple_SLHA_INPUT) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(THDMIISNMSSMBCsimple_REFERENCES) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(THDMIISNMSSMBCsimple_GNUPLOT) $(THDMIISNMSSMBCsimple_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBTHDMIISNMSSMBCsimple_DEP)
		$(Q)-rm -f $(EXETHDMIISNMSSMBCsimple_DEP)
		$(Q)-rm -f $(LLTHDMIISNMSSMBCsimple_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBTHDMIISNMSSMBCsimple)
		$(Q)-rm -f $(LLTHDMIISNMSSMBCsimple_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBTHDMIISNMSSMBCsimple_OBJ)
		$(Q)-rm -f $(EXETHDMIISNMSSMBCsimple_OBJ)
		$(Q)-rm -f $(LLTHDMIISNMSSMBCsimple_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBTHDMIISNMSSMBCsimple_SRC)
		$(Q)-rm -f $(LIBTHDMIISNMSSMBCsimple_HDR)
		$(Q)-rm -f $(LIBTHDMIISNMSSMBCsimple_CXXQFT_HDR)
		$(Q)-rm -f $(EXETHDMIISNMSSMBCsimple_SRC)
		$(Q)-rm -f $(LLTHDMIISNMSSMBCsimple_SRC)
		$(Q)-rm -f $(LLTHDMIISNMSSMBCsimple_MMA)
		$(Q)-rm -f $(METACODE_STAMP_THDMIISNMSSMBCsimple)
		$(Q)-rm -f $(THDMIISNMSSMBCsimple_INCLUDE_MK)
		$(Q)-rm -f $(THDMIISNMSSMBCsimple_SLHA_INPUT)
		$(Q)-rm -f $(THDMIISNMSSMBCsimple_REFERENCES)
		$(Q)-rm -f $(THDMIISNMSSMBCsimple_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXETHDMIISNMSSMBCsimple_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(THDMIISNMSSMBCsimple_TARBALL) \
		$(LIBTHDMIISNMSSMBCsimple_SRC) $(LIBTHDMIISNMSSMBCsimple_HDR) $(LIBTHDMIISNMSSMBCsimple_CXXQFT_HDR) \
		$(EXETHDMIISNMSSMBCsimple_SRC) \
		$(LLTHDMIISNMSSMBCsimple_SRC) $(LLTHDMIISNMSSMBCsimple_MMA) \
		$(THDMIISNMSSMBCsimple_MK) $(THDMIISNMSSMBCsimple_INCLUDE_MK) \
		$(THDMIISNMSSMBCsimple_SLHA_INPUT) $(THDMIISNMSSMBCsimple_REFERENCES) \
		$(THDMIISNMSSMBCsimple_GNUPLOT)

$(LIBTHDMIISNMSSMBCsimple_SRC) $(LIBTHDMIISNMSSMBCsimple_HDR) $(LIBTHDMIISNMSSMBCsimple_CXXQFT_HDR) $(EXETHDMIISNMSSMBCsimple_SRC) $(LLTHDMIISNMSSMBCsimple_SRC) $(LLTHDMIISNMSSMBCsimple_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_THDMIISNMSSMBCsimple)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_THDMIISNMSSMBCsimple): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_THDMIISNMSSMBCsimple)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_THDMIISNMSSMBCsimple)"
		@echo "Note: to regenerate THDMIISNMSSMBCsimple source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_THDMIISNMSSMBCsimple)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_THDMIISNMSSMBCsimple):
		@true
endif

$(LIBTHDMIISNMSSMBCsimple_DEP) $(EXETHDMIISNMSSMBCsimple_DEP) $(LLTHDMIISNMSSMBCsimple_DEP) $(LIBTHDMIISNMSSMBCsimple_OBJ) $(EXETHDMIISNMSSMBCsimple_OBJ) $(LLTHDMIISNMSSMBCsimple_OBJ) $(LLTHDMIISNMSSMBCsimple_LIB): \
	CPPFLAGS += $(MODTHDMIISNMSSMBCsimple_SUBMOD_INC) $(MODTHDMIISNMSSMBCsimple_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBTHDMIISNMSSMBCsimple_DEP) $(EXETHDMIISNMSSMBCsimple_DEP) $(LLTHDMIISNMSSMBCsimple_DEP) $(LIBTHDMIISNMSSMBCsimple_OBJ) $(EXETHDMIISNMSSMBCsimple_OBJ) $(LLTHDMIISNMSSMBCsimple_OBJ) $(LLTHDMIISNMSSMBCsimple_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLTHDMIISNMSSMBCsimple_OBJ) $(LLTHDMIISNMSSMBCsimple_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBTHDMIISNMSSMBCsimple): $(LIBTHDMIISNMSSMBCsimple_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBTHDMIISNMSSMBCsimple) $(MODTHDMIISNMSSMBCsimple_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLTHDMIISNMSSMBCsimple_LIB): $(LLTHDMIISNMSSMBCsimple_OBJ) $(LIBTHDMIISNMSSMBCsimple) $(MODTHDMIISNMSSMBCsimple_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBTHDMIISNMSSMBCsimple_DEP) $(EXETHDMIISNMSSMBCsimple_DEP)
ALLSRC += $(LIBTHDMIISNMSSMBCsimple_SRC) $(EXETHDMIISNMSSMBCsimple_SRC)
ALLLIB += $(LIBTHDMIISNMSSMBCsimple)
ALLEXE += $(EXETHDMIISNMSSMBCsimple_EXE)
ALLMODDEP += $(MODTHDMIISNMSSMBCsimple_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLTHDMIISNMSSMBCsimple_DEP)
ALLSRC += $(LLTHDMIISNMSSMBCsimple_SRC)
ALLLL  += $(LLTHDMIISNMSSMBCsimple_LIB)
endif
