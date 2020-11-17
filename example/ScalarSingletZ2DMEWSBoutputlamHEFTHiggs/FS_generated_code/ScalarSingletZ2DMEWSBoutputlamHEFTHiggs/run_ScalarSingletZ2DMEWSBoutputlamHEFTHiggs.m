Get["models/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs/ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_librarylink.m"];

handle = FSScalarSingletZ2DMEWSBoutputlamHEFTHiggsOpenHandle[
    fsSettings -> {
        precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
        maxIterations -> 0,                (* FlexibleSUSY[1] *)
        solver -> 1,     (* FlexibleSUSY[2] *)
        calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
        poleMassLoopOrder -> 4,            (* FlexibleSUSY[4] *)
        ewsbLoopOrder -> 4,                (* FlexibleSUSY[5] *)
        betaFunctionLoopOrder -> 4,        (* FlexibleSUSY[6] *)
        thresholdCorrectionsLoopOrder -> 4,(* FlexibleSUSY[7] *)
        higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
        higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
        higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
        higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
        forceOutput -> 0,                  (* FlexibleSUSY[12] *)
        topPoleQCDCorrections -> 3,        (* FlexibleSUSY[13] *)
        betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
        forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
        poleMassScale -> 0,                (* FlexibleSUSY[17] *)
        eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
        eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
        eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
        eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
        eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
        calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
        thresholdCorrections -> 124111421, (* FlexibleSUSY[24] *)
        higgs3loopCorrectionRenScheme -> 0,(* FlexibleSUSY[25] *)
        higgs3loopCorrectionAtAsAs -> 1,   (* FlexibleSUSY[26] *)
        higgs3loopCorrectionAbAsAs -> 1,   (* FlexibleSUSY[27] *)
        higgs3loopCorrectionAtAtAs -> 1,   (* FlexibleSUSY[28] *)
        higgs3loopCorrectionAtAtAt -> 1,   (* FlexibleSUSY[29] *)
        higgs4loopCorrectionAtAsAsAs -> 1, (* FlexibleSUSY[30] *)
        parameterOutputScale -> 0          (* MODSEL[12] *)
    },
    fsSMParameters -> {
        alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
        GF -> 1.166378700*^-5,  (* SMINPUTS[2] *)
        alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
        MZ -> 91.1876,          (* SMINPUTS[4] *)
        mbmb -> 4.18,           (* SMINPUTS[5] *)
        Mt -> 173.34,           (* SMINPUTS[6] *)
        Mtau -> 1.77699,        (* SMINPUTS[7] *)
        Mv3 -> 0,               (* SMINPUTS[8] *)
        MW -> 80.385,           (* SMINPUTS[9] *)
        Me -> 0.000510998902,   (* SMINPUTS[11] *)
        Mv1 -> 0,               (* SMINPUTS[12] *)
        Mm -> 0.1056583715,     (* SMINPUTS[13] *)
        Mv2 -> 0,               (* SMINPUTS[14] *)
        md2GeV -> 0.00475,      (* SMINPUTS[21] *)
        mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
        ms2GeV -> 0.104,        (* SMINPUTS[23] *)
        mcmc -> 1.27,           (* SMINPUTS[24] *)
        CKMTheta12 -> 0,
        CKMTheta13 -> 0,
        CKMTheta23 -> 0,
        CKMDelta -> 0,
        PMNSTheta12 -> 0,
        PMNSTheta13 -> 0,
        PMNSTheta23 -> 0,
        PMNSDelta -> 0,
        PMNSAlpha1 -> 0,
        PMNSAlpha2 -> 0,
        alphaEm0 -> 1/137.035999074,
        Mh -> 125.09
    },
    fsModelParameters -> {
        muH2Input -> 0,
        LamSHInput -> 0,
        LamSInput -> 0,
        muS2Input -> 0,
        QEWSB -> 0,
        Qin -> 0
    }
];

spectrum    = FSScalarSingletZ2DMEWSBoutputlamHEFTHiggsCalculateSpectrum[handle];
observables = FSScalarSingletZ2DMEWSBoutputlamHEFTHiggsCalculateObservables[handle];
FSScalarSingletZ2DMEWSBoutputlamHEFTHiggsCloseHandle[handle];

Print[spectrum];
