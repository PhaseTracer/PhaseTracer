
FSModelName = "THDMIISNMSSMBCsimple";
FSEigenstates = SARAH`EWSB;
FSDefaultSARAHModel = "THDMS";

OnlyLowEnergyFlexibleSUSY = False;

(* input parameters *)

MINPAR = {
    {1, lambdaNMSSM},
    {2, kappaNMSSM},
    {3, AlambdaNMSSM},
    {4, AkappaNMSSM},
    {5, AtopNMSSM},
    {6, mstopL},
    {7, mstopR},
    {8, MEWSB},
    {10, TanBeta},  
    {13, vSIN}
};

EXTPAR = {
};

HighScale = Sqrt[mstopL*mstopR];

HighScaleFirstGuess = Sqrt[mstopL*mstopR];

deltaLambda2 = ( 3*Yu[3,3]^2*AtopNMSSM^2 / (8*Pi^2*mstopL*mstopR) ) * ( 1 - AtopNMSSM^2 / (12*mstopL*mstopR)  );  

(* BCs from 1407.4134 -- note: lambdas are not all the same *)
(* Lambda3 = Lambda3_1407.4134 + Lambda4_1407.4134 *)
(* Lambda4 = Lambda4_1407.4134 *)
HighScaleInput = {
    {Lambda1, 0.25*(0.6*g1^2 +g2^2)},
    {Lambda2, 0.25*(0.6*g1^2 +g2^2) + deltaLambda2},
    {Lambda3, 0.25*(g2^2 - 0.6*g1^2) + lambdaNMSSM^2 - 0.5* g2^2},
    {Lambda4, lambdaNMSSM^2 - 0.5* g2^2},
    {Lambda5, lambdaNMSSM^2},
    {Lambda6, lambdaNMSSM^2},
    {Lambda7, -lambdaNMSSM*kappaNMSSM},
    {Lambda8, kappaNMSSM^2},
    {vS, vSIN},
    {M123, lambdaNMSSM * AlambdaNMSSM},
    {M5, -kappaNMSSM*AkappaNMSSM}
};

EWSBOutputParameters = { M112, M222, M332};

SUSYScale = MEWSB;

SUSYScaleFirstGuess = MEWSB;

SUSYScaleInput = {
};

LowScale = LowEnergyConstant[MZ];

LowScaleFirstGuess = LowEnergyConstant[MZ];

LowScaleInput = {
   {v1, 2 MZMSbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Cos[ArcTan[TanBeta]]},
   {v2, 2 MZMSbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Sin[ArcTan[TanBeta]]},
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic}
};

InitialGuessAtLowScale = {
   {v1, LowEnergyConstant[vev] Cos[ArcTan[TanBeta]]},
   {v2, LowEnergyConstant[vev] Sin[ArcTan[TanBeta]]},
   {vS, vSIN},
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic}
};

DefaultPoleMassPrecision = MediumPrecision;
HighPoleMassPrecision    = {hh, Ah, Hpm};
MediumPoleMassPrecision  = {};
LowPoleMassPrecision     = {};

ExtraSLHAOutputBlocks = {
   {EFFHIGGSCOUPLINGS, NoScale,
           {{1, FlexibleSUSYObservable`CpHiggsPhotonPhoton},
            {2, FlexibleSUSYObservable`CpHiggsGluonGluon},
            {3, FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton},
            {4, FlexibleSUSYObservable`CpPseudoScalarGluonGluon} } }
};
