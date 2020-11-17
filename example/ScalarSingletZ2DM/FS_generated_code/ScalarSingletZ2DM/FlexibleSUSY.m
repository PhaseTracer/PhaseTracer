FSModelName = "ScalarSingletZ2DM";
FSEigenstates = SARAH`EWSB;
FSDefaultSARAHModel = ScalarSingletZ2DM;
FlexibleEFTHiggs = False;

(* input parameters *)

MINPAR = {
    {1, LamHInput},
    {2, LamSHInput},
    {3, LamSInput},
    {4, muS2Input}
};

EXTPAR = {
    {0, QEWSB},
    {1, Qin}
};    

HighScaleFirstGuess=Qin;

HighScale=Qin;

HighScaleInput = {
   {LamH, LamHInput},
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS2, muS2Input}
};

SUSYScale = QEWSB;

SUSYScaleFirstGuess = QEWSB;

EWSBOutputParameters = {muH2};

SUSYScaleInput = {};

InitialGuessAtSUSYScale = {};

LowScale = LowEnergyConstant[MZ];

LowScaleFirstGuess = LowEnergyConstant[MZ];

LowScaleInput = {
   {v, 2 MZMSbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2]},
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic}
};

InitialGuessAtLowScale = {
   {v, LowEnergyConstant[vev]},
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic}
};

EffectiveMu = muH;

SMParticles = {
    Electron, TopQuark, BottomQuark,
    VectorP, VectorZ, VectorG, VectorW, Neutrino,
    Hp, Ah (* goldstones *)
};

OnlyLowEnergyFlexibleSUSY = False;

DefaultPoleMassPrecision = HighPrecision;
HighPoleMassPrecision    = {};
MediumPoleMassPrecision  = {};
LowPoleMassPrecision     = {};
