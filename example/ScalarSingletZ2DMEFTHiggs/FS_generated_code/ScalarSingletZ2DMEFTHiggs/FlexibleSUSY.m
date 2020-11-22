FSModelName = "ScalarSingletZ2DMEFTHiggs";
FSEigenstates = SARAH`EWSB;
FSDefaultSARAHModel = ScalarSingletZ2DM;
FlexibleEFTHiggs = True;

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

EWSBOutputParameters = {muH2};

SUSYScale = QEWSB;

SUSYScaleFirstGuess = QEWSB;

LowScale = LowEnergyConstant[MZ];

LowScaleFirstGuess = LowEnergyConstant[MZ];

MatchingScaleInput = {
    {v, VEV}
    };

SUSYScaleInput = {
   {LamH, LamHInput},
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS2, muS2Input}
};

InitialGuessAtSUSYScale = {
   {v, LowEnergyConstant[vev]},
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS2, muS2Input},
   {LamH, LamHInput}
};

EffectiveMu = muH2;

HighScaleFirstGuess=Qin;

HighScale=Qin;

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
FSLoopLibraries = { FSSOFTSUSY };
FSFeynArtsAvailable = False;
FSFormCalcAvailable = False;
