FSModelName = "ScalarSingletZ2DM";
FSEigenstates = SARAH`EWSB;
FSDefaultSARAHModel = ScalarSingletZ2DM;
FlexibleEFTHiggs = True;

(* input parameters *)

MINPAR = {
    {1, muH2Input},
    {2, LamSHInput},
    {3, LamSInput},
    {4, muS2Input}
};

EXTPAR = {
    {0, QEWSB},
    {1, Qin}
};    

EWSBOutputParameters = {LamH};

SUSYScale = QEWSB;

SUSYScaleFirstGuess = QEWSB;


LowScale = LowEnergyConstant[MZ];

LowScaleFirstGuess = LowEnergyConstant[MZ];


SUSYScaleMatching = {
    {v, VEV}
    };


SUSYScaleInput = {
   {muH2, muH2Input},
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS2, muS2Input}
};

InitialGuessAtSUSYScale = {
   {v, LowEnergyConstant[vev]},
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS2, muS2Input},
   {muH2, muH2Input}
};

EffectiveMu = muH;



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
