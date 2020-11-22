FSModelName = "ScalarSingletZ2DMMhInput";
FSEigenstates = SARAH`EWSB;
FSDefaultSARAHModel = ScalarSingletZ2DM;
FlexibleEFTHiggs = False;

(* input parameters *)

MINPAR = {
    {1, MhInput},
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
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS2, muS2Input}
};

SUSYScale = QEWSB;

SUSYScaleFirstGuess = QEWSB;

EWSBOutputParameters = {muH2};

SUSYScaleInput = {
  FSFindRoot[ { LamH }, { MhInput - Pole[M[hh]] } ]
};

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
   {Ye, Automatic},
   { LamH, MhInput^2 / LowEnergyConstant[vev]^2 }
(* If you do not know analytic expression can use: 
 FSFindRoot[ { LamH }, { MhInput - M[hh] } ] *)  
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
FSLoopLibraries = { FSSOFTSUSY };
FSFeynArtsAvailable = False;
FSFormCalcAvailable = False;
