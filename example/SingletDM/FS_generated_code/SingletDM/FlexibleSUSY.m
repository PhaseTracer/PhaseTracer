FSModelName = "SingletDM";
FSEigenstates = SARAH`EWSB;
FSDefaultSARAHModel = SingletDM;
FlexibleEFTHiggs = True;

(* input parameters *)

MINPAR = {
    {1, HiggsIN},
    {2, LamSHInput},
    {3, LamSInput},
    {4, muSInput}
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
   {muH, HiggsIN},
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS, muSInput}
};

InitialGuessAtSUSYScale = {
   {v, LowEnergyConstant[vev]},
   {LamSH, LamSHInput},
   {LamS, LamSInput},
   {muS, muSInput},
   {muH, HiggsIN}
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
