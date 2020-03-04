(* ::Package:: *)

Off[General::spell]

Model`Name = "THDMS";
Model`NameLaTeX ="Two Higgs Doublet Model plus Singlet";
Model`Authors = "F.Staub, G.Pozzo, P.Athron";
Model`Date = "2017-04-28";

(* 2017-04-28: Added ScalarField[[3]], S, phiS, sigmaS, modified the lagrangian*)
(* 2014-11-06: Changed sign in Lagrangian *)
(* 2015-11-16: added conj[H1].H2 term *)

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False};
Gauge[[2]]={WB, SU[2], left,        g2,True};
Gauge[[3]]={G,  SU[3], color,       g3,False};


(* Chiral Superfields *)

FermionFields[[1]] = {q, 3, {uL,dL},     1/6, 2, 3};  
FermionFields[[2]] = {l, 3, {vL,eL},    -1/2, 2, 1};
FermionFields[[3]] = {d, 3, conj[dR],   1/3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR],  -2/3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR],     1, 1,  1};

ScalarFields[[1]] =  {H1, 1, {H10, H1m},     -1/2, 2,  1};
ScalarFields[[2]] =  {H2, 1, {H2p, H20},     1/2, 2,  1};
ScalarFields[[3]] =  {S, 1, Sing,              0, 1,  1}; 


     


(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][Additional]= {
	{LagHC, { AddHC->True}},
	{LagNoHC,{ AddHC->False}}
};


LagNoHC = -(M112 conj[H1].H1 + M222 conj[H2].H2 + 1/2 Lambda1 conj[H1].H1.conj[H1].H1 + \
		 1/2 Lambda2 conj[H2].H2.conj[H2].H2 + Lambda3 conj[H2].H2.conj[H1].H1 - Lambda4 (conj[H2].H1).(conj[H1].H2) 
		+ M332 conj[S].S + Lambda5 conj[S].S.conj[H1].H1 + Lambda6 conj[S].S.conj[H2].H2 + Lambda8 conj[S].S.conj[S].S);


LagHC = -(-Yd H1.d.q -Ye H1.e.l +Yu H2.u.q
 - M123 H1.H2.S -1/3 M5 S.S.S + Lambda7 conj[S].conj[S].H1.H2);

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW}
};     
        
        
          	

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs]= 
{    {H10, {v1, 1/Sqrt[2]}, {sigma1, I/Sqrt[2]},{phi1, 1/Sqrt[2]}},
     {H20, {v2, 1/Sqrt[2]}, {sigma2, I/Sqrt[2]},{phi2, 1/Sqrt[2]}},
     {Sing  , {vS, 1/Sqrt[2]}, {sigmaS, I/Sqrt[2]},{phiS, 1/Sqrt[2]}}     }; 
 

DEFINITION[EWSB][MatterSector]=   
    { {{phi1, phi2, phiS}, {hh, ZH}}, 
      {{sigma1, sigma2, sigmaS}, {Ah, ZA}}, 
      {{H1m, conj[H2p]},{Hm,ZP}}, 
      {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
      {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
      {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}       };  


(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  DL, conj[DR]},
 Fe ->{  EL, conj[ER]},
 Fu ->{  UL, conj[UR]},
 Fv ->{  vL, 0}};

DEFINITION[EWSB][GaugeES]={
 Fd1 ->{  FdL, 0},
 Fd2 ->{  0, FdR},
 Fu1 ->{  Fu1, 0},
 Fu2 ->{  0, Fu2},
 Fe1 ->{  Fe1, 0},
 Fe2 ->{  0, Fe2}};



