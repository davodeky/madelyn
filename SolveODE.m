
function[] = SolveODE()
%This function solves the ordinary differential equations
global DEF %global DEF defines model parameters
global SYSTEM %global SYSTEM defines sytem parameters
global DRUG %global DRUG defines drug parameters
global DDI %global DDI enhances drug parameters for DDI prediction
global STUDY %global STUDY defines study design parameters
global MODEL %global MODEL defines parameters important for modeling
global RES %global RES saves the results for post-processing
%__Numbering variables for ODE solver___________________________________________
CompNo = SYSTEM.CompNo; %number of organs / tissues + blood
CoSeNo = SYSTEM.CompNo + SYSTEM.SegNo; %compartments + intestinal segments
SubNo = SYSTEM.SubNo; %number of subcompartments
IndNo = STUDY.IndNo; %number of virtual individuals
DrugNo = DRUG.DrugNo; %number of drugs
DDINo = 1; %number of DDI simulations
ODENo = MODEL.ODENo; %number of ODE equations
%%% Model structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LU = DEF.lung; AD = DEF.adipose; BO = DEF.bone;
BR = DEF.brain; GO = DEF.gonads; HE = DEF.heart;
KI = DEF.kidney; MU = DEF.muscle; SK = DEF.skin;
TH = DEF.thymus; GU = DEF.gut; SP = DEF.spleen;
PA = DEF.pancreas; LI = DEF.liver; LN = DEF.lymphnode;
RE = DEF.remaining; VB = DEF.plasma; AB = DEF.RBC;
ST = CompNo + DEF.stomach; DU = CompNo + DEF.duodenum;
JE = CompNo + DEF.jejunum; IL = CompNo + DEF.ileum;
CN = CompNo + DEF.colon; FS = CompNo + DEF.faeces;
C2D6 = CoSeNo + DEF.CYP2D6; C3A4 = CoSeNo + DEF.CYP3A4;
C3A5 = CoSeNo + DEF.CYP3A5; C2J2 = CoSeNo + DEF.CYP2J2;
%%% No of equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DRUG.DrugNo == 1
 NumEquations = ODENo * DrugNo;
else
 NumEquations = ODENo * DDINo;
end
%%% intestinal segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sto = DEF.stomach; duo = DEF.duodenum; jej = DEF.jejunum;
ile = DEF.ileum; col = DEF.colon;
%%% CYP enzymes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CYP2D6 = DEF.CYP2D6; CYP3A4 = DEF.CYP3A4;
CYP3A5 = DEF.CYP3A5; CYP2J2 = DEF.CYP2J2;
%__System data__________________________________________________________________
%%% Volumes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vvb = SYSTEM.Vvein; %volume for the venous blood pool
Vab = SYSTEM.Vartery; %volume for the arterial blood pool
Vvas = SYSTEM.Vvas; %vascular volume for each compartment
Vint = SYSTEM.Vint; %interstitial volume for each compartment
Vcel = SYSTEM.Vcel; %intracellular volume for each compartment
Vlum = SYSTEM.VlumCAT; %volume of the intestinal lumen
Vent = SYSTEM.VentCAT; %volume of enterocytes in each segment of the intestine
%%% Blood and lymph flows %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CO = SYSTEM.Qorg(:, LU); %cardiac output
Qorg = SYSTEM.Qorg; %regional blood flows
QHA = SYSTEM.QHA; %hepatic arterial blodo flow
QBY = SYSTEM.QBY; %blood flow of the liver bypass
Porg = Qorg.*(1-SYSTEM.HCT); %plasma flow

Ltot = SYSTEM.TotLymphFlow ; %total lymph flow
Lorg = SYSTEM.Lorg; %regional lymph flows
QL = Qorg - Lorg; %substract regional blood from lymph flows
PL = Porg - Lorg; %substract regionalplasma from lymph flows
%%% GI tract %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ITT = 1./SYSTEM.TransitT; %intestinal transit time
AbCYPin = SYSTEM.CYPseg_AB .* 10^3; %abundance of intestinal enzymes in [pmol]
kdCYPin = SYSTEM.CYPin_kdeg; %degradation rate of intestinal enzymes
%%% Liver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WLI = SYSTEM.Worg(:,LI); %liver weight
MPPGL = SYSTEM.MPPGL; %microsomal protein per gram liver (MPPGL)
AbCYPhe = SYSTEM.CYPhe_AB; %hepatic CYP enzyme abundance
kdCYPhe = SYSTEM.CYPhe_kdeg; %degradation rate of hepatic CYP enzymes
%__Drug data____________________________________________________________________
%%% PhysChem properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MW = DDI.MolW; %molecular weight
fuine = DDI.fuine; %fraction unbound in the interstitial space
fucel = DDI.fucel; %fraction unbound in the intracellular space
BP = DDI.BP; %blood-to-plasma ratio
%%% Absorption %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLab = DDI.CLab; %absoprtion flux from the lumen into the enterocytes
LagR = DDI.LagRate; %lag rate to delay Cmax/Tmax - parameter is artificial
%%% Distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jin = DDI.Jin; %flux from the interstitial to the intracellular space
Jout = DDI.Jout; %flux from the intracellular to the interstitial space
%%% Metabolism and Elimination %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLint = DDI.CLint; %intrinsic clearance of an unspecified enzyme pathway
CLbi = DDI.CLbi; %biliary clearance
CLre = DDI.CLre; %renal clearance
CLad = DDI.CLad; %additional plasma clearance
VmCYP = DDI.Vmax_CYP; %Vmax for CYP enzymes
KmCYP = DDI.Km_CYP; %KM for CYP enzymes
CiCYP = DDI.CLint_CYP; %intrinisc clearance for CYP enzymes
%%% DDIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fco = STUDY.DDIMat{2,1}; %factor for concentration
 %1 = automatic concentration is used
%0 = concentration of the perpetrator is used
%factor for concentration of the perpetrator in the relevant compartment
FcelLI = STUDY.DDIMat{2,2}; %intracellular concentration in the liver
FentDU = STUDY.DDIMat{2,3}; %enterocytic concentration in the duodenum
FentJE = STUDY.DDIMat{2,4}; %enterocytic concentration in the jejunum
FentIL = STUDY.DDIMat{2,5}; %enterocytic concentration in the ileum
%fraction unbound ofr the perpetrator
fuLIcel = STUDY.DDIMat{2,6}; %fraction unbound in the liver
fuGUcel = STUDY.DDIMat{2,7}; %fraction unbound in the enterocytes
FKiCYP = 0; %factor for competitive inhibition
 %0 = no competitive inhibition
%1 = competitive inhibition is considered
KiCYP = 0; %inhibition constant
kinactCYP = 0; %maximum inactivation rate constant
KappCYP = 0; %apparent enzyme inhibition constant
IndMaxCYP = 0; %maximum fold of induction
IC50CYP = 0; %half-maximum induction
%__Study design_________________________________________________________________
Dose = reshape(STUDY.DoseEventMat(:, 3, :), STUDY.NoEvents, DRUG.DrugNo);
AdminRoute = reshape(STUDY.DoseEventMat(:, 4, :), STUDY.NoEvents, DRUG.DrugNo);
%rep = DDI.DDINo / DRUG.DrugNo;
rep = 1;
Dose = repmat(Dose, 1, rep);
AdminRoute = repmat(AdminRoute, 1, rep);
StartT = STUDY.DoseEventMat(:, 1, 1); %start time for drug administration
EndT = STUDY.DoseEventMat(:, 2, 1); %end time for drug administration
NP = STUDY.DoseEventMat(:, 5, 1); %resolution for one dosing event
NumPoints = sum(STUDY.DoseEventMat(:, 5, 1)); %resolution for entire simulation
%__Set up initial conditions____________________________________________________
%predefine matrices for the time and the concentration output
Conc = zeros(NumPoints, NumEquations*IndNo);
%set up the right concentration for multiple dosing
MultConc = [STUDY.DoseEventMat(:, 6, 1) + 1,...
 STUDY.DoseEventMat(:, 5, 1) + STUDY.DoseEventMat(:, 6, 1)];
disp('Start ODE');
for ind = 1:IndNo
 %save time and concentration for each subject
 TInd = zeros(NumPoints, 1);
 CInd = zeros(NumPoints, NumEquations);
 %command to load a function to estimate the initial concentration C0
 C0 = Initialise_Conc();
 for ne = 1:STUDY.NoEvents
 %initial concentration for multiple dosing M0
 M0 = Initialise_Mult(C0);
%__Solve equations______________________________________________________________
 %use a stiff solver - alternatively ode45 might be used
 sol = ode15s(@rhs_function, [StartT(ne,1),EndT(ne,1)], M0);

%__Save solution to a vector____________________________________________________
 %time vector for each multiple dosing step is generated
 T = linspace(StartT(ne, 1), EndT(ne, 1), NP(ne, 1))';

 %evaluate the solution sol at each timepoint T
 C = deval(sol, T)';

 TInd(MultConc(ne, 1) : MultConc(ne,2), 1) = T;
 CInd(MultConc(ne, 1) : MultConc(ne,2), :) = C;

 %combine concentration for each individual
 for eq = 1:NumEquations
 Conc(:, (ind-1) * ODENo * DDINo + eq) = CInd(:, eq);
 end
 fprintf('No of event: = %g\n',ne);
 end
 fprintf('No of subject: = %g\n',ind);
end
%__Save solution of the ODE solver globally for post-processing_________________
RES.Time = TInd; %time in h
RES.Conc = Conc; %concentration in microM
MODEL.NP = NP;
MODEL.NumPoints = NumPoints;
disp('End ODE');
%===============================================================================
%%% USED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
%__Initialise conditions for single dose / first dose___________________________
function C0 = Initialise_Conc()

 %it is assumed that each concentration is 0
 C0 = zeros(1, NumEquations);

 for cdi = 1:DDINo
 %%% initial values for CYP abundance (necessary for MBI & induction) %%%
 c2D6 = ODENo*(cdi-1) + C2D6 + SubNo(C2D6); %CYP2D6
 c3A4 = ODENo*(cdi-1) + C3A4 + SubNo(C3A4); %CYP3A4
 c3A5 = ODENo*(cdi-1) + C3A5 + SubNo(C3A5); %CYP3A5
 c2J2 = ODENo*(cdi-1) + C2J2 + SubNo(C2J2); %CYP2J2

 %initial values for hepatic and intestinal CYP2D6
 C0(1, c2D6) = AbCYPhe(ind, CYP2D6);
 C0(1, c2D6+1) = AbCYPin(ind, CYP2D6, duo);
 C0(1, c2D6+2) = AbCYPin(ind, CYP2D6, jej);
 C0(1, c2D6+3) = AbCYPin(ind, CYP2D6, ile);

 %initial values for hepatic and intestinal CYP3A4
 C0(1, c3A4) = AbCYPhe(ind, CYP3A4);
 C0(1, c3A4+1) = AbCYPin(ind, CYP3A4, duo);
 C0(1, c3A4+2) = AbCYPin(ind, CYP3A4, jej);
 C0(1, c3A4+3) = AbCYPin(ind, CYP3A4, ile);
 %initial values for hepatic and intestinal CYP3A5
 C0(1, c3A5) = AbCYPhe(ind, CYP3A5);
 C0(1,c3A5+1) = AbCYPin(ind, CYP3A5, duo);
 C0(1,c3A5+2) = AbCYPin(ind, CYP3A5, jej);
 C0(1,c3A5+3) = AbCYPin(ind, CYP3A5, ile);

 %initial value for hepatic CYP2J2
 C0(1, c2J2) = AbCYPhe(ind, CYP2J2);
 %%% initial values for drug concentration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 switch AdminRoute(1, cdi)

 %venous concentration in the case of intravenous administration
 case DEF.iv
 vb = ODENo*(cdi-1) + VB + SubNo(VB);
 C0(1, vb) = ((Dose(1, cdi) * 1000) / MW(cdi)) / Vvb(ind);

 %stomach concentration in the case of oral drug administration
 case DEF.oral
 st = ODENo*(cdi-1) + ST + SubNo(ST);
 C0(1, st) = ((Dose(1,cdi) * 1000) / MW(cdi)) / Vlum(ind, sto);
 end
 end
end
%__Initialise conditions for multiple dosing____________________________________
function M0 = Initialise_Mult(C0)

 %single or first dose
 if ne == 1
 %C0 has already been defined for a single / first dose
 M0 = C0;

 else
 %for all multiple doses, the concentration is defined by the solution sol
 M0 = CInd(MultConc(ne-1, 2), :);

 for mdi = 1:DDINo
 switch AdminRoute(ne,mdi)

 %venous concentration in the case of intravenous administration
 case DEF.iv
 vb = ODENo*(mdi-1) + VB + SubNo(VB);
 M0(1, vb) = (((Dose(ne, mdi) * 1000) / MW(mdi)) / ...
 Vvb(ind)) + CInd(MultConc(ne-1, 2), vb);

 %stomach concentration in the case of oral drug administration
 case DEF.oral
 st = ODENo*(mdi-1) + ST + SubNo(ST);
 M0(1, st) = (((Dose(ne, mdi) * 1000) / MW(mdi)) / ...
 Vlum(ind, sto)) + CInd(MultConc(ne-1, 2), st);
 end
 end
 end
end
%__Define the right hand site of the equations__________________________________
function dtdy = rhs_function(~, y)

 %prepare the output as a column vector
 dtdy = zeros(NumEquations, 1);

 for d = 1:DDINo
 %define index for each compartment
 lu = ODENo*(d-1) + LU + SubNo(LU); %lung
 ad = ODENo*(d-1) + AD + SubNo(AD); %adipose
 bo = ODENo*(d-1) + BO + SubNo(BO); %bone
 br = ODENo*(d-1) + BR + SubNo(BR); %brain
 go = ODENo*(d-1) + GO + SubNo(GO); %gonads
 he = ODENo*(d-1) + HE + SubNo(HE); %heart
 ki = ODENo*(d-1) + KI + SubNo(KI); %kidney
 mu = ODENo*(d-1) + MU + SubNo(MU); %muscle
 sk = ODENo*(d-1) + SK + SubNo(SK); %skin
 th = ODENo*(d-1) + TH + SubNo(TH); %thymus
 gu = ODENo*(d-1) + GU + SubNo(GU); %gut
 sp = ODENo*(d-1) + SP + SubNo(SP); %spleen
 pa = ODENo*(d-1) + PA + SubNo(PA); %pancreas
 li = ODENo*(d-1) + LI + SubNo(LI); %liver
 ln = ODENo*(d-1) + LN + SubNo(LN); %lymphnode
 re = ODENo*(d-1) + RE + SubNo(RE); %remaining
 vb = ODENo*(d-1) + VB + SubNo(VB); %venous
 ab = ODENo*(d-1) + AB + SubNo(AB); %arterial

 st = ODENo*(d-1) + ST + SubNo(ST); %stomach
 du = ODENo*(d-1) + DU + SubNo(DU); %duodenum
 je = ODENo*(d-1) + JE + SubNo(JE); %jejunum
 il = ODENo*(d-1) + IL + SubNo(IL); %ileum
 cn = ODENo*(d-1) + CN + SubNo(CN); %colon
 fs = ODENo*(d-1) + FS + SubNo(FS); %faeces

 c2D6 = ODENo*(d-1) + C2D6 + SubNo(C2D6); %CYP2D6
 c3A4 = ODENo*(d-1) + C3A4 + SubNo(C3A4); %CYP3A4
 c3A5 = ODENo*(d-1) + C3A5 + SubNo(C3A5); %CYP3A5
 c2J2 = ODENo*(d-1) + C2J2 + SubNo(C2J2); %CYP2J2

 %define mechanism-based inhibition (MBI)
 MechLI = sum((kinactCYP(:, d, :) .* y((li+2)*Fco(:, d)+FcelLI(:, d)) .* ...
fuLIcel(:, d, ind)) ./ (KappCYP(:, d, :) + y((li+2)*Fco(:, d)+FcelLI(:, d)) .* ...
fuLIcel(:, d, ind)), 1);

 MechDU = sum((kinactCYP(:, d, :) .* y((du+2)*Fco(:, d)+FentDU(:, d)) .* ...
fuGUcel(:, d, ind)) ./ (KappCYP(:, d, :) + y((du+2)*Fco(:, d)+FentDU(:, d)) .* ...
fuGUcel(:, d, ind)), 1);
 MechJE = sum((kinactCYP(:, d, :) .* y((je+2)*Fco(:, d)+FentJE(:, d)) .* ...
fuGUcel(:, d, ind)) ./ (KappCYP(:, d, :) + y((je+2)*Fco(:, d)+FentJE(:, d)) .* ...
fuGUcel(:, d, ind)), 1);
 MechIL = sum((kinactCYP(:, d, :) .* y((il+2)*Fco(:, d)+FentIL(:, d)) .* ...
fuGUcel(:, d, ind)) ./ (KappCYP(:, d, :) + y((il+2)*Fco(:, d)+FentIL(:, d)) .* ...
fuGUcel(:, d, ind)), 1);

 %define induction
 InduLI = sum((IndMaxCYP(:, d, :) .* y((li+2)*Fco(:, d)+FcelLI(:, d)) .* ...
fuLIcel(:, d, ind)) ./ (IC50CYP(:, d, :) + y((li+2)*Fco(:, d)+FcelLI(:, d)) .* ...
fuLIcel(:, d, ind)), 1);

 InduDU = sum((IndMaxCYP(:, d, :) .* y((du+2)*Fco(:, d)+FentDU(:, d)) .* ... 
fuGUcel(:, d, ind)) ./ (IC50CYP(:, d, :) + y((du+2)*Fco(:, d)+FentDU(:, d)) .* ...
fuGUcel(:, d, ind)), 1);
 InduJE = sum((IndMaxCYP(:, d, :) .* y((je+2)*Fco(:, d)+FentJE(:, d)) .* ...
fuGUcel(:, d, ind)) ./ (IC50CYP(:, d, :) + y((je+2)*Fco(:, d)+FentJE(:, d)) .* ...
fuGUcel(:, d, ind)), 1);
 InduIL = sum((IndMaxCYP(:, d, :) .* y((il+2)*Fco(:, d)+FentIL(:, d)) .* ...
fuGUcel(:, d, ind)) ./ (IC50CYP(:, d, :) + y((il+2)*Fco(:, d)+FentIL(:, d)) .* ...
fuGUcel(:, d, ind)), 1);

 %calculate the dynamic CYP abundance
 %CYP abundance should be zero in this case???
 dtdy(c2D6) = kdCYPhe(ind, CYP2D6) * AbCYPhe(ind, CYP2D6) * (1 + InduLI ...
(CYP2D6)) - (kdCYPhe(ind, CYP2D6) + MechLI(CYP2D6)) * y(c2D6);
 dtdy(c2D6+1) = kdCYPin(ind, CYP2D6) * AbCYPin(ind, CYP2D6, duo) * (1 + ...
InduDU(CYP2D6)) - (kdCYPin(ind, CYP2D6) + MechDU(CYP2D6)) * y(c2D6+1);
 dtdy(c2D6+2) = kdCYPin(ind, CYP2D6) * AbCYPin(ind, CYP2D6, jej) * (1 + ...
InduJE(CYP2D6)) - (kdCYPin(ind, CYP2D6) + MechJE(CYP2D6)) * y(c2D6+2);
 dtdy(c2D6+3) = kdCYPin(ind, CYP2D6) * AbCYPin(ind, CYP2D6, ile) * (1 + ...
InduIL(CYP2D6)) - (kdCYPin(ind, CYP2D6) + MechIL(CYP2D6)) * y(c2D6+3);

 dtdy(c3A4) = kdCYPhe(ind, CYP3A4) * AbCYPhe(ind, CYP3A4) * (1 + InduLI ...
(CYP3A4)) - (kdCYPhe(ind, CYP3A4) + MechLI(CYP3A4)) * y(c3A4);
 dtdy(c3A4+1) = kdCYPin(ind, CYP3A4) * AbCYPin(ind, CYP3A4, duo) * (1 + ...
InduDU(CYP3A4)) - (kdCYPin(ind, CYP3A4) + MechDU(CYP3A4)) * y(c3A4+1);
 dtdy(c3A4+2) = kdCYPin(ind, CYP3A4) * AbCYPin(ind, CYP3A4, jej) * (1 + ...
InduJE(CYP3A4)) - (kdCYPin(ind, CYP3A4) + MechJE(CYP3A4)) * y(c3A4+2);
 dtdy(c3A4+3) = kdCYPin(ind, CYP3A4) * AbCYPin(ind, CYP3A4, ile) * (1 + ...
InduIL(CYP3A4)) - (kdCYPin(ind, CYP3A4) + MechIL(CYP3A4)) * y(c3A4+3);

 dtdy(c3A5) = kdCYPhe(ind, CYP3A5) * AbCYPhe(ind, CYP3A5) * (1 + InduLI ...
(CYP3A5)) - (kdCYPhe(ind, CYP3A5) + MechLI(CYP3A5)) * y(c3A5);
 dtdy(c3A5+1) = kdCYPin(ind, CYP3A5) * AbCYPin(ind, CYP3A5, duo) * (1 + ...
InduDU(CYP3A5)) - (kdCYPin(ind, CYP3A5) + MechDU(CYP3A5)) * y(c3A5+1);
 dtdy(c3A5+2) = kdCYPin(ind, CYP3A5) * AbCYPin(ind, CYP3A5, jej) * (1 + ...
InduJE(CYP3A5)) - (kdCYPin(ind, CYP3A5) + MechJE(CYP3A5)) * y(c3A5+2);
 dtdy(c3A5+3) = kdCYPin(ind, CYP3A5) * AbCYPin(ind, CYP3A5, ile) * (1 + ...
InduIL(CYP3A5)) - (kdCYPin(ind, CYP3A5) + MechIL(CYP3A5)) * y(c3A5+3);

 dtdy(c2J2) = kdCYPhe(ind, CYP2J2) * AbCYPhe(ind, CYP2J2) * (1 + InduLI ...
(CYP2J2)) - (kdCYPhe(ind, CYP2J2) + MechLI(CYP2J2)) * y(c2J2);

 %define competitive inhibition
 CYPComLI = reshape(1 + sum((((y((li+2)*Fco(:, d)+FcelLI(:, d)) .* fuLIcel...
(:, d, ind) .* FKiCYP(:, d, :)) ./ KiCYP(:, d, :))), 1), SYSTEM.CYPliNo, 1);

 CYPComDU = reshape(1 + sum((((y((du+2)*Fco(:, d)+FentDU(:, d)) .* fuGUcel...
(:, d, ind) .* FKiCYP(:, d, :)) ./ KiCYP(:, d, :))), 1), SYSTEM.CYPliNo, 1);
 CYPComJE = reshape(1 + sum((((y((je+2)*Fco(:, d)+FentJE(:, d)) .* fuGUcel...
(:, d, ind) .* FKiCYP(:, d, :)) ./ KiCYP(:, d, :))), 1), SYSTEM.CYPliNo, 1);
 CYPComIL = reshape(1 + sum((((y((il+2)*Fco(:, d)+FentIL(:, d)) .* fuGUcel...
(:, d, ind) .* FKiCYP(:, d, :)) ./ KiCYP(:, d, :))), 1), SYSTEM.CYPliNo, 1);

 %define enzymatic metbaolism
 CYPMetLI = (VmCYP(:, d) ./ ((KmCYP(:, d) .* CYPComLI) + y(li+2)*fucel(ind,...
LI, d))) +...
 (CiCYP(:, d) ./ CYPComLI);

 CYPMetDU = (VmCYP(:, d) ./ ((KmCYP(:, d) .* CYPComDU) + y(du+2)*fucel(ind,...
GU, d))) +...
 (CiCYP(:, d) ./ CYPComDU);
 CYPMetJE = (VmCYP(:, d) ./ ((KmCYP(:, d) .* CYPComJE) + y(je+2)*fucel(ind,...
GU, d))) +...
 (CiCYP(:, d) ./ CYPComJE);
 CYPMetIL = (VmCYP(:, d) ./ ((KmCYP(:, d) .* CYPComIL) + y(il+2)*fucel(ind,...
GU, d))) +...
 (CiCYP(:, d) ./ CYPComIL);

 %vascular space of compartments
 dtdy(lu) = (1/Vvas(ind,LU)) * (Qorg(ind,LU)*y(vb) - QL(ind,LU)*y(lu) -...
Porg(ind,LU)*(y(lu)/BP(d)) + PL(ind,LU)*y(lu+1));
 dtdy(ad) = (1/Vvas(ind,AD)) * (Qorg(ind,AD)*y(ab) - QL(ind,AD)*y(ad) -...
Porg(ind,AD)*(y(ad)/BP(d)) + PL(ind,AD)*y(ad+1));
 dtdy(bo) = (1/Vvas(ind,BO)) * (Qorg(ind,BO)*y(ab) - QL(ind,BO)*y(bo) -...
Porg(ind,BO)*(y(bo)/BP(d)) + PL(ind,BO)*y(bo+1));
 dtdy(br) = (1/Vvas(ind,BR)) * (Qorg(ind,BR)*y(ab) - QL(ind,BR)*y(br) -...
Porg(ind,BR)*(y(br)/BP(d)) + PL(ind,BR)*y(br+1));
 dtdy(go) = (1/Vvas(ind,GO)) * (Qorg(ind,GO)*y(ab) - QL(ind,GO)*y(go) -...
Porg(ind,GO)*(y(go)/BP(d)) + PL(ind,GO)*y(go+1));
 dtdy(he) = (1/Vvas(ind,HE)) * (Qorg(ind,HE)*y(ab) - QL(ind,HE)*y(he) -...
Porg(ind,HE)*(y(he)/BP(d)) + PL(ind,HE)*y(he+1));
 dtdy(ki) = (1/Vvas(ind,KI)) * (Qorg(ind,KI)*y(ab) - QL(ind,KI)*y(ki) -...
Porg(ind,KI)*(y(ki)/BP(d)) + PL(ind,KI)*y(ki+1) - CLre(ind,d)*y(ki));
 dtdy(mu) = (1/Vvas(ind,MU)) * (Qorg(ind,MU)*y(ab) - QL(ind,MU)*y(mu) -...
Porg(ind,MU)*(y(mu)/BP(d)) + PL(ind,MU)*y(mu+1));
 dtdy(sk) = (1/Vvas(ind,SK)) * (Qorg(ind,SK)*y(ab) - QL(ind,SK)*y(sk) -...
Porg(ind,SK)*(y(sk)/BP(d)) + PL(ind,SK)*y(sk+1));
 dtdy(th) = (1/Vvas(ind,TH)) * (Qorg(ind,TH)*y(ab) - QL(ind,TH)*y(th) -...
Porg(ind,TH)*(y(th)/BP(d)) + PL(ind,TH)*y(th+1));
 dtdy(gu) = (1/Vvas(ind,GU)) * (Qorg(ind,GU)*y(ab) - QL(ind,GU)*y(gu) -...
Porg(ind,GU)*(y(gu)/BP(d)) + PL(ind,GU)*y(gu+1));
 dtdy(sp) = (1/Vvas(ind,SP)) * (Qorg(ind,SP)*y(ab) - QL(ind,SP)*y(sp) -...
Porg(ind,SP)*(y(sp)/BP(d)) + PL(ind,SP)*y(sp+1));
 dtdy(pa) = (1/Vvas(ind,PA)) * (Qorg(ind,PA)*y(ab) - QL(ind,PA)*y(pa) -...
Porg(ind,PA)*(y(pa)/BP(d)) + PL(ind,PA)*y(pa+1));
 dtdy(li) = (1/Vvas(ind,LI)) * (QHA(ind)*y(ab) + QBY(ind)*y(ab) + QL(ind,...
GU)*y(gu) + QL(ind,SP)*y(sp) + QL(ind,PA)*y(pa) - QL(ind,LI)*y(li) - Porg(ind,LI)*...
(y(li)/BP(d)) + PL(ind,LI)*y(li+1));
 dtdy(ln) = (1/Vvas(ind,LN)) * (Qorg(ind,LN)*y(ab) - Qorg(ind,LN)*y(ln) +...
Ltot(ind)*y(ln+1) - Ltot(ind)*y(ln));
 dtdy(re) = (1/Vvas(ind,RE)) * (Qorg(ind,RE)*y(ab) - QL(ind,RE)*y(re) -...
Porg(ind,RE)*(y(re)/BP(d)) + PL(ind,RE)*y(re+1));

 %interstitial space of compartments
 dtdy(lu+1) = (1/Vint(ind,LU)) * (Porg(ind,LU)*(y(lu)/BP(d)) - PL(ind,LU)*y...
(lu+1) - Lorg(ind,LU)*y(lu+1) - Jin(ind,LU,d)*y(lu+1)*fuine(ind,LU,d) + Jout(ind,...
LU,d)*y(lu+2)*fucel(ind,LU,d));
 dtdy(ad+1) = (1/Vint(ind,AD)) * (Porg(ind,AD)*(y(ad)/BP(d)) - PL(ind,AD)*y...
(ad+1) - Lorg(ind,AD)*y(ad+1) - Jin(ind,AD,d)*y(ad+1)*fuine(ind,AD,d) + Jout(ind,...
AD,d)*y(ad+2)*fucel(ind,AD,d));
 dtdy(bo+1) = (1/Vint(ind,BO)) * (Porg(ind,BO)*(y(bo)/BP(d)) - PL(ind,BO)*y...
(bo+1) - Lorg(ind,BO)*y(bo+1) - Jin(ind,BO,d)*y(bo+1)*fuine(ind,BO,d) + Jout(ind,...
BO,d)*y(bo+2)*fucel(ind,BO,d));
 dtdy(br+1) = (1/Vint(ind,BR)) * (Porg(ind,BR)*(y(br)/BP(d)) - PL(ind,BR)*y...
(br+1) - Lorg(ind,BR)*y(br+1) - Jin(ind,BR,d)*y(br+1)*fuine(ind,BR,d) + Jout(ind,...
BR,d)*y(br+2)*fucel(ind,BR,d));
 dtdy(go+1) = (1/Vint(ind,GO)) * (Porg(ind,GO)*(y(go)/BP(d)) - PL(ind,GO)*y...
(go+1) - Lorg(ind,GO)*y(go+1) - Jin(ind,GO,d)*y(go+1)*fuine(ind,GO,d) + Jout(ind,...
GO,d)*y(go+2)*fucel(ind,GO,d));
 dtdy(he+1) = (1/Vint(ind,HE)) * (Porg(ind,HE)*(y(he)/BP(d)) - PL(ind,HE)*y...
(he+1) - Lorg(ind,HE)*y(he+1) - Jin(ind,HE,d)*y(he+1)*fuine(ind,HE,d) + Jout(ind,...
HE,d)*y(he+2)*fucel(ind,HE,d));
 dtdy(ki+1) = (1/Vint(ind,KI)) * (Porg(ind,KI)*(y(ki)/BP(d)) - PL(ind,KI)*y...
(ki+1) - Lorg(ind,KI)*y(ki+1) - Jin(ind,KI,d)*y(ki+1)*fuine(ind,KI,d) + Jout(ind,...
KI,d)*y(ki+2)*fucel(ind,KI,d));
 dtdy(mu+1) = (1/Vint(ind,MU)) * (Porg(ind,MU)*(y(mu)/BP(d)) - PL(ind,MU)*y...
(mu+1) - Lorg(ind,MU)*y(mu+1) - Jin(ind,MU,d)*y(mu+1)*fuine(ind,MU,d) + Jout(ind,...
MU,d)*y(mu+2)*fucel(ind,MU,d));
 dtdy(sk+1) = (1/Vint(ind,SK)) * (Porg(ind,SK)*(y(sk)/BP(d)) - PL(ind,SK)*y...
(sk+1) - Lorg(ind,SK)*y(sk+1) - Jin(ind,SK,d)*y(sk+1)*fuine(ind,SK,d) + Jout(ind,...
SK,d)*y(sk+2)*fucel(ind,SK,d));
 dtdy(th+1) = (1/Vint(ind,TH)) * (Porg(ind,TH)*(y(th)/BP(d)) - PL(ind,TH)*y...
(th+1) - Lorg(ind,TH)*y(th+1) - Jin(ind,TH,d)*y(th+1)*fuine(ind,TH,d) + Jout(ind,...
TH,d)*y(th+2)*fucel(ind,TH,d));
 dtdy(gu+1) = (1/Vint(ind,GU)) * (Porg(ind,GU)*(y(gu)/BP(d)) - PL(ind,GU)*y...
(gu+1) - Lorg(ind,GU)*y(gu+1)+...
 Jin(ind,GU,d)*y(du+2)*fucel(ind,GU,d)+ ...
 Jin(ind,GU,d)*y(je+2)*fucel(ind,GU,d)+ ...
 Jin(ind,GU,d)*y(il+2)*fucel(ind,GU,d)+ ...
 Jin(ind,GU,d)*y(cn+2)*fucel(ind,GU,d));
 dtdy(sp+1) = (1/Vint(ind,SP)) * (Porg(ind,SP)*(y(sp)/BP(d)) - PL(ind,SP)*y...
(sp+1) - Lorg(ind,SP)*y(sp+1) - Jin(ind,SP,d)*y(sp+1)*fuine(ind,SP,d) + Jout(ind,...
SP,d)*y(sp+2)*fucel(ind,SP,d));
 dtdy(pa+1) = (1/Vint(ind,PA)) * (Porg(ind,PA)*(y(pa)/BP(d)) - PL(ind,PA)*y...
(pa+1) - Lorg(ind,PA)*y(pa+1) - Jin(ind,PA,d)*y(pa+1)*fuine(ind,PA,d) + Jout(ind,...
PA,d)*y(pa+2)*fucel(ind,PA,d));
 dtdy(li+1) = (1/Vint(ind,LI)) * (Porg(ind,LI)*(y(li)/BP(d)) - PL(ind,LI)*y...
(li+1) - Lorg(ind,LI)*y(li+1) - Jin(ind,LI,d)*y(li+1)*fuine(ind,LI,d) + Jout(ind,...
LI,d)*y(li+2)*fucel(ind,LI,d));
 dtdy(ln+1) = (1/Vint(ind,LN)) * (Lorg(ind,AD)*y(ad+1) + Lorg(ind,BO)*y...
(bo+1) + Lorg(ind,BR)*y(br+1) +...
 Lorg(ind,GO)*y(go+1) + Lorg(ind,HE)*y...
(he+1) + Lorg(ind,KI)*y(ki+1) +...
 Lorg(ind,MU)*y(mu+1) + Lorg(ind,SK)*y...
(sk+1) + Lorg(ind,TH)*y(th+1) +...
 Lorg(ind,GU)*y(gu+1) + Lorg(ind,SP)*y...
(sp+1) + Lorg(ind,PA)*y(pa+1) +...
 Lorg(ind,LI)*y(li+1) + Lorg(ind,RE)*y...
(re+1) + Lorg(ind,LU)*y(lu+1) - Ltot(ind)*y(ln+1)-...
 Jin(ind,LN,d)*y(ln+1)*fuine(ind,LN,d) +...
Jout(ind,LN,d)*y(ln+2)*fucel(ind,LN,d));
 dtdy(re+1) = (1/Vint(ind,RE)) * (Porg(ind,RE)*(y(re)/BP(d)) - PL(ind,RE)*y...
(re+1) - Lorg(ind,RE)*y(re+1) - Jin(ind,RE,d)*y(re+1)*fuine(ind,RE,d) + Jout(ind,...
RE,d)*y(re+2)*fucel(ind,RE,d));

 %intracellular space of compartments
 dtdy(lu+2) = (1/Vcel(ind,LU)) * (Jin(ind,LU,d)*y(lu+1)*fuine(ind,LU,d) -...
Jout(ind,LU,d)*y(lu+2)*fucel(ind,LU,d));
 dtdy(ad+2) = (1/Vcel(ind,AD)) * (Jin(ind,AD,d)*y(ad+1)*fuine(ind,AD,d) -...
Jout(ind,AD,d)*y(ad+2)*fucel(ind,AD,d));
 dtdy(bo+2) = (1/Vcel(ind,BO)) * (Jin(ind,BO,d)*y(bo+1)*fuine(ind,BO,d) -...
Jout(ind,BO,d)*y(bo+2)*fucel(ind,BO,d));
 dtdy(br+2) = (1/Vcel(ind,BR)) * (Jin(ind,BR,d)*y(br+1)*fuine(ind,BR,d) -...
Jout(ind,BR,d)*y(br+2)*fucel(ind,BR,d));
 dtdy(go+2) = (1/Vcel(ind,GO)) * (Jin(ind,GO,d)*y(go+1)*fuine(ind,GO,d) -...
Jout(ind,GO,d)*y(go+2)*fucel(ind,GO,d));
 dtdy(he+2) = (1/Vcel(ind,HE)) * (Jin(ind,HE,d)*y(he+1)*fuine(ind,HE,d) -...
Jout(ind,HE,d)*y(he+2)*fucel(ind,HE,d));
 dtdy(ki+2) = (1/Vcel(ind,KI)) * (Jin(ind,KI,d)*y(ki+1)*fuine(ind,KI,d) -...
Jout(ind,KI,d)*y(ki+2)*fucel(ind,KI,d));
 dtdy(mu+2) = (1/Vcel(ind,MU)) * (Jin(ind,MU,d)*y(mu+1)*fuine(ind,MU,d) -...
Jout(ind,MU,d)*y(mu+2)*fucel(ind,MU,d));
 dtdy(sk+2) = (1/Vcel(ind,SK)) * (Jin(ind,SK,d)*y(sk+1)*fuine(ind,SK,d) -...
Jout(ind,SK,d)*y(sk+2)*fucel(ind,SK,d));
 dtdy(th+2) = (1/Vcel(ind,TH)) * (Jin(ind,TH,d)*y(th+1)*fuine(ind,TH,d) -...
Jout(ind,TH,d)*y(th+2)*fucel(ind,TH,d));
 dtdy(sp+2) = (1/Vcel(ind,SP)) * (Jin(ind,SP,d)*y(sp+1)*fuine(ind,SP,d) -...
Jout(ind,SP,d)*y(sp+2)*fucel(ind,SP,d));
 dtdy(pa+2) = (1/Vcel(ind,PA)) * (Jin(ind,PA,d)*y(pa+1)*fuine(ind,PA,d) -...
Jout(ind,PA,d)*y(pa+2)*fucel(ind,PA,d));
 dtdy(li+2) = (1/Vcel(ind,LI)) * (Jin(ind,LI,d)*y(li+1)*fuine(ind,LI,d) -...
Jout(ind,LI,d)*y(li+2)*fucel(ind,LI,d) -...
 (((CYPMetLI(CYP2D6)*y(c2D6)*MPPGL(ind)*WLI...
(ind)*1000) +...
 (CYPMetLI(CYP3A4)*y(c3A4)*MPPGL(ind)*WLI...
(ind)*1000) +...
 (CYPMetLI(CYP3A5)*y(c3A5)*MPPGL(ind)*WLI...
(ind)*1000) +...
 (CYPMetLI(CYP2J2)*y(c2J2)*MPPGL(ind)*WLI...
(ind)*1000) +...
 (CLint(d)*MPPGL(ind)*WLI(ind)*1000) +...
CLbi(ind,d))*y(li+2)*fucel(ind,LI,d)));
 dtdy(ln+2) = (1/Vcel(ind,LN)) * (Jin(ind,LN,d)*y(ln+1)*fuine(ind,LN,d) -...
Jout(ind,LN,d)*y(ln+2)*fucel(ind,LN,d));
 dtdy(re+2) = (1/Vcel(ind,RE)) * (Jin(ind,RE,d)*y(re+1)*fuine(ind,RE,d) -...
Jout(ind,RE,d)*y(re+2)*fucel(ind,RE,d));

 %blood pools
 dtdy(vb) = (1/Vvb(ind)) * (QL(ind,AD)*y(ad) + QL(ind,BO)*y(bo) + QL(ind,...
BR)*y(br) + QL(ind,GO)*y(go) +...
 QL(ind,HE)*y(he) + QL(ind,KI)*y(ki) + QL(ind,...
MU)*y(mu) + QL(ind,SK)*y(sk) +...
 QL(ind,TH)*y(th) + QL(ind,LI)*y(li) + Qorg...
(ind,LN)*y(ln) + QL(ind,RE)*y(re) +...
 Ltot(ind)*y(ln) - CO(ind)*y(vb));

 dtdy(ab) = (1/Vab(ind)) * (QL(ind,LU)*y(lu) - Qorg(ind,AD)*y(ab) - Qorg...
(ind,BO)*y(ab) - Qorg(ind,BR)*y(ab) -...
 Qorg(ind,GO)*y(ab) - Qorg(ind,HE)*y(ab) - Qorg...
(ind,KI)*y(ab) - Qorg(ind,MU)*y(ab) -...
 Qorg(ind,SK)*y(ab) - Qorg(ind,TH)*y(ab) - Qorg...
(ind,GU)*y(ab) - Qorg(ind,SP)*y(ab) -...
 Qorg(ind,PA)*y(ab) - QHA(ind)*y(ab) - QBY(ind)...
*y(ab) - Qorg(ind,LN)*y(ab) - Qorg(ind,RE)*y(ab) -...
 CLad(d)*y(ab));
 %intestinal lumen
 dtdy(st) = (1/Vlum(ind,sto)) * (-ITT(ind,sto)*y(st)*Vlum(ind,sto));

 dtdy(du) = (1/Vlum(ind,duo)) * (ITT(ind,sto)*y(st)*Vlum(ind,sto) - ITT...
(ind,duo)*y(du)*Vlum(ind,duo) - CLab(ind,duo,d)*y(du));
 dtdy(je) = (1/Vlum(ind,jej)) * (ITT(ind,duo)*y(du)*Vlum(ind,duo) - ITT...
(ind,jej)*y(je)*Vlum(ind,jej) - CLab(ind,jej,d)*y(je));
 dtdy(il) = (1/Vlum(ind,ile)) * (ITT(ind,jej)*y(je)*Vlum(ind,jej) - ITT...
(ind,ile)*y(il)*Vlum(ind,ile) - CLab(ind,ile,d)*y(il));
 dtdy(cn) = (1/Vlum(ind,col)) * (ITT(ind,ile)*y(il)*Vlum(ind,ile) - ITT...
(ind,col)*y(cn)*Vlum(ind,col) - CLab(ind,col,d)*y(cn));
 dtdy(fs) = ITT(ind,col)*y(cn)*Vlum(ind,col);

 %artificial uptake when a lag rate is necessary to delay Cmax / Tmax
 dtdy(du+1) = ((CLab(ind,duo,d)*y(du)) / Vlum(ind,duo) - LagR(d)*y(du+1));
 dtdy(je+1) = ((CLab(ind,jej,d)*y(je)) / Vlum(ind,jej) - LagR(d)*y(je+1));
 dtdy(il+1) = ((CLab(ind,ile,d)*y(il)) / Vlum(ind,ile) - LagR(d)*y(il+1));
 dtdy(cn+1) = ((CLab(ind,col,d)*y(cn)) / Vlum(ind,col) - LagR(d)*y(cn+1));
 %enetrocytes
 dtdy(du+2) = (1/Vent(ind,duo)) * (LagR(d)*y(du+1)*Vlum(ind,duo) - Jin(ind,...
GU,d)*y(du+2)*fucel(ind,GU,d) -...
 ((CYPMetDU(CYP2D6)*y(c2D6+1) + CYPMetDU...
(CYP3A4)*y(c3A4+1) + CYPMetDU(CYP3A5)*y(c3A5+1)) * y(du+2)*fucel(ind,GU,d)));
 dtdy(je+2) = (1/Vent(ind,jej)) * (LagR(d)*y(je+1)*Vlum(ind,jej) - Jin(ind,...
GU,d)*y(je+2)*fucel(ind,GU,d) -...
 ((CYPMetJE(CYP2D6)*y(c2D6+2) + CYPMetJE...
(CYP3A4)*y(c3A4+2) + CYPMetJE(CYP3A5)*y(c3A5+2)) * y(je+2)*fucel(ind,GU,d)));
 dtdy(il+2) = (1/Vent(ind,ile)) * (LagR(d)*y(il+1)*Vlum(ind,ile) - Jin(ind,...
GU,d)*y(il+2)*fucel(ind,GU,d) -...
 ((CYPMetIL(CYP2D6)*y(c2D6+3) + CYPMetIL...
(CYP3A4)*y(c3A4+3) + CYPMetIL(CYP3A5)*y(c3A5+3)) * y(il+2)*fucel(ind,GU,d)));
 dtdy(cn+2) = (1/Vent(ind,col)) * (LagR(d)*y(cn+1)*Vlum(ind,col) - Jin(ind,...
GU,d)*y(cn+2)*fucel(ind,GU,d));
 end
end
end