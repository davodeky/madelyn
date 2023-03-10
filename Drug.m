function[] = Drug()
%This function loads the revelant drug files from the drug file library and
%performs the in vitro-to-in vivo extrapolation
global DEF %global DEF defines model parameters
global SYSTEM %global SYSTEM defines sytem parameters
global DRUG %global DRUG defines drug parameters
global DDI %global DDI enhances drug parameters for DDI prediction
global STUDY %global STUDY defines study design parameters
global MODEL %global MODEL defines parameters important for modeling
%define variables used in the script
dru = 5;
DDINo = 1;
%__PreProcessing________________________________________________________________
%prepare all drug parameters for each drug file in the library
%%% physchem characteristics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DRUG.MolW = zeros(1, DEF.DrugLibNo); %molecular weight in [g/mol]
DRUG.logP = zeros(1, DEF.DrugLibNo); %octanol water partition coef
DRUG.pka = zeros(1, DEF.DrugLibNo);
DRUG.BP = zeros(1, DEF.DrugLibNo); %blood-to-plasma ratio
DRUG.fu = zeros(1, DEF.DrugLibNo); %fraction unbound
DRUG.protein = zeros(1, DEF.DrugLibNo); %main binding protein
%%% absorption %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%apparent permeability in 10^-6 cm/sec
DRUG.Papp = zeros(1, DEF.DrugLibNo);
%a rate constant can be introduced to match the observed Tmax
DRUG.LagRate = zeros(1, DEF.DrugLibNo);
%%% distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%universal KP Scalar for all compartments at once
DRUG.KpScalarAll = zeros(1, DEF.DrugLibNo);
%Kp scalar for single compartments
DRUG.KpScalar = zeros(SYSTEM.CompNo, DEF.DrugLibNo);
%%% metabolism & elimination %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DRUG.Vmax_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
DRUG.Km_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
DRUG.CLint_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
DRUG.CLint = zeros(1, DEF.DrugLibNo);
DRUG.CLrenal = zeros(1, DEF.DrugLibNo);
DRUG.CLrenalCV = zeros(1, DEF.DrugLibNo);
DRUG.CLbile = zeros(1, DEF.DrugLibNo);
DRUG.CLbileCV = zeros(1, DEF.DrugLibNo);
DRUG.CLadditional = zeros(1, DEF.DrugLibNo);
DRUG.CLadditionalCV = zeros(1, DEF.DrugLibNo);
%%% drug-drug interactions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DRUG.Ki_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
DRUG.kinact_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
DRUG.Kapp_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
DRUG.IndMax_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
DRUG.IC50_CYP = zeros(SYSTEM.CYPliNo, DEF.DrugLibNo);
%__load the drug models_________________________________________________________
%load the user-defined drug model
% addpath('acetaminophen');
%convert the string of the user chosen drug model to the function in the library
for d = 1:DRUG.DrugNo
 switch d
 case 1
 drug1 = str2func(DEF.name{dru, DRUG.DrugName(d)});
 drug1();

 case 2
 drug2 = str2func(DEF.name{dru, DRUG.DrugName(d)});
 drug2();

 case 3
 drug3 = str2func(DEF.name{dru, DRUG.DrugName(d)});
 drug3();
 end
end
%__PostProcessing_______________________________________________________________
%Extract only non-zero values for used drugs in the simulation based on MW
ParGen = [DRUG.MolW; DRUG.logP; DRUG.pka; DRUG.BP; DRUG.fu; DRUG.protein; ...
 DRUG.Papp; DRUG.LagRate; ...
 DRUG.KpScalarAll; ...
 DRUG.CLint; DRUG.CLrenal; DRUG.CLrenalCV;...
 DRUG.CLbile; DRUG.CLbileCV; ...
 DRUG.CLadditional; DRUG.CLadditionalCV];
ParSystem = [DRUG.KpScalar; DRUG.MolW];
ParCYP = [DRUG.Vmax_CYP; DRUG.Km_CYP; DRUG.CLint_CYP;...
 DRUG.Ki_CYP; DRUG.kinact_CYP; DRUG.Kapp_CYP;...
 DRUG.IndMax_CYP; DRUG.IC50_CYP; DRUG.MolW; ];
%delete zero values
ParGen(:, all (~any (ParGen), 1)) = [];
ParSystem(:, all (~any (ParSystem), 1)) = [];
ParCYP(:, all (~any (ParCYP), 1)) = [];
%factor for DDI studies
if DRUG.DrugNo == 1
 %in the case of only one drug, only auto-induction and inhibition is used
 DDI.DDINo = 1;

else
 DDI.DDINo = 2*DRUG.DrugNo;
end
FacDDI = DDI.DDINo / DRUG.DrugNo;
% FacDDI = 1; %only one drug

DDI.MolW = repmat(ParGen(1, :), 1, FacDDI);
DDI.logP = repmat(ParGen(2, :), 1, FacDDI);
DDI.pka = repmat(ParGen(3, :), 1, FacDDI);
DDI.BP = repmat(ParGen(4, :), 1, FacDDI);
DDI.fu = repmat(ParGen(5, :), 1, FacDDI);
DDI.protein = repmat(ParGen(6, :), 1, FacDDI);
DDI.Papp = repmat(ParGen(7, :), 1, FacDDI);
DDI.LagRate = repmat(ParGen(8, :), 1, FacDDI);
%if there is no delay, the lag rate is set to a very high value
DDI.LagRate(DDI.LagRate == 0) = 1000;
DDI.KpScalarAll = repmat(ParGen(9, :), 1, FacDDI);
DDI.KpScalar = repmat(ParSystem(1:SYSTEM.CompNo, :), 1, FacDDI);
%Kp Scalar cannot be 0, because they are multiplied to Kpu
DDI.KpScalarAll(DDI.KpScalarAll == 0) = 1;
DDI.KpScalar(DDI.KpScalar == 0) = 1;
%Vmax is converted from pmol/min/pmol to micromol/h/pmol
DDI.Vmax_CYP = repmat(ParCYP(1:SYSTEM.CYPliNo, :), 1, FacDDI) .*10^-6 .* 60;
%if Km is 0, it is converted to 1 to prevent dividing by 0 in the code
DDI.Km_CYP = repmat(ParCYP(SYSTEM.CYPliNo+1:2*SYSTEM.CYPliNo, :), 1, FacDDI);
DDI.Km_CYP(DDI.Km_CYP == 0) = 1;
%CLint (CYP) is converted from microL/min/pmol to L/h/pmol
DDI.CLint_CYP = repmat(ParCYP(2*SYSTEM.CYPliNo+1:3*SYSTEM.CYPliNo, :), ...
 1, FacDDI) .*10^-6 .* 60;
%CLint (tot hep) is converted from microL/min/mg to L/h/mg
DDI.CLint = repmat(ParGen(10, :), 1, FacDDI) .*10^-6 .* 60;
DDI.CLrenal = repmat(ParGen(11, :), 1, FacDDI);
DDI.CLrenalCV = repmat(ParGen(12, :), 1, FacDDI);
DDI.CLbile = repmat(ParGen(13, :), 1, FacDDI);
DDI.CLbileCV = repmat(ParGen(14, :), 1, FacDDI);
DDI.CLadd = repmat(ParGen(15, :), 1, FacDDI);
DDI.CLaddCV = repmat(ParGen(16, :), 1, FacDDI);
DDI.Ki_CYP = repmat(ParCYP(3*SYSTEM.CYPliNo+1:4*SYSTEM.CYPliNo, :), 1,FacDDI);
DDI.kinact_CYP = repmat(ParCYP(4*SYSTEM.CYPliNo+1:5*SYSTEM.CYPliNo, :), 1,FacDDI);
DDI.Kapp_CYP = repmat(ParCYP(5*SYSTEM.CYPliNo+1:6*SYSTEM.CYPliNo, :), 1,FacDDI);
DDI.IndMax_CYP = repmat(ParCYP(6*SYSTEM.CYPliNo+1:7*SYSTEM.CYPliNo, :), 1,FacDDI);
DDI.IC50_CYP = repmat(ParCYP(7*SYSTEM.CYPliNo+1:8*SYSTEM.CYPliNo, :), 1,FacDDI);
%if Ki / Kapp / IC50 is 0, it is converted to 1 to prevent dividing by 0 in the code
DDI.Ki_CYP(DDI.Ki_CYP == 0) = 1;
DDI.Kapp_CYP(DDI.Kapp_CYP == 0) = 1;
DDI.IC50_CYP(DDI.IC50_CYP == 0) = 1;
%%% name of drugs for outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DrugName = cell(DRUG.DrugNo, 1); 
DDIName {d} = DrugName{d};
%save the name of the drug into a new vector
for d = 1:DRUG.DrugNo
 DrugName{d} = DEF.name{5, DRUG.DrugName(d)};
end
%DDI predictions need a combined name
Help = ' + ';
for d = 1:DRUG.DrugNo
 switch d
 case 1
 DDIName{d} = DrugName{d};
 case 2
 DDIName{d} = DrugName{d};
 DDIName{d+1} = [DrugName{d-1}, Help, DrugName{d}];
 DDIName{d+2} = [DrugName{d}, Help, DrugName{d-1}];
 case 3
 DDIName{d} = DrugName{d};
 DDIName{d+1} = [DrugName{d-2}, Help, DrugName{d-1} , Help,DrugName{d}];
 DDIName{d+2} = [DrugName{d-1}, Help, DrugName{d-2} , Help,DrugName{d}];
 DDIName{d+3} = [DrugName{d}, Help, DrugName{d-2} ,Help, DrugName{d-1}];
 end
end
%save the name for DDI predictions globally
DDI.Name = DDIName;
DDI.No = DDINo;
%__drug absorption______________________________________________________________
%%% Effective permeability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in vitro - in vivo relationship comes from Sun et al., 2002
PeffMen = 10.^ (0.6795 .* log10(DDI.Papp) - 0.3355);
PeffMen = permute( repmat(PeffMen, STUDY.IndNo, 1, SYSTEM.SegNo), [1, 3, 2]);
%PeffAll for each of the segments will be based on the difference in length
%following the approach by Darwich et al., 2010
PeffMen = PeffMen .* (SYSTEM.LengthGU ./ sum(SYSTEM.LengthGU, 2));
%absorption clearance, which will be used in the model
DDI.CLab = repmat(SYSTEM.PSA, 1, 1, 1) .* PeffMen .* 3600 .* 0.001;
%__drug distribution____________________________________________________________
%%% Define plasma concentration of the main binding protein
%concentration of the plasma-binding protein in [g/L]
ProtConc = zeros(STUDY.IndNo, 1);
%reference for the calculation of age-dependency (reference is 30 years)
ProtRef = zeros(STUDY.IndNo, 1);
%partition coefficient for plasma-binding proteins into the tissue
KpPR = zeros(STUDY.IndNo, SYSTEM.CompNo, 1);
for d = 1:1
 if DDI.protein(d) == DEF.albumin
 ProtConc(:, d) = SYSTEM.Albumin;
 ProtRef(:, d) = 45.6;
 KpPR(:, :, d) = SYSTEM.KpHSA;

 elseif DDI.protein(d) == DEF.AAG
 ProtConc(:, d) = SYSTEM.AAG;
 ProtRef(:, d) = 0.798;
 KpPR(:, :, d) = SYSTEM.KpHSA;
 end
end
%%% pH-dependent parameters for distribution (Rodgers & Rowland) %%%%%%%%%%%%%%%
%this are the X,Y and Z parameter from Rodgers & Rowland
DDI.KRio = zeros(SYSTEM.CompNo, 1);
for com = 1:SYSTEM.CompNo
 for d = 1:1
 DDI.KRio(com, d) = 1 + 10.^(DDI.pka(d) - SYSTEM.pH(com));
 end
end
%%% vegetable oil:water partition coefficient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prediction according to Poulin & Theil (2002), because it is assumed that
%it predicts partitioning into the adipose tissue better
DDI.logD = 1.115.*(abs(DDI.logP)) - 1.35 - log10(DDI.KRio(DEF.plasma, :));
%%% fraction unbound in plasma - age-dependent changes %%%%%%%%%%%%%%%%%%%%%%%%%
DDI.fup = 1 ./ (1 + ((((1./(repmat(DDI.fu, STUDY.IndNo, 1))) - 1)./ ...
 ProtRef).*ProtConc));
%%% fraction unbound in the interstitial space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fuint is based on the semiphysiological Bmax model
%it is only calculated for albumin, because AAG doesn't distribute to Vint
DDI.fuine = ones(STUDY.IndNo, SYSTEM.CompNo, 1);
for d = 1:1
 if DDI.protein(d) == DEF.albumin
 DDI.fuine(:, :, d) = 1 ./ (((SYSTEM.KpHSA./SYSTEM.FraEW) .* ...
 ((1./DDI.fup(:,d)) - 1)) + 1);
 end
end
%%% fraction unbound in the cellular space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fucel is based on Rodgers & Rowland
%it is only calculated for albumin, because AAG does not need to be considered
DDI.fucel = ones(STUDY.IndNo, SYSTEM.CompNo, 1);
for d = 1:1
 if DDI.protein(d) == DEF.albumin
 DDI.fucel(:, :, d) = 1 ./ (1 + (((DDI.logP(d).*SYSTEM.FraNL + ...
 (0.3.*DDI.logP(d) + 0.7).*SYSTEM.FraNP)./ ...
 DDI.KRio(DEF.plasma, d)) +...
 SYSTEM.KpHSA.*ProtConc(:, d)));
 end
end
%%% tissue-partition coefficient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DDI.KaPR = zeros(STUDY.IndNo, 1);
DDI.Kpu = zeros(STUDY.IndNo, SYSTEM.CompNo, 1);
for d = 1:1

 for ind = 1:STUDY.IndNo
 DDI.KaPR(ind, d) = ((1 ./ DDI.fup(ind, d)) - 1 -...
 ((DDI.logP(d) .* SYSTEM.FraNL(ind, DEF.plasma) + ((0.3.*DDI.logP(d) + 0.7) .* SYSTEM.FraNP(ind, DEF.plasma))) ./ DDI.KRio(DEF.plasma, d))).*...
 (1 ./ ProtConc(ind, d));

 for com = 1:SYSTEM.CompNo
 DDI.Kpu(ind, com, d) = (((DDI.KRio(com, d) .* SYSTEM.FraIW(ind, com)) ./ DDI.KRio(DEF.plasma, d)) +...SYSTEM.FraEW(ind, com) +...
 ((DDI.logP(d) .* SYSTEM.FraNL(ind, com) + (0.3.* DDI.logP(d) + 0.7) .* SYSTEM.FraNP(ind, com)) ./ DDI.KRio(DEF.plasma, d)) + ...
 (DDI.KaPR(ind, d) .* KpPR(ind, com, d) *ProtConc(ind, d))) .*...
 DDI.KpScalarAll(d) .* DDI.KpScalar(com, d);
 end

 DRUG.Kpu(ind, DEF.adipose, d) = (((DDI.KRio(DEF.adipose, d) .* SYSTEM.FraIW(ind, DEF.adipose)) ./ DDI.KRio(DEF.plasma, d)) +...
 SYSTEM.FraEW(ind, DEF.adipose) + ...
 ((DDI.logD(d) .* SYSTEM.FraNL(ind, DEF.adipose) + (0.3 .* DDI.logD(d) + 0.7) .* SYSTEM.FraNP(ind, DEF.adipose)) ./ DDI.KRio(DEF.plasma, d)) +...
 (DDI.KaPR(ind, d) .* KpPR(ind, DEF.adipose, d) * ProtConc(ind, d))) .*...
 DDI.KpScalarAll(d) .* DDI.KpScalar(DEF.adipose, d);
 end
end
%%% flux through the membrane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DDI.Jout = repmat(SYSTEM.Qorg, 1, 1, 1);
%prepare variable for faster calculations
BP = permute (repmat (DDI.BP, STUDY.IndNo, 1, SYSTEM.OrgNo), [1, 3, 2]);
DDI.Jin = (abs((DDI.Kpu(:, 1:16, :) - (1./BP)) .* (DDI.fucel(:, 1:16, :) ./ ...
 DDI.fuine(:, 1:16, :))) .* DDI.Jout);
%%% Volume of distribution [L/kg] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DDI.Vss = zeros(STUDY.IndNo, 1);
for d = 1:1
 DDI.Vss(:, d) = ((SYSTEM.Vorg(:, DEF.plasma) ./ DDI.fup(:, d)) + ...
 (sum(SYSTEM.Vorg(:, 1:SYSTEM.OrgNo) .*...
 DDI.Kpu(:, 1:SYSTEM.OrgNo, d), 2)))./SYSTEM.Weight;
end
%__drug elimination_____________________________________________________________
%prepare empyt vectors
DDI.CLre = zeros(STUDY.IndNo, 1);
DDI.CLbi = zeros(STUDY.IndNo, 1);
DDI.CLad = zeros(STUDY.IndNo, 1);
for d = 1:1
 %renal clearance in [L/h]
 DDI.CLre(:, d) = normrnd(DDI.CLrenal(d),...
 ((DDI.CLrenalCV(d)/100) .* DDI.CLrenal(d)), ...
 STUDY.IndNo, 1);

 %biliary clearance in [L/h]
 DDI.CLbi(:, d) = normrnd(DDI.CLbile(d),...
 ((DDI.CLbileCV(d)./100) .* DDI.CLbile(d)), ...
 STUDY.IndNo, 1);

 %additional plasma clearance in [L/h]
 DDI.CLad(:, d) = normrnd(DDI.CLadd (d),...
 ((DDI.CLaddCV(d)./100) .* DDI.CLadd (d)), ...
 STUDY.IndNo, 1);
end

%link the renal clearance to the GFR of each individual
for d = 1
 DDI.CLre(:, d) = DDI.CLre(:, d) .*...
 (SYSTEM.GFR ./ (130 - 10.*SYSTEM.Sex)) .* ...
 (DDI.fup(:, d)./DDI.fu(d));
end
%ritonavir has an impact on the renal clearance of rivaroxaban, which is not
%mechanistically in the model, but will be considered
% global Index
% for dr = 1:DRUG.DrugNo
%  if DRUG.DrugName(dr) == DEF.rivaroxaban
%  if DRUG.DrugName(dr-1) == DEF.ritonavir
%  Index = find(DRUG.DrugName == DEF.rivaroxaban);
%  DDI.CLre(:, Index + DRUG.DrugNo) = 0.5 .* DDI.CLre(:, Index + DRUG.DrugNo);
%  end
%  end
% end
%__prepare DDI matrix___________________________________________________________
%%% Variables for DDI matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define three different drugs A-c
drugA = 1; drugB = 2; drugC = 3;
%six different DDIs are possible
ddi_1 = 1; %drug A alone
ddi_2 = 2; %drug B alone
ddi_3 = 3; %drug A with drug B / drug C alone
ddi_4 = 4; %drug B with drug A / drug A with drug B and C
ddi_5 = 5; %drug B with drug A and C
ddi_6 = 6; %drug C with drug A and B
%abbreviation for organs to define concentration for the DDI matrix
LI = DEF.liver; DU = SYSTEM.CompNo + DEF.duodenum;
JE = SYSTEM.CompNo + DEF.jejunum; IL = SYSTEM.CompNo + DEF.ileum;
%organs can be divided into three different subcompartments
cel = 2;
%intestinal segments can be divided into fluid, transit uptake and enterocytes
ent = 2;
%__DDI Matrix___________________________________________________________________
%%% Generate interaction parameters for all drugs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ki for CYP enzymes
Ki_CYP = permute( repmat( DDI.Ki_CYP(:, 1:DRUG.DrugNo), 1, 1, 1), ...
 [2, 3, 1]);
switch DRUG.DrugNo
 case 2
 %ddi_1 and ddi_2 are drug A and drug B alone
 Ki_CYP(drugB, ddi_1, :) = 1; Ki_CYP(drugA, ddi_2, :) = 1;

 case 3
 %ddi_1, ddi_2 and ddi_3 are drug A, drug B and drug C alone
 Ki_CYP(drugB, ddi_1, :) = 1; Ki_CYP(drugC, ddi_1, :) = 1;
 Ki_CYP(drugA, ddi_2, :) = 1; Ki_CYP(drugC, ddi_2, :) = 1;
 Ki_CYP(drugA, ddi_3, :) = 1; Ki_CYP(drugB, ddi_3, :) = 1;
end
%generate a factor to prevent dividing by 0 when Ki = 0
Fki_CYP = ones(DRUG.DrugNo, 1, SYSTEM.CYPliNo);
Fki_CYP(Ki_CYP == 1) = 0;
%MBI for CYP enzymes
kinact_CYP = permute( repmat( DDI.kinact_CYP(:, 1:DRUG.DrugNo), ...
 1, 1, 1), [2, 3, 1]);
Kapp_CYP = permute( repmat( DDI.Kapp_CYP(:, 1:DRUG.DrugNo),...
 1, 1, 1), [2, 3, 1]);
switch DRUG.DrugNo
 case 2
 kinact_CYP(drugB, ddi_1, :) = 0; kinact_CYP(drugA, ddi_2, :) = 0;

 Kapp_CYP(drugB, ddi_1, :) = 1; Kapp_CYP(drugA, ddi_2, :) = 1;
 case 3
 kinact_CYP(drugB, ddi_1, :) = 0; kinact_CYP(drugC, ddi_1, :) = 0;
 kinact_CYP(drugA, ddi_2, :) = 0; kinact_CYP(drugC, ddi_2, :) = 0;
 kinact_CYP(drugA, ddi_3, :) = 0; kinact_CYP(drugB, ddi_3, :) = 0;

 Kapp_CYP(drugB, ddi_1, :) = 1; Kapp_CYP(drugC, ddi_1, :) = 1;
 Kapp_CYP(drugA, ddi_2, :) = 1; Kapp_CYP(drugC, ddi_2, :) = 1;
 Kapp_CYP(drugA, ddi_3, :) = 1; Kapp_CYP(drugB, ddi_3, :) = 1;
end
% %Induction for CYP enzymes
IndMax_CYP = permute( repmat( DDI.IndMax_CYP(:, 1:DRUG.DrugNo), ...
 1, 1, 1), [2, 3, 1]);
IC50_CYP = permute( repmat( DDI.IC50_CYP(:, 1:DRUG.DrugNo),...
 1, 1, 1), [2, 3, 1]);
switch DRUG.DrugNo
 case 2
 IndMax_CYP(drugB, ddi_1, :) = 0; IndMax_CYP(drugA, ddi_2, :) = 0;

 IC50_CYP(drugB, ddi_1, :) = 1; IC50_CYP(drugA, ddi_2, :) = 1;
 case 3
 IndMax_CYP(drugB, ddi_1, :) = 0; IndMax_CYP(drugC, ddi_1, :) = 0;
 IndMax_CYP(drugA, ddi_2, :) = 0; IndMax_CYP(drugC, ddi_2, :) = 0;
 IndMax_CYP(drugA, ddi_3, :) = 0; IndMax_CYP(drugB, ddi_3, :) = 0;

 IC50_CYP(drugB, ddi_1, :) = 1; IC50_CYP(drugC, ddi_1, :) = 1;
 IC50_CYP(drugA, ddi_2, :) = 1; IC50_CYP(drugC, ddi_2, :) = 1;
 IC50_CYP(drugA, ddi_3, :) = 1; IC50_CYP(drugB, ddi_3, :) = 1;
end
%%% choose the correct concentration and fu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%multiplying factor for the concentration
Fco = ones(DRUG.DrugNo, 1);
switch DRUG.DrugNo
 case 2
 Fco(drugA, ddi_4) = 0;
 Fco(drugB, ddi_3) = 0;

 case 3
 Fco(drugA, ddi_5) = 0; Fco(drugA, ddi_6) = 0;
 Fco(drugB, ddi_4) = 0; Fco(drugB, ddi_6) = 0;
 Fco(drugC, ddi_4) = 0; Fco(drugC, ddi_5) = 0;
end
%additive factor for the concentration depending on the compartment
FcelLI = zeros(DRUG.DrugNo, 1);
FentDU = zeros(DRUG.DrugNo, 1);
FentJE = zeros(DRUG.DrugNo, 1);
FentIL = zeros(DRUG.DrugNo, 1);
%fraction unbound for DDI predictions
fuLIcel = permute (repmat (DDI.fucel(:, DEF.liver, 1:1), 1, 1,1),...
    [3, 2, 1]);
fuGUent = permute (repmat (DDI.fucel(:, DEF.gut, :), 1, 1, 1),...
    [2, 3, 1]);
switch DRUG.DrugNo
 case 2
 FcelLI(drugA, ddi_4) = LI + SYSTEM.SubNo(LI) + cel;
 FcelLI(drugB, ddi_3) = MODEL.ODENo + LI + SYSTEM.SubNo(LI) + cel;

 fuLIcel(drugA, ddi_2, :) = zeros(STUDY.IndNo, 1);
 fuLIcel(drugB, ddi_1, :) = zeros(STUDY.IndNo, 1);

 FentDU(drugA, ddi_4) = DU + SYSTEM.SubNo(DU) + ent;
 FentDU(drugB, ddi_3) = MODEL.ODENo + DU + SYSTEM.SubNo(DU) + ent;

 FentJE(drugA, ddi_4) = JE + SYSTEM.SubNo(JE) + ent;
 FentJE(drugB, ddi_3) = MODEL.ODENo + JE + SYSTEM.SubNo(JE) + ent;

 FentIL(drugA, ddi_4) = IL + SYSTEM.SubNo(IL) + ent;
 FentIL(drugB, ddi_3) = MODEL.ODENo + IL + SYSTEM.SubNo(IL) + ent;

 fuGUent(drugA, ddi_2, :) = zeros(STUDY.IndNo, 1);
 fuGUent(drugB, ddi_1, :) = zeros(STUDY.IndNo, 1);

 case 3
 FcelLI(drugA, [ddi_5, ddi_6]) = LI + SYSTEM.SubNo(LI) + cel;
 FcelLI(drugB, [ddi_4, ddi_6]) = MODEL.ODENo + LI + SYSTEM.SubNo(LI) + cel;
 FcelLI(drugC, [ddi_4, ddi_5]) = 2*MODEL.ODENo + LI + SYSTEM.SubNo(LI) +cel;

 fuLIcel(drugA, [ddi_2, ddi_3], :) = zeros(1, 2, STUDY.IndNo);
 fuLIcel(drugB, [ddi_1, ddi_3], :) = zeros(1, 2, STUDY.IndNo);
 fuLIcel(drugC, [ddi_1, ddi_2], :) = zeros(1, 2, STUDY.IndNo);

 FentDU(drugA, [ddi_5, ddi_6]) = DU + SYSTEM.SubNo(DU) + ent;
 FentDU(drugB, [ddi_4, ddi_6]) = MODEL.ODENo + DU + SYSTEM.SubNo(DU) + ent;
 FentDU(drugC, [ddi_4, ddi_5]) = 2*MODEL.ODENo + DU + SYSTEM.SubNo(DU) +ent;

 FentJE(drugA, [ddi_5, ddi_6]) = JE + SYSTEM.SubNo(JE) + ent;
 FentJE(drugB, [ddi_4, ddi_6]) = MODEL.ODENo + JE + SYSTEM.SubNo(JE) + ent;
 FentJE(drugC, [ddi_4, ddi_5]) = 2*MODEL.ODENo + JE + SYSTEM.SubNo(JE) +ent;

 FentIL(drugA, [ddi_5, ddi_6]) = IL + SYSTEM.SubNo(IL) + ent;
 FentIL(drugB, [ddi_4, ddi_6]) = MODEL.ODENo + IL + SYSTEM.SubNo(IL) + ent;
 FentIL(drugC, [ddi_4, ddi_5]) = 2*MODEL.ODENo + IL + SYSTEM.SubNo(IL) +ent;

 fuGUent(drugA, [ddi_2, ddi_3], :) = zeros(1, 2, STUDY.IndNo);
 fuGUent(drugB, [ddi_1, ddi_3], :) = zeros(1, 2, STUDY.IndNo);
 fuGUent(drugC, [ddi_1, ddi_2], :) = zeros(1, 2, STUDY.IndNo);
end
%%% Combine all data into one matrix to inform the ODE system %%%%%%%%%%%%%%%%%%
STUDY.DDIMat = {'Fco', 'FcelLI', 'FentDU', 'FentJE', 'FentIL', 'fuLIcel','fuGUent';...
 Fco, FcelLI, FentDU, FentJE, FentIL, fuLIcel, fuGUent};

STUDY.DDIMat_CYP = {'Fki', 'Ki', 'kinact', 'Kapp', 'IndMax', 'IC50';...
 Fki_CYP, Ki_CYP, kinact_CYP, Kapp_CYP, IndMax_CYP, IC50_CYP};
end
