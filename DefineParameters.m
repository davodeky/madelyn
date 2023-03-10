%Define parameters

function[] = DefineParameters()
%This function defines the population, model structure (PBPK compartments),
%the drugs, the virtual trial design and the simulation settings to be used
global DEF %global DEF defines model parameters
global SYSTEM %global SYSTEM defines sytem parameters
global DRUG %global DRUG defines drug parameters
global STUDY %global STUDY defines study design parameters
global MODEL %global MODEL defines parameters important for modeling
global OBS %global OBS saves observed parameters for the output
global NAME

%define parameters used in the script
pop = 1; com = 2; seg = 3; cyp = 4; dru = 5; pro = 6; adm = 7;

%try to fix the NAME issue
NAME.acetaminophen = 1;

%__Define Model parameters and the structure____________________________________
%%% Population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEF.AgingCaucasian = 1; DEF.name{pop, DEF.AgingCaucasian} = 'AgingCaucasian';
%%% Compartments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBC = red blood cells
DEF.lung = 1; DEF.name{com, DEF.lung} = 'lung';
DEF.adipose = 2; DEF.name{com, DEF.adipose} = 'adipose';
DEF.bone = 3; DEF.name{com, DEF.bone} = 'bone';
DEF.brain = 4; DEF.name{com, DEF.brain} = 'brain';
DEF.gonads = 5; DEF.name{com, DEF.gonads} = 'gonads';
DEF.heart = 6; DEF.name{com, DEF.heart} = 'heart';
DEF.kidney = 7; DEF.name{com, DEF.kidney} = 'kidney';
DEF.muscle = 8; DEF.name{com, DEF.muscle} = 'muscle';
DEF.skin = 9; DEF.name{com, DEF.skin} = 'skin';
DEF.thymus = 10; DEF.name{com, DEF.thymus} = 'thymus';
DEF.gut = 11; DEF.name{com, DEF.gut} = 'gut';
DEF.spleen = 12; DEF.name{com, DEF.spleen} = 'spleen';
DEF.pancreas = 13; DEF.name{com, DEF.pancreas} = 'pancreas';
DEF.liver = 14; DEF.name{com, DEF.liver} = 'liver';
DEF.lymphnode = 15; DEF.name{com, DEF.lymphnode} = 'lymphnode';
DEF.remaining = 16; DEF.name{com, DEF.remaining} = 'remaining';
DEF.plasma = 17; DEF.name{com, DEF.plasma} = 'plasma';
DEF.RBC = 18; DEF.name{com, DEF.RBC} = 'RBC';
Comp = [DEF.lung DEF.adipose DEF.bone DEF.brain DEF.gonads DEF.heart ...
 DEF.kidney DEF.muscle DEF.skin DEF.thymus DEF.gut DEF.spleen ...
 DEF.pancreas DEF.liver DEF.lymphnode DEF.remaining DEF.plasma ...
 DEF.RBC];
Org = [DEF.lung DEF.adipose DEF.bone DEF.brain DEF.gonads DEF.heart ...
 DEF.kidney DEF.muscle DEF.skin DEF.thymus DEF.gut DEF.spleen ...
 DEF.pancreas DEF.liver DEF.lymphnode DEF.remaining];
%how many compartments are used with and without the blood?
CompNo = length(Comp); SYSTEM.CompNo = CompNo;
OrgNo = length(Org); SYSTEM.OrgNo = OrgNo;
%%% intestinal segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEF.stomach = 1; DEF.name{seg, DEF.stomach} = 'stomach';
DEF.duodenum = 2; DEF.name{seg, DEF.duodenum} = 'duodenum';
DEF.jejunum = 3; DEF.name{seg, DEF.jejunum} = 'jejunum';
DEF.ileum = 4; DEF.name{seg, DEF.ileum} = 'ileum';
DEF.colon = 5; DEF.name{seg, DEF.colon} = 'colon';
DEF.faeces = 6; DEF.name{seg, DEF.faeces} = 'faeces';
Seg = [DEF.stomach DEF.duodenum DEF.jejunum DEF.ileum DEF.colon DEF.faeces];
%how many intestinal segments are modelled?
SegNo = length(Seg); SYSTEM.SegNo = SegNo;
%%% CYP enzymes (for dynamic abundance for MBI and induction) %%%%%%%%%%%%%%%%%%
DEF.CYP2D6 = 1; DEF.name{cyp, DEF.CYP2D6} = 'CYP2D6';
DEF.CYP3A4 = 2; DEF.name{cyp, DEF.CYP3A4} = 'CYP3A4';
DEF.CYP3A5 = 3; DEF.name{cyp, DEF.CYP3A5} = 'CYP3A5';
DEF.CYP2J2 = 4; DEF.name{cyp, DEF.CYP2J2} = 'CYP2J2';
CYPli = [DEF.CYP2D6 DEF.CYP3A4 DEF.CYP3A5 DEF.CYP2J2];
CYPin = [DEF.CYP2D6 DEF.CYP3A4 DEF.CYP3A5];
%how many CYP enzymes are in the liver (li) or intestine (in)?
CYPliNo = length(CYPli); SYSTEM.CYPliNo = CYPliNo;
CYPinNo = length(CYPin); SYSTEM.CYPinNo = CYPinNo;
%%% Subcompartments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the number of subcompartments are also the number of equations
Sub(DEF.lung) = 2; Sub(DEF.adipose) = 2; Sub(DEF.bone) = 2;
Sub(DEF.brain) = 2; Sub(DEF.gonads) = 2; Sub(DEF.heart) = 2;
Sub(DEF.kidney) = 2; Sub(DEF.muscle) = 2; Sub(DEF.skin) = 2;
Sub(DEF.thymus) = 2; Sub(DEF.gut) = 1; Sub(DEF.spleen) = 2;
Sub(DEF.pancreas) = 2; Sub(DEF.liver) = 2; Sub(DEF.lymphnode) = 2;
Sub(DEF.remaining) = 2; Sub(DEF.plasma) = 0; Sub(DEF.RBC) = 0;
Sub(CompNo + DEF.stomach) = 0; Sub(CompNo + DEF.duodenum) = 2;
Sub(CompNo + DEF.jejunum) = 2; Sub(CompNo + DEF.ileum) = 2;
Sub(CompNo + DEF.colon) = 2; Sub(CompNo + DEF.faeces) = 0;
Sub(CompNo + SegNo + DEF.CYP2D6) = 3; Sub(CompNo + SegNo + DEF.CYP3A4) = 3;
Sub(CompNo + SegNo + DEF.CYP3A5) = 3; Sub(CompNo + SegNo + DEF.CYP2J2) = 0;
%define the number of subcompartments
SubNo = zeros(1, CompNo + SegNo + CYPliNo);
for tot = 1:(CompNo + SegNo + CYPliNo)
 if tot == 0
 SubNo(tot) = 0;

 else
 SubNo(tot) = sum(Sub(1:tot-1));
 end
end
SYSTEM.SubNo = SubNo;
%number of ODEs to solve
ODENo = CompNo + SegNo + CYPliNo + sum(Sub); MODEL.ODENo = ODENo;
%%% Drugs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEF.acetaminophen = 1; DEF.name{dru, DEF.acetaminophen} = 'acetaminophen';
DEF.ritonavir = 2; DEF.name{dru, DEF.ritonavir} = 'ritonavir';
DEF.rivaroxaban = 3; DEF.name{dru, DEF.rivaroxaban} = 'rivaroxaban';
DrugLib = [DEF.acetaminophen DEF.ritonavir DEF.rivaroxaban];
%How many drug files are in the libraray?
DEF.DrugLibNo = length(DrugLib);
%%% Main binding protein %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AAG = alpha acidic glycoprotein
DEF.albumin = 1; DEF.name{pro, DEF.albumin} = 'albumin';
DEF.AAG = 2; DEF.name{pro, DEF.AAG} = 'AAG';
%%% Route of administration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEF.iv = 1; DEF.name{adm, DEF.iv} = 'iv';
DEF.oral = 2; DEF.name{adm, DEF.oral} = 'oral';
%===============================================================================
%__The user chooses the simulation settings_____________________________________
%%% Drug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%how many drugs should be simulated in parallel?
DRUG.DrugNo = 1;
%enter the name of each simulated drug
DRUG.DrugName = zeros(1, DRUG.DrugNo);
for d = 1:DRUG.DrugNo
 switch d
 case 1 %drug 1
 DRUG.DrugName(d) = DEF.acetaminophen;

 case 2 %drug 2
 DRUG.DrugName(d) = DEF.ritonavir;

 case 3 %drug 3
 DRUG.DrugName(d) = DEF.rivaroxaban;

 end
end
%%% Virtual study design %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STUDY.TrialNo = 10; %number of trials
STUDY.IndTrialNo = 12; %number of individuals per trial
%total number of subjects to be simulated
STUDY.IndNo = STUDY.TrialNo * STUDY.IndTrialNo;
STUDY.PropFem = 0; %proportion of women
STUDY.AgeMin = 20; %minimal age in the simulation in [years]
STUDY.AgeMax = 50; %maximum age in the simulation in [years]
STUDY.Resolution = 10; %resolution of each time unit
%prepare empty vectors for dosing regimen
AdminRoute = zeros(1, DRUG.DrugNo); %route of administration
NoDoses = zeros(1, DRUG.DrugNo); %number of doses
Dose = zeros(1, DRUG.DrugNo); %dose in [mg]
StartDose = zeros(1, DRUG.DrugNo); %when is the drug given? [h]
DoseIntervall = zeros(1, DRUG.DrugNo); %intervall between doses in [h]
LastTime = zeros(1, DRUG.DrugNo); %prolongation in [h]
for d = 1:DRUG.DrugNo
 switch d
 case 1
 AdminRoute(d) = DEF.oral;
 NoDoses(d) = 7;
 Dose(d) = 800;
 StartDose(d) = 0;
 DoseIntervall(d) = 24;
 LastTime(d) = 0;

 case 2
 AdminRoute(d) = DEF.oral;
 NoDoses(d) = 7;
 Dose(d) = 100;
 StartDose(d) = 0;
 DoseIntervall(d) = 24;
 LastTime(d) = 0;

 case 3
 AdminRoute(d) = DEF.oral;
 NoDoses(d) = 1;
 Dose(d) = 10;
 StartDose(d) = 144;
 DoseIntervall(d) = 24;
 LastTime(d) = 24;
 end
end
%%% Enter observed data to compare to the simulated outcome %%%%%%%%%%%%%%%%%%%%
OBS.Time_Drug1 = [];
OBS.Conc_Drug1 = [];
OBS.SD_Drug1 = [];

OBS.Time_Drug2 = [];
OBS.Conc_Drug2 = [];
OBS.SD_Drug2 = [];
OBS.Time_Drug3 = [];
OBS.Conc_Drug3 = [];
OBS.SD_Drug3 = [];
OBS.Time_DDI1 = [];
OBS.Conc_DDI1 = [];
OBS.SD_DDI1 = [];
OBS.Time_DDI2 = [];
OBS.Conc_DDI2 = [];
OBS.SD_DDI2 = [];
OBS.Time_DDI3 = [];
OBS.Conc_DDI3 = [];
OBS.SD_DDI3 = [];
%===============================================================================
%__set up a dose event matrix___________________________________________________
%define the name of columns for the dose event matrix
start = 1; ende = 2; dose = 3; admin = 4; res = 5; mmr = 6;
NoColDoseEvent = length([start, ende, dose, admin, res, mmr]);
Regimen = zeros( max(NoDoses), NoColDoseEvent, DRUG.DrugNo);
%%% Combine all dosing events in one matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define the start and end of each dosing event for each drug
%dose and route of administration are assigned to each dosing event
for d = 1:DRUG.DrugNo
 for m = 1:NoDoses(d)
 %a case is defined for single doses
 if m == 1
 Regimen(1, start, d) = 0 + StartDose(d);
 Regimen(1, ende, d) = Regimen(1, start, d) + DoseIntervall(d) + ...
 LastTime(d);

 %a case is defined for multiple doses
 else
 Regimen(1, start, d) = 0 + StartDose(d);
 Regimen(1, ende, d) = Regimen(1, start, d) + DoseIntervall(d);

 Regimen(m, start, d) = Regimen(m-1, ende, d);
 Regimen(m, ende, d) = Regimen(m, start, d) + DoseIntervall(d);
 %a case is defined for the last dose to prolong the elimination phase
 if m == NoDoses(d)
 Regimen(m, ende, d) = Regimen(m, start, d) +...
 DoseIntervall(d) +...
 LastTime(d);
 end
 end

 Regimen(m, dose, d) = Dose(d);
 Regimen(m, admin, d) = AdminRoute(d);
 Regimen(m, res, d) = (Regimen(m, ende, d) - Regimen(m, start, d)) .* ...
 STUDY.Resolution;

 if m == 1
 Regimen(1, mmr, d) = 0;

 else
 Regimen(m, mmr, d) = Regimen(m-1, mmr, d) +...
 (Regimen(m-1, ende, d) - Regimen(m-1, start, d)) .*...
 STUDY.Resolution;
 end
 end
end
%%% find unique dosing events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract dosing events for each drug and delete zero values
for d = 1:DRUG.DrugNo
 switch d
 case 1
 Drug1Regimen = sortrows(Regimen(:, :, d), start);
 Drug1Regimen( ~any(Drug1Regimen, 2), :) = [];

 case 2
 Drug2Regimen = sortrows(Regimen(:, :, d), start);
 Drug2Regimen( ~any(Drug2Regimen, 2), :) = [];

 case 3
 Drug3Regimen = sortrows(Regimen(:, :, d), start);
 Drug3Regimen( ~any(Drug3Regimen, 2), :) = [];
 end
end
%combine dosing regimens of the different drugs and look for unique events
switch DRUG.DrugNo
 case 1
 Drug1RegMat = Drug1Regimen;
 Drug1RegFun = RegimenFun(Drug1RegMat, NoDoses(1));

 case 2
 Drug1RegMat = [Drug1Regimen; Drug2Regimen];
 Drug1RegFun = RegimenFun(Drug1RegMat, NoDoses(1));

 Drug2RegMat = [Drug2Regimen; Drug1Regimen];
 Drug2RegFun = RegimenFun(Drug2RegMat, NoDoses(2));

 case 3
 Drug1RegMat = [Drug1Regimen; Drug2Regimen; Drug3Regimen];
 Drug1RegFun = RegimenFun(Drug1RegMat, NoDoses(1));

 Drug2RegMat = [Drug2Regimen; Drug1Regimen; Drug3Regimen];
 Drug2RegFun = RegimenFun(Drug2RegMat, NoDoses(2));

 Drug3RegMat = [Drug3Regimen; Drug1Regimen; Drug2Regimen];
 Drug3RegFun = RegimenFun(Drug3RegMat, NoDoses(3));
end
%prolongation of the terminal time needs to be considered for all drugs
switch DRUG.DrugNo
 case 1
 StartTime = Drug1RegFun(:, start);
 EndTime = Drug1RegFun(:, ende);
 case 2
 StartTime = [Drug1RegFun(:, start); Drug2RegFun(:, start)];
 EndTime = [Drug1RegFun(:, ende); Drug2RegFun(:, ende)];

 case 3
 StartTime = [Drug1RegFun(:, start); Drug2RegFun(:, start); Drug3RegFun(:,start)];
 EndTime = [Drug1RegFun(:, ende); Drug2RegFun(:, ende); Drug3RegFun(:,ende)];
end
%find unique start and end time and calculate the resolution
UniqueStartT = unique(StartTime);
UniqueEndT = unique(EndTime);
for d = 1:DRUG.DrugNo
 if NoDoses(d) ~= 0
 if LastTime(d) ~= 0
 UniqueEndT(end-1) = UniqueEndT(end);
 UniqueEndT(end) = [];
 end
 end
end
ResolutionEvent = (UniqueEndT - UniqueStartT) .* STUDY.Resolution;
ResolutionAll = zeros(length(UniqueStartT), 1);
for a = 1:length(UniqueStartT)
 if a == 1
 ResolutionAll(a) = 0;

 else
 ResolutionAll(a) = ResolutionAll(a-1) +...
 (UniqueEndT(a-1) - UniqueStartT(a-1)) .* STUDY.Resolution;
 end
end
%each drug becomes the same start / end time and resolution
switch DRUG.DrugNo
 case 1
 Drug1RegFun(:, start) = UniqueStartT;
 Drug1RegFun(:, ende) = UniqueEndT;
 Drug1RegFun(:, res) = ResolutionEvent;
 Drug1RegFun(:, mmr) = ResolutionAll;

 case 2
 Drug1RegFun(:, start) = UniqueStartT;
 Drug1RegFun(:, ende) = UniqueEndT;
 Drug1RegFun(:, res) = ResolutionEvent;
 Drug1RegFun(:, mmr) = ResolutionAll;

 Drug2RegFun(:, start) = UniqueStartT;
 Drug2RegFun(:, ende) = UniqueEndT;
 Drug2RegFun(:, res) = ResolutionEvent;
 Drug2RegFun(:, mmr) = ResolutionAll;

 case 3
 Drug1RegFun(:, start) = UniqueStartT;
 Drug1RegFun(:, ende) = UniqueEndT;
 Drug1RegFun(:, res) = ResolutionEvent;
 Drug1RegFun(:, mmr) = ResolutionAll;

 Drug2RegFun(:, start) = UniqueStartT;
 Drug2RegFun(:, ende) = UniqueEndT;
 Drug2RegFun(:, res) = ResolutionEvent;
 Drug2RegFun(:, mmr) = ResolutionAll;

 Drug3RegFun(:, start) = UniqueStartT;
 Drug3RegFun(:, ende) = UniqueEndT;
 Drug3RegFun(:, res) = ResolutionEvent;
 Drug3RegFun(:, mmr) = ResolutionAll;
end
%what is the number of total events in the simulation?
NoEvents = length(Drug1RegFun(:, start)); STUDY.NoEvents = NoEvents;
%comine the dosing regimen for all drugs
DoseEventMat = zeros(NoEvents, NoColDoseEvent, DRUG.DrugNo);
switch DRUG.DrugNo
 case 1
 DoseEventMat(:, :, 1) = Drug1RegFun;

 case 2
 DoseEventMat(:, :, 1) = Drug1RegFun;
 DoseEventMat(:, :, 2) = Drug2RegFun;

 case 3
 DoseEventMat(:, :, 1) = Drug1RegFun;
 DoseEventMat(:, :, 2) = Drug2RegFun;
 DoseEventMat(:, :, 3) = Drug3RegFun;
end
%save the dose event matrix globally
STUDY.DoseEventMat = DoseEventMat;
STUDY.Dose = Dose;
STUDY.NoDoses = NoDoses;
%===============================================================================
%%% USED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
function DrugReg = RegimenFun(DrugMat, NoDoses)
 %dosing events from other drugs are set to zero
 DrugMat(NoDoses + 1:end, [dose, admin, res, mmr]) = 0;

 %drug events are sorted based on the start time
 DrugSort = sortrows(DrugMat, start);

 %find unique time points
 [~, idx] = unique(DrugSort(:, start));

 %ouput results of the used function
 DrugReg = DrugSort(idx, :);
end
end