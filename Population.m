
function[] = Population()
%This function generates the virtual population based on the user settings
%Attention: the normrnd command needs to Statistics and Machine Learning Toolbox
global DEF %global DEF defines model parameters
global SYSTEM %global SYSTEM defines sytem parameters
global STUDY %global STUDY defines study design parameters
%__Age distribution_____________________________________________________________
%Weibull distribution was fitted to data from Eurostat (Stader et al., 2018)
Age = round(61.73.*(-log(1 - rand(STUDY.IndNo, 1))).^(1/1.55));
%Weibull distribution are infinity and need to be truncated
for ind = 1:STUDY.IndNo
 while Age(ind) < STUDY.AgeMin || Age(ind) > STUDY.AgeMax
 Age(ind) = round(61.73.*(-log(1 - rand(1, 1))).^(1/1.55));
 end
end
%save age globally
SYSTEM.Age = Age;
%__Sex distribution_____________________________________________________________
%number of women in the simulation - round ensure the number to be an integer
FemNo = round(STUDY.PropFem .* STUDY.IndNo);
%generate a matrix with random numbers
IndexSex = randperm(STUDY.IndNo);
%assign randomly a 1 to women and a 0 to men
Sex = zeros(STUDY.IndNo, 1);
Sex(IndexSex(1 : FemNo)) = 1;
%save sex globally
SYSTEM.Sex = Sex;
%__Anthropometric parameters____________________________________________________
%Equations are published in Stader et al., 2018
%%% Body height in [cm] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Height_Mean = -0.0039.*Age.^2 + 0.238.*Age - 12.5.*Sex + 176;
Height_CV = 3.8;
Height_Min = [150, 140];
Height_Max = [200, 180];
SYSTEM.Height = Calc_SysPar(Height_Mean, Height_CV, Height_Min, Height_Max);
%%% Body weight in [kg] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Weight_Mean = -0.0039.*Age.^2 + 1.12.*SYSTEM.Height + 0.611.*Age - ...
 0.424.*Sex - 137;
Weight_CV = 15.2;
Weight_Min = [50, 40];
Weight_Max = [110, 90];
Weight_NoVa = -0.0039.*Age.^2 + 1.12.*Height_Mean + 0.611.*Age - 0.424.*Sex - 137;
SYSTEM.Weight = Calc_SysPar(Weight_Mean, Weight_CV, Weight_Min, Weight_Max);
%%% Body Surface Area (BSA) according to DuBois & DuBois %%%%%%%%%%%%%%%%%%%%%%%
SYSTEM.BSA = 0.007184 .* SYSTEM.Height.^0.725 .* SYSTEM.Weight.^0.425;
BSA_NoVa = 0.007184 .* Height_Mean.^0.725 .* Weight_NoVa.^0.425;
%__CYP enzymes phenotypes_______________________________________________________
%define variables
em = 1; %EM = enhanced metaboliseres (2 allels) - assign a random 1
pm = 2; %PM = poor metaboliseres (0 allels) - assign a random 0
um = 3; %UM = ultrarapid metaboliseres (4 allels) - assign a random 2
%we set up the frequencies for CYP enzymes
FreqCYP = zeros(SYSTEM.CYPliNo, length([em, pm, um]));
%enter the frequencies when available
FreqCYP(DEF.CYP2D6, pm) = 0.08;
FreqCYP(DEF.CYP2D6, um) = 0.05;
FreqCYP(DEF.CYP3A5, pm) = 0.83;
%how many subjects per phenotype?
PhenNoCYP = round(FreqCYP .* STUDY.IndNo);
PhenNoCYP(:, em) = STUDY.IndNo - PhenNoCYP(:, pm) - PhenNoCYP(:, um);
%prepare vector for loop
PheCYP = ones(STUDY.IndNo, SYSTEM.CYPliNo);
CYP_EM = num2cell(zeros(SYSTEM.CYPliNo), 1);
CYP_PM = num2cell(zeros(SYSTEM.CYPliNo), 1);
CYP_UM = num2cell(zeros(SYSTEM.CYPliNo), 1);
CYP_All = num2cell(zeros(SYSTEM.CYPliNo), 1);
for cyp = 1:SYSTEM.CYPliNo
 CYP_EM{cyp} = ones(PhenNoCYP(cyp, em), 1); %EMs get a 1
 CYP_PM{cyp} = zeros(PhenNoCYP(cyp, pm), 1); %PMs get a 0
 CYP_UM{cyp} = 2*ones(PhenNoCYP(cyp, um), 1); %UMs get a 3

 %combine all phenotypes
 CYP_All{cyp} = [CYP_EM{cyp}; CYP_PM{cyp}; CYP_UM{cyp}];

 %randomly assign 0 (PM) / 1 (EM) / 2 (UM) based on frequencies
 PheCYP(:, cyp) = CYP_All {cyp} (randperm(length(CYP_All {cyp} (:))));
end
%__Blood parameters_____________________________________________________________
%all parameters are from Stader et al., 2018
%%% haematocrit (HCT) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HCT_Mean = 0.443 - 0.033.*Sex;
HCT_CV = 14.4;
HCT_Min = [0.3, 0.3];
HCT_Max = [0.5, 0.5];
SYSTEM.HCT = Calc_SysPar(HCT_Mean, HCT_CV, HCT_Min, HCT_Max);
%%% albumin (HSA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Albumin_Mean = - 0.0709.*Age + 47.7;
Albumin_CV = 7.9;
Albumin_Min = [35.8, 35.8];
Albumin_Max = [50.2, 50.2];
SYSTEM.Albumin = Calc_SysPar(Albumin_Mean, Albumin_CV, Albumin_Min, Albumin_Max);
%%% alpha-acidic glycoprotein (AAG) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AAG_Mean = 0.798*ones(STUDY.IndNo, 1);
AAG_CV = 24.3;
AAG_Min = [0.476, 0.476];
AAG_Max = [1.22, 1.22];
SYSTEM.AAG = Calc_SysPar(AAG_Mean, AAG_CV, AAG_Min, AAG_Max);
%__Organ weights (Worg) in [kg]_________________________________________________
%all parameters are from Stader et al., 2018
SYSTEM.Worg = zeros(STUDY.IndNo, SYSTEM.CompNo);
%%% Lungs (LU) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WLU_Mean = exp(0.00771.*Age + 0.0279.*SYSTEM.Height - 5.58);
WLU_CV = 0;
WLU_Min = [0.00453.*SYSTEM.Weight, 0.00453.*SYSTEM.Weight,];
WLU_Max = [0.0122.*SYSTEM.Weight, 0.0122.*SYSTEM.Weight,];
WLU_NoVa = exp(0.00771.*Age + 0.0279.*Height_Mean - 5.58);
SYSTEM.Worg(:, DEF.lung) = Calc_SysPar(WLU_Mean, WLU_CV, WLU_Min, WLU_Max);
%%% Adipose (AD) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WAD_Mean = abs(0.68.*SYSTEM.Weight - 0.56.*SYSTEM.Height + 6.1.*Sex + 65);
WAD_CV = 29.6;
WAD_Min = [0.10.*SYSTEM.Weight, 0.21.*SYSTEM.Weight];
WAD_Max = [0.42.*SYSTEM.Weight, 0.59.*SYSTEM.Weight];
WAD_NoVa = abs(0.68.*Weight_NoVa - 0.56.*Height_Mean + 6.1.*Sex + 65);
SYSTEM.Worg(:, DEF.adipose) = Calc_SysPar(WAD_Mean, WAD_CV, WAD_Min, WAD_Max);
%%% Bone (BO) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WBO_Mean = exp(0.024.*SYSTEM.Height - 1.9);
WBO_CV = 13.2;
WBO_Min = [0.079.*SYSTEM.Weight, 0.079.*SYSTEM.Weight];
WBO_Max = [0.13.*SYSTEM.Weight, 0.13.*SYSTEM.Weight];
WBO_NoVa = exp(0.024.*Height_Mean - 1.9);
SYSTEM.Worg(:, DEF.bone) = Calc_SysPar(WBO_Mean, WBO_CV, WBO_Min, WBO_Max);
%%% Brain (BR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WBR_Mean = exp(-0.00075.*Age + 0.00778.*SYSTEM.Height - 0.97);
WBR_CV = 9.0;
WBR_Min = [0.0168.*SYSTEM.Weight, 0.0168.*SYSTEM.Weight];
WBR_Max = [0.0295.*SYSTEM.Weight, 0.0295.*SYSTEM.Weight];
WBR_NoVa = exp(-0.00075.*Age + 0.00778.*Height_Mean - 0.97);
SYSTEM.Worg(:, DEF.brain) = Calc_SysPar(WBR_Mean, WBR_CV, WBR_Min, WBR_Max);
%%% Gonads (GO) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WGO_Mean = -0.00022.*Age - 0.00034.*SYSTEM.Weight - 0.030.*Sex + 0.082;
WGO_CV = 34.8;
WGO_Min = [0.00032.*SYSTEM.Weight, 0.00008.*SYSTEM.Weight];
WGO_Max = [0.00066.*SYSTEM.Weight, 0.00017.*SYSTEM.Weight];
WGO_NoVa = -0.00022.*Age - 0.00034.*Weight_NoVa - 0.030.*Sex + 0.082;
SYSTEM.Worg(:, DEF.gonads) = Calc_SysPar(WGO_Mean, WGO_CV, WGO_Min, WGO_Max);
%%% Heart (HE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WHE_Mean = 0.34.*SYSTEM.BSA + 0.0018.*Age - 0.36;
WHE_CV = 19.9;
WHE_Min = [0.0018.*SYSTEM.Weight, 0.0018.*SYSTEM.Weight];
WHE_Max = [0.0076.*SYSTEM.Weight, 0.0076.*SYSTEM.Weight];
WHE_NoVa = 0.34.*BSA_NoVa + 0.0018.*Age - 0.36;
SYSTEM.Worg(:, DEF.heart) = Calc_SysPar(WHE_Mean, WHE_CV, WHE_Min, WHE_Max);
%%% Kidney (KI) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WKI_Mean = -0.00038.*Age - 0.056.*Sex + 0.33;
WKI_CV = 21.2;
WKI_Min = [0, 0];
WKI_Max = [1, 1];
WKI_NoVa = -0.00038.*Age - 0.056.*Sex + 0.33;
SYSTEM.Worg(:, DEF.kidney) = Calc_SysPar(WKI_Mean, WKI_CV, WKI_Min, WKI_Max);
%%% Muscle (MU) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WMU_Mean = 17.9.*SYSTEM.BSA - 0.0667.*Age - 5.68.*Sex - 1.22;
WMU_CV = 11.8;
WMU_Min = [0.310.*SYSTEM.Weight, 0.199.*SYSTEM.Weight];
WMU_Max = [0.459.*SYSTEM.Weight, 0.388.*SYSTEM.Weight];
WMU_NoVa = 17.9.*BSA_NoVa - 0.0667.*Age - 5.68.*Sex - 1.22;
SYSTEM.Worg(:, DEF.muscle) = Calc_SysPar(WMU_Mean, WMU_CV, WMU_Min, WMU_Max);
%%% Skin (SK) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WSK_Mean = exp(-0.0058.*Age - 0.37.*Sex + 1.13);
WSK_CV = 8.3;
WSK_Min = [0.0094.*SYSTEM.Weight, 0.0054.*SYSTEM.Weight];
WSK_Max = [0.0479.*SYSTEM.Weight, 0.0411.*SYSTEM.Weight];
WSK_NoVa = exp(-0.0058.*Age - 0.37.*Sex + 1.13);
SYSTEM.Worg(:, DEF.skin) = Calc_SysPar(WSK_Mean, WSK_CV, WSK_Min, WSK_Max);
%%% Thymus (TH) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WTH_Mean = 0.0221*ones(STUDY.IndNo, 1);
WTH_CV = 44.8;
WTH_Min = [0.00016.*SYSTEM.Weight, 0.00016.*SYSTEM.Weight];
WTH_Max = [0.00046.*SYSTEM.Weight, 0.00046.*SYSTEM.Weight];
WTH_NoVa = 0.0221*ones(STUDY.IndNo, 1);
SYSTEM.Worg(:, DEF.thymus) = Calc_SysPar(WTH_Mean, WTH_CV, WTH_Min, WTH_Max);
%%% Gut (GU) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WGU_Mean = 3E-06.*SYSTEM.Height.^2.49;
WGU_CV = 7.3;
WGU_Min = [0, 0];
WGU_Max = [2, 2];
WGU_NoVa = 3E-06.*Height_Mean.^2.49;
SYSTEM.Worg(:, DEF.gut) = Calc_SysPar(WGU_Mean, WGU_CV, WGU_Min, WGU_Max);
%%% Spleen (SP) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WSP_Mean = exp(1.13.*SYSTEM.BSA - 3.93);
WSP_CV = 51.7;
WSP_Min = [0.00098.*SYSTEM.Weight, 0.00098.*SYSTEM.Weight];
WSP_Max = [0.00321.*SYSTEM.Weight, 0.00321.*SYSTEM.Weight];
WSP_NoVa = exp(1.13.*BSA_NoVa - 3.93);
SYSTEM.Worg(:, DEF.spleen) = Calc_SysPar(WSP_Mean, WSP_CV, WSP_Min, WSP_Max);
%%% Pancreas (PA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WPA_Mean = 0.103*ones(STUDY.IndNo, 1);
WPA_CV = 27.8;
WPA_Min = [0, 0];
WPA_Max = [1, 1];
WPA_NoVa = 0.103*ones(STUDY.IndNo, 1);
SYSTEM.Worg(:, DEF.pancreas) = Calc_SysPar(WPA_Mean, WPA_CV, WPA_Min, WPA_Max);
%%% Liver (LI) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WLI_Mean = exp(0.87.*SYSTEM.BSA - 0.014.*Age - 1.0);
WLI_CV = 23.7;
WLI_Min = [0.0147.*SYSTEM.Weight, 0.0147.*SYSTEM.Weight];
WLI_Max = [0.0332.*SYSTEM.Weight, 0.0332.*SYSTEM.Weight];
WLI_NoVa = exp(0.87.*BSA_NoVa - 0.014.*Age - 1.0);
SYSTEM.Worg(:, DEF.liver) = Calc_SysPar(WLI_Mean, WLI_CV, WLI_Min, WLI_Max);
%%% Lymphnode (LN) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data are from Gill et al., 2016
WLN_Mean = 0.00386.*SYSTEM.Weight;
WLN_CV = 30;
WLN_Min = [0, 0];
WLN_Max = [1, 1];
WLN_NoVa = 0.00386.*Weight_NoVa;
SYSTEM.Worg(:, DEF.lymphnode) = Calc_SysPar(WLN_Mean, WLN_CV, WLN_Min, WLN_Max);
%%% Blood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%weight of the entire blood
WBL = -0.0098.*Age - 1.89.*Sex + 6.06;
WorgBlood = normrnd(WBL, (10.4/100).*WBL);
%blood weight is split into plasma and red blood cells
SYSTEM.Worg(:, DEF.plasma) = WorgBlood.*(1 - 0.91.*SYSTEM.HCT);
SYSTEM.Worg(:, DEF.RBC) = WorgBlood - SYSTEM.Worg(:, DEF.plasma);
%weight of venous (VB) and arterial blood (AB)
SYSTEM.Wvein = (2/3).*WorgBlood;
SYSTEM.Wartery = (1/3).*WorgBlood;
%%% Remaining (RE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SYSTEM.Worg(:, DEF.remaining) = SYSTEM.Weight - ...
 SYSTEM.Worg(:, DEF.lung) -...
 SYSTEM.Worg(:, DEF.adipose) -...
 SYSTEM.Worg(:, DEF.bone) -...
 SYSTEM.Worg(:, DEF.brain) -...
 SYSTEM.Worg(:, DEF.gonads) -...
 SYSTEM.Worg(:, DEF.heart) -...
 SYSTEM.Worg(:, DEF.kidney) -...
 SYSTEM.Worg(:, DEF.muscle) -...
 SYSTEM.Worg(:, DEF.skin) -...
 SYSTEM.Worg(:, DEF.thymus) -...
 SYSTEM.Worg(:, DEF.gut) -...
 SYSTEM.Worg(:, DEF.spleen) -...
 SYSTEM.Worg(:, DEF.pancreas) -...
 SYSTEM.Worg(:, DEF.liver) -...
 SYSTEM.Worg(:, DEF.lymphnode) -...
 SYSTEM.Worg(:, DEF.plasma) -...
 SYSTEM.Worg(:, DEF.RBC);
%because of the random variability, the remaining organ weight can be negative
%in these cases, the additional random variability is removed
for ind = 1:STUDY.IndNo
 if SYSTEM.Worg(ind, DEF.remaining) < 0
 SYSTEM.Height(ind) = Height_Mean(ind);
 SYSTEM.Weight(ind) = Weight_NoVa(ind);
 SYSTEM.BSA(ind) = BSA_NoVa(ind);

 SYSTEM.Worg(ind, DEF.lung) = WLU_NoVa(ind);
 SYSTEM.Worg(ind, DEF.adipose) = WAD_NoVa(ind);
 SYSTEM.Worg(ind, DEF.bone) = WBO_NoVa(ind);
 SYSTEM.Worg(ind, DEF.brain) = WBR_NoVa(ind);
 SYSTEM.Worg(ind, DEF.gonads) = WGO_NoVa(ind);
 SYSTEM.Worg(ind, DEF.heart) = WHE_NoVa(ind);
 SYSTEM.Worg(ind, DEF.kidney) = WKI_NoVa(ind);
 SYSTEM.Worg(ind, DEF.muscle) = WMU_NoVa(ind);
 SYSTEM.Worg(ind, DEF.skin) = WSK_NoVa(ind);
 SYSTEM.Worg(ind, DEF.thymus) = WTH_NoVa(ind);
 SYSTEM.Worg(ind, DEF.gut) = WGU_NoVa(ind);
 SYSTEM.Worg(ind, DEF.spleen) = WSP_NoVa(ind);
 SYSTEM.Worg(ind, DEF.pancreas) = WPA_NoVa(ind);
 SYSTEM.Worg(ind, DEF.liver) = WLI_NoVa(ind);
 SYSTEM.Worg(ind, DEF.lymphnode) = WLN_NoVa(ind);
 SYSTEM.Worg(ind, DEF.plasma) = WBL(ind) * (1 - 0.91.*HCT_Mean(ind));
 SYSTEM.Worg(ind, DEF.RBC) = WBL(ind) - SYSTEM.Worg(ind, DEF.plasma);
 end
end
SYSTEM.Worg(:, DEF.remaining) = SYSTEM.Weight - ...
 SYSTEM.Worg(:, DEF.lung) -...
 SYSTEM.Worg(:, DEF.adipose) -...
 SYSTEM.Worg(:, DEF.bone) -...
 SYSTEM.Worg(:, DEF.brain) -...
 SYSTEM.Worg(:, DEF.gonads) -...
 SYSTEM.Worg(:, DEF.heart) -...
 SYSTEM.Worg(:, DEF.kidney) -...
 SYSTEM.Worg(:, DEF.muscle) -...
 SYSTEM.Worg(:, DEF.skin) -...
 SYSTEM.Worg(:, DEF.thymus) -...
 SYSTEM.Worg(:, DEF.gut) -...
 SYSTEM.Worg(:, DEF.spleen) -...
 SYSTEM.Worg(:, DEF.pancreas) -...
 SYSTEM.Worg(:, DEF.liver) -...
 SYSTEM.Worg(:, DEF.lymphnode) -...
 SYSTEM.Worg(:, DEF.plasma) -...
 SYSTEM.Worg(:, DEF.RBC);
%__Organ density (OrgDen) [kg/L]________________________________________________
%Organ densities are from the ICRP 1975 (Snyder) and 2002 (Valentin)
SYSTEM.OrgDen = ones(1, SYSTEM.CompNo);
SYSTEM.OrgDen(DEF.lung) = 1.0; %free of blood and air
SYSTEM.OrgDen(DEF.adipose) = 0.916;
SYSTEM.OrgDen(DEF.bone) = 1.9;
SYSTEM.OrgDen(DEF.brain) = 1.04;
SYSTEM.OrgDen(DEF.gonads) = 1.045; %combined from males and females
SYSTEM.OrgDen(DEF.heart) = 1.03;
SYSTEM.OrgDen(DEF.kidney) = 1.05;
SYSTEM.OrgDen(DEF.muscle) = 1.041;
SYSTEM.OrgDen(DEF.skin) = 1.1;
SYSTEM.OrgDen(DEF.thymus) = 1.025;
SYSTEM.OrgDen(DEF.gut) = 1.042;
SYSTEM.OrgDen(DEF.spleen) = 1.06;
SYSTEM.OrgDen(DEF.pancreas) = 1.045;
SYSTEM.OrgDen(DEF.liver) = 1.08; %Heineman et al. (1999)
SYSTEM.OrgDen(DEF.lymphnode) = 1.0; %no data available; assume 1.0
SYSTEM.OrgDen(DEF.plasma) = 1.027;
SYSTEM.OrgDen(DEF.RBC) = 1.09;
%use the same organ density for each subject
SYSTEM.OrgDen = repmat(SYSTEM.OrgDen, STUDY.IndNo, 1);
%use the weighted mean of all used tissues for the remaining organ
SYSTEM.OrgDen(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.OrgDen, 2)./ ...
 sum(SYSTEM.Worg, 2);
%__Organ volume (Vorg) [L]______________________________________________________
SYSTEM.Vorg = SYSTEM.Worg./SYSTEM.OrgDen;
SYSTEM.Vvein = SYSTEM.Wvein./1.06;
SYSTEM.Vartery = SYSTEM.Wartery./1.06;
%__Tissue composition___________________________________________________________
%Values published by Gill et al., 2016 and Jamei et al., 2009 are used
%thymus data are from rat (Rodgers & Rowland, 2005)
%gonad data are from Pierson et al., 1978, Bieri & Privali, 1965 and Diagne et al.,1983
%lymphnode data are from Zhu et al., 1996
SYSTEM.FraEW = ones(1, SYSTEM.CompNo); %fraction of extracellular water
SYSTEM.FraIW = ones(1, SYSTEM.CompNo); %fraction of intracellular water
SYSTEM.FraNL = ones(1, SYSTEM.CompNo); %fraction of neutral lipids
SYSTEM.FraNP = ones(1, SYSTEM.CompNo); %fraction of phospholipids
SYSTEM.AP = ones(1, SYSTEM.CompNo); %acidic phospholipids in [mg/g]
SYSTEM.KpHSA = ones(1, SYSTEM.CompNo); %partition coefficient of albumin
SYSTEM.FraEW(DEF.lung) = 0.348; SYSTEM.FraIW(DEF.lung) = 0.463;
SYSTEM.FraEW(DEF.adipose) = 0.141; SYSTEM.FraIW(DEF.adipose) = 0.039;
SYSTEM.FraEW(DEF.bone) = 0.098; SYSTEM.FraIW(DEF.bone) = 0.341;
SYSTEM.FraEW(DEF.brain) = 0.092; SYSTEM.FraIW(DEF.brain) = 0.678;
SYSTEM.FraEW(DEF.gonads) = 0.239; SYSTEM.FraIW(DEF.gonads) = 0.561;
SYSTEM.FraEW(DEF.heart) = 0.313; SYSTEM.FraIW(DEF.heart) = 0.445;
SYSTEM.FraEW(DEF.kidney) = 0.283; SYSTEM.FraIW(DEF.kidney) = 0.5;
SYSTEM.FraEW(DEF.muscle) = 0.091; SYSTEM.FraIW(DEF.muscle) = 0.669;
SYSTEM.FraEW(DEF.skin) = 0.623; SYSTEM.FraIW(DEF.skin) = 0.0947;
SYSTEM.FraEW(DEF.thymus) = 0.150; SYSTEM.FraIW(DEF.thymus) = 0.626;
SYSTEM.FraEW(DEF.gut) = 0.267; SYSTEM.FraIW(DEF.gut) = 0.451;
SYSTEM.FraEW(DEF.spleen) = 0.208; SYSTEM.FraIW(DEF.spleen) = 0.58;
SYSTEM.FraEW(DEF.pancreas) = 0.12; SYSTEM.FraIW(DEF.pancreas) = 0.664;
SYSTEM.FraEW(DEF.liver) = 0.165; SYSTEM.FraIW(DEF.liver) = 0.586;
SYSTEM.FraEW(DEF.lymphnode) = 0.208; SYSTEM.FraIW(DEF.lymphnode) = 0.58;
SYSTEM.FraEW(DEF.plasma) = 0.945; SYSTEM.FraIW(DEF.plasma) = 0;
SYSTEM.FraEW(DEF.RBC) = 0; SYSTEM.FraIW(DEF.RBC) = 0.666;
SYSTEM.FraNL(DEF.lung) = 0.003; SYSTEM.FraNP(DEF.lung) = 0.009;
SYSTEM.FraNL(DEF.adipose) = 0.79; SYSTEM.FraNP(DEF.adipose) = 0.002;
SYSTEM.FraNL(DEF.bone) = 0.074; SYSTEM.FraNP(DEF.bone) = 0.0011;
SYSTEM.FraNL(DEF.brain) = 0.051; SYSTEM.FraNP(DEF.brain) = 0.0565;
SYSTEM.FraNL(DEF.gonads) = 0.007; SYSTEM.FraNP(DEF.gonads) = 0.0077;
SYSTEM.FraNL(DEF.heart) = 0.015; SYSTEM.FraNP(DEF.heart) = 0.0166;
SYSTEM.FraNL(DEF.kidney) = 0.0207; SYSTEM.FraNP(DEF.kidney) = 0.0162;
SYSTEM.FraNL(DEF.muscle) = 0.0238; SYSTEM.FraNP(DEF.muscle) = 0.0072;
SYSTEM.FraNL(DEF.skin) = 0.0248; SYSTEM.FraNP(DEF.skin) = 0.0111;
SYSTEM.FraNL(DEF.thymus) = 0.017; SYSTEM.FraNP(DEF.thymus) = 0.0092;
SYSTEM.FraNL(DEF.gut) = 0.0487; SYSTEM.FraNP(DEF.gut) = 0.0163;
SYSTEM.FraNL(DEF.spleen) = 0.0201; SYSTEM.FraNP(DEF.spleen) = 0.0198;
SYSTEM.FraNL(DEF.pancreas) = 0.041; SYSTEM.FraNP(DEF.pancreas) = 0.0093;
SYSTEM.FraNL(DEF.liver) = 0.0348; SYSTEM.FraNP(DEF.liver) = 0.0252;
SYSTEM.FraNL(DEF.lymphnode) = 0.0201; SYSTEM.FraNP(DEF.lymphnode) = 0.0198;
SYSTEM.FraNL(DEF.plasma) = 0.35; SYSTEM.FraNP(DEF.plasma) = 0.225;
SYSTEM.FraNL(DEF.RBC) = 0.17; SYSTEM.FraNP(DEF.RBC) = 0.29;
SYSTEM.AP(DEF.lung) = 0.5; SYSTEM.KpHSA(DEF.lung) = 0.212;
SYSTEM.AP(DEF.adipose) = 0.4; SYSTEM.KpHSA(DEF.adipose) = 0.021;
SYSTEM.AP(DEF.bone) = 0.67; SYSTEM.KpHSA(DEF.bone) = 0.1;
SYSTEM.AP(DEF.brain) = 0.4; SYSTEM.KpHSA(DEF.brain) = 0.048;
SYSTEM.AP(DEF.gonads) = 1.23; SYSTEM.KpHSA(DEF.gonads) = 0.048;
SYSTEM.AP(DEF.heart) = 3.07; SYSTEM.KpHSA(DEF.heart) = 0.157;
SYSTEM.AP(DEF.kidney) = 2.48; SYSTEM.KpHSA(DEF.kidney) = 0.13;
SYSTEM.AP(DEF.muscle) = 2.49; SYSTEM.KpHSA(DEF.muscle) = 0.025;
SYSTEM.AP(DEF.skin) = 1.32; SYSTEM.KpHSA(DEF.skin) = 0.277;
SYSTEM.AP(DEF.thymus) = 2.3; SYSTEM.KpHSA(DEF.thymus) = 0.075;
SYSTEM.AP(DEF.gut) = 2.84; SYSTEM.KpHSA(DEF.gut) = 0.158;
SYSTEM.AP(DEF.spleen) = 2.81; SYSTEM.KpHSA(DEF.spleen) = 0.097;
SYSTEM.AP(DEF.pancreas) = 1.67; SYSTEM.KpHSA(DEF.pancreas) = 0.06;
SYSTEM.AP(DEF.liver) = 5.09; SYSTEM.KpHSA(DEF.liver) = 0.086;
SYSTEM.AP(DEF.lymphnode) = 2.81; SYSTEM.KpHSA(DEF.lymphnode) = 0.097;
SYSTEM.AP(DEF.plasma) = 0.04; SYSTEM.KpHSA(DEF.plasma) = 1.0;
SYSTEM.AP(DEF.RBC) = 0.44; SYSTEM.KpHSA(DEF.RBC) = 0.0;
%KpALB does not change with age (Yan et al., 1968)
%there is no variability of tissue composition parameters
SYSTEM.FraEW = repmat(SYSTEM.FraEW, STUDY.IndNo, 1);
SYSTEM.FraIW = repmat(SYSTEM.FraIW, STUDY.IndNo, 1);
SYSTEM.FraNL = repmat(SYSTEM.FraNL, STUDY.IndNo, 1);
SYSTEM.FraNP = repmat(SYSTEM.FraNP, STUDY.IndNo, 1);
SYSTEM.AP = repmat(SYSTEM.AP, STUDY.IndNo, 1);
SYSTEM.KpHSA = repmat(SYSTEM.KpHSA, STUDY.IndNo, 1);
%remaining organ will always be the weighted mean
SYSTEM.FraEW(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.FraEW, 2)./sum(SYSTEM.Worg, 2);
SYSTEM.FraIW(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.FraIW, 2)./sum(SYSTEM.Worg, 2);
SYSTEM.FraNL(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.FraNL, 2)./sum(SYSTEM.Worg, 2);
SYSTEM.FraNP(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.FraNP, 2)./sum(SYSTEM.Worg, 2);
SYSTEM.AP(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.AP, 2)./sum(SYSTEM.Worg,2);
SYSTEM.KpHSA(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.KpHSA, 2)./sum(SYSTEM.Worg, 2);
%__Subcompartment volume [L]____________________________________________________
%%% fraction of vascular space of each tissue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%values are from Gill et al., 2016
%thymus is from Shah and Betts, 2012
%gonads is calculated from the Vint (Valentin, 2002)
SYSTEM.FraVas = ones(STUDY.IndNo, SYSTEM.CompNo);
SYSTEM.FraVas(:, DEF.lung) = 0.185*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.adipose) = 0.031*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.bone) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.brain) = -0.000545 .* Age + 0.056;
SYSTEM.FraVas(:, DEF.gonads) = 0.069*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.heart) = 0.042*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.kidney) = 0.07*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.muscle) = 0.027*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.skin) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.thymus) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.gut) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.spleen) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.pancreas) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.liver) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.lymphnode) = 0.05*ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.plasma) = ones(STUDY.IndNo, 1);
SYSTEM.FraVas(:, DEF.RBC) = ones(STUDY.IndNo, 1);
%weighted mean for the remaining organ
SYSTEM.FraVas(:, DEF.remaining) = sum(SYSTEM.Worg.*SYSTEM.FraVas, 2)./ ...
 sum(SYSTEM.Worg, 2);
SYSTEM.Vvas = SYSTEM.FraVas.*SYSTEM.Vorg;
SYSTEM.Vint = (SYSTEM.Vorg.*SYSTEM.FraEW) - (SYSTEM.Vvas.*(1 - SYSTEM.HCT));
SYSTEM.Vcel = SYSTEM.Vorg - SYSTEM.Vvas - SYSTEM.Vint;
%__pH of the each organ_________________________________________________________
%values are from Schmitt, 2008
SYSTEM.pH = zeros(1, SYSTEM.OrgNo);
SYSTEM.pH(DEF.lung) = 6.6;
SYSTEM.pH(DEF.adipose) = 7.1;
SYSTEM.pH(DEF.bone) = 7.0;
SYSTEM.pH(DEF.brain) = 7.1;
SYSTEM.pH(DEF.gonads) = 7.0;
SYSTEM.pH(DEF.heart) = 7.1;
SYSTEM.pH(DEF.kidney) = 7.22;
SYSTEM.pH(DEF.muscle) = 7.0;
SYSTEM.pH(DEF.skin) = 7.0;
SYSTEM.pH(DEF.thymus) = 7.0; %no data; take global value of Rodgers % Rowland,2005
SYSTEM.pH(DEF.gut) = 7.0;
SYSTEM.pH(DEF.spleen) = 7.0;
SYSTEM.pH(DEF.pancreas) = 7.0; %no data; take global value of Rodgers % Rowland,2005
SYSTEM.pH(DEF.liver) = 7.23;
SYSTEM.pH(DEF.lymphnode) = 7.0;
SYSTEM.pH(DEF.plasma) = 7.4; %Valentin, 2002
SYSTEM.pH(DEF.RBC) = 7.21; %Waddell, 1969
SYSTEM.pH(DEF.remaining) = 7.0; %no data; take global value of Rodgers % Rowland, 05
%__blood flows [L/h]____________________________________________________________
%data are from Stader et al., 2018
FraQorg = zeros(STUDY.IndNo, SYSTEM.OrgNo);
SYSTEM.Qorg = zeros(STUDY.IndNo, SYSTEM.OrgNo);
%%% Cardiac output (CO) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CO_Mean = 159.*SYSTEM.BSA - 1.56.*Age + 114;
CO = normrnd(CO_Mean, (21.1/100).*CO_Mean);
%%% regional blood flows - fraction of CO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FraQorg(:, DEF.lung) = 100*ones(STUDY.IndNo, 1);
FraQorg(:, DEF.adipose) = (0.044 + 0.027.*Sex) .* Age + 2.4.*Sex + 3.9;

FraQorg(:, DEF.bone) = 5*ones(STUDY.IndNo, 1);
FraQorg(:, DEF.brain) = exp(-0.48.*SYSTEM.BSA + 0.04.*Sex + 3.5);
FraQorg(:, DEF.gonads) = -0.03.*Sex + 0.05;
FraQorg(:, DEF.heart) = 4 + 1.*Sex;
FraQorg(:, DEF.kidney) = -8.7.*SYSTEM.BSA + 0.29.*SYSTEM.Height - ...
 0.081.*Age - 13;

FraQorg(:, DEF.muscle) = -6.4.*Sex + 17.5;
FraQorg(:, DEF.skin) = 5*ones(STUDY.IndNo, 1);
FraQorg(:, DEF.thymus) = 1.5*ones(STUDY.IndNo, 1);
FraQorg(:, DEF.liver) = -0.108.*Age + 1.04.*Sex + 27.9;
%hepatic arterial blood flow appears to be independent of age
FraQHA = 6.5*ones(STUDY.IndNo,1);
FraQPV = FraQorg(:, DEF.liver) - FraQHA;
%values of gut, spleen and pancreas are scaled via hepatic blood flow from ICRP
FraQorg(:, DEF.gut) = ((2.*Sex + 14).*FraQPV)./(1.5.*Sex + 19);
FraQorg(:, DEF.spleen) = (3.*FraQPV)./(1.5.*Sex + 19);
FraQorg(:, DEF.pancreas) = (1.*FraQPV)./(1.5.*Sex + 19);
%blood flow bypassing the portal vein organs gut, spleen and pancreas
FraQBY = FraQPV -...
 FraQorg(:, DEF.gut) -...
 FraQorg(:, DEF.spleen) -...
 FraQorg(:, DEF.pancreas);
FraQorg(:, DEF.lymphnode) = 1.65*ones(STUDY.IndNo, 1);
FraQorg(:, DEF.remaining) = FraQorg(:, DEF.lung) -...
 FraQorg(:, DEF.adipose) -...
 FraQorg(:, DEF.bone) -...
 FraQorg(:, DEF.gonads) -...
 FraQorg(:, DEF.heart) -...
 FraQorg(:, DEF.kidney) -...
 FraQorg(:, DEF.muscle) -...
 FraQorg(:, DEF.skin) -...
 FraQorg(:, DEF.thymus) -...
 FraQorg(:, DEF.liver) -...
 FraQorg(:, DEF.lymphnode);
%divide the fraction of blood flows by 100
FraQorg = FraQorg./100;
%%% regional blood flows %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SYSTEM.Qorg = FraQorg.*CO;
SYSTEM.QHA = (FraQHA./100).*CO;
SYSTEM.QBY = (FraQBY./100).*CO;
%__lymph flows in [L/h]_________________________________________________________
%data are from Gill et al., 2016
%%% total lymph flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LT = 0.00386.*SYSTEM.Weight;
SYSTEM.TotLymphFlow = normrnd(LT, (30/100).*LT);
%%% fraction of regional lymph flows %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FraLorg = zeros(1, SYSTEM.OrgNo);
FraLorg(DEF.lung) = 0.03;
FraLorg(DEF.adipose) = 0.128;
FraLorg(DEF.bone) = 0;
FraLorg(DEF.brain) = 0.0105;
FraLorg(DEF.gonads) = 0.013;
FraLorg(DEF.heart) = 0.01;
FraLorg(DEF.kidney) = 0.085;
FraLorg(DEF.muscle) = 0.16;
FraLorg(DEF.skin) = 0.073;
FraLorg(DEF.thymus) = 0.011;
FraLorg(DEF.gut) = 0.12;
FraLorg(DEF.spleen) = 0;
FraLorg(DEF.pancreas) = 0.003;
FraLorg(DEF.liver) = 0.33;
FraLorg(DEF.lymphnode) = 0;
%it is assumed that the fraction of lymph flow is similar for all individuals
FraLorg = repmat(FraLorg, STUDY.IndNo, 1);
%the remainig organ get the rest of the lymph flow
FraLorg(:, DEF.remaining) = 1 -...
 FraLorg(:, DEF.lung) -...
 FraLorg(:, DEF.adipose) -...
 FraLorg(:, DEF.bone) -...
 FraLorg(:, DEF.brain) -...
 FraLorg(:, DEF.gonads) -...
 FraLorg(:, DEF.heart) -...
 FraLorg(:, DEF.kidney) -...
 FraLorg(:, DEF.muscle) -...
 FraLorg(:, DEF.skin) -...
 FraLorg(:, DEF.thymus) -...
 FraLorg(:, DEF.gut) -...
 FraLorg(:, DEF.spleen) -...
 FraLorg(:, DEF.pancreas) -...
 FraLorg(:, DEF.liver);

%calculate the lymph flow
SYSTEM.Lorg = SYSTEM.TotLymphFlow.*FraLorg;
%__liver parameters_____________________________________________________________
%%% MPPGL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MPPGL = microsomal protein per gram liver in [mg / g liver]
%Barter et al., 2008
MPPGL_Mean = 10.^(0.0000024.*Age.^3 - 0.00038.*Age.^2 + 0.0158.*Age + 1.407);
MPPGL_CV = 46;
MPPGL_Min = [10, 10];
MPPGL_Max = [110, 110];
SYSTEM.MPPGL = Calc_SysPar(MPPGL_Mean, MPPGL_CV, MPPGL_Min, MPPGL_Max);
%__hepatic CYP enzyme abundance [pmol/mg]_______________________________________
%hepatic CYP abundance is from Achour et al., 2014
CYP2D6_Mean = 12.6*ones(STUDY.IndNo, 1);
CYP2D6_CV = 74;
CYP2D6_Min = [4.2, 4.2];
CYP2D6_Max = [38, 38];
CYP2J2_Mean = 1.2*ones(STUDY.IndNo, 1);
CYP2J2_CV = 58;
CYP2J2_Min = [0.4, 0.4];
CYP2J2_Max = [3.6, 3.6];
CYP3A4_Mean = 93.0*ones(STUDY.IndNo, 1);
CYP3A4_CV = 81;
CYP3A4_Min = [18.6, 18.6];
CYP3A4_Max = [601 601];
%AhC = Abundance of hepatic CYP enzymes
AhC(:, DEF.CYP2D6) = Calc_SysPar(CYP2D6_Mean, CYP2D6_CV, CYP2D6_Min, CYP2D6_Max);
AhC(:, DEF.CYP3A4) = Calc_SysPar(CYP3A4_Mean, CYP3A4_CV, CYP3A4_Min, CYP3A4_Max);
AhC(:, DEF.CYP3A5) = 0.41.*AhC(:, DEF.CYP3A4) + 56.1;
AhC(:, DEF.CYP2J2) = Calc_SysPar(CYP2J2_Mean, CYP2J2_CV, CYP2J2_Min, CYP2J2_Max);
SYSTEM.CYPhe_AB = AhC.*PheCYP;
%degradation rate of hepatic CYP enzymes in [1/h]
SYSTEM.CYPhe_kdeg = zeros(1, SYSTEM.CYPliNo);
SYSTEM.CYPhe_kdeg(DEF.CYP2D6) = 0.0143;
SYSTEM.CYPhe_kdeg(DEF.CYP3A4) = 0.0077;
SYSTEM.CYPhe_kdeg(DEF.CYP3A5) = 0.0193;
SYSTEM.CYPhe_kdeg(DEF.CYP2J2) = 0.01;
SYSTEM.CYPhe_kdeg = repmat(SYSTEM.CYPhe_kdeg, STUDY.IndNo, 1);
%__kidney parameters____________________________________________________________
%%% glomerular filtration rate (GFR) in [mL/min] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GFR = exp(-0.0079.*Age + 0.5.*SYSTEM.BSA + 4.2);
SYSTEM.GFR = normrnd(GFR, (14.7/100).*GFR);
%__parameters of the GI tract___________________________________________________
%%% Volumes of intestinal segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%total volume of intestinal segments in [L]
%division of the intestinal volume to different segments is given in the ICRP
SYSTEM.VsegCAT = zeros(STUDY.IndNo, SYSTEM.SegNo);
SYSTEM.VsegCAT(:, DEF.stomach) = 0.153 .* SYSTEM.Vorg(:, DEF.gut);
SYSTEM.VsegCAT(:, DEF.duodenum) = 0.060 .* SYSTEM.Vorg(:, DEF.gut);
SYSTEM.VsegCAT(:, DEF.jejunum) = 0.279 .* SYSTEM.Vorg(:, DEF.gut);
SYSTEM.VsegCAT(:, DEF.ileum) = 0.316 .* SYSTEM.Vorg(:, DEF.gut);
SYSTEM.VsegCAT(:, DEF.colon) = 0.192 .* SYSTEM.Vorg(:, DEF.gut);
%vascular space of the gut in [L]
%Gill et al. (2016) report a fraction of 0.05 for the gut
SYSTEM.VvasCAT = zeros(STUDY.IndNo, SYSTEM.SegNo);
SYSTEM.VvasCAT(:, DEF.duodenum) = 0.05 .* SYSTEM.VsegCAT(:, DEF.duodenum);
SYSTEM.VvasCAT(:, DEF.jejunum) = 0.05 .* SYSTEM.VsegCAT(:, DEF.jejunum);
SYSTEM.VvasCAT(:, DEF.ileum) = 0.05 .* SYSTEM.VsegCAT(:, DEF.ileum);
SYSTEM.VvasCAT(:, DEF.colon) = 0.05 .* SYSTEM.VsegCAT(:, DEF.colon);
SYSTEM.Vvas(:, DEF.gut) = sum(SYSTEM.VvasCAT, 2);
%interstitial space of the gut in [L]
SYSTEM.VintCAT = SYSTEM.FraEW(:, DEF.gut) .* SYSTEM.VsegCAT - SYSTEM.VvasCAT;
%the stomach is not modelled with interstitial space
SYSTEM.VintCAT(:, DEF.stomach) = 0;
SYSTEM.Vint(:, DEF.gut) = sum(SYSTEM.VintCAT, 2);
%%% Length of the intestine in [cm] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%length are from the ICRP (Valentin, 2002)
SYSTEM.LengthGU = zeros(STUDY.IndNo, SYSTEM.SegNo);
SYSTEM.LengthGU(:, DEF.duodenum) = 0.091 .* (1.6.*SYSTEM.Height);
SYSTEM.LengthGU(:, DEF.jejunum) = 0.426 .* (1.6.*SYSTEM.Height);
SYSTEM.LengthGU(:, DEF.ileum) = 0.483 .* (1.6.*SYSTEM.Height);
SYSTEM.LengthGU(:, DEF.colon) = 0.52 .* SYSTEM.Height + 18.5;
%%% Radius of the intestine in [cm] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%it is assumed that the intestine is a cylinder
%1000 converts L to cm³
SYSTEM.RadiusGU = sqrt((SYSTEM.VsegCAT .* 1000) ./ (pi .* SYSTEM.LengthGU));
%%% Surface of the intestine in [cm²] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SurfaceGU = 2 .* pi .* SYSTEM.RadiusGU .* SYSTEM.LengthGU;
%the intestinal surface is enlarged by plicae circulares, villi and microvilli
%data from Helander & Fändriks, 2014
FPC = zeros(1, SYSTEM.SegNo); %factor for plicae circulares
FVI = zeros(1, SYSTEM.SegNo); %factor for villi
FMV = zeros(1, SYSTEM.SegNo); %factor for microvilli
FPC(DEF.duodenum) = 1.6; FVI(DEF.duodenum) = 6.5; FMV(DEF.duodenum) = 14.6;
FPC(DEF.jejunum) = 1.6; FVI(DEF.jejunum) = 8.6; FMV(DEF.jejunum) = 9.2;
FPC(DEF.ileum) = 1.6; FVI(DEF.ileum) = 4.5; FMV(DEF.ileum) = 15.7;
FPC(DEF.colon) = 1.0; FVI(DEF.colon) = 6.5; FMV(DEF.colon) = 1.0;
FPC = repmat(FPC, STUDY.IndNo, 1);
FVI = repmat(FVI, STUDY.IndNo, 1);
FMV = repmat(FMV, STUDY.IndNo, 1);
SYSTEM.PSA = SurfaceGU.*FPC.*FVI.*FMV;
%%% enterocyte space of the gut in [L] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SYSTEM.VentCAT = (3040 .* (SurfaceGU ./ 0.0001)) .* 770 .* 1E-15;
%stomach does not have enterocytes
SYSTEM.VentCAT(:, DEF.stomach) = 0;
%luminal space of the gut in [L]
SYSTEM.VlumCAT = SYSTEM.VsegCAT - SYSTEM.VvasCAT - SYSTEM.VintCAT - ...
 SYSTEM.VentCAT;

%%% intestinal transit times [h] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%uniform distribution between 0.25 (Yu & Amidon 1999) and 0.4 (Jamei et al. 2009)
%for the gastric emptying time (GET) - for liquids only
GET = 0.25 + (0.40 - 0.25)*rand(STUDY.IndNo, 1);
SIT_Mean = 3.4*ones(STUDY.IndNo, 1);
SIT_CV = 40.2;
SIT_Min = [0.5, 0.5];
SIT_Max = [9.5, 9.5];
SIT = Calc_SysPar(SIT_Mean, SIT_CV, SIT_Min, SIT_Max);
CNT_Mean = 17.2*ones(STUDY.IndNo, 1);
CNT_CV = 20.0;
CNT_Min = [17.2/5, 17.2/5];
CNT_Max = [17.2*5, 17.2*5];
CNT = Calc_SysPar(CNT_Mean, CNT_CV, CNT_Min, CNT_Max);
%separation into different segments is based on length (Darwich et al., 2010)
SYSTEM.TransitT = zeros(STUDY.IndNo, SYSTEM.SegNo);
SYSTEM.TransitT(:, DEF.stomach) = GET;
SYSTEM.TransitT(:, DEF.duodenum) = SIT.*0.091;
SYSTEM.TransitT(:, DEF.jejunum) = SIT.*0.426;
SYSTEM.TransitT(:, DEF.ileum) = SIT.*0.483;
SYSTEM.TransitT(:, DEF.colon) = CNT;
%__intestinal CYP enzymes_______________________________________________________
%minimum / maximum is an arbitrary 3-fold difference
CYP2D6_Mean = 0.8*ones(STUDY.IndNo, 1);
CYP2D6_CV = 60;
CYP2D6_Min = [0.8/3, 0.8/3];
CYP2D6_Max = [0.8*3, 0.8*3];
CYP3A4_Mean = 66.2*ones(STUDY.IndNo, 1);
CYP3A4_CV = 60;
CYP3A4_Min = [66.3/3, 66.3/3];
CYP3A4_Max = [66.3*3, 66.3*3];
CYP3A5_Mean = 24.6*ones(STUDY.IndNo, 1);
CYP3A5_CV = 60;
CYP3A5_Min = [24.6/3, 24.6/3];
CYP3A5_Max = [24.6*3, 24.6*3];
%AiC = Abundance of intestinal CYP enzymes
AiC(:, DEF.CYP2D6) = Calc_SysPar(CYP2D6_Mean, CYP2D6_CV, CYP2D6_Min, CYP2D6_Max);
AiC(:, DEF.CYP3A4) = Calc_SysPar(CYP3A4_Mean, CYP3A4_CV, CYP3A4_Min, CYP3A4_Max);
AiC(:, DEF.CYP3A5) = Calc_SysPar(CYP3A5_Mean, CYP3A5_CV, CYP3A5_Min, CYP3A5_Max);
%intestinal CYP abundance [nmol]
CYPin_AB = AiC .* PheCYP(:, 1:SYSTEM.CYPinNo);
%separation into different segments is done according to Paine et al., 1997
SYSTEM.CYPseg_AB = zeros(STUDY.IndNo, SYSTEM.CYPinNo, SYSTEM.SegNo);
SYSTEM.CYPseg_AB(:, :, DEF.duodenum) = 0.136.*CYPin_AB;
SYSTEM.CYPseg_AB(:, :, DEF.jejunum) = 0.544.*CYPin_AB;
SYSTEM.CYPseg_AB(:, :, DEF.ileum) = 0.320.*CYPin_AB;
%degradation rate of intestinal CYP enzymes in [1/h]
SYSTEM.CYPin_kdeg = 0.03.*ones(STUDY.IndNo, SYSTEM.CYPliNo);
%===============================================================================
%%% USED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
function SysPar = Calc_SysPar(Mean, CV, Min, Max)

 %read gender differences
 MinMale = Min(1); MinFemale = Min(2);
 MaxMale = Max(1); MaxFemale = Max(2);

 %generate parameter
 SysPar = normrnd(Mean, (CV/100).*Mean);

 %truncate parameter at minimum and maximum observed values
 for sub = 1:STUDY.IndNo
 if SYSTEM.Sex(sub) == 1
 if SysPar(sub) < MinFemale || SysPar(sub) > MaxFemale
 SysPar(sub) = normrnd(Mean(sub), (CV/1000).*Mean(sub));
 end

 else
 if SysPar(sub) < MinMale || SysPar(sub) > MaxMale
 SysPar(sub) = normrnd(Mean(sub), (CV/1000).*Mean(sub));
 end
 end
 end
end
end