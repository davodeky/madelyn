
function[] = PostProcessing()
%This function processes the data from the ODE solution and outputs the results
%Attention: Some statistical calculations (geomean, prctile) require
%the Statistical and Machine Learning toolbox
global DEF %global DEF defines model parameters
global SYSTEM %global SYSTEM defines sytem parameters
global DRUG %global DRUG defines drug parameters
global DDI %global DDI enhances drug parameters for DDI prediction
global STUDY %global STUDY defines study design parameters
global MODEL %global MODEL defines parameters important for modeling
global OBS %global OBS saves observed parameters for the output
global RES %global RES saves the results for post-processing
global MEAN %global MEAN saves the mean of PK parameters
global STDV %global STDV saves the standard deviation of PK parameters
global GEOM %global GEOM saves the geometric mean of PK parameters
global PERC %global PERC saves the percentiles of PK parameters
%post-processing is exemplarily shown for venous blood concentration
%variables used in this script
InDoEv = reshape(STUDY.DoseEventMat(:, 3, :), STUDY.NoEvents, DRUG.DrugNo);
IFD = zeros(1, DRUG.DrugNo); %index for the first dosing event
ILD = zeros(1, DRUG.DrugNo); %index for the last dosing event
for drug = 1:DRUG.DrugNo
 IFD(drug) = find(InDoEv(:, drug), 1, 'first');
 ILD(drug) = find(InDoEv(:, drug), 1, 'last');
end
IFD = repmat(IFD, 1, DDI.DDINo/DRUG.DrugNo);
ILD = repmat(ILD, 1, DDI.DDINo/DRUG.DrugNo);
Dose = repmat(STUDY.Dose, STUDY.IndNo, DDI.DDINo/DRUG.DrugNo);
%__Extract the concentration from the ODE solution______________________________
%extract the venous blood concentration from the solution
VenousConc = zeros(MODEL.NumPoints, STUDY.IndNo, DDI.DDINo);
for d = 1:DDI.DDINo
 for ind = 1:STUDY.IndNo
 venous = (ind-1)*MODEL.ODENo*DDI.DDINo + MODEL.ODENo*(d-1) + ...
 DEF.plasma + SYSTEM.SubNo(DEF.plasma);
 VenousConc(:, ind, d) = RES.Conc(:, venous) .* DDI.MolW(d);
 end
end
%calculate statistics for venous blood concentration
VenousConcMean = mean(VenousConc, 2);
VenousConcPerc = prctile(VenousConc, [5, 95], 2);
%__Calculate PK parameters______________________________________________________
%Cmax / Tmax / AUCt
Venous_PKparaT = Calc_PKparaT(VenousConc, RES.Time);
%extract AUCt for the last dose event for the extrapolation of the AUCinf
VB_AUCtLast = zeros(STUDY.IndNo, DDI.DDINo);
for d = 1:DDI.DDINo
 VB_AUCtLast(:, d) = Venous_PKparaT{3} (:, ILD(d), d);
end
%%% Elimination rate and extrapolation of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
Venous_PKparaINF = Calc_PKparaINF(VenousConc, RES.Time, VB_AUCtLast);
%%% Statistics for PK parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Venous_PKparaSTAT = Calc_PKparaSTAT(Venous_PKparaT, Venous_PKparaINF);
for dn = 1:DDI.DDINo
 MEAN.Venous_CmaxFirst(:, dn) = Venous_PKparaSTAT{1,1} (:, IFD(dn), dn);
 MEAN.Venous_TmaxFirst(:, dn) = Venous_PKparaSTAT{2,1} (:, IFD(dn), dn);
 MEAN.Venous_AUCtFirst(:, dn) = Venous_PKparaSTAT{3,1} (:, IFD(dn), dn);

 MEAN.Venous_CmaxLast(:, dn) = Venous_PKparaSTAT{1,1} (:, ILD(dn), dn);
 MEAN.Venous_TmaxLast(:, dn) = Venous_PKparaSTAT{2,1} (:, ILD(dn), dn);
 MEAN.Venous_AUCtLast(:, dn) = Venous_PKparaSTAT{3,1} (:, ILD(dn), dn);

 STDV.Venous_CmaxFirst(:, dn) = Venous_PKparaSTAT{1,2} (:, IFD(dn), dn);
 STDV.Venous_TmaxFirst(:, dn) = Venous_PKparaSTAT{2,2} (:, IFD(dn), dn);
 STDV.Venous_AUCtFirst(:, dn) = Venous_PKparaSTAT{3,2} (:, IFD(dn), dn);

 STDV.Venous_CmaxLast(:, dn) = Venous_PKparaSTAT{1,2} (:, ILD(dn), dn);
 STDV.Venous_TmaxLast(:, dn) = Venous_PKparaSTAT{2,2} (:, ILD(dn), dn);
 STDV.Venous_AUCtLast(:, dn) = Venous_PKparaSTAT{3,2} (:, ILD(dn), dn);

 GEOM.Venous_CmaxFirst(:, dn) = Venous_PKparaSTAT{1,3} (:, IFD(dn), dn);
 GEOM.Venous_TmaxFirst(:, dn) = Venous_PKparaSTAT{2,3} (:, IFD(dn), dn);
 GEOM.Venous_AUCtFirst(:, dn) = Venous_PKparaSTAT{3,3} (:, IFD(dn), dn);

 GEOM.Venous_CmaxLast(:, dn) = Venous_PKparaSTAT{1,3} (:, ILD(dn), dn);
 GEOM.Venous_TmaxLast(:, dn) = Venous_PKparaSTAT{2,3} (:, ILD(dn), dn);
 GEOM.Venous_AUCtLast(:, dn) = Venous_PKparaSTAT{3,3} (:, ILD(dn), dn);

 PERC.Venous_CmaxFirst(:, dn) = Venous_PKparaSTAT{1,4} (:, IFD(dn), dn);
 PERC.Venous_TmaxFirst(:, dn) = Venous_PKparaSTAT{2,4} (:, IFD(dn), dn);
 PERC.Venous_AUCtFirst(:, dn) = Venous_PKparaSTAT{3,4} (:, IFD(dn), dn);

 PERC.Venous_CmaxLast(:, dn) = Venous_PKparaSTAT{1,4} (:, ILD(dn), dn);
 PERC.Venous_TmaxLast(:, dn) = Venous_PKparaSTAT{2,4} (:, ILD(dn), dn);
 PERC.Venous_AUCtLast(:, dn) = Venous_PKparaSTAT{3,4} (:, ILD(dn), dn);
end
MEAN.Venous_Thalf = Venous_PKparaSTAT{4,1};
STDV.Venous_Thalf = Venous_PKparaSTAT{4,2};
GEOM.Venous_Thalf = Venous_PKparaSTAT{4,3};
PERC.Venous_Thalf = Venous_PKparaSTAT{4,4};
MEAN.Venous_AUCinf = Venous_PKparaSTAT{5,1};
STDV.Venous_AUCinf = Venous_PKparaSTAT{5,2};
GEOM.Venous_AUCinf = Venous_PKparaSTAT{5,3};
PERC.Venous_AUCinf = Venous_PKparaSTAT{5,4};
MEAN.Venous_CLF = Venous_PKparaSTAT{6,1};
STDV.Venous_CLF = Venous_PKparaSTAT{6,2};
GEOM.Venous_CLF = Venous_PKparaSTAT{6,3};
PERC.Venous_CLF = Venous_PKparaSTAT{6,4};
MEAN.Venous_VDF = Venous_PKparaSTAT{7,1};
STDV.Venous_VDF = Venous_PKparaSTAT{7,2};
GEOM.Venous_VDF = Venous_PKparaSTAT{7,3};
PERC.Venous_VDF = Venous_PKparaSTAT{7,4};
MEAN.Venous_AUCratio = Venous_PKparaSTAT{8,1};
STDV.Venous_AUCratio = Venous_PKparaSTAT{8,2};
GEOM.Venous_AUCratio = Venous_PKparaSTAT{8,3};
PERC.Venous_AUCratio = Venous_PKparaSTAT{8,4};
%__output results_______________________________________________________________
%figure for concentration output
Create_Figure(RES.Time, VenousConcMean, VenousConcPerc);
%output parameters on screen
switch DRUG.DrugNo
 case 1
 d = 1;
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.Venous_CLF
(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n\n', [MEAN.Venous_VDF
(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);

 case 2
 for d = 1:DDI.DDINo
 switch d
 case 1
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 case 2
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);

 case 3
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 fprintf('AUC ratio %g %g %g\n\n', [MEAN.
Venous_AUCratio(1), PERC.Venous_AUCratio(1,1), PERC.Venous_AUCratio(3,1)]);

 case 4
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 fprintf('AUC ratio %g %g %g\n\n', [MEAN.
Venous_AUCratio(2), PERC.Venous_AUCratio(1,2), PERC.Venous_AUCratio(3,2)]);
 end
 end

 case 3
 for d = 1:DDI.DDINo
 switch d
 case 1
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 case 2
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 case 3
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 case 4
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 fprintf('AUC ratio %g %g %g\n\n', [MEAN.
Venous_AUCratio(1), PERC.Venous_AUCratio(1,1), PERC.Venous_AUCratio(3,1)]);

 case 5
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 fprintf('AUC ratio %g %g %g\n\n', [MEAN.
Venous_AUCratio(2), PERC.Venous_AUCratio(1,2), PERC.Venous_AUCratio(3,2)]);

 case 6
 fprintf('Calculated PK parameters for %s:\n', DDI.Name{d});
 fprintf('Cmax [ng/mL]: %g %g %g\n', [MEAN.
Venous_CmaxFirst(d), PERC.Venous_CmaxFirst(1,d), PERC.Venous_CmaxFirst(3,d)]);
 fprintf('Tmax [h]: %g %g %g\n', [MEAN.
Venous_TmaxFirst(d), PERC.Venous_TmaxFirst(1,d), PERC.Venous_TmaxFirst(3,d)]);
 fprintf('AUCt [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCtLast(d), PERC.Venous_AUCtLast(1,d), PERC.Venous_AUCtLast(3,d)]);
 fprintf('AUCinf [ng*mL/h]: %g %g %g\n', [MEAN.
Venous_AUCinf(d), PERC.Venous_AUCinf(1,d), PERC.Venous_AUCinf(3,d)]);
 fprintf('CL [L/h]: %g %g %g\n', [MEAN.
Venous_CLF(d), PERC.Venous_CLF(1,d), PERC.Venous_CLF(3,d)]);
 fprintf('VD [L/kg]: %g %g %g\n', [MEAN.
Venous_VDF(d), PERC.Venous_VDF(1,d), PERC.Venous_VDF(3,d)]);
 fprintf('AUC ratio %g %g %g\n\n', [MEAN.
Venous_AUCratio(3), PERC.Venous_AUCratio(1,3), PERC.Venous_AUCratio(3,3)]);
 end
 end
end
%===============================================================================
%%% USED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
%__Calculate PK parameters______________________________________________________
%%% Cmax/Tmax/AUCt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PKparaT = Calc_PKparaT(Conc, Time)

 %prepare matrices for all calculated parameters
 Cmax = zeros(STUDY.IndNo, STUDY.NoEvents, DDI.DDINo);
 Tmax = zeros(STUDY.IndNo, STUDY.NoEvents, DDI.DDINo);
 AUCt = zeros(STUDY.IndNo, STUDY.NoEvents, DDI.DDINo);
 AUC = zeros(MODEL.NumPoints, STUDY.IndNo, DDI.DDINo);

 %Time needs to be available for each individual and drug
 Time = repmat(Time, 1, STUDY.IndNo, DDI.DDINo);

 %point in time matrix to start and end each dosing event
 TimePointST = reshape(STUDY.DoseEventMat(:, 1, :) .* STUDY.Resolution, ...
 STUDY.NoEvents, DRUG.DrugNo);
 TimePointEN = reshape(STUDY.DoseEventMat(:, 2, :) .* STUDY.Resolution, ...
 STUDY.NoEvents, DRUG.DrugNo);

 %the first index for the start point cannot be 0, but needs to be 1
 TimePointST(1, :, :) = 1;

 TimePointST = repmat(TimePointST, 1, DDI.DDINo / DRUG.DrugNo);
 TimePointEN = repmat(TimePointEN, 1, DDI.DDINo / DRUG.DrugNo);
 %Extract Cmax and Tmax for each dosing event and calculate AUCt
 for sub = 1:STUDY.IndNo

 for dr = 1:DDI.DDINo

 for ev = 1:STUDY.NoEvents

 %Calculation of Cmax / Tmax are only done if a dose is given
 if TimePointST(ev, dr) ~= 0
 for t = TimePointST(ev, dr) : TimePointEN(ev, dr)
 if Conc(t, sub, dr) > Cmax(sub, ev, dr)
 Cmax(sub, ev, dr) = Conc(t, sub, dr);
 Tmax(sub, ev, dr) = Time(t, sub, dr);
 end

%trapezoidal method to calculate the AUC
if t == TimePointST(ev, dr)
 AUC(t, sub, dr) = 0;
 else
 AUC(t, sub, dr) = ((Conc(t, sub, dr) + ...
 Conc(t-1, sub, dr)) .* ...
 (Time(t, sub, dr) - ...
 Time(t-1, sub, dr)) ./ 2) + ...
 AUC(t-1, sub, dr);
 end
 end
 end

 if TimePointEN(ev, dr) ~= 0
 AUCt(sub, ev, dr) = AUC(TimePointEN(ev, dr), sub, dr);
 end

 end

 end

 end

 %save all parameters in one cell array
 PKparaT = {Cmax, Tmax, AUCt};
end
%%% extrapolate parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PKparaINF = Calc_PKparaINF(Conc, Time, AUCt)

 %take the log of the concentration
 LogConc = log10(Conc);

 %calculate the slope and beta for the last 5 time points
 slope = zeros(STUDY.IndNo, 2, DDI.DDINo);
 beta = zeros(STUDY.IndNo, DDI.DDINo);

 for dr = 1:DDI.DDINo
 for sub = 1:STUDY.IndNo
 slope(sub, :, dr) = polyfit(Time(end-4:end),...
 LogConc(end-4:end, sub, dr), 1);
 beta(sub, dr) = abs(slope(sub, 1, dr));
 end
 end

 %calculate the half-life and extrapolate the AUC to infitity
 Thalf = log(2) ./ beta;
 AUCinf = AUCt + (reshape(Conc(end, :, :), STUDY.IndNo, DDI.DDINo) ./ beta);

 %calculate the AUC ratio for DDI prediction
 AUCratio = ones(STUDY.IndNo, DRUG.DrugNo);
 if DRUG.DrugNo == 1
 AUCratio(:, 1) = 1.0;

 else
 for de = 1:DRUG.DrugNo
 AUCratio(:, de) = AUCinf(:, de + DRUG.DrugNo) ./ AUCinf(:, de);
 end
 end
 CLF = (Dose.*1000) ./ AUCinf;
 VDF = (CLF ./ beta) ./ SYSTEM.Weight;
 %save all results in one cell array
 PKparaINF = {Thalf, AUCinf, CLF, VDF, AUCratio};
end
%%% calculate statistics of PK parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PKparaSTAT = Calc_PKparaSTAT(ParTIME, ParINF)

 per = [5, 50, 95];
 Mean_Cmax = mean(ParTIME{1}, 1); Geom_Cmax = geomean(ParTIME{1}, 1);
 SD_Cmax = std(ParTIME{1}, '', 1); Perc_Cmax = prctile(ParTIME{1}, per,
1);

 Mean_Tmax = mean(ParTIME{2}, 1); Geom_Tmax = geomean(ParTIME{2}, 1);
 SD_Tmax = std(ParTIME{2}, '', 1); Perc_Tmax = prctile(ParTIME{2}, per,
1);

 Mean_AUCt = mean(ParTIME{3}, 1); Geom_AUCt = geomean(ParTIME{3}, 1);
 SD_AUCt = std(ParTIME{3}, '', 1); Perc_AUCt = prctile(ParTIME{3}, per,
1);

 Mean_Thalf = mean(ParINF{1}, 1); Geom_Thalf = geomean(ParINF{1}, 1);
 SD_Thalf = std(ParINF{1}, '', 1); Perc_Thalf = prctile(ParINF{1}, per, 1);

 Mean_AUCi = mean(ParINF{2}, 1); Geom_AUCi = geomean(ParINF{2}, 1);
 SD_AUCi = std(ParINF{2}, '', 1); Perc_AUCi = prctile(ParINF{2}, per, 1);

 Mean_CLF = mean(ParINF{3}, 1); Geom_CLF = geomean(ParINF{3}, 1);
 SD_CLF = std(ParINF{3}, '', 1); Perc_CLF = prctile(ParINF{3}, per, 1);

 Mean_VDF = mean(ParINF{4}, 1); Geom_VDF = geomean(ParINF{4}, 1);
 SD_VDF = std(ParINF{4}, '', 1); Perc_VDF = prctile(ParINF{4}, per, 1);

 Mean_AUCra = mean(ParINF{5}, 1); Geom_AUCra = geomean(ParINF{5}, 1);
 SD_AUCra = std(ParINF{5}, '', 1); Perc_AUCra = prctile(ParINF{5}, per, 1);

 %save results
 PKparaSTAT = {Mean_Cmax, SD_Cmax, Geom_Cmax, Perc_Cmax;...
 Mean_Tmax, SD_Tmax, Geom_Tmax, Perc_Tmax;...
 Mean_AUCt, SD_AUCt, Geom_AUCt, Perc_AUCt;...
 Mean_Thalf, SD_Thalf, Geom_Thalf, Perc_Thalf; ...
 Mean_AUCi, SD_AUCi, Geom_AUCi, Perc_AUCi;...
 Mean_CLF, SD_CLF, Geom_CLF, Perc_CLF;...
 Mean_VDF, SD_VDF, Geom_VDF, Perc_VDF;...
 Mean_AUCra, SD_AUCra, Geom_AUCra, Perc_AUCra};
end
%__Output plasma concentration__________________________________________________
function FigPlot = Create_Figure(Time, Mean, Perc)

 %prepare simulation time for the visualisation of the 95% CI
 DDIplotT = [Time', fliplr(Time')];
 DDIplotL = [Time(2:MODEL.NumPoints)', fliplr(Time(2:MODEL.NumPoints)')];
 switch DRUG.DrugNo

 case 1
 Plot1 = figure;

 %first subplot shows concentration
 subplot(1,2,1); hold on;
 %draw area between percentiles
 Plot1t = area(Time, [Perc(:,1,1), Perc(:,2,1) - Perc(:,1,1)],
'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot1t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot1t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time, Mean(:,1), '-k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1, 'or',
'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Venous Conc. [ng/mL]', 'fontweight', 'bold', 'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;
 %draw area between percentiles
 Plot1t = area(Time, [Perc(:,1,1), Perc(:,2,1) - Perc(:,1,1)],
'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot1t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot1t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,1), 'k',
'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1, 'or',
'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold', 'fontsize', 12);
 text(-0.3,1.05, DDI.Name{1}, 'Units', 'normalized', 'fontweight',
'bold' , 'fontsize', 14);
 set(Plot1, 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.6]);

 case 2
 for dr = 1:DDI.DDINo
 switch dr

case 1
 Plot1 = figure;

 %first subplot shows concentration
 subplot(1,2,1); hold on;
 %draw area between percentiles
 Plot1t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot1t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot1t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time, Mean(:,dr), '-k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Venous Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;
 %draw area between percentiles
 Plot1t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot1t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot1t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,dr),
'k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot1, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 2
 Plot2 = figure;

 %first subplot shows concentration
 subplot(1,2,1); hold on;
 %draw area between percentiles
 Plot2t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot2t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot2t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time, Mean(:,dr), '-k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Venous Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;
 %draw area between percentiles
 Plot2t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot2t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot2t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,dr),
'k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot2, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 3
 Plot3 = figure;

%first subplot shows concentration
 subplot(1,2,1); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(:,1,1)', fliplr(Perc(:,2,1)')];
 DDI1_area_Perp = [Perc(:,1,3)', fliplr(Perc(:,2,3)')];
 DDI1_area_SuPe = [Perc(:,2,1)', fliplr(Perc(:,1,3)')];

%fill different areas with different colours
 fill(DDIplotT, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotT, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 h = fill(DDIplotT, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(h,'facealpha',0.1)

 %plot the concentrations
 plot(Time, Mean(:,1), 'k', 'LineWidth', 1.5);
 plot(Time, Mean(:,3), '--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI1, OBS.Conc_DDI1, OBS.SD_DDI1, '+r',
'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot is in log scale
 subplot(1,2,2); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(2:MODEL.NumPoints,1,1)', fliplr
(Perc(2:MODEL.NumPoints,2,1)')];
 DDI1_area_Perp = [Perc(2:MODEL.NumPoints,1,3)', fliplr
(Perc(2:MODEL.NumPoints,2,3)')];
 DDI1_area_SuPe = [Perc(2:MODEL.NumPoints,2,1)', fliplr
(Perc(2:MODEL.NumPoints,1,3)')];

%fill different areas with different colours
 fill(DDIplotL, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotL, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 hlog = fill(DDIplotL, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(hlog,'facealpha',0.1)

%plot the concentrations
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,1),
'k', 'LineWidth', 1.5);
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,3),
'--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI1, OBS.Conc_DDI1, OBS.SD_DDI1, '+r',
'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot3, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 4
 Plot4 = figure;

%first subplot shows concentration
 subplot(1,2,1); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(:,1,2)', fliplr(Perc(:,2,2)')];
 DDI1_area_Perp = [Perc(:,1,4)', fliplr(Perc(:,2,4)')];
 DDI1_area_SuPe = [Perc(:,2,2)', fliplr(Perc(:,1,4)')];

%fill different areas with different colours
 fill(DDIplotT, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotT, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 h = fill(DDIplotT, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(h,'facealpha',0.1)

%plot the concentrations
 plot(Time, Mean(:,2), 'k', 'LineWidth', 1.5);
 plot(Time, Mean(:,4), '--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI2, OBS.Conc_DDI2, OBS.SD_DDI2, '+r',
'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);
 %second subplot is in log scale
 subplot(1,2,2); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(2:MODEL.NumPoints,1,2)', fliplr
(Perc(2:MODEL.NumPoints,2,2)')];
 DDI1_area_Perp = [Perc(2:MODEL.NumPoints,1,4)', fliplr
(Perc(2:MODEL.NumPoints,2,4)')];
 DDI1_area_SuPe = [Perc(2:MODEL.NumPoints,2,2)', fliplr
(Perc(2:MODEL.NumPoints,1,4)')];

%fill different areas with different colours
 fill(DDIplotL, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotL, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 hlog = fill(DDIplotL, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(hlog,'facealpha',0.1)

%plot the concentrations
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,2),
'k', 'LineWidth', 1.5);
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,4),
'--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI2, OBS.Conc_DDI2, OBS.SD_DDI2, '+r',
'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot4, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);
 end
 end

 case 3
 for dr = 1:DDI.DDINo
 switch dr

case 1
 Plot1 = figure;

 %first subplot shows concentration
 subplot(1,2,1); hold on;
 %draw area between percentiles
 Plot1t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot1t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot1t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time, Mean(:,dr), '-k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Venous Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;
 %draw area between percentiles
 Plot1t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot1t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot1t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,dr),
'k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot1, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 2
 Plot2 = figure;

 %first subplot shows concentration
 subplot(1,2,1); hold on;
 %draw area between percentiles
 Plot2t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot2t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot2t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time, Mean(:,dr), '-k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Venous Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;
 %draw area between percentiles
 Plot2t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot2t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot2t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,dr),
'k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot2, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 3
 Plot3 = figure;

 %first subplot shows concentration
 subplot(1,2,1); hold on;
 %draw area between percentiles
 Plot3t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot3t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot3t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time, Mean(:,dr), '-k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug3, OBS.Conc_Drug3, OBS.SD_Drug3,
'or', 'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Venous Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;
 %draw area between percentiles
 Plot3t = area(Time, [Perc(:,1,dr), Perc(:,2,dr) - Perc(:,1,
dr)], 'LineStyle', 'none');
 %the area between x.axis and 5% CI should be white
 Plot3t(1).FaceColor = [1 1 1];
 %the area between 5 and 95% CI gets a light green
 Plot3t(2).FaceColor = [0.88 0.94 0.85];
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,dr),
'k', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug3, OBS.Conc_Drug3, OBS.SD_Drug3,
'or', 'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot3, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 4
 Plot4 = figure;

%first subplot shows concentration
 subplot(1,2,1); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(:,1,1)', fliplr(Perc(:,2,1)')];
 DDI1_area_Perp = [Perc(:,1,4)', fliplr(Perc(:,2,4)')];
 DDI1_area_SuPe = [Perc(:,2,1)', fliplr(Perc(:,1,4)')];

%fill the area between percentiles
 fill(DDIplotT, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotT, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 h = fill(DDIplotT, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(h,'facealpha',0.1)

%plot the concentration
 plot(Time, Mean(:,1), 'k', 'LineWidth', 1.5);
 plot(Time, Mean(:,4), '--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI1, OBS.Conc_DDI1, OBS.SD_DDI1, '+r',
'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(2:MODEL.NumPoints,1,1)', fliplr
(Perc(2:MODEL.NumPoints,2,1)')];
 DDI1_area_Perp = [Perc(2:MODEL.NumPoints,1,4)', fliplr
(Perc(2:MODEL.NumPoints,2,4)')];
 DDI1_area_SuPe = [Perc(2:MODEL.NumPoints,2,1)', fliplr
(Perc(2:MODEL.NumPoints,1,4)')];

%fill the area between percentiles
 fill(DDIplotL, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotL, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 hlog = fill(DDIplotL, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(hlog,'facealpha',0.1)

%plot the concentration
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,1),
'k', 'LineWidth', 1.5);
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,4),
'--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug1, OBS.Conc_Drug1, OBS.SD_Drug1,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI1, OBS.Conc_DDI1, OBS.SD_DDI1, '+r',
'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot4, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 5
 Plot5 = figure;

%first subplot shows concentration
 subplot(1,2,1); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(:,1,2)', fliplr(Perc(:,2,2)')];
 DDI1_area_Perp = [Perc(:,1,5)', fliplr(Perc(:,2,5)')];
 DDI1_area_SuPe = [Perc(:,2,2)', fliplr(Perc(:,1,5)')];

%fill the area between percentiles
 fill(DDIplotT, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotT, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 h = fill(DDIplotT, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(h,'facealpha',0.1)

%plot the concentration
 plot(Time, Mean(:,2), 'b', 'LineWidth', 1.5);
 plot(Time, Mean(:,5), '--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI2, OBS.Conc_DDI2, OBS.SD_DDI2, '+r',
'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(2:MODEL.NumPoints,1,2)', fliplr
(Perc(2:MODEL.NumPoints,2,2)')];
 DDI1_area_Perp = [Perc(2:MODEL.NumPoints,1,5)', fliplr
(Perc(2:MODEL.NumPoints,2,5)')];
 DDI1_area_SuPe = [Perc(2:MODEL.NumPoints,2,2)', fliplr
(Perc(2:MODEL.NumPoints,1,5)')];

%fill the area between percentiles
 fill(DDIplotL, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotL, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 hlog = fill(DDIplotL, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(hlog,'facealpha',0.1)

%plot the concentration
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,2),
'b', 'LineWidth', 1.5);
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,5),
'--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug2, OBS.Conc_Drug2, OBS.SD_Drug2,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI2, OBS.Conc_DDI2, OBS.SD_DDI2, '+r',
'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Plot5, 'Units', 'normalized', 'Position', [0.2 0.2 0.6
0.6]);

 case 6
 Venous6 = figure;

%first subplot shows concentration
 subplot(1,2,1); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(:,1,3)', fliplr(Perc(:,2,3)')];
 DDI1_area_Perp = [Perc(:,1,6)', fliplr(Perc(:,2,6)')];
 DDI1_area_SuPe = [Perc(:,2,3)', fliplr(Perc(:,1,6)')];

 %fill the area between percentiles
 fill(DDIplotT, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotT, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 h = fill(DDIplotT, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(h,'facealpha',0.1)

%plot the concentration
 plot(Time, Mean(:,3), 'b', 'LineWidth', 1.5);
 plot(Time, Mean(:,6), '--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug3, OBS.Conc_Drug3, OBS.SD_Drug3,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI3, OBS.Conc_DDI3, OBS.SD_DDI3, '+r',
'LineWidth', 1.5);
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 set(gca, 'fontsize', 12);

 %second subplot shows concentration on a log scale
 subplot(1,2,2); hold on;

%draw area between percentiles
 DDI1_area_Subst = [Perc(2:MODEL.NumPoints,1,3)', fliplr
(Perc(2:MODEL.NumPoints,2,3)')];
 DDI1_area_Perp = [Perc(2:MODEL.NumPoints,1,6)', fliplr
(Perc(2:MODEL.NumPoints,2,3)')];
 DDI1_area_SuPe = [Perc(2:MODEL.NumPoints,2,3)', fliplr
(Perc(2:MODEL.NumPoints,1,6)')];

%fill the area between percentiles
 fill(DDIplotL, DDI1_area_Subst, [0.88 0.94 0.85],
'Edgecolor', 'none');
 fill(DDIplotL, DDI1_area_Perp, [0.81 0.92 1], 'Edgecolor',
'none');
 hlog = fill(DDIplotL, DDI1_area_SuPe, [0.35 0.63 0],
'Edgecolor', 'none');
 set(hlog,'facealpha',0.1)

%plot the concentration
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,3),
'b', 'LineWidth', 1.5);
 plot(Time(2:MODEL.NumPoints), Mean(2:MODEL.NumPoints,6),
'--b', 'LineWidth', 1.5);
 errorbar(OBS.Time_Drug3, OBS.Conc_Drug3, OBS.SD_Drug3,
'or', 'LineWidth', 1.5);
 errorbar(OBS.Time_DDI3, OBS.Conc_DDI3, OBS.SD_DDI3, '+r',
'LineWidth', 1.5);
 set(gca, 'fontsize', 12);
 set(gca, 'Yscale', 'log');
 set(gca, 'ycolor', 'k');
 set(gca, 'xcolor', 'k');
 xlabel('Time [h]', 'fontweight', 'bold', 'fontsize', 12);
 ylabel('Plasma Conc. [ng/mL]', 'fontweight', 'bold',
'fontsize', 12);
 text(-0.3,1.05, DDI.Name{dr}, 'Units', 'normalized',
'fontweight', 'bold' , 'fontsize', 14);
 set(Venous6, 'Units', 'normalized', 'Position', [0.2 0.2
0.6 0.6]);
 end
 end
 end
FigPlot = 1;
end
end