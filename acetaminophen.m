
function[] = acetaminophen()

%name of the function needs to be the drug name in small letters
%otherwise the code will not work

%IMPORTANT: THE USER IS RESPONSIBLE FOR SENSFUL VALUES - NO CHECK IS PERFORMED!

global DRUG NAME DEF

NAME = 1;

%__PhysChem_____________________________________________________________________
%molecular weight in g/mol
%DRUG.MolW(1)                          = 151.163;
DRUG.MolW(1)                          = 151.163;

%octanol-water partition coefficient
DRUG.logP(1)                          = 0.49;                        

%decide if the drug is neutral (neutral), a monoprotic base (mono_base),
%a diprotic base /di_base), a monoprotic acid (mono_acid), 
%a diprotic acid (di_acid) or a zwitterion (zwitterion)
%DRUG.type(1)                          = NAME.neutral;             %http://websites.umich.edu/~chemstu/content_weeks/F_06_Week10/pictures_and_notes/drugsp855.pdf

%pKa of the strongest acid / base (important for distribution and dissolution)
%pka = 0 if not used (e.g. for neutral or monoprotic compounds)
%if the substance has more than 2 pKa, use the strongest one
DRUG.pka1(1)                          = 9.38;                        %https://pubchem.ncbi.nlm.nih.gov/compound/Acetaminophen#section=Environmental-Bioconcentration
DRUG.pka2(1)                          = 0.0;                        %not sure how it's neutral and has this pKa but references claim both?

%blood-plasma ratio
DRUG.BP(1)                            = 1.37;                      %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3991461/

%fraction unbound in plasma
DRUG.fu(1)                            = 0.759;                      %https://pubmed.ncbi.nlm.nih.gov/7832576/
%main binding protein (albumin or AAG)
DRUG.protein(1)                       = NAME.albumin;                   %https://pubmed.ncbi.nlm.nih.gov/7832576/

%__Absorption___________________________________________________________________
%parameters for the first order absorption model
%if "predicted" is switched on (default), the values will be predicted
%index for the prediction vector
fa = 1;    ka = 2;    Qgut = 3;    fugut = 4;
%if fa, ka and Qgut are user-defined and not predicted, a CV can be attached
%ka is in [1/h] and Qgut in [L/h]

%DRUG.predicted(fa, 1)                 = DRUG.ON;
DRUG.fa(1)                            = 1.0;                        %reference
DRUG.faCV(1)                          = 30.0;                       %reference

%DRUG.predicted(ka, 1)                 = DRUG.ON;
DRUG.ka(1)                            = 1.0;                        %reference
DRUG.kaCV(1)                          = 30.0;                       %reference

%DRUG.predicted(Qgut, 1)               = DRUG.ON;
DRUG.Qgut(1)                          = 1.0;                        %reference
DRUG.QgutCV(1)                        = 30.0;                       %reference

%DRUG.predicted(fugut, 1)              = DRUG.ON;
DRUG.fugut(1)                         = 1.0;                        %reference

%apparent permeability from Caco-2 cells in 10^-6 cm/s
%IMPORTANT: in the model, this is assumed to be passive only
DRUG.Papp(1)                          = 263;                        %https://www.sciencedirect.com/science/article/pii/S2214750019305542

%rate of permeability to enhance / slow down intestinal permeability in 1/h
DRUG.kPerUP(1)                        = 0.3675;                        %check units: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6824060/

%lag time to delay absorption in [h]
DRUG.LagTime(1)                       = 0.0;                        %https://aasldpubs.onlinelibrary.wiley.com/doi/10.1002/hep4.1406 ,assumed to be instantaneous uptake in most models

%factor if absorption in the colon occurs (=1) or not (=0)
DRUG.FabsColon(1)                     = 1.0;                        %https://pubmed.ncbi.nlm.nih.gov/34689/

%Factor for the permeability from the enterocyte into the systemic circulation
DRUG.GUPermScalar(1)                  = 1.0;                   %reference

%unspecific clearance in the intestine (e.g. spontaneous hydrolysis, UGT in L/h)
DRUG.CLun(1)                          = 0.0;                        %what is this? 0 because assuming very fast absorption?

%__Distribution_________________________________________________________________
%scalar to alter the predicted partition coefficient for all compartments or 
%specific compartments only
DRUG.KpScalarAll(1)                   = 1.0;                        %should this be altered based on 85% clearance in 24h?
% DRUG.KpScalar(NAME.adipose, 1)        = 1.0;                        %reference

%create a ratio between influx and efflux into the intracellular space
DRUG.JinScalarAll(1)                  = 1.0;                        %https://pubmed.ncbi.nlm.nih.gov/7039926/
% DRUG.JinScalar(NAME.adipose, 1)       = 1.0;                        %rapid and even distribution widely cited
% DRUG.JinScalar(NAME.muscle, 1)        = 1.0;
% DRUG.JinScalar(NAME.liver, 1)         = 1.0;

%restrict diffusion into the intracellular space to better capture VdF
DRUG.FinScalarAll(1)                  = 1.0;                        %reference
% DRUG.FinScalar(NAME.adipose, 1)       = 1.0;                        %reference;

%__Elimination__________________________________________________________________
%either CLint or Vmax and KM need to be entered per isoform enzyme per pathway
%4 different pathways can be assigned to each enzymes
%name should be: DRUG.Clint_name of enzyme_number of pathway(NAME.name)
%CLint is in microL/min/pmol / Vmax is in pmol/min/pmol / KM is in microM
%reference for all: https://onlinelibrary.wiley.com/doi/10.1111/j.1440-1681.2008.05029.x

%formation of paracetamol glucuronide
DRUG.Vmax_mG_1(1)       = 6.4*10^10;                      
DRUG.KM_mG_1(1)         = 7000;                        

%formation of paracetamol sulphate
DRUG.Vmax_mS_1(2)       = 7.4*10^11;                       
DRUG.KM_mS_1(2)         = 97;                      

%formation of paracetamol cysteine and mercapturate conjugates
DRUG.Vmax_mCM_2(3)       = 2.5*10^8;                        
DRUG.KM_mCM_2(3)         = 300;                                         

%hepatic intrinsic clearance not assigned to an enzyme in microL/min/mg
DRUG.CLint(1)                         = 0.0;                        %reference

%renal clearance in L/h
DRUG.CLrenal(1)                       = 0.942;                    	 %only 2-5% excreted in unchanged form through urine, https://link.springer.com/article/10.1007/BF00558162 using healthy patient values
DRUG.CLrenalCV(1)                     = 15.0;                       %idk what this is??

%biliary clearance in L/h
DRUG.CLbile(1)                        = 3.88;                        %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1401122/?page=1
DRUG.CLbileCV(1)                      = 0.0;                       %reference

%additional plasma clearance; can be used to asign an in vivo clearance
DRUG.CLadditional(1)                  = 20.22;                        %https://pubmed.ncbi.nlm.nih.gov/284720/
DRUG.CLadditionalCV(1)                = 0.0;                       %reference

%__Transporter__________________________________________________________________
%CLpd is the passive diffusion in L/h/Mio cells
DRUG.CLpd_LI(1)                       = 0.0;                        %reference
DRUG.CLpd_GU(1)                       = 0.0;                        %reference

%naming includes the location
%ha = hepatic and in = intestinal
%CLint is in microl/min/pmol / Vmax is in pmol/min/pmol / KM is in microM
DRUG.CLint_TRA_LI(1, 1)    = 0.0;                        %reference

% DRUG.Vmax_TRA_GU(NAME.MDR1, 1)        = 0.0;                        %reference
DRUG.Vmax_TRA_GU(1, 1)        = 0.0;
DRUG.KM_TRA_GU(NAME.MDR1, 1)          = 0.0;                        %reference

%secretion rate into the gut if biliary transporters are used in [L/h]
DRUG.Kbile(1)                         = 0.0;                        %reference

%__DDI__________________________________________________________________________
%this part is for perpetrator drug models

%%% competitive inhibition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ki = inhibition constant in microM
DRUG.Ki_CYP(name.CYP2C9, 1)           = 0.0;                     %reference
DRUG.Ki_CYP(NAME.CYP2C19, 1)          = 0.0;                       %reference
DRUG.Ki_CYP(NAME.CYP2D6, 1)           = 0.0;                    %reference
DRUG.Ki_CYP(NAME.CYP3A4, 1)           = 0.0;                 %reference

DRUG.Ki_TRA_LI(NAME.MDR1, 1)          = 0.0;                  %reference
DRUG.Ki_TRA_LI(NAME.OATP1B1, 1)       = 0.0;                      %reference
DRUG.Ki_TRA_LI(NAME.OATP1B3, 1)       = 0.0;                      %reference
DRUG.Ki_TRA_LI(NAME.BCRP, 1)          = 0.0;                       %reference
% 
DRUG.Ki_TRA_GU(NAME.MDR1, 1)          = 0.0;                  %reference
DRUG.Ki_TRA_GU(NAME.BCRP, 1)          = 0.0;                       %reference

%%% mechanism-based inhibition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%assume no inhibition
%kinact = maximum inactivation rate constant in 1/h
%Kapp   = apparent enzyme inhibition constant for MBI, where kinact/2 in microM
% DRUG.kinact_CYP(NAME.CYP3A4, 1)       = 0.0;              %reference
% DRUG.Kapp_CYP(NAME.CYP3A4, 1)         = 0.0;                   %reference
% 
% DRUG.kinact_CYP(NAME.CYP3A5, 1)       = 0.0;             %reference
% DRUG.Kapp_CYP(NAME.CYP3A5, 1)         = 0.0;                  %reference
% 
% DRUG.kinact_CYP(NAME.CYP2J2, 1)       = 0.0;            %reference
% DRUG.Kapp_CYP(NAME.CYP2J2, 1)         = 0.0;                  %reference
% 
% %%% induction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %IndMax is the maximum fold of induction
% %IC50 is the half maximum inhibitory concentration
% DRUG.IndMax_CYP(NAME.CYP3A4, 1)       = 0.0;                       %reference
% DRUG.IC50_CYP(NAME.CYP3A4, 1)         = 0.0;                       %reference
% 
% DRUG.IndMax_UGT(NAME.UGT1A1, 1)       = 0.0;                        %reference
% DRUG.IC50_UGT(NAME.UGT1A1, 1)         = 0.0;                       %reference

end









