function [H, f, F, AF, SF, Long, VCTR, GRDVctr, midLinegrd, indGrdReg] = simulate_acoustics(...
    cnt, noms_ss_cnt, ind_deb_ss_cnt, ind_fin_ss_cnt, ind_LmrksCnt, name_LmrksCnt, sex, constrNaso, constrMin, nbGrd)

% function [H, f, F, AF, SF, Long, VCTR, GRDVctr, midLinegrd, indGrdReg] = simulate_acoustics(...
%    cnt, namesCnts, indStartCnts, indEndCnts, indLmrksCnt, nameLmrksCnt, sex,...
%    {constrNaso, constrMin, nbGrd});
% 
% Given the 2D contours of an articulation, returns the acoustic transfer
% function, the formants, the area function and the sagittal function.
%
% The main steps of the function are:
%   1. Calculation of the lower and upper vocal tract
%   2. Calculation of sagittal distance
%   3. Calculation of the area function
%   4. Calcuation of the transfer function and the formants
% 
% Inputs
%     cnt(nbPts,iXY)       : Articulator contours, split by NaN (no NaN at the end)
%     namesCnts(nbOrg)   : Name of the articulators contours (in cells, ex: noms_ss_cnt={'jaw','tongue'})
%                            See below in the help the required contours and their names
%     indStartCnts(nbOrg): Start indices of each articulator in 'cnt'
%     indEndCnts(nbOrg): End indices of each articulator in 'cnt'
%     indLmrksCnt(NbLmrk) : Indices of the landmarks in 'cnt'
%                            See below in the help the required landmarks and their names
%     nameLmrksCnt(NbLmrk): Names of the landmarks
%                            See below in the help the required landmarks and their names
%     sex(1)               : Sex of the subject: 'm', 'f' or 'u' (undisclosed, we take average data between m and f)
%     constrNaso(1)        : Size of the opening under which we close the nasopharyngeal port [default = 0.5]
%     constrMin(1)         : Imposed minimal constriction in cm2 (set for the area function and the subsequent statistics) [default = 0]
%     nbGrd(1)             : Initialisation of the number of grids (note. the final result might be slightly different) [default = 100]
% 
% Outputs
%     H(nbFreq)            : Amplitude of the transfer function (plot: '20*log10(abs(H))')
%                            Output of the function 'vtn2frm_ftr'
%     f(nbFreq)            : Frequencies
%     F(nbForm)            : Formants
%                            Output of the function 'vtn2frm_ftr'
%     AF(nbGrid)           : Area function
%     SF(nbGrid)           : Sagittal function (sagittal distances, before the alpha-beta model)
%     Long(nbGrid)         : Length of each tube
%     VCTR(struct)         : Contours of the vocal tract
%            VCTR.cntMainInt(nbPts,iXY): Main contour of the internal vocal tract, from larynx toi lips
%            VCTR.cntMainExt(nbPts,iXY): Main contour of the external vocal tract, from larynx toi lips
%            VCTR.cntIslandsInt(nbPts,iXY): Islands of the internal vocal tract
%            VCTR.cntIslandsExt(nbPts,iXY): Islands of the external vocal tract
%     GRDVctr(nbGrid)      : Structure with information for each grid
%                            Output of the function 'crossarea_2_tract_2D'
%            GRDVctr(iGrd).X(nbPts): X values of the intersections of the vocal tract and the line grid
%            GRDVctr(iGrd).Y(nbPts): Y values of the intersections of the vocal tract and the line grid
%            GRDVctr(iGrd).X2D_grd_Int(1): X values of ONE point of the intersections of the internal vocal tract and the line grid
%            GRDVctr(iGrd).Y2D_grd_Int(1): Y values of ONE point of the intersections of the internal vocal tract and the line grid
%            GRDVctr(iGrd).X2D_grd_Ext(1): X values of ONE point of the intersections of the external vocal tract and the line grid
%            GRDVctr(iGrd).Y2D_grd_Ext(1): Y values of ONE point of the intersections of the external vocal tract and the line grid
%            GRDVctr(iGrd).D(1): Sagittal distance of the grid (before correction with cosalp)
%            GRDVctr(iGrd).Xc(1): X value of the vocal tract center of this grid
%            GRDVctr(iGrd).Yc(1): Y value of the vocal tract center of this grid
%            GRDVctr(iGrd).cosalp(1): Correction parameter to compensate oblique grid lines
%     midLinegrd(nbPts,iXY): Estimated vocal tract midline used to build the grid
%     indGrdReg(nbGrid)    : Detected region of the vocal tract for each grid line
%                            as briefly defined in Soquet et al. 2002:
%                                 1 -> Larynx
%                                 2 -> Low pharynx
%                                 3 -> Mid pharynx
%                                 4 -> Oropharynx
%                                 5 -> Velum
%                                 6 -> Hard palate
%                                 7 -> Alveolar region
%                                 8 -> Labial region
% 
% 
% *******************************************************************
% Required contours
% 
% 'palateFull'
% 'pharynx'
% 'jaw'
% 'epiglottis'
% 'lowerLip'
% 'upperLip'
% 'tongue'
% 'backLarynx'
% 
% *******************************************************************
% Required landmarks
% 
% Palate:
% 'PUT' = Palate Upper Teeth
% 'PSS' = Palate Start Soft
% 'PBEL' = Palate Bone End Low
% 'PVT' = Palate Velum Tip
% 'PBEU' = Palate Bone End Up
% 
% Pharynx:
% 'PhL' = Pharynx Low
% 'PhU' = Pharynx Up
% 
% Epiglottis:
% 'ET' = Epiglottis Tip
% 
% Lower lip:
% 'LLVO' = Lower Lip Vermillon Out
% 'LLVI' = Lower Lip Vermillon In
% 
% Upper lip:
% 'ULVO' = Upper Lip Vermillon Out
% 'ULVI' = Upper Lip Vermillon In
% 
% Tongue:
% 'TJUB' = Tongue Jaw Up Back
% 'TJLB' = Tongue Jaw Low Back
% 'TS' = Tongue Sub
% 'TE' = Tongue Epiglottis
% 'THTL' = Tongue Hyoid Tip Low
% 
% Glottis:
% 'GF' = Glottis Front
% 
% BackLarynx:
% 'BLG' = BackLarynx Glottis
% 'BLPh' = BackLarynx Pharynx
% 
% *******************************************************************
% Example plots of the results
% 
% -----------------------------------------
% 1: Plot of the vocal tract with the grid
% 
% figure
% FIG
% plotf(cnt);
% plotf(VCTR.cntMainInt, 'r', 'LineWidth', 2);
% plotf(VCTR.cntMainExt, 'm', 'LineWidth', 2);
% for iGrd = 1:length(GRDVctr)
%     if mod(length(GRDVctr(iGrd).X),2) == 0
%         for ii = 1:2:length(GRDVctr(iGrd).X)
%             plot(GRDVctr(iGrd).X(ii+[0,1]), GRDVctr(iGrd).Y(ii+[0,1]), 'k');
%         end  % for ii = 1:length(GRD(iGrd).X)/2
%     end  % if mod(length(GRD(iGrd).X),2) == 0
% end  % for iGrd = 1:nbGrdCircle
% 
% -----------------------------------------
% 2: Plot of the area function
% 
% figure;
% FIG2
% plot_FA_oral_nasal(AF',Long',0,'haut','b',2,1);
% axis([0 sum(Long)+1 0 max(AF)+1])
% 
% -----------------------------------------
% 3: Plot of the acoustic transfer function
% 
% figure
% FIG2
% plot(f, 20*log10(abs(H)))
% 
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 09/10/2015



%% ======================================================================
% Optional inputs

constrNaso_default = 0.5;
constrMin_default = 0;
nbGrd_default = 100;

if nargin < 8
    constrNaso = constrNaso_default;
end  % if nargin < 8
if isnan(constrNaso) || isinf(constrNaso)
    constrNaso = constrNaso_default;
end  % if isnan(constrNaso) || isinf(constrNaso)

if nargin < 9
    constrMin = constrMin_default;
end  % if nargin < 9
if isnan(constrMin) || isinf(constrMin)
    constrMin = constrMin_default;
end  % if isnan(constrMin) || isinf(constrMin)

if nargin < 10
    nbGrd = nbGrd_default;
end  % if nargin < 10
if isnan(nbGrd) || isinf(nbGrd)
    nbGrd = nbGrd_default;
end  % if isnan(constrMin) || isinf(constrMin)



%% ======================================================================
% Renaming of the inputs

cntAlign = cnt;
noms_ss_cntAlign = noms_ss_cnt;
ind_deb_ss_cntAlign = ind_deb_ss_cnt;
ind_fin_ss_cntAlign = ind_fin_ss_cnt;
ind_LmrksCntAlign = ind_LmrksCnt;
name_LmrksCntAlign = name_LmrksCnt;


%% ======================================================================
% Constants

% Dilatation of the contours to fusion Lower Lip - Jaw and Upper Lip - Palate
lDilat = 0.025; % in cm

% Percentage of the external velum contour (starting at PVT) that we remove
% to calculate the nasopharyngeal constriction
pVelRemove = 0.2;

% Number of points of each of the internal and external vocal tract
% contours to estimate the first rough midline
nbPtsVctrMidlineInit = 500;

% Number of grid lines on the first grid based on the circle
nbGrdCircle = 20;

% Percentage of allowed distance when forcing the addition of the first and
% last grid lines
pDistLine = 0.33;

% Frequencies for the acoustic simulation
nbfreq = 256; Fmax = 5000; Fmin = Fmax / nbfreq;
f = linspace(Fmin, Fmax, nbfreq);


% Alpha-Beta Model ---------------------------------------------
% From Soquet 2002
if strcmp(sex,'m')
    % Valeurs de alpha-beta en fonction de la zone du conduit vocal pour une femme
    alp = [0.78, 1.86, 1.74, 1.99, 1.84, 1.82, 2.67, 2.42];
    bet = [0.7, 1.22, 1.23, 0.81, 0.93, 1.43, 1.48, 1.67];
elseif strcmp(sex,'f')  % if strcmp(listSex{iSubject},'m')
    % Valeurs de alpha-beta en fonction de la zone du conduit vocal pour un homme
    alp = [1.11, 1.79, 1.34, 0.73, 1.39, 1.34, 1.92, 4.72];
    bet = [2.35, 1.38, 1.62, 1.81, 1.08, 1.51, 1.20, 2.48];
elseif strcmp(sex,'u')  % if strcmp(listSex{iSubject},'m')
    % Average between male and females
    alp_m = [0.78, 1.86, 1.74, 1.99, 1.84, 1.82, 2.67, 2.42];
    bet_m = [0.7, 1.22, 1.23, 0.81, 0.93, 1.43, 1.48, 1.67];
    alp_f = [1.11, 1.79, 1.34, 0.73, 1.39, 1.34, 1.92, 4.72];
    bet_f = [2.35, 1.38, 1.62, 1.81, 1.08, 1.51, 1.20, 2.48];
    alp = (alp_m + alp_f) / 2;
    bet = (bet_m + bet_f) / 2;
end  % if strcmp(listSex{iSubject},'m')
% Zones
iLar = 1; % Larynx
iLowPha = 2; % Low pharynx
iMidPha = 3; % Mid pharynx
iOroPha = 4; % Oropharynx
iVel = 5; % Velum
iHardPal = 6; % Hard palate
iAlv = 7; % Alveolar region
iLab = 8; % Labial region
% ---------------------------------------------


%% ======================================================================
% Landmarks and contours used thoughout the function

% Landmark indices in the contour
iULVO = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'ULVO')));
iULVI = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'ULVI')));
iLLVO = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'LLVO')));
iLLVI = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'LLVI')));
iPhL = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'PhL')));
iPhU = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'PhU')));
iPVT = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'PVT')));
iPBEU = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'PBEU')));
iTE = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'TE')));
iTJUB = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'TJUB')));
iBLG = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'BLG')));
iBLPh = ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'BLPh')));

% Landmark values
PhL = cntAlign(iPhL,:);
TE = cntAlign(iTE,:);
TJUB = cntAlign(iTJUB,:);
PVT = cntAlign(iPVT,:);
TS = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'TS'))),:);
ET = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'ET'))),:);
PBEL = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'PBEL'))),:);
PUT = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'PUT'))),:);
PSS = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'PSS'))),:);
TJLB = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'TJLB'))),:);
THTL = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'THTL'))),:);
GF = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'GF'))),:);
BLPh = cntAlign(ind_LmrksCntAlign(find(strcmp(name_LmrksCntAlign,'BLPh'))),:);

% Middle of TS and ET
TS_ET = (TS+ET)./2;

% Epiglottis contour
iOrg = find(strcmp(noms_ss_cntAlign, 'epiglottis'));
indEpi = ind_deb_ss_cntAlign(iOrg):ind_fin_ss_cntAlign(iOrg);
cntEpi = cntAlign(indEpi,:);

% Lower Lip contour
iOrg = find(strcmp(noms_ss_cntAlign, 'lowerLip'));
indLL = ind_deb_ss_cntAlign(iOrg):ind_fin_ss_cntAlign(iOrg);
cntLL = cntAlign(indLL,:);

% Jaw contour
iOrg = find(strcmp(noms_ss_cntAlign, 'jaw'));
indJ = ind_deb_ss_cntAlign(iOrg):ind_fin_ss_cntAlign(iOrg);
cntJ = cntAlign(indJ,:);
% Other naming
indJaw = indJ;
cntJaw = cntJ;

% Upper Lip contour
iOrg = find(strcmp(noms_ss_cntAlign, 'upperLip'));
indUL = ind_deb_ss_cntAlign(iOrg):ind_fin_ss_cntAlign(iOrg);
cntUL = cntAlign(indUL,:);

% Palate contour
iOrg = find(strcmp(noms_ss_cntAlign, 'palateFull'));
indP = ind_deb_ss_cntAlign(iOrg):ind_fin_ss_cntAlign(iOrg);
cntP = cntAlign(indP,:);

% Pharynx contour
iOrg = find(strcmp(noms_ss_cntAlign, 'pharynx'));
indPha = ind_deb_ss_cntAlign(iOrg):ind_fin_ss_cntAlign(iOrg);
cntPha = cntAlign(indPha,:);


%% ======================================================================
% 1. Calculation of the lower and upper vocal tract

%---------------------------------------------------------------------
% 'Vertical' tangent plane to the lips to close the vocal tract at the lips

cntLocLipUp = cntAlign(iULVO:((iULVO<iULVI)*2-1):iULVI,:);
cntLocLipLow = cntAlign(iLLVO:((iLLVO<iLLVI)*2-1):iLLVI,:);

% We simply try all the combinations
thresh = 1e-3;
indCouples = [];
for iPtUp = 1:size(cntLocLipUp,1)
    for iPtLow = 1:size(cntLocLipLow,1)
        ptUp = cntLocLipUp(iPtUp,:);
        ptLow = cntLocLipLow(iPtLow,:);
        Droite = [ptLow; ptUp - ptLow];
        ptsInterUp = remDouble(intersection_line_contour(Droite, cntLocLipUp([1:end, 1],:)));
        ptsInterLow = remDouble(intersection_line_contour(Droite, cntLocLipLow([1:end, 1],:)));
        ptsInterUp = remDouble([ptsInterUp; ptUp], thresh);
        ptsInterLow = remDouble([ptsInterLow; ptLow], thresh);
        if (size(ptsInterUp,1) == 1) && (size(ptsInterLow,1) == 1)
            indCouples = [indCouples; iPtUp, iPtLow];
        end  % if (size(ptsInterUp,1) == 1) && (size(ptsInterLow,1) == 1)
    end  % for iPtLow = 1:size(cntLocLipLow,1)
end  % for iPtUp = 1:size(cntLocLipUp,1)

% We keep the solution that leaves the 2 contours on its right
nbUp = size(cntLocLipUp,1);
nbLow = size(cntLocLipLow,1);
for iCouple = 1:size(indCouples,1)
    ptUp = cntLocLipUp(indCouples(iCouple,1),:);
    ptLow = cntLocLipLow(indCouples(iCouple,2),:);
    vectLowRight = ptUp - ptLow;
    crossUp = cross(repmat([vectLowRight, 0],nbUp,1), [bsxfun(@minus,cntLocLipUp,ptLow), zeros(nbUp,1)]);
    crossLow = cross(repmat([vectLowRight, 0],nbLow,1), [bsxfun(@minus,cntLocLipLow,ptLow), zeros(nbLow,1)]);
    if all(crossUp(:,3)<=0) & all(crossLow(:,3)<=0)
        iCoupleOK = iCouple;
    end  % if all(crossUp(:,3)<=0) & all(crossLow(:,3)<=0)
end  % for iCouple = 1:size(indCouples,1)

% Indices in the global contours
iULClose = iULVO + ((iULVO<iULVI)*2-1) * (indCouples(iCoupleOK,1) - 1);
iLLClose = iLLVO + ((iLLVO<iLLVI)*2-1) * (indCouples(iCoupleOK,2) - 1);

%---------------------------------------------------------------------
% Fusion Lower Lip - Jaw and Upper Lip - Palate

% Highest point of intersection between the dilated contours for LL-Jaw
l = 0;
ptsInter = [];
% Loop until we are sure at least one intersection exists
while isempty(ptsInter)
    l = l + lDilat;
    % Dilatation of the contours
    cntLLDilat = cntDilateErode(cntLL, l);
    cntJDilat = cntDilateErode(cntJ, l);
    % Intersection
    ptsInter = remDouble(intersection_cnt_cnt_2D(cntLLDilat(1:end-1,:), cntJDilat));
end  % while isempty(ptsInter)
% Highest point
ptInterLow = ptsInter(argmax(ptsInter(:,2)),:);

% Lowest point of intersection between the dilated contours for UL-palate
l = 0;
ptsInter = [];
% Loop until we are sure at least one intersection exists
while isempty(ptsInter)
    l = l + lDilat;
    % Dilatation of the contours
    cntULDilat = cntDilateErode(cntUL, l);
    cntPDilat = cntDilateErode(cntP, l);
    % Intersection
    ptsInter = remDouble(intersection_cnt_cnt_2D(cntULDilat(1:end-1,:), cntPDilat));
end  % while isempty(ptsInter)
% Lowest point
ptInterUp = ptsInter(argmin(ptsInter(:,2)),:);

% Closest point of each contour to the intersection point
iLLFusion = argmin(dist_C(cntLL, ptInterLow')) + indLL(1) - 1;
iJFusion = argmin(dist_C(cntJ, ptInterLow')) + indJ(1) - 1;
iULFusion = argmin(dist_C(cntUL, ptInterUp')) + indUL(1) - 1;
iPFusion = argmin(dist_C(cntP, ptInterUp')) + indP(1) - 1;

%---------------------------------------------------------------------
% Closest point of the glottis point 'GF' on the epiglottis

iGFEpi = argmin(dist_C(cntEpi, GF')) + indEpi(1) - 1;

%---------------------------------------------------------------------
% Closest point of the backLarynx point 'BLPh' on the pharynx

iBLPhPh = argmin(dist_C(cntPha, BLPh')) + indPha(1) - 1;

%---------------------------------------------------------------------
% Closest point of the tongue point 'TE' on the epiglottis

iTEEpi = argmin(dist_C(cntEpi, TE')) + indEpi(1) - 1;

%---------------------------------------------------------------------
% Closest point of the tongue point 'TJUB' on the jaw

iTJUBJaw = argmin(dist_C(cntJaw, TJUB')) + indJaw(1) - 1;


%---------------------------------------------------------------------
% Potential nasopharyngeal closure

% Section of the palate on which we calculate the constriction
ind = iPVT:((iPVT<iPBEU)*2-1):iPBEU;
indVelConstr = ind(round(pVelRemove*length(ind)):end);
cntVel = cntAlign(indVelConstr,:);

% Constriction
D = dist_C(cntVel, cntPha');
[constrNaso, iMin] = min(D(:));
[iMinVel,iMinPha] = ind2sub(size(D),iMin);

% Indices of the constriction in the global contours
iNasoConstrVel = indVelConstr(1) + ((indVelConstr(1)<indVelConstr(end))*2-1) * (iMinVel - 1);
iNasoConstrPha = iMinPha + indPha(1) - 1;

%---------------------------------------------------------------------
% Lower vocal tract contour from larynx to lips

% GFEpi -> TEEpi -> TE -> TJUB -> TJUBJaw -> JFusion -> LLFusion -> LLClose
ind = [iGFEpi:((iGFEpi<iTEEpi)*2-1):iTEEpi,...
    iTE:((iTE<iTJUB)*2-1):iTJUB,...
    iTJUBJaw:((iTJUBJaw<iJFusion)*2-1):iJFusion,...
    iLLFusion:((iLLFusion<iLLClose)*2-1):iLLClose];

cntVctrLow = cntAlign(ind,:);

%---------------------------------------------------------------------
% Upper vocal tract contour from larynx to lips

% If the Velum-Pharynx constriction is bigger than a threshold (no closure)
% BLG -> BLPh -> iBLPhPh -> PhU -> PBEU -> PFusion -> ULFusion -> ULClose

% If the Velum-Pharynx constriction is smaller than a threshold (closure)
% BLG -> BLPh -> iBLPhPh -> NasoConstrPha -> NasoConstrVel -> PFusion -> ULFusion -> ULClose

if constrNaso > constrNaso
    ind = [iBLG:((iBLG<iBLPh)*2-1):iBLPh,...
        iBLPhPh:((iBLPhPh<iPhU)*2-1):iPhU,...
        iPBEU:((iPBEU<iPFusion)*2-1):iPFusion,...
        iULFusion:((iULFusion<iULClose)*2-1):iULClose];
else  % if constrNaso > constrNaso
    ind = [iBLG:((iBLG<iBLPh)*2-1):iBLPh,...
        iBLPhPh:((iBLPhPh<iNasoConstrPha)*2-1):iNasoConstrPha,...
        iNasoConstrVel:((iNasoConstrVel<iPFusion)*2-1):iPFusion,...
        iULFusion:((iULFusion<iULClose)*2-1):iULClose];
end  % if constrNaso > constrNaso

cntVctrUp = cntAlign(ind,:);

%---------------------------------------------------------------------
% Cleaning the contours from the self-intersections

for iCnt = 1:2
    
    switch iCnt
        case 1  % switch iCnt
            cnt = cntVctrLow;
        case 2  % switch iCnt
            cnt = cntVctrUp;
    end  % switch iCnt
    
    % Resampling the contour to be sure there is consecute point too close
    L = sum(distdiff(cnt));
    gapMin = 1e-2;
    nbPtsMax = round(L / gapMin);
    nbPtsResamp = min([round(size(cnt,1)*2), nbPtsMax]);
    [Xspl, Yspl] = bez2spl(cnt(:,1)', cnt(:,2)', nbPtsResamp-1);
    cnt = [Xspl', Yspl'];
    
    % Initialisation of the output contour
    cntOut = cnt(1,:);
    % Initialisation of segment to consider
    iPt = 1;
    iPtPlus1 = 2;
    pt = cnt(1,:);
    % Initialisation of the flag to count or discard the last point of the
    % segment in the final contour
    countPt = 1;
    while dist_C(pt, cnt(end,:)') > 1e-2

        % Segment
        seg = [pt; cnt(iPtPlus1,:)];
        % Rest of the contour
        cntRest1 = cnt(1:iPt,:); % Before the considered segment
        cntRest2 = cnt(iPtPlus1:end,:); % After the considered segment
        % Intersection between the segment and the first rest of the contour
        % First rough quick intersection
        [ptsInter1, ou1] = intersection_courbe(seg, cntRest1);
        % Slower but more precise intersection
        if ~isempty(ptsInter1)
            ptsInter1 = intersection_cnt_cnt_2D(seg, cntRest1);
        end  % if ~isempty(ptsInter)
        % We remove the possible intersection at the connection
        if ~isempty(ptsInter1)
            indKO = find(dist_C(ptsInter1, cntRest1(end,:)')<1e-8);
            ptsInter1(indKO,:) = [];
        end  % if ~isempty(ptsInter)
        % Intersection between the segment and the second rest of the contour
        % First rough quick intersection
        [ptsInter2, ou2] = intersection_courbe(seg, cntRest2);
        % Slower but more precise intersection
        if ~isempty(ptsInter2)
            ptsInter2 = intersection_cnt_cnt_2D(seg, cntRest2);
        end  % if ~isempty(ptsInter)
        % We remove the possible intersection at the connection
        if ~isempty(ptsInter2)
            indKO = find(dist_C(ptsInter2, cntRest2(1,:)')<1e-8);
            ptsInter2(indKO,:) = [];
        end  % if ~isempty(ptsInter)
        % Concatenation of the intersection results
        ptsInter = [ptsInter1; ptsInter2];
        % We finally get rid off intersections equal the first point of the
        % segment (already counted)
        if ~isempty(ptsInter)
            % indKO = find(dist_C(ptsInter, seg(1,:)')<1e-4);
            indKO = find(dist_C(ptsInter, seg(1,:)')<1e-8);
            ptsInter(indKO,:) = [];
        end  % if ~isempty(ptsInter)
        
        % There is no intersection
        if isempty(ptsInter)
            % The second segment point is added to the output contour only if flag is 1
            if countPt
                cntOut = [cntOut; seg(2,:)];
            end  % if countPt
            % Next segment to consider
            pt = seg(2,:);
            iPt = iPt + 1;
            iPtPlus1 = iPtPlus1 + 1;
            
            % There at least one intersection
        else  % if isempty(ptsInter)
            % We consider only the closest intersection to the segment start
            if size(ptsInter,1) > 1
                ptsInter = ptsInter(argmin(dist_C(ptsInter, seg(1,:)')),:);
            end  % if size(ptsInter,1) > 1
            % We add the point to the output contour
            cntOut = [cntOut; ptsInter];
            % The next segment will start at this intersection
            pt = ptsInter;
            % We change the flag from 0 to 1 or 1 to 0
            countPt = mod(countPt+1,2);
            % Adding of a separation in the output contour if there is a break
            if ~countPt
                cntOut = [cntOut; NaN, NaN];
            end  % if ~countPt
            
        end  % if isempty(ptsInter)
        
    end  % while iPt < size(cnt,1)
    
    if isnan(cntOut(end,1))
        cntOut = cntOut(1:end-1,:);
    end  % if isnan(cntOut(end,1))
    
    % Linking of the segments
    
    % Splitting in segments
    [nbZones, indStart, indEnd] = nb_contig_zones(cntOut(:,1)');
    
    % Starting with the first segment
    cntLink = cntOut(indStart(1):indEnd(1),:);
    % Remaining segments to link
    indSegRest = setdiff(1:nbZones,1);
    % Loop on the remaining segments
    while ~isempty(indSegRest)
        % Is the end of the already linked contour connected with the
        % start of a remaining segment?
        D = dist_C(cntLink(end,:), cntOut(indStart(indSegRest),:)');
        indDOK = find(D<1e-3);
        iZone = [];
        if ~isempty(indDOK)
            iZone = indSegRest(indDOK(argmin(D(indDOK))));
        end  % if ~isempty(indDOK)
                
        if ~isempty(iZone)
            % There is indeed a segment to connect
            cntLink = [cntLink; cntOut(indStart(iZone):indEnd(iZone),:)];
            indSegRest = setdiff(indSegRest, iZone);
            
        else  % if ~isempty(iZone)
            % There is no segment to connect, we add the next segment after
            % a NaN
            cntLink = [cntLink; NaN, NaN; cntOut(indStart(indSegRest(1)):indEnd(indSegRest(1)),:)];
            indSegRest = indSegRest(2:end);
        end  % if ~isempty(iZone)
        
    end  % while ~isempty(indSegRest)
    
    switch iCnt
        case 1  % switch iCnt
            cntVctrLowClean = cntLink;
        case 2  % switch iCnt
            cntVctrUpClean = cntLink;
    end  % switch iCnt
end  % for iCnt = 1:2


%% ======================================================================
% 2. Calculation of sagittal distance

% Splitting of each contour into main tract and islands
for iCnt = 1:2
    switch iCnt
        case 1  % switch iCnt
            cnt = cntVctrLowClean;
        case 2  % switch iCnt
            cnt = cntVctrUpClean;
    end  % switch iCnt
    
    [nbZones, indStart, indEnd] = nb_contig_zones(cnt(:,1)');

    % Init
    cntMain = [];
    cntIslands = [];
    for iZone = 1:nbZones
        cntZone = cnt(indStart(iZone):indEnd(iZone),:);
        if isferme(cntZone(:,1)',cntZone(:,2)')
            cntIslands = [cntIslands; cntZone; NaN, NaN];
        else  % if isferme(cnt(indStart(iZone):indEnd(iZone),1)',cnt(indStart(iZone):indEnd(iZone),2)')
            cntMain = cntZone;
        end  % if isferme(cnt(indStart(iZone):indEnd(iZone),1)',cnt(indStart(iZone):indEnd(iZone),2)')
    end  % for iZone = 1:nbZones
    if ~isempty(cntIslands); if isnan(cntIslands(end,1));
        cntIslands = cntIslands(1:end-1,:);
	end; end;  % if ~isempty(cntIslands); if isnan(cntIslands(end,1));
    switch iCnt
        case 1  % switch iCnt
            cntVctrLowMain = cntMain;
            cntVctrLowIslands = cntIslands;
        case 2  % switch iCnt
            cntVctrUpMain = cntMain;
            cntVctrUpIslands = cntIslands;
    end  % switch iCnt
end  % for iCnt = 1:2

% Save
clear VCTR
VCTR.cntMainInt = cntVctrLowMain;
VCTR.cntMainExt = cntVctrUpMain;
VCTR.cntIslandsInt = cntVctrLowIslands;
VCTR.cntIslandsExt = cntVctrUpIslands;

% Rough estimation of the midline: approximation by a circle
try
    cntLow = perform_curve_resampling(cntVctrLowMain, nbPtsVctrMidlineInit, 'nbpts')';
    cntUp = perform_curve_resampling(cntVctrUpMain, nbPtsVctrMidlineInit, 'nbpts')';
    cntMid = (cntLow + cntUp) / 2;
    Circle = circle_least_square(cntMid) ;
catch
    Circle = circle_least_square([cntVctrLowMain; cntVctrUpMain]) ;
end  %  try

% Definition of a resclicing grid based on the circle
vectLow = cntVctrLowMain(1,:) - Circle.centre;
vectUp = cntVctrUpMain(1,:) - Circle.centre;
angleStart = max([atan2(vectLow(2),vectLow(1)), atan2(vectUp(2),vectUp(1))]);
vectLow = cntVctrLowMain(end,:) - Circle.centre;
vectUp = cntVctrUpMain(end,:) - Circle.centre;
angleEnd = min([atan2(vectLow(2),vectLow(1)), atan2(vectUp(2),vectUp(1))]);
% End angles between -pi/2 and 3*pi/2 to ensure continuity
if angleEnd<=-pi/2
    angleEnd = angleEnd + 2*pi;
end  % if angleEnd<=-pi/2
anglesGrdCircle = linspace(angleStart+0.01, angleEnd-0.01, nbGrdCircle);

% Amelioration of the grid at the larynx and the lips based on the circle
% and these known points
% Draw a line starting from the glottis, encompassing the circle and
% finishing at the lips
% Connection point on the circle at the larynx side
ptLar = (cntVctrLowMain(1,:)+cntVctrUpMain(1,:)) / 2;
vectLar = ptLar - Circle.centre;
aLar = atan2(vectLar(2),vectLar(1));
% We add the glottis line only if the glottis is outside the circle
angleConnectLar = NaN;
if dist(Circle.centre, ptLar') > Circle.radius
    aLarAdd = acos(Circle.radius / dist(Circle.centre, ptLar'));
    angleConnectLar = wrapTo2Pi(aLar + aLarAdd);
    if angleConnectLar >= 3*pi/2
        % Angles between -pi/2 and 3*pi/2
        angleConnectLar = angleConnectLar - 2*pi;
    end  % if angleConnectLar >= 3*pi/2
end  % if dist(Circle.Centre, ptLar') > Circle.radius
% Connection point on the circle at the lips side
ptLips = (cntVctrLowMain(end,:)+cntVctrUpMain(end,:)) / 2;
vectLips = ptLips - Circle.centre;
aLips = atan2(vectLips(2),vectLips(1));
if aLips<=-pi/2
    % Angles between -pi/2 and 3*pi/2
    aLips = aLips + 2*pi;
end  % if angleEnd<=-pi/2
% We add the lips line only if the lip point is outside the circle
angleConnectLips = NaN;
if dist(Circle.centre, ptLips') > Circle.radius
    aLipsRem = acos(Circle.radius / dist(Circle.centre, ptLips'));
    angleConnectLips = wrapTo2Pi(aLips - aLipsRem);
    if angleConnectLips >= 3*pi/2
        % Angles between -pi/2 and 3*pi/2
        angleConnectLips = angleConnectLips - 2*pi;
    end  % if angleConnectLar >= 3*pi/2
end  % if dist(Circle.Centre, ptLar') > Circle.radius
% Keep the circle partion located between angleConnectLar and angleConnectLips
if isnan(angleConnectLar)
    angleLow = angleStart+0.01;
else  % if isnan(angleConnectLar)
    angleLow = angleConnectLar;
end  % if isnan(angleConnectLar)
if isnan(angleConnectLips)
    angleUp = angleEnd-0.01;
else  % if isnan(angleConnectLar)
    angleUp = angleConnectLips;
end  % if isnan(angleConnectLar)
t = linspace(angleLow, angleUp, 200);
cirleArc = bsxfun(@plus, Circle.centre, Circle.radius * [cos(t)', sin(t)']);
% Create an estimate of the midline based on the circle and larynx and lips
if isnan(angleConnectLar)
    pt1 = [];
else  % if isnan(angleConnectLar)
    pt1 = ptLar + 0.01 * norme_vecteur(cirleArc(1,:)-ptLar);
end  % if isnan(angleConnectLar)
if isnan(angleConnectLips)
    ptEnd = [];
else  % if isnan(angleConnectLar)
    ptEnd = ptLips + 0.01 * norme_vecteur(cirleArc(end,:)-ptLips);
end  % if isnan(angleConnectLar)
midLineRough = [pt1; cirleArc; ptEnd];
% Resample the grid in equidistant points
midLineCircle = perform_curve_resampling(midLineRough, nbGrdCircle, 'nbpts')';

% Complete the grid with associated angles
vectDiff = midLineCircle(2:end,:) - midLineCircle(1:end-1,:);
vectDiffPerp = [vectDiff(:,2), -vectDiff(:,1)];
vectDiffPerpMid = (vectDiffPerp(1:end-1,:) + vectDiffPerp(2:end,:))/2;
anglesMidLineCircle = [...
    atan2(cntVctrUpMain(1,2)-cntVctrLowMain(1,2), cntVctrUpMain(1,1)-cntVctrLowMain(1,1));...
    atan2(vectDiffPerpMid(:,2),vectDiffPerpMid(:,1));...
    atan2(cntVctrUpMain(end,2)-cntVctrLowMain(end,2), cntVctrUpMain(end,1)-cntVctrLowMain(end,1))];


% Loop on the grid lines
clear GRDCircle
for iGrd = 1:nbGrdCircle
    
    % Grid line
    lineGrdInit = [midLineCircle(iGrd,:); cos(anglesMidLineCircle(iGrd)), sin(anglesMidLineCircle(iGrd))];
    % Projection of the circle center on this grid line
    CGrd = projection_point_droite(Circle.centre,droiteMatrix2Struct(lineGrdInit));
    lineGrd = [CGrd; norme_vecteur(midLineCircle(iGrd,:)-CGrd)];
    
    % Intersections line - vctr
    ptsLowMain = []; ptsUpMain = []; ptsLowIslands = []; ptsUpIslands = [];
    ptsLowMain = remDouble(intersection_line_contour(lineGrd, cntVctrLowMain),1e-6);
    ptsUpMain = remDouble(intersection_line_contour(lineGrd, cntVctrUpMain),1e-6);
    if ~isempty(cntVctrLowIslands)
        ptsLowIslands = remDouble(intersection_line_contour(lineGrd, cntVctrLowIslands));
    end  % if ~isempty(cntVctrLowIslands)
    if ~isempty(cntVctrUpIslands)
        ptsUpIslands = remDouble(intersection_line_contour(lineGrd, cntVctrUpIslands));
    end  % if ~isempty(cntVctrLowIslands)
    ptsLow = [ptsLowMain; ptsLowIslands];
    ptsUp = [ptsUpMain; ptsUpIslands];
    
    % Abcisse of the points on the line starting from grid centre (only 1 dimension)
    dotLow = [];
    dotUp = [];
    if ~isempty(ptsLow)
        dotLow = dot(bsxfun(@minus, ptsLow, CGrd)', repmat(lineGrd(2,:),size(ptsLow,1),1)');
    end  % if ~isempty(ptsLow)
    if ~isempty(ptsUp)
        dotUp = dot(bsxfun(@minus, ptsUp, CGrd)', repmat(lineGrd(2,:),size(ptsUp,1),1)');
    end  % if ~isempty(ptsUp)
        
    % We keep only the points with a positive abcisse on this line
    indLowOK = find(dotLow>0);
    indUpOK = find(dotUp>0);
    
    % Points ordering
    [~, indSortLow] = sort([dotLow(indLowOK)]);
    [~, indSortUp] = sort([dotUp(indUpOK)]);
    ptsLowSort = ptsLow(indLowOK(indSortLow),:);
    ptsUpSort = ptsUp(indUpOK(indSortUp),:);
    if (size(ptsUpSort,1) > 1) & ~mod(size(ptsUpSort,1),2)
        % In the rare case where the line grid intersects an even time the
        % upper vocal tract, we remove the last intersection because it
        % means the vocal tract for this plane id not closed otherwise
        ptsUpSort = ptsUpSort(1:end-1,:);
    end  % if (size(ptsUpSort,1) > 1) & ~mod(size(ptsUpSort,1),2)
    ptsSort = [ptsLowSort; ptsUpSort];
    
    if isempty(ptsLowSort) || isempty(ptsUpSort)
        GRDCircle(iGrd).X = [NaN, NaN];
        GRDCircle(iGrd).Y = [NaN, NaN];
        GRDCircle(iGrd).X2D_grd_Int = NaN;
        GRDCircle(iGrd).Y2D_grd_Int = NaN;
        GRDCircle(iGrd).X2D_grd_Ext = NaN;
        GRDCircle(iGrd).Y2D_grd_Ext = NaN;
        
    else  % if isempty(ptsLow) || isempty(ptsUp)
        
        % Save the grid values
        GRDCircle(iGrd).X = ptsSort(:,1)';
        GRDCircle(iGrd).Y = ptsSort(:,2)';
        GRDCircle(iGrd).X2D_grd_Int = ptsLowSort(1,1)';
        GRDCircle(iGrd).Y2D_grd_Int = ptsLowSort(1,2)';
        GRDCircle(iGrd).X2D_grd_Ext = ptsUpSort(1,1)';
        GRDCircle(iGrd).Y2D_grd_Ext = ptsUpSort(1,2)';
    end  % if isempty(ptsLow) || isempty(ptsUp)
    
end  % for iGrd = 1:nbGrdCircle

% Extension of the grid (with the center point of each grid)
[~,~,GRDCircle] = crossarea_2_tract_2D(GRDCircle);

% First estimation of the midline according to the grid centres
X = [GRDCircle.Xc]';
Y = [GRDCircle.Yc]';
indOK = find(~isnan(X));
midLineEst1 = [X(indOK), Y(indOK)];
% Addition of the first and last grid lines
midLineEst1 = [(cntVctrLowMain(1,:)+cntVctrUpMain(1,:))/2;...
    midLineEst1;...
    (cntVctrLowMain(end,:)+cntVctrUpMain(end,:))/2];

% Smooth
midLineEst = smoothcurv(midLineEst1,2);

% Definition of the grid based on the estimated midLine
% Points of the midline
midLinegrd = perform_curve_resampling(midLineEst, nbGrd, 'nbpts')';
% Angle of each midline segment
diffMidLine = diff(midLinegrd);
angleSegs = atan2(diffMidLine(:,2),diffMidLine(:,1));
% All the angles between -pi/2 and 3*pi/2 to ensure continuity
indAnglesModif = find(angleSegs<=-pi/2);
angleSegs(indAnglesModif) = angleSegs(indAnglesModif) + 2*pi;
% Angle of the normal of the segments
anglesNorms = angleSegs - pi/2;
% Angle of the normal at each point of the midline
anglesGrd = [anglesNorms(1);...
    (anglesNorms(1:end-1) + anglesNorms(2:end)) / 2;...
    anglesNorms(end)];

% Safety contour in order to avoid intersections
% occuring in the oppsite side of the vocal tract
% Safety contour
cntSafety = [cntVctrLowMain(1,:); TE; cntAlign(round((iTE+iTJUB)/2),:); TJLB; cntAlign(iJFusion,:); cntVctrLowMain(end,:)];
% Direction of movement of the safety points
v = norme_vecteur(cntSafety(2:end,:) - cntSafety(1:end-1,:));
vPerp = [-v(:,2), v(:,1)];
dirSafety = [vPerp(1,:); (vPerp(1:end-1,:)+vPerp(2:end,:))/2; vPerp(end,:)];
% Gap to move the points
gapSafety = 0.1;
% As long as there are some points inside the safety contour, we reduce the
% safety contour
cntSafetyClose = [cntSafety; cntSafety(end,1), min(cntSafety(:,2))-1; cntSafety(1,1), min(cntSafety(:,2))-1; cntSafety(1,:)];
IN = find(inpolygon(cntVctrLowClean(:,1),cntVctrLowClean(:,2),cntSafetyClose(:,1),cntSafetyClose(:,2)));
while ~isempty(IN)
    cntSafety = cntSafety + gapSafety * dirSafety;
    cntSafetyClose = [cntSafety; cntSafety(end,1), min(cntSafety(:,2))-1; cntSafety(1,1), min(cntSafety(:,2))-1; cntSafety(1,:)];
    IN = find(inpolygon(cntVctrLowClean(:,1),cntVctrLowClean(:,2),cntSafetyClose(:,1),cntSafetyClose(:,2)));
end  % while ~isempty(IN)
% One more time to me sure
cntSafety = cntSafety + gapSafety * dirSafety;
% Lenghtening at the beginning and at the end
cntSafety(1,:) = cntSafety(1,:) + 5 * norme_vecteur(cntSafety(1,:)-cntSafety(2,:));
cntSafety(end,:) = cntSafety(end,:) + 5 * norme_vecteur(cntSafety(end,:)-cntSafety(end-1,:));

% Loop on the grid lines
clear GRD
for iGrd = 1:nbGrd
    
    % Grid line
    lineGrd = [midLinegrd(iGrd,:); cos(anglesGrd(iGrd)), sin(anglesGrd(iGrd))];
    GRD(iGrd).angleGrd = anglesGrd(iGrd);
    GRD(iGrd).C = midLinegrd(iGrd,:);
    GRD(iGrd).lineGrd = lineGrd;
    
    % Intersections line - vctr
    ptsLowMain = []; ptsUpMain = []; ptsLowIslands = []; ptsUpIslands = [];
    ptsLowMain = remDouble(intersection_line_contour(lineGrd, cntVctrLowMain));
    ptsUpMain = remDouble(intersection_line_contour(lineGrd, cntVctrUpMain));
    if ~isempty(cntVctrLowIslands)
        ptsLowIslands = remDouble(intersection_line_contour(lineGrd, cntVctrLowIslands));
    end  % if ~isempty(cntVctrLowIslands)
    if ~isempty(cntVctrUpIslands)
        ptsUpIslands = remDouble(intersection_line_contour(lineGrd, cntVctrUpIslands));
    end  % if ~isempty(cntVctrLowIslands)
    ptsLow = [ptsLowMain; ptsLowIslands];
    ptsUp = [ptsUpMain; ptsUpIslands];
    
    discardGrd = 0;
    if isempty(ptsLow) || isempty(ptsUp)
        discardGrd = 1; % We just discard the grid if there is no intersection withz the up/low vocal tract
        
    else  % if isempty(ptsLow) || isempty(ptsUp)
        
        % Intersections line - safety contour
        ptRef = remDouble(intersection_line_contour(lineGrd, cntSafety));
        % We keep the closest point
        ptRef = ptRef(argmin(dist_C(ptRef, lineGrd(1,:)')),:);
        
        % Abcisse of the points on the line starting from ref point (only 1 dimension)
        vectRef = lineGrd(2,:);
        dotLow = dot(bsxfun(@minus, ptsLow, ptRef)', repmat(vectRef,size(ptsLow,1),1)');
        dotUp = dot(bsxfun(@minus, ptsUp, ptRef)', repmat(vectRef,size(ptsUp,1),1)');
        
        % We keep only the points with a positive abcisse on this line
        indLowOK = find(dotLow>0);
        indUpOK = find(dotUp>0);
        
        % We test again with the remaining points and we discard the grid
        % if there is no intersection withz the up/low vocal tract
        if isempty(indLowOK) || isempty(indUpOK)
            discardGrd = 1;
        else  % if isempty(indLowOK) || isempty(indUpOK)
            % Points ordering
            [~, indSortLow] = sort([dotLow(indLowOK)]);
            [~, indSortUp] = sort([dotUp(indUpOK)]);
            ptsLowSort = ptsLow(indLowOK(indSortLow),:);
            ptsUpSort = ptsUp(indUpOK(indSortUp),:);
            if (size(ptsUpSort,1) > 1) & ~mod(size(ptsUpSort,1),2)
                % In the rare case where the line grid intersects an even time the
                % upper vocal tract, we remove the last intersection because it
                % means the vocal tract for this plane is not closed otherwise
                ptsUp(indUpOK(indSortUp(end)),:) = [];
                % Update the values that have consequently changed
                dotUp = dot(bsxfun(@minus, ptsUp, ptRef)', repmat(vectRef,size(ptsUp,1),1)');
                indUpOK = find(dotUp>0);
                [~, indSortUp] = sort([dotUp(indUpOK)]);
                ptsUpSort = ptsUp(indUpOK(indSortUp),:);
            end  % if (size(ptsUpSort,1) > 1) & ~mod(size(ptsUpSort,1),2)
            ptsSort = [ptsLowSort; ptsUpSort];
            
            % If one point of the upper vocal tract is "lower" than any  point of the lower
            % vocal tract, there is intersection of the 2 contours. In that
            % case we define the grid as a single point , in the middle of
            % these points (will be counted later as a zero area in the
            % 'crossarea_2_tract_2D' functionl
            if dotLow(indLowOK(indSortLow(end))) >= dotUp(indUpOK(indSortUp(1)))
                ptsLowSort = ptsLow(indLowOK(indSortLow(end)),:);
                ptsUpSort = ptsUp(indUpOK(indSortUp(1)),:);
                ptsSort = (ptsLowSort + ptsUpSort) / 2;
            end  % if dotLow(indLowOK(indSortLow(end))) >= dotUp(indUpOK(indSortUp(1)))
            
            if mod(size(ptsSort,1),2) & (size(ptsSort,1) > 1)
                discardGrd = 1;
            end  % if mod(size(ptsSort,1),2) & (size(ptsSort,1) > 1)
            
            % Save the grid values
            GRD(iGrd).X = ptsSort(:,1)';
            GRD(iGrd).Y = ptsSort(:,2)';
            GRD(iGrd).X2D_grd_Int = ptsLowSort(1,1)';
            GRD(iGrd).Y2D_grd_Int = ptsLowSort(1,2)';
            GRD(iGrd).X2D_grd_Ext = ptsUpSort(1,1)';
            GRD(iGrd).Y2D_grd_Ext = ptsUpSort(1,2)';
            
        end  % if isempty(indLowOK) || isempty(indUpOK)
        
    end  % if isempty(ptsLow) || isempty(ptsUp)
    
    if discardGrd
        GRD(iGrd).X = [NaN, NaN];
        GRD(iGrd).Y = [NaN, NaN];
        GRD(iGrd).X2D_grd_Int = NaN;
        GRD(iGrd).Y2D_grd_Int = NaN;
        GRD(iGrd).X2D_grd_Ext = NaN;
        GRD(iGrd).Y2D_grd_Ext = NaN;
    end  % if discardGrd
    
end  % for iGrd = 1:nbGrd

% We force the first and last grid lines
clear GRD1 GRDend
% First
GRD1.angleGrd = atan2(cntVctrUpMain(1,2)-cntVctrLowMain(1,2), cntVctrUpMain(1,1)-cntVctrLowMain(1,1));
GRD1.C = mean([cntVctrLowMain(1,:); cntVctrUpMain(1,:)]);
GRD1.lineGrd = [GRD1.C; cos(GRD1.angleGrd), sin(GRD1.angleGrd)];
GRD1.X = [cntVctrLowMain(1,1), cntVctrUpMain(1,1)];
GRD1.Y = [cntVctrLowMain(1,2), cntVctrUpMain(1,2)];
GRD1.X2D_grd_Int = GRD1.X(1);
GRD1.Y2D_grd_Int = GRD1.Y(1);
GRD1.X2D_grd_Ext = GRD1.X(2);
GRD1.Y2D_grd_Ext = GRD1.Y(2);
% Last
GRDend.angleGrd = atan2(cntVctrUpMain(end,2)-cntVctrLowMain(end,2), cntVctrUpMain(end,1)-cntVctrLowMain(end,1));
GRDend.C = mean([cntVctrLowMain(end,:); cntVctrUpMain(end,:)]);
GRDend.lineGrd = [GRDend.C; cos(GRDend.angleGrd), sin(GRDend.angleGrd)];
GRDend.X = [cntVctrLowMain(end,1), cntVctrUpMain(end,1)];
GRDend.Y = [cntVctrLowMain(end,2), cntVctrUpMain(end,2)];
GRDend.X2D_grd_Int = GRDend.X(1);
GRDend.Y2D_grd_Int = GRDend.Y(1);
GRDend.X2D_grd_Ext = GRDend.X(2);
GRDend.Y2D_grd_Ext = GRDend.Y(2);

% We verify that the first and last grid lines are not too close from the
% next/previous ones
% Min observed distance
mD = 1000;
for iGRD = 1:length(GRD)-1
    mD = min([mD, distance_points_droite(GRD(iGRD+1).lineGrd, GRD(iGRD).C)]);
end  % for ii = 1:length(GRD)-1
is1TooClose = 0;
if distance_points_droite(GRD1.lineGrd, GRD(1).C) < (mD * pDistLine)
    is1TooClose = 1;
end  % if distance_points_droite(GRD1.lineGrd, GRD(1).C)
isEndTooClose = 0;
if distance_points_droite(GRDend.lineGrd, GRD(end).C) < (mD * pDistLine)
    isEndTooClose = 1;
end  % if distance_points_droite(GRD1.lineGrd, GRD(1).C)
    
% Fusion
GRDtmp = GRD;
clear GRD
if ~is1TooClose
    GRD(1) = GRD1;
    GRD(2:length(GRDtmp)+1) = GRDtmp;
else  % if ~is1TooClose
    GRD(1) = GRD1;
    GRD(2:length(GRDtmp)) = GRDtmp(2:end);
end  % if ~is1TooClose
if ~isEndTooClose
    GRD(length(GRD)+1) = GRDend;
else  % if ~isEndTooClose
    GRD(length(GRD)) = GRDend;
end  % if ~isEndTooClose

% Function of sagittal distance
indGrdOK = find(~isnan([GRD.X2D_grd_Int]));
[DistSag, Long, GRDVctr] = crossarea_2_tract_2D(GRD(indGrdOK));
% NaN sagittal distances are in fact equal to zero (occlusion)
DistSag(find(isnan(DistSag))) = 0;
% Sagittal function
SF = DistSag;

% We close potential Warning messages boxes opened in 'crossarea_2_tract_2D'
allHandle = allchild(0);
allTag = get(allHandle, 'Tag');
isMsgbox = strncmp(allTag, 'Msgbox_', 7);
delete(allHandle(isMsgbox));


%% ======================================================================
% 3. Calcul of the area function according to the alpha-beta model
% From Soquet 2002

%-----------------------------------------------------------
% Split of the vocal tract into regions

% Larynx region
% We project PhL on the line grid and see when the direction of projection
% change
ptRef = PhL;
iGrd = 0;
valTest = 1;
while valTest > 0
    iGrd = iGrd + 1;
    lineGrd = GRDVctr(iGrd).lineGrd;
    ptRefproj = projection_point_droite(ptRef,droiteMatrix2Struct(lineGrd)) ;
    vectCross = cross([lineGrd(2,:),0], [ptRef-ptRefproj,0]);
    valTest = vectCross(3);
end  % while valTest > 0
iGrdLarMax = max([1, iGrd-1]);

% Low pharynx region
% We project ET on the line grid and see when the direction of projection
% change
ptRef = ET;
iGrd = iGrdLarMax;
valTest = 1;
while valTest > 0
    iGrd = iGrd + 1;
    lineGrd = GRDVctr(iGrd).lineGrd;
    ptRefproj = projection_point_droite(ptRef,droiteMatrix2Struct(lineGrd)) ;
    vectCross = cross([lineGrd(2,:),0], [ptRef-ptRefproj,0]);
    valTest = vectCross(3);
end  % while valTest > 0
iGrdLowPhaMax = iGrd - 1;

% Oro-pharynx region
% We project PVT on the line grid and see when the direction of projection
% change
ptRef = PVT;
iGrd = iGrdLowPhaMax;
valTest = 1;
while valTest > 0
    iGrd = iGrd + 1;
    lineGrd = GRDVctr(iGrd).lineGrd;
    ptRefproj = projection_point_droite(ptRef,droiteMatrix2Struct(lineGrd)) ;
    vectCross = cross([lineGrd(2,:),0], [ptRef-ptRefproj,0]);
    valTest = vectCross(3);
end  % while valTest > 0
iGrdOroPhaMax = iGrd - 1;

% Mid-pharynx region
% Middle between Low pharynx and oro pharynx
indGrdPha = iGrdLowPhaMax+1:iGrdOroPhaMax;
linePharynx = [struct2vect(GRDVctr, indGrdPha, 'C', 1)', struct2vect(GRDVctr, indGrdPha, 'C', 2)'];
DD = distdiff(linePharynx);
iGrdMidPhaMax = min(find((cumsum(DD)-(sum(DD)/2))>=0)) + iGrdLowPhaMax;

% Velum Region
% We project PBEL on the line grid and see when the direction of projection
% change
ptRef = PBEL;
iGrd = iGrdOroPhaMax;
valTest = 1;
while valTest > 0
    iGrd = iGrd + 1;
    lineGrd = GRDVctr(iGrd).lineGrd;
    ptRefproj = projection_point_droite(ptRef,droiteMatrix2Struct(lineGrd)) ;
    vectCross = cross([lineGrd(2,:),0], [ptRef-ptRefproj,0]);
    valTest = vectCross(3);
end  % while valTest > 0
iGrdVelMax = iGrd - 1;

% Hard palate region
% We project the middle of PUT and PSS on the line grid and see when the direction of projection
% change
ptRef = (PUT+PSS)/2;
iGrd = iGrdVelMax;
valTest = 1;
while valTest > 0
    iGrd = iGrd + 1;
    lineGrd = GRDVctr(iGrd).lineGrd;
    ptRefproj = projection_point_droite(ptRef,droiteMatrix2Struct(lineGrd)) ;
    vectCross = cross([lineGrd(2,:),0], [ptRef-ptRefproj,0]);
    valTest = vectCross(3);
end  % while valTest > 0
iGrdHardPalMax = iGrd - 1;

% Alevolar region
% We project PUT on the line grid and see when the direction of projection
% change
ptRef = PUT;
iGrd = iGrdHardPalMax;
valTest = 1;
while valTest > 0
    iGrd = iGrd + 1;
    lineGrd = GRDVctr(iGrd).lineGrd;
    ptRefproj = projection_point_droite(ptRef,droiteMatrix2Struct(lineGrd)) ;
    vectCross = cross([lineGrd(2,:),0], [ptRef-ptRefproj,0]);
    valTest = vectCross(3);
end  % while valTest > 0
iGrdAlvMax = iGrd - 1;


% Indices of the grids for each region
lastGrd = 0;
indGrdLar = 1:iGrdLarMax;
lastGrd = max([lastGrd, indGrdLar]);
indGrdLowPha = lastGrd+1:iGrdLowPhaMax;
lastGrd = max([lastGrd, indGrdLowPha]);
indGrdMidPha = lastGrd+1:iGrdMidPhaMax;
lastGrd = max([lastGrd, indGrdMidPha]);
indGrdOroPha = lastGrd+1:iGrdOroPhaMax;
lastGrd = max([lastGrd, indGrdOroPha]);
indGrdVel = lastGrd+1:iGrdVelMax;
lastGrd = max([lastGrd, indGrdVel]);
indGrdHardPal = lastGrd+1:iGrdHardPalMax;
lastGrd = max([lastGrd, indGrdHardPal]);
indGrdAlv = lastGrd+1:iGrdAlvMax;
lastGrd = max([lastGrd, indGrdAlv]);
indGrdLab = lastGrd+1:length(GRDVctr);
indGrdReg = [...
    iLar*ones(length(indGrdLar),1);...
    iLowPha*ones(length(indGrdLowPha),1);...
    iMidPha*ones(length(indGrdMidPha),1);...
    iOroPha*ones(length(indGrdOroPha),1);...
    iVel*ones(length(indGrdVel),1);...
    iHardPal*ones(length(indGrdHardPal),1);...
    iAlv*ones(length(indGrdAlv),1);...
    iLab*ones(length(indGrdLab),1)];

% Corrected distance of each grid
distGrdCorr = [GRDVctr.D] .* [GRDVctr.cosalp];
% NaN distances are in fact equal to zero (occlusion)
indGrdOccl = find(isnan(distGrdCorr));
distGrdCorr(indGrdOccl) = 0;

% Estimated area for each grid with alpha-beta model
vectAlpha = [...
    repmat(alp(iLar),1,length(find(indGrdReg==iLar))),...
    repmat(alp(iLowPha),1,length(find(indGrdReg==iLowPha))),...
    repmat(alp(iMidPha),1,length(find(indGrdReg==iMidPha))),...
    repmat(alp(iOroPha),1,length(find(indGrdReg==iOroPha))),...
    repmat(alp(iVel),1,length(find(indGrdReg==iVel))),...
    repmat(alp(iHardPal),1,length(find(indGrdReg==iHardPal))),...
    repmat(alp(iAlv),1,length(find(indGrdReg==iAlv))),...
    repmat(alp(iLab),1,length(find(indGrdReg==iLab)))];
vectBeta = [...
    repmat(bet(iLar),1,length(find(indGrdReg==iLar))),...
    repmat(bet(iLowPha),1,length(find(indGrdReg==iLowPha))),...
    repmat(bet(iMidPha),1,length(find(indGrdReg==iMidPha))),...
    repmat(bet(iOroPha),1,length(find(indGrdReg==iOroPha))),...
    repmat(bet(iVel),1,length(find(indGrdReg==iVel))),...
    repmat(bet(iHardPal),1,length(find(indGrdReg==iHardPal))),...
    repmat(bet(iAlv),1,length(find(indGrdReg==iAlv))),...
    repmat(bet(iLab),1,length(find(indGrdReg==iLab)))];
areaGrd = vectAlpha .* (distGrdCorr.^vectBeta);

% Area function
AF = mean([areaGrd(1:end-1); areaGrd(2:end)]);
% An occlusion must stay an occlusion and not lost by the averaging
indAFOccl = unique([indGrdOccl'-1, indGrdOccl']);
indAFOccl(find(indAFOccl==0)) = [];
indAFOccl(find(indAFOccl==length(areaGrd))) = [];
AF(indAFOccl) = 0;
% Imposed minimal constriction
AF = max([AF; constrMin*ones(size(AF))]);

%% ======================================================================
% 4. 4. Calcuation of the transfer function and the formants

% Removing the nasal part
clear AF_NASAL_no
AF_NASAL_no.A_mid = NaN;
AF_NASAL_no.A_left = NaN;
AF_NASAL_no.A_right = NaN;

% Clean area function
[AFclean, Lengthclean] = clean_area_function(AF, Long);
clear ORAL_cavites;
ORAL_cavites.A = AFclean;
ORAL_cavites.L = Lengthclean;

% Compute acoustics
[H, F] = vtn2frm_ftr_wall_vibration(ORAL_cavites, AF_NASAL_no);

return



