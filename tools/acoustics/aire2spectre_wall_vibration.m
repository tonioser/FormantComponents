function [H_eq, F_form, H_partiel, Z_nasal_in, Z_bucal_in] = aire2spectre_wall_vibration(ORAL, NASAL, nbfreq, Fmax, Fmin);

% function [H] = aire2spectre_wall_vibration(ORAL, NASAL, nbFreq, Fmax, Fmin);
%
% Returns the acoustic transfer function given the area functions of the
% oral and nasal tracts in input
%
%
% Inputs
%   ORAL      : Structure containing the parameters of the oral area function
%       ORAL.A(nbOralTubes)         : Areas of the tubes in cm2 (glottis -> lips)
%       ORAL.L(nbOralTubes)	        : Lengths of the tubes in cm (glottis -> lips)
%
%   NASAL     : Structure containing the parameters of the nasal area function
%       NASAL.A_mid(nbNasophaTubes) : Areas of the tubes of the nasopharyngeal port in cm2
%       NASAL.L_mid(nbNasophaTubes) : Lengths of the tubes of the nasopharyngeal port in cm
%       NASAL.ind_connect_mid       : Index in the oral area function of the connection of the nasopharyngeal port
%
%       NASAL.A_left(nbNasLeft)     : Areas of the tubes of the left nasal cavity in cm2
%       NASAL.L_left(nbNasLeft)     : Lengths of the tubes of the left nasal cavity in cm
%       NASAL.A_right(nbNasRight)   : Areas of the tubes of the right nasal cavity in cm2
%       NASAL.L_right(nbNasRight)   : Lengths of the tubes of the right nasal cavity in cm
%
%       NASAL.A_sinleft(nbSinLeft)  : Areas of the tubes of the left maxillary sinus in cm2
%       NASAL.L_sinleft(nbSinLeft)  : Lengths of the tubes of the left maxillary sinus in cm
%       NASAL.A_sinright(nbSinRight): Areas of the tubes of the right maxillary sinus in cm2
%       NASAL.L_sinright(nbSinRight): Lengths of the tubes of the right maxillary sinus in cm
%       NASAL.ind_connect_sinleft   : Index in the left nasal cavity area function of the connection of the left maxillary sinus
%       NASAL.ind_connect_sinright  : Index in the left right cavity area function of the connection of the right maxillary sinus
%   nbFreq(1) : Number of frequency points calculated (optional, [default = 256])
%   Fmax(1)   : Maximal frequency calculated in Hz (optional, [default = 5000])
%   Fmin(1)   : Minimal frequency calculated in Hz (optional, [default = Fmax / nbfreq])
%
%
% Outputs
%     H(nbFreq) : Transfer function (amplitude plot: '20*log10(abs(H))')
%
% Author: Antoine Serrurier
% Date: 19/12/2005


%===================================================================
% Inputs and constants

if nargin == 2
    nbfreq = 256;
    Fmax = 5000;
end;
if nargin == 3
    Fmax = 5000;
end;

global CONST_DAT

% Mrayati constants
c = 35100 ; % cm/s
rho = 1.14e-3 ; % g/cm^3
lambda = 5.5e-5 ;
eta = 1.4 ;
mu = 1.86e-4 ;
cp = 0.24 ;
bp = 1600 ; % dynes.s/cm
mp = 1.4 ; % g/cm²

CONST_DAT=[c ; rho ; lambda ; eta ; mu ; cp ; bp ; mp];

% Frequencies
f = linspace(Fmax/nbfreq, Fmax, nbfreq);
if nargin == 5
    f = linspace(Fmin, Fmax, nbfreq);
end  % if nargin == 5
w=2*pi*f;

%===================================================================
% Left nasal cavity transfer function

% Radiation impedance
Zr_nasal_left = rho/(2*pi*c)*(w.^2)+j*8*rho/(3*pi.*sqrt( pi*NASAL.A_left(end))) * w;

% Transfer function
if any(isnan(NASAL.A_left))
    H_nasal_left = NaN * ones(size(f)); Y_nasal_left = NaN * ones(size(f));
elseif any(~NASAL.A_left)
    H_nasal_left = zeros(size(f)); Y_nasal_left = zeros(size(f));
else
    % With sinus
    if isfield(NASAL,'A_sinleft')
        if all(~isnan(NASAL.A_sinleft))
            % Transfer function of the cavity from the sinus to the end
            ind_OK = NASAL.ind_connect_sinleft:length(NASAL.A_left);
            [H_nasal_outleft, Y_nasal_outleft] = spectrelec_wall_vibration(w, NASAL.A_left(ind_OK)', Zr_nasal_left, NASAL.L_left(ind_OK)',1);
            % Transfer function of the sinus
            [H_sinleft, Y_sinleft] = spectrelec_wall_vibration(w, NASAL.A_sinleft', +Inf*ones(size(Zr_nasal_left)), NASAL.L_sinleft',1);
            if all(isnan(Y_sinleft)); Y_sinleft = 0; end
            % Resulting impedance at the junction
            Zr_nasal_out_sin_left = 1 ./ (Y_nasal_outleft + Y_sinleft);
            % Transfer function from the cavum to the junction
            ind_OK = 1:NASAL.ind_connect_sinleft-1;
            [H_nasal_leftin, Y_nasal_left] = spectrelec_wall_vibration(w, NASAL.A_left(ind_OK)', Zr_nasal_out_sin_left, NASAL.L_left(ind_OK)',1);
            % Overall transfer function
            H_nasal_left = H_nasal_leftin .* Y_nasal_outleft ./ (Y_nasal_outleft + Y_sinleft) .* H_nasal_outleft;
        else  % if all(~isnan(NASAL.A_sinleft))
            [H_nasal_left, Y_nasal_left] = spectrelec_wall_vibration(w, NASAL.A_left', Zr_nasal_left, NASAL.L_left',1);
        end  % if all(~isnan(NASAL.A_sinleft))
    	% Without sinus
    else  % if isfield(NASAL,'A_sinleft')
        [H_nasal_left, Y_nasal_left] = spectrelec_wall_vibration(w, NASAL.A_left', Zr_nasal_left, NASAL.L_left',1);
    end  % if isfield(NASAL,'A_sinleft')
end  % if any(~A_nasal_mid(ind_OK))

%===================================================================
% Right nasal cavity transfer function

% Radiation impedance
Zr_nasal_right = rho/(2*pi*c)*(w.^2)+j*8*rho/(3*pi.*sqrt( pi*NASAL.A_right(end))) * w;

% Transfer function
if any(isnan(NASAL.A_right))
    H_nasal_right = NaN * ones(size(f)); Y_nasal_right = NaN * ones(size(f));
elseif any(~NASAL.A_right)
    H_nasal_right = zeros(size(f)); Y_nasal_right = zeros(size(f));
else
    % With sinus
    if isfield(NASAL,'A_sinright')
        if all(~isnan(NASAL.A_sinright))
            % Transfer function of the cavity from the sinus to the end
            ind_OK = NASAL.ind_connect_sinright:length(NASAL.A_right);
            [H_nasal_outright, Y_nasal_outright] = spectrelec_wall_vibration(w, NASAL.A_right(ind_OK)', Zr_nasal_right, NASAL.L_right(ind_OK)',1);
            % Transfer function of the sinus
            [H_sinright, Y_sinright] = spectrelec_wall_vibration(w, NASAL.A_sinright', +Inf*ones(size(Zr_nasal_right)), NASAL.L_sinright',1);
            if all(isnan(Y_sinright)); Y_sinright = 0; end
            % Resulting impedance at the junction
            Zr_nasal_out_sin_right = 1 ./ (Y_nasal_outright + Y_sinright);
            % Transfer function from the cavum to the junction
            ind_OK = 1:NASAL.ind_connect_sinright-1;
            [H_nasal_rightin, Y_nasal_right] = spectrelec_wall_vibration(w, NASAL.A_right(ind_OK)', Zr_nasal_out_sin_right, NASAL.L_right(ind_OK)',1);
            % Overall transfer function
            H_nasal_right = H_nasal_rightin .* Y_nasal_outright ./ (Y_nasal_outright + Y_sinright) .* H_nasal_outright;
        else  % if all(~isnan(NASAL.A_sinright))
            [H_nasal_right, Y_nasal_right] = spectrelec_wall_vibration(w, NASAL.A_right', Zr_nasal_right, NASAL.L_right',1);
        end  % if all(~isnan(NASAL.A_sinright))
    	% Without sinus
    else  % if isfield(NASAL,'A_sinright')
        [H_nasal_right, Y_nasal_right] = spectrelec_wall_vibration(w, NASAL.A_right', Zr_nasal_right, NASAL.L_right',1);
    end  % if isfield(NASAL,'A_sinright')
end  % if any(~A_nasal_mid(ind_OK))

%===================================================================
% Nasoparyngeal port transfer function

% Radiation impedance
switch num2str([any(isnan(Y_nasal_left)); any(isnan(Y_nasal_right))])'
    case '00'
        Zr_nasal_mid =  1 ./ (Y_nasal_left + Y_nasal_right);
    case '01'
        Zr_nasal_mid =  1 ./ Y_nasal_left;
    case '10'
        Zr_nasal_mid =  1 ./ Y_nasal_right;
    case '11'
        Zr_nasal_mid = rho/(2*pi*c)*(w.^2)+j*8*rho/(3*pi.*sqrt( pi*NASAL.A_mid(end))) * w;
end  % switch num2str([any(isnan(Y_nasal_left)); any(isnan(Y_nasal_right))])'

% Transfer function
if any(isnan(NASAL.A_mid))
    H_nasal_mid = NaN * ones(size(f)); Y_nasal_mid = NaN * ones(size(f));
elseif any(~NASAL.A_mid)
    H_nasal_mid = zeros(size(f)); Y_nasal_mid = zeros(size(f));
else
    [H_nasal_mid, Y_nasal_mid] = spectrelec_wall_vibration(w, NASAL.A_mid', Zr_nasal_mid, NASAL.L_mid',1);
end  % if any(~A_nasal_mid(ind_OK))

%===================================================================
% Mouth transfer function (= dfrom the junction with the nasal tract to the lips)

% Case if there is no junction with the nasal tract
if ~isfield(NASAL, 'ind_connect_mid') % Si pas de connection
    if isfield(ORAL, 'ind_connect_sinleft') & isfield(ORAL, 'ind_connect_sinright')
        NASAL.ind_connect_mid = max([ORAL.ind_connect_sinleft, ORAL.ind_connect_sinright]) + 1;
    elseif isfield(ORAL, 'ind_connect_sinleft') & ~isfield(ORAL, 'ind_connect_sinright')
        NASAL.ind_connect_mid = ORAL.ind_connect_sinleft + 1;
    elseif ~isfield(ORAL, 'ind_connect_sinleft') & isfield(ORAL, 'ind_connect_sinright')
        NASAL.ind_connect_mid = ORAL.ind_connect_sinright + 1;
    else
        NASAL.ind_connect_mid = 1;
    end  % if isfield(ORAL, 'ind_connect_sinleft') & isfield(ORAL, 'ind_connect_sinright')
end  % if ~isfield(NASAL, 'ind_connect_mid')

% Indices from the junction with the nasal tract to the lips
ind_tuyau_bucal_tmp = NASAL.ind_connect_mid:length(ORAL.A);
% Non-NaN tubes only
ind_tuyau_bucal_OK = ind_tuyau_bucal_tmp(find(~isnan(ORAL.A(ind_tuyau_bucal_tmp))));

% Transfer function
if any(isnan(ORAL.A(ind_tuyau_bucal_OK))) || isempty(ind_tuyau_bucal_OK)
    H_bucal = NaN * ones(size(f)); Y_bucal = NaN * ones(size(f));
elseif any(~ORAL.A(ind_tuyau_bucal_OK))
    H_bucal = zeros(size(f)); Y_bucal = zeros(size(f));
else
    % Radiation impedance
    Zr_bucal =  rho/(2*pi*c)*(w.^2)+j*8*rho/(3*pi.*sqrt(pi*ORAL.A(ind_tuyau_bucal_OK(end))))*w;
    [H_bucal, Y_bucal] = spectrelec_wall_vibration(w, ORAL.A(ind_tuyau_bucal_OK)', Zr_bucal, ORAL.L(ind_tuyau_bucal_OK)');
end  % if any(~A_nasal_mid(ind_OK))

%===================================================================
% Pharyngeal cavity transfer function

% Indices of the pharynx
ind_tuyau_pha_OK = [];
if ~isempty(ind_tuyau_bucal_OK)
    ind_tuyau_pha_OK = find(~isnan(ORAL.A(1:ind_tuyau_bucal_OK(1)-1)));
end  % if ~isempty(ind_tuyau_bucal_OK)

% Transfer function
if any(isnan(ORAL.A(ind_tuyau_pha_OK))) || isempty(ind_tuyau_pha_OK)
    H_pha = NaN * ones(size(f)); Y_pha = NaN * ones(size(f));
elseif any(~ORAL.A(ind_tuyau_pha_OK))
    H_pha = zeros(size(f)); Y_pha = zeros(size(f));
else
    % Radiation impedance
    switch num2str([any(isnan(Y_nasal_mid)); any(isnan(Y_bucal))])'
        case '00'
            Zr_pha =  1 ./ (Y_nasal_mid + Y_bucal);
        case '01'
            Zr_pha =  1 ./ Y_nasal_mid;
        case '10'
            Zr_pha =  1 ./ Y_bucal;
        case '11'
            Zr_pha = rho/(2*pi*c)*(w.^2)+j*8*rho/(3*pi.*sqrt( pi*ORAL.A(ind_tuyau_pha_OK(end)))) * w;
    end  % switch num2str([any(isnan(Y_nasal_left)); any(isnan(Y_nasal_right))])'

    %-----------------------------------------------------------
    % Piriformis sinus

    % Create them artificially if they do not exist
    if ~isfield(ORAL,'A_sinleft')
        ORAL.A_sinleft = 0; ORAL.L_sinleft = 1; ORAL.ind_connect_sinleft = 1;
    end  % if ~isfield(ORAL,'A_sinleft')
    if ~isfield(ORAL,'A_sinright')
        ORAL.A_sinright = 0; ORAL.L_sinright = 1; ORAL.ind_connect_sinright = 1;
    end  % if ~isfield(ORAL,'A_sinright')

    % If there is a NaN, put it to 0 and close the sinus at the input
    if any(isnan(ORAL.A_sinleft))
        A_sinleft = repl_val(ORAL.A_sinleft, NaN, 0); ORAL.A_sinleft = A_sinleft;
        L_sinleft = repl_val(ORAL.L_sinleft, NaN, 1); ORAL.L_sinleft = L_sinleft;
        ORAL.A_sinleft(1) = 0;
    end  % if any(isnan(ORAL.A_sinleft))
    if any(isnan(ORAL.A_sinright))
        A_sinright = repl_val(ORAL.A_sinright, NaN, 0); ORAL.A_sinright = A_sinright;
        L_sinright = repl_val(ORAL.L_sinright, NaN, 1); ORAL.L_sinright = L_sinright;
        ORAL.A_sinright(1) = 0;
    end  % if any(isnan(ORAL.A_sinright))

    % Sorting of the sinus upstream/downstream
    [ind_connect_amont, ind_LR_1st] = max([ORAL.ind_connect_sinleft, ORAL.ind_connect_sinright]);
    ind_connect_aval = min([ORAL.ind_connect_sinleft, ORAL.ind_connect_sinright]);

    % Transfer function of the pharynx from the upstream sinus to junction with nasal tract
    ind_OK = ind_connect_amont:length(ind_tuyau_pha_OK);
    [H_pha_outamont, Y_pha_outamont] = spectrelec_wall_vibration(w, ORAL.A(ind_OK)', Zr_pha, ORAL.L(ind_OK)');

    % Transfer function of the upstream sinus
    if ind_LR_1st == 1;
        A_sinamont = ORAL.A_sinleft; L_sinamont = ORAL.L_sinleft;
        A_sinaval = ORAL.A_sinright; L_sinaval = ORAL.L_sinright;
    else  % if ind_LR_1st == 1;
        A_sinamont = ORAL.A_sinright; L_sinamont = ORAL.L_sinright;
        A_sinaval = ORAL.A_sinleft; L_sinaval = ORAL.L_sinleft;
    end  % if ind_LR_1st == 1;
    [H_sinamont, Y_sinamont] = spectrelec_wall_vibration(w, A_sinamont', +Inf*ones(size(Zr_pha)), L_sinamont');
    if all(isnan(Y_sinamont)); Y_sinamont = 0; end

    % Impedance at the junction
    Zr_pha_sinamont = 1 ./ (Y_pha_outamont + Y_sinamont);
    % Transfer function
    H_pha_sinamont = H_pha_outamont .* Y_pha_outamont ./ (Y_pha_outamont + Y_sinamont);

    % Transfer function of the pharynx from the downstream sinus to the upstream sinus
    ind_OK = ind_connect_aval:ind_connect_amont-1;
    % Check if 2 tubes connected at the same location
    if isempty(ind_OK)
        Y_pha_amontaval = Y_pha_outamont + Y_sinamont;
        H_pha_amontaval_loc = 1;
    else  % if isempty(ind_OK)
        [H_pha_amontaval_loc, Y_pha_amontaval] = spectrelec_wall_vibration(w, ORAL.A(ind_OK)', Zr_pha_sinamont, ORAL.L(ind_OK)');
    end  % if isempty(ind_OK)

    % Transfer function of the downstream sinus
    [H_sinaval, Y_sinaval] = spectrelec_wall_vibration(w, A_sinaval', +Inf*ones(size(Zr_pha)), L_sinaval');
    if all(isnan(Y_sinaval)); Y_sinaval = 0; end

    % Impedance at the junction
    Zr_pha_sinaval = 1 ./ (Y_pha_amontaval + Y_sinaval);
    % Transfer function
    H_pha_sinaval = H_pha_sinamont .* H_pha_amontaval_loc .* Y_pha_amontaval ./ (Y_pha_amontaval + Y_sinaval);

    % Transfer function of the pharynx from the glottis to the downstream sinus
    ind_OK = 1:ind_connect_aval-1;
    % Check if there is nothing downstream
    if isempty(ind_OK)
        Y_pha = 1 ./ Zr_pha_sinaval;
        H_pha = H_pha_sinaval;
    else  % if isempty(ind_OK)
        [H_pha_glotte, Y_pha] = spectrelec_wall_vibration(w, ORAL.A(ind_OK)', Zr_pha_sinaval, ORAL.L(ind_OK)');
        H_pha = H_pha_glotte .* H_pha_sinaval;
    end  % if isempty(ind_OK)
    %-----------------------------------------------------------

end  % if any(~A_nasal_mid(ind_OK))

%***************************************
while 0
    figure
    clf; hold on; grid on;
    plot(f, 20*log10(abs(H_pha)))
end  % while 0
%***************************************

%===================================================================
% Overall transfer function

% Depends on which sections exist
switch num2str([any(isnan(H_pha)); any(isnan(H_bucal)); any(isnan(H_nasal_mid)); any(isnan(H_nasal_left)); any(isnan(H_nasal_right))])'
    case '00000'
        H_eq = H_pha .*...
            (H_bucal .* Y_bucal + H_nasal_mid .*...
            ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left)) .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '00001'
        H_eq = H_pha .*...
            (H_bucal .* Y_bucal + H_nasal_mid .* H_nasal_left .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '00010'
        H_eq = H_pha .*...
            (H_bucal .* Y_bucal + H_nasal_mid .* H_nasal_right .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '00011'
        H_eq = H_pha .*...
            (H_bucal .* Y_bucal + H_nasal_mid .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '00100'
        H_eq = H_pha .*...
            (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left + Y_bucal);
    case '00101'
        H_eq = H_pha .*...
            (H_bucal .* Y_bucal + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_left + Y_bucal);
    case '00110'
        H_eq = H_pha .*...
            (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right) ./ (Y_nasal_right + Y_bucal);
    case '00111'
        H_eq = H_pha .* H_bucal;
    case '01000'
        H_eq = H_pha .* H_nasal_mid .*...
            ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left));
    case '01001'
        H_eq = H_pha .* H_nasal_mid .* H_nasal_left;
    case '01010'
        H_eq = H_pha .* H_nasal_mid .* H_nasal_right;
    case '01011'
        H_eq = H_pha .* H_nasal_mid;
    case '01100'
        H_eq = H_pha .*...
            (H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left);
    case '01101'
        H_eq = H_pha .* H_nasal_left;
    case '01110'
        H_eq = H_pha .* H_nasal_right;
    case '01111'
        H_eq = H_pha;
    case '10000'
        H_eq = (H_bucal .* Y_bucal + H_nasal_mid .*...
            ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left)) .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '10001'
        H_eq = (H_bucal .* Y_bucal + H_nasal_mid .* H_nasal_left .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '10010'
        H_eq = (H_bucal .* Y_bucal + H_nasal_mid .* H_nasal_right .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '10011'
        H_eq = (H_bucal .* Y_bucal + H_nasal_mid .* Y_nasal_mid) ./...
            (Y_bucal + Y_nasal_mid);
    case '10100'
        H_eq = (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left + Y_bucal);
    case '10101'
        H_eq = (H_bucal .* Y_bucal + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_left + Y_bucal);
    case '10110'
        H_eq = (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right) ./ (Y_nasal_right + Y_bucal);
    case '10111'
        H_eq = H_bucal;
    case '11000'
        H_eq = H_nasal_mid .*...
            ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left));
    case '11001'
        H_eq = H_nasal_mid .* H_nasal_left;
    case '11010'
        H_eq = H_nasal_mid .* H_nasal_right;
    case '11011'
        H_eq = H_nasal_mid;
    case '11100'
        H_eq = (H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left);
    case '11101'
        H_eq = H_nasal_left;
    case '11110'
        H_eq = H_nasal_right;
    case '11111'
        H_eq = NaN;
end

%***************************************
while 0
    figure
    clf; hold on; grid on;
    plot(f, 20*log10(abs(H_eq)))
end  % while 0
%***************************************

%===================================================================
% No formant detection in this function

F_form = [];

%===================================================================
% Output the partial transfer function if required

if nargout > 2
    % Depends on which sections exist
    switch num2str([any(isnan(H_pha)); any(isnan(H_bucal)); any(isnan(H_nasal_mid)); any(isnan(H_nasal_left)); any(isnan(H_nasal_right))])'
    	case '00000'
    		H_partiel.H_oral  = H_pha .* H_bucal .* Y_bucal ./ (Y_bucal + Y_nasal_mid);
    		H_partiel.H_nasal = H_pha .* H_nasal_mid .* ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left)) .* Y_nasal_mid ./...
                (Y_bucal + Y_nasal_mid);
    	case '00001'
    		H_partiel.H_oral  = H_pha .* H_bucal .* Y_bucal ./ (Y_bucal + Y_nasal_mid);
    		H_partiel.H_nasal = H_pha .* H_nasal_mid .* H_nasal_left .* Y_nasal_mid ./ (Y_bucal + Y_nasal_mid);
    	case '00010'
    		H_partiel.H_oral  = H_pha .* H_bucal .* Y_bucal ./ (Y_bucal + Y_nasal_mid);
    		H_partiel.H_nasal = H_pha .* H_nasal_mid .* H_nasal_right .* Y_nasal_mid ./ (Y_bucal + Y_nasal_mid);
    	case '00011'
    		H_partiel.H_oral  = H_pha .* H_bucal .* Y_bucal ./ (Y_bucal + Y_nasal_mid);
    		H_partiel.H_nasal = H_pha .* H_nasal_mid .* Y_nasal_mid ./ (Y_bucal + Y_nasal_mid);
    	case '00100'
            H_partiel = 'To do!'
            errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .*...
    		%       (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left + Y_bucal);
    	case '00101'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .*...
    		%       (H_bucal .* Y_bucal + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_left + Y_bucal);
    	case '00110'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .*...
    		%       (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right) ./ (Y_nasal_right + Y_bucal);
    	case '00111'
    		H_partiel = 'To do!'
    		% errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .* H_bucal;
    	case '01000'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .* H_nasal_mid .*...
    		%       ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left));
    	case '01001'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .* H_nasal_mid .* H_nasal_left;
    	case '01010'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .* H_nasal_mid .* H_nasal_right;
    	case '01011'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .* H_nasal_mid;
    	case '01100'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .*...
    		%       (H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left);
    	case '01101'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .* H_nasal_left;
    	case '01110'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha .* H_nasal_right;
    	case '01111'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_pha;
    	case '10000'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_bucal .* Y_bucal + H_nasal_mid .*...
    		%       ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left)) .* Y_nasal_mid) ./...
    		%       (Y_bucal + Y_nasal_mid);
    	case '10001'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_bucal .* Y_bucal + H_nasal_mid .* H_nasal_left .* Y_nasal_mid) ./...
    		%       (Y_bucal + Y_nasal_mid);
    	case '10010'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_bucal .* Y_bucal + H_nasal_mid .* H_nasal_right .* Y_nasal_mid) ./...
    		%       (Y_bucal + Y_nasal_mid);
    	case '10011'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_bucal .* Y_bucal + H_nasal_mid .* Y_nasal_mid) ./...
    		%       (Y_bucal + Y_nasal_mid);
    	case '10100'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left + Y_bucal);
    	case '10101'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_bucal .* Y_bucal + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_left + Y_bucal);
    	case '10110'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_bucal .* Y_bucal + H_nasal_right .* Y_nasal_right) ./ (Y_nasal_right + Y_bucal);
    	case '10111'
    		H_partiel = 'To do!'
    		% errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_bucal;
    	case '11000'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_nasal_mid .*...
    		%       ((H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left));
    	case '11001'
    		H_partiel.H_oral  = NaN;
    		H_partiel.H_nasal = H_nasal_mid .* H_nasal_left;
    	case '11010'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_nasal_mid .* H_nasal_right;
    	case '11011'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_nasal_mid;
    	case '11100'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = (H_nasal_right .* Y_nasal_right + H_nasal_left .* Y_nasal_left) ./ (Y_nasal_right + Y_nasal_left);
    	case '11101'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_nasal_left;
    	case '11110'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = H_nasal_right;
    	case '11111'
    		H_partiel = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    		%H_eq = NaN;
    end
end  % if nargout > 2

%===================================================================
% Output the nasal input impedance if required

if nargout > 3
    % Depends on which sections exist
    switch num2str([any(isnan(H_pha)); any(isnan(H_bucal)); any(isnan(H_nasal_mid)); any(isnan(H_nasal_left)); any(isnan(H_nasal_right))])'
    	case '00000'
    		Z_nasal_in = 1 ./ Y_nasal_mid;
    		%Z_nasal_in = 1 ./ Y_bucal;
    	case '00001'
    		Z_nasal_in = 1 ./ Y_nasal_mid;
    		%Z_nasal_in = 1 ./ Y_bucal;
    	case '00010'
    		Z_nasal_in = 1 ./ Y_nasal_mid;
    		%Z_nasal_in = 1 ./ Y_bucal;
    	case '00011'
    		Z_nasal_in = 1 ./ Y_nasal_mid;
    	case '00100'
    		Z_nasal_in = 1 ./ (Y_nasal_left + Y_nasal_right);
    	case '00101'
    		Z_nasal_in = 1 ./ Y_nasal_left;
    	case '00110'
    		Z_nasal_in = 1 ./ Y_nasal_right;
    	case '00111'
    		Z_nasal_in = NaN;
    	case '01000'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01001'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01010'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01011'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01100'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01101'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01110'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01111'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10000'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10001'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10010'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10011'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10100'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10101'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10110'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10111'
    		Z_nasal_in = 'To do!'
    		% errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11000'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11001'
    		Z_nasal_in = 1 ./ Y_nasal_mid;
    	case '11010'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11011'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11100'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11101'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11110'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11111'
    		Z_nasal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    end
end  % if nargout > 2

%===================================================================
% Output the buccal input impedance if required

if nargout > 4
    % Depends on which sections exist
    switch num2str([any(isnan(H_pha)); any(isnan(H_bucal)); any(isnan(H_nasal_mid)); any(isnan(H_nasal_left)); any(isnan(H_nasal_right))])'
    	case '00000'
    		Z_bucal_in = 1 ./ Y_bucal;
    		%Z_bucal_in = 1 ./ Y_bucal;
    	case '00001'
    		Z_bucal_in = 1 ./ Y_bucal;
    		%Z_bucal_in = 1 ./ Y_bucal;
    	case '00010'
    		Z_bucal_in = 1 ./ Y_bucal;
    		%Z_bucal_in = 1 ./ Y_bucal;
    	case '00011'
    		Z_bucal_in = 1 ./ Y_nasal_mid;
    	case '00100'
    		Z_bucal_in = 1 ./ Y_bucal;
    	case '00101'
    		Z_bucal_in = 1 ./ Y_bucal;
    	case '00110'
    		Z_bucal_in = 1 ./ Y_bucal;
    	case '00111'
    		Z_bucal_in = 1 ./ Y_bucal;
    	case '01000'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
   	case '01001'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01010'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01011'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01100'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01101'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01110'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '01111'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10000'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10001'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10010'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10011'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10100'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10101'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10110'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '10111'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11000'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11001'
    		Z_bucal_in = 1 ./ Y_nasal_mid;
    	case '11010'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11011'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11100'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11101'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11110'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    	case '11111'
    		Z_bucal_in = 'To do!'
    		errordlg('"aire2spectre_wall_vibration": Case to treat!')
    end
end  % if nargout > 2

clear global CONST_DAT;

return

