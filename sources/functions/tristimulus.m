function [Tristimulus, LMS_to_RGB, RGB_to_XYZ, Color_Matching_EE] = tristimulus(L_cone_E, M_cone_E, S_cone_E, V_Lambda)

WL=(400:1:700);

fundamentals = [L_cone_E; M_cone_E; S_cone_E];
RGB_to_LMS = [L_cone_E(WL==700), L_cone_E(WL==545), L_cone_E(WL==435);...
    M_cone_E(WL==700), M_cone_E(WL==545), M_cone_E(WL==435);...
    S_cone_E(WL==700), S_cone_E(WL==545), S_cone_E(WL==435)];

LMS_to_RGB = inv(RGB_to_LMS);

Color_Matching = LMS_to_RGB*fundamentals;

% Normalize color matching functions to sum to 1 across the visible
% spectrum

Color_Matching_Sum = sum(Color_Matching,2);

Color_Matching_EE = Color_Matching./Color_Matching_Sum;

Color_Matching_EE = Color_Matching_EE./max(Color_Matching_EE(2,:));

LMS_to_RGB(1,:) = LMS_to_RGB(1,:)/Color_Matching_Sum(1);
LMS_to_RGB(2,:) = LMS_to_RGB(2,:)/Color_Matching_Sum(2);
LMS_to_RGB(3,:) = LMS_to_RGB(3,:)/Color_Matching_Sum(3);

% Transform color matching functions into RGB space using usual equations
% R/R+G+B, G/R+G+B and B/R+G+B

Red = zeros(1,length(WL));
Green = zeros(1,length(WL));
Blue = zeros(1,length(WL));

for x = 1:length(Color_Matching_EE)
    Red(x) = Color_Matching_EE(1,x)/sum(Color_Matching_EE(1:3,x));
    Green(x) = Color_Matching_EE(2,x)/sum(Color_Matching_EE(1:3,x));
    Blue(x) = Color_Matching_EE(3,x)/sum(Color_Matching_EE(1:3,x));
end

% Calculate RGB to XYZ Conversion Matrix

Color_Matching_EE = Color_Matching_EE';

opts = optimset('Display','off');

x0=[1,1,1];

Best_Fit = lsqnonlin(@(x)(Color_Matching_EE(:, 1).*x(1)+Color_Matching_EE(:,2).*x(2)+Color_Matching_EE(:,3).*x(3)) - V_Lambda', x0, [], [], opts);

% Derive Equation for Alychne line of zero luminance

Alychne_Var = [Best_Fit(1)-Best_Fit(3), Best_Fit(2)-Best_Fit(3), Best_Fit(3)];

Alychne_slope = -Alychne_Var(1)/Alychne_Var(2);

Alychne_yint = -Alychne_Var(3)/Alychne_Var(2);

% Derive equation for Cr to Cg line

CrCg_slope = -100/99;
CrCg_yint = -CrCg_slope;

% Derive rg chromaticity coordinates for red primary

Cr = [0, 0, 0];

Cr(1) = (CrCg_yint-Alychne_yint)/(Alychne_slope-CrCg_slope);
Cr(2) = CrCg_slope*Cr(1) + CrCg_yint;
Cr(3) = 1-Cr(1)-Cr(2);

% Choose best chromaticity coordinates for blue and green primaries

Cb = [0,0,0];
Cg = [0,0,0];

CbCg_xtangent = Red(WL==500)-0.025;
CbCg_ytangent = Green(Red == Red(WL==500));
CbCg_slope = -3.25;
CbCg_yint = CbCg_ytangent - (CbCg_slope*CbCg_xtangent);

Cb(1) = (Alychne_yint-CbCg_yint)/(CbCg_slope-Alychne_slope);
Cb(2) = Alychne_slope*Cb(1)+Alychne_yint;
Cb(3) = 1-Cb(1)-Cb(2);

Cg(1) = (CrCg_yint-CbCg_yint)/(CbCg_slope-CrCg_slope);
Cg(2) = CbCg_slope*Cg(1)+CbCg_yint;
Cg(3) = 1-Cg(1)-Cg(2);

% Derive RGB to XYZ transformation

RGB_to_XYZ = [Cr(1), Cr(2), Cr(3);...
    Cg(1), Cg(2), Cg(3);...
    Cb(1), Cb(2), Cb(3)];

RGB_to_XYZ = inv(RGB_to_XYZ);
RGB_to_XYZ = RGB_to_XYZ';
RGB_to_XYZ(1,:) = RGB_to_XYZ(1,:)/sum(RGB_to_XYZ(1,:))*100;
RGB_to_XYZ(2,:) = RGB_to_XYZ(2,:)/sum(RGB_to_XYZ(2,:))*100;
RGB_to_XYZ(3,:) = RGB_to_XYZ(3,:)/sum(RGB_to_XYZ(3,:))*100;

RGB_to_XYZ = round(RGB_to_XYZ,7);

%Convert Spectral Locus to XYZ

Color_Matching_EE = Color_Matching_EE';

Tristimulus = zeros(3,61);

for x=1:length(Color_Matching_EE)
    Tristimulus(:,x) = RGB_to_XYZ*Color_Matching_EE(:,x);
end

Tristimulus=Tristimulus./max(Tristimulus(2,:));

Tristimulus(2,:) = V_Lambda;


end
