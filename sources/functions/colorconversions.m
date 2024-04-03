function [XYZ_to_LMS, sRGB_to_XYZ, sRGB_to_LMS, D65_XYZ, iPhoneWhite_XYZ] = colorconversions(Tristimulus, LMS_to_RGB, RGB_to_XYZ, iPhone_Primaries, D65_Spectra)

% Calculate sRGB_to_XYZ matrix using iPhone Primaries and Neitz Color Matching Functions

Tristimulus = Tristimulus(:,1:5:end)';

iPhoneRed_XYZ = 683.*sum(iPhone_Primaries(:,1).*Tristimulus);
iPhoneGreen_XYZ = 683.*sum(iPhone_Primaries(:,2).*Tristimulus);
iPhoneBlue_XYZ = 683.*sum(iPhone_Primaries(:,3).*Tristimulus);
iPhoneWhite_XYZ = iPhoneRed_XYZ+ iPhoneGreen_XYZ + iPhoneBlue_XYZ;

D65_XYZ = sum(D65_Spectra.*Tristimulus);

iPhoneRed_xy = iPhoneRed_XYZ/sum(iPhoneRed_XYZ);
iPhoneGreen_xy = iPhoneGreen_XYZ/sum(iPhoneGreen_XYZ);
iPhoneBlue_xy = iPhoneBlue_XYZ/sum(iPhoneBlue_XYZ);
iPhoneWhite_xy = iPhoneWhite_XYZ/sum(iPhoneWhite_XYZ);

D65_xy = D65_XYZ./sum(D65_XYZ);

% Derive sRGB to XYZ tristimulus matrix using D65 as the target white point

chromaticityMatrix = [iPhoneRed_xy(1)/iPhoneRed_xy(2), iPhoneGreen_xy(1)/iPhoneGreen_xy(2), iPhoneBlue_xy(1)/iPhoneBlue_xy(2);...
    1,1,1;...
    iPhoneRed_xy(3)/iPhoneRed_xy(2), iPhoneGreen_xy(3)/iPhoneGreen_xy(2), iPhoneBlue_xy(3)/iPhoneBlue_xy(2)];

whitePointMatrix = [D65_xy(1)/D65_xy(2);1;D65_xy(3)/D65_xy(2)];

max_RGB = inv(chromaticityMatrix)*whitePointMatrix;

max_RGBM = [max_RGB(1), 0, 0;...
    0, max_RGB(2), 0;...
    0, 0, max_RGB(3)];

sRGB_to_XYZ = chromaticityMatrix*max_RGBM;

iPhoneWhite_XYZ = sRGB_to_XYZ*[1;1;1];

XYZ_to_sRGB = inv(sRGB_to_XYZ);

% Calculate XYZ_to_LMS transformation matrix

LMS_to_XYZ = RGB_to_XYZ*LMS_to_RGB;
XYZ_to_LMS = inv(LMS_to_XYZ);

% Calculate RGB to LMS matrix

sRGB_to_LMS = XYZ_to_LMS*sRGB_to_XYZ;

end