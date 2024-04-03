function [RG_Match, Y_Match, Rayleigh_Range, loss] = anomaloscope(normal_L, normal_M, anomalous_L, anomalous_M)

% Import spectral power data for anomaloscope primaries

Py = readmatrix('Nagel.xlsx', 'Range', 'C12:C312')';
Pg = readmatrix('Nagel.xlsx', 'Range', 'D12:D312')';
Pr = readmatrix('Nagel.xlsx', 'Range', 'B12:B312')';

Py=Py.*1000;
Pg=Pg.*1000;
Pr=Pr.*1000;

%Calculate cone responses and yellow setting for color normal

Ey_L_Norm = sum(normal_L.*Py);
Eg_L_Norm = sum(normal_L.*Pg);
Er_L_Norm = sum(normal_L.*Pr);

Y_GL_Norm = (Eg_L_Norm*90)/Ey_L_Norm;
Y_RL_Norm = (Er_L_Norm*90)/Ey_L_Norm;

Ey_M_Norm = sum(normal_M.*Py);
Eg_M_Norm = sum(normal_M.*Pg);
Er_M_Norm = sum(normal_M.*Pr);

Y_GM_Norm = (Eg_M_Norm*90)/Ey_M_Norm;
Y_RM_Norm = (Er_M_Norm*90)/Ey_M_Norm;

Y_L_Norm = zeros(1,74);
Y_M_Norm = zeros(1,74);

for x=0:73
    Y_L_Norm(x+1)=(Y_RL_Norm*(x/73))+(Y_GL_Norm*(1-(x/73)));
    Y_M_Norm(x+1)=(Y_RM_Norm*(x/73))+(Y_GM_Norm*(1-(x/73)));
end

%Calculate cone responses and yellow setting for anomalous

Ey_L = sum(anomalous_L.*Py);
Eg_L = sum(anomalous_L.*Pg);
Er_L = sum(anomalous_L.*Pr);

Y_GL = (Eg_L*90)/Ey_L;
Y_RL = (Er_L*90)/Ey_L;

Ey_M = sum(anomalous_M.*Py);
Eg_M = sum(anomalous_M.*Pg);
Er_M = sum(anomalous_M.*Pr);

Y_GM = (Eg_M*90)/Ey_M;
Y_RM = (Er_M*90)/Ey_M;

%Plot Linear Response of Each Cone as R/R+G Value

Y_L=zeros(1,74);
Y_Match=zeros(1,74);

for x=0:73
    Y_L(x+1)=(Y_RL*(x/73))+(Y_GL*(1-(x/73)));
    Y_Match(x+1)=(Y_RM*(x/73))+(Y_GM*(1-(x/73)));
end

%Calculate Rayleigh Match for both color normal and color anomalous

x0=40;
fun_norm = @(x) ((Y_RL_Norm*(x/73))+(Y_GL_Norm*(1-(x/73))))-((Y_RM_Norm*(x/73))+(Y_GM_Norm*(1-(x/73))));

options = optimset('Display', 'off');

RG_Norm = fzero(fun_norm,x0, options);
Y_Norm = Y_RL_Norm*(RG_Norm/73)+Y_GL_Norm*(1-(RG_Norm/73));

x1=40;
fun = @(x) ((Y_RL*(x/73))+(Y_GL*(1-(x/73))))-((Y_RM*(x/73))+(Y_GM*(1-(x/73))));

RG_Match=fzero(fun,x1, options);

Y_Match=Y_RL*(RG_Match/73)+Y_GL*(1-(RG_Match/73));

%Calculate Match Range vs Color Normal Range

normalrange=0.345;
maximum = 73;

Normal_Threshold = RG_Norm+(normalrange/2);

Normal_Contrast = ((Y_RL_Norm*(Normal_Threshold/73))+(Y_GL_Norm*(1-(Normal_Threshold/73))))...
    /((Y_RM_Norm*(Normal_Threshold/73))+(Y_GM_Norm*(1-(Normal_Threshold/73))));

x2 = 1;
fun = @(x) ((Y_RL*(x/73))+(Y_GL*(1-(x/73))))/((Y_RM*(x/73))+(Y_GM*(1-(x/73)))) - Normal_Contrast;

RG_Threshold = fzero(fun,x2, options);

Threshold=RG_Threshold-RG_Match;

Rayleigh_Range=double(Threshold*2);

%Calculate Percent of Color Loss

loss = ((maximum/normalrange)-(maximum/Rayleigh_Range))/(maximum/normalrange);

end