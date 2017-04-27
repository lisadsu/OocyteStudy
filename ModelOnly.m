%model only
%load depthdata, tries to use the model (kelvinFit3 and Kelvinmodel2)
clear all

basePath = 'Z:\Data\IVF\OocyteClinicalStudy\RawData\DepthData';

cd(basePath);
[fileName, pathName] = uigetfile('*.mat');
nameToSave = strcat(fileName,'_parameters');

load(fileName)
pressureApplied = 0.345; %psi
pipSize = 70; %micron
Fin = pressureApplied * 6894.744 * pi* (pipSize/2*(10^-6))^2 %Pressure*conv to N/m^2 * area in m^2
convFactor = pipRefOpeningPixels/70; %Find conversion factor from pipette opening size (pixel/micron)

A = (aspiration_depth - offsetVal) * 10^-6 / convFactor; % convert from pixels to meters
%addzero point to beginning
A(2:end+1)=A;
A(1)=0; 



tauTryList = .02:.02:.2;
fValList = zeros(1,length(tauTryList));

for kk = 1:length(tauTryList)
    
    start_params(1) = .2; % k0
    start_params(2) = .2; % k1
    start_params(3) = tauTryList(kk);%tauStart; % tau
    start_params(4) = 5; % n1_inv (slope of creep)
    
    [xfine, yfit, k0, k1, n0, n1, F0, tau, fval] = KelvinFit3(t, A, Fin, 1, start_params);
    fValList(kk) = fval;
end

figure();
clf;
start_params(3) = tauTryList(fValList == min(fValList));
[xfine, yfit, k0, k1, n0, n1, F0, tau, fval] = KelvinFit3(t, A, Fin, 1, start_params);
fval
[k1 n1 tau k0]


save([basePath '/ExtractedParameters/' nameToSave '.mat'], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A', 'offsetVal');
