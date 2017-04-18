% depth detection measurement for GUI
% Livia Zarnescu 8-31-15
%

function [paramFit, extraFig] = MeasureEmbryoAspiration(ROIframes, t, params, embryoNum, ...
    manualMeasure, filePathProc, procFileName, handles, extraFig, startFrame)

framesToGet = round(params.frameRate*params.frameMult);
if nargin < 10
    startFrame = round(params.frameStartMult*params.frameRate);
end

if ~ishandle(extraFig)
    extraFig = figure;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. extract aspiration depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if manualMeasure
    
    aspiration_vals = zeros(1, length(t));
        
    % depth detection
    manual_depth = GetAspirationDepthManual(framesToGet, ...
        ROIframes(:,:,startFrame:end), ...
        extraFig);
    
    % adjust for varying vector lengths
    % chop off end if there is mismatch
    minLength = min(size(manual_depth,2), size(t,2));
    t = t(1:minLength);
    minLength = min(size(aspiration_vals,2),size(manual_depth,2));
    aspiration_vals = aspiration_vals(1:minLength) + manual_depth(1:minLength);
    aspiration_depth = aspiration_vals;
    
else
    aspiration_depth = ...
        GetAspirationDepthAuto(framesToGet, ...
        ROIframes(:,:,startFrame:end), ...
        params.adaptHist, extraFig);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Fit to Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = t(t < params.cropVal);
aspiration_depth = aspiration_depth(1:length(t));

if isfield(params, 'offsetVal')
    offsetVal = params.offsetVal;
else
    offsetVal = 6;
end

A = (aspiration_depth - offsetVal) * 40 * 10^-6 / params.convFactor; % convert from pixels to meters

tauTryList = .02:.02:.2;
fValList = [];

for kk = 1:length(tauTryList)
    
    start_params(1) = .2; % k0
    start_params(2) = .2; % k1
    start_params(3) = tauTryList(kk);%tauStart; % tau
    start_params(4) = 5; % n1_inv (slope of creep)
    Fin = params.Fin;
    
    [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 0, start_params);
    fValList = [fValList fval];
end

if(isfield(handles, 'PlotAxes'))
    axes(handles.PlotAxes);
    cla;
end

start_params(3) = tauTryList(fValList == min(fValList));
[xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 1, start_params);

paramFit = struct;
paramFit.xfine = xfine;
paramFit.yfit = yfit;
paramFit.k0 = k0;
paramFit.k1 = k1;
paramFit.n0 = n0;
paramFit.n1 = n1;
paramFit.F0 = F0;
paramFit.tau = tau;
paramFit.fval = fval;
paramFit.t = t;
paramFit.A = A;
paramFit.aspiration_depth = aspiration_depth;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Save all data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if manualMeasure
    
    save([filePathProc '\' procFileName num2str(embryoNum) ...
        '.mat'], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A', 'offsetVal');

else
    
    % save params if good
    if ~exist([filePathProc '\AutoMeasure\'], 'dir')
        mkdir([filePathProc], '\AutoMeasure');
    end
    
%     save([filePathProc '\AutoMeasure\' procFileName num2str(embryoNum) ...
%         '.mat'], 'xfine', 'yfit', ...
%         'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
%         'aspiration_depth', 'A');
   
    save([filePathProc '\' procFileName num2str(embryoNum) ...
        '.mat'], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A', 'offsetVal');

end

if(isfield(handles, 'MeasPipBtn'))
    set(handles.MeasPipBtn, 'Enable', 'on');
end
