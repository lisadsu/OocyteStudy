% recalc mech params
% taking into account curvature of cell edge into pipette

%% 1. bulk params

clearvars;

% input date
type = 'Mouse Oocyte';
date = '10-8-15';

baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\';
dataDir = [baseDir date ' analysis\AutoMeasure\'];

% find all .mat files in directory
allFiles = dir([dataDir 'aspiration_data_' strrep(date, '-', '_') '*.mat']);

curvatureAdjusted = 1;
convFactorAdjusted = 1;
curvatureOffset = 6*curvatureAdjusted;
convFactor = 108 + 20*convFactorAdjusted;
    
for i = 1:length(allFiles)
    
    i
    currFile = allFiles(i);
    load([dataDir currFile.name]);
    Fin = F0;
    
    display('Old params:');
    display([k1 n1 tau k0]);
    
    % note that curvature and convFactor have been adjusted
    % so we have convFactor = 40 um / 128 pixels
    % aspiration depth should be 6 pixels smaller to account for curvature
    % of cell edge
    % aspiration depth is raw data in units of pixels, A is in meters
    A = (aspiration_depth - curvatureOffset) * 40 * 10^-6 / convFactor; % convert from pixels to meters
    
%     figure(5);
%     clf;
    tauTryList = .02:.02:.2;
    fValList = [];
    
    for kk = 1:length(tauTryList)
        
        start_params(1) = .2; % k0
        start_params(2) = .2; % k1
        start_params(3) = tauTryList(kk);%tauStart; % tau
        start_params(4) = 5; % n1_inv (slope of creep)
        
        [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 0, start_params);
        fValList = [fValList fval];
    end
    
%     figure(5);
%     clf;
    start_params(3) = tauTryList(fValList == min(fValList));
    [xfine yfit k0 k1 n0 n1 F0 tau fval] = KelvinFit3(t, A, Fin, 0, start_params);
    
    display('New params:');
    display([k1 n1 tau k0]);
    display(' ');
    display(' ');

    save([dataDir currFile.name], 'xfine', 'yfit', ...
        'k0', 'k1', 'n0', 'n1', 'tau', 'F0', 'fval', 't', ...
        'aspiration_depth', 'A', 'curvatureAdjusted', 'convFactorAdjusted');
        
end


%% 2. Zona only
% 10-2-15, 10-8-15 and 10-21-15 are DONE ... don't do them again!

% clearvars
% date = '10-2-15';
% numEmbryos = 5;
% embryoNums = 1:5;
% extraString = '';
% baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\';
% dataDir = [baseDir date ' analysis\Zona\'];
% 
% load([dataDir 'allZona.mat']);
% cellSize = cellSize*108/128;
% 
% for i = 1:numEmbryos
%     
%     i
%     
%     % individual files
%     currFileName = [dataDir 'E' num2str(embryoNums(i)) extraString '_zona.mat'];
%     
%     if exist(currFileName, 'file')
%         
%         % depthList is currently in microns, keep but adjust
%         load(currFileName);
%         depthList = depthList*108/128;
%         cellCurr = cellSize(i);
%         save(currFileName, 'depthList', 'pressureList', 'cellCurr');
%         
%     end
%     
%     depthAll{i} = depthAll{i}*108/128;
%  
% end
% 
% save([dataDir 'allZona.mat'], 'depthAll', 'pressureAll', 'cellSize');













