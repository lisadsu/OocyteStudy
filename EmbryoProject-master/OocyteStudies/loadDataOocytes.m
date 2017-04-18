% load mouse oocyte data

function outputStruct = loadDataOocytes()

outputStruct = struct;
numExperiments = 12;
experimentType = zeros(1,numExperiments); % 0 is non-fert, 1 is fert
dateList = cell(1,numExperiments);
numOocytes = zeros(1,numExperiments); % total num in experimental group
oocyteNums = cell(1,numExperiments); % actual labels of meas oocytes
fertInfo = cell(1,numExperiments);
blastForm = cell(1,numExperiments);
hatchInfo = cell(1,numExperiments);
maturationEnv = cell(1,numExperiments); % 0 = in vivo, 1 = KSOM, 2 = MM
morphologyInfo = cell(1,numExperiments); % 0 = GV, 1 = M1, 2 = MII, -1 = fragmented
measHour = cell(1,numExperiments); % for IVF experiments, time at which mech was measured
fileNameApp = cell(1,numExperiments); % string to append to E[n] filename

% how much to add to k1 relative to what it would be if 0.2 psi were used
% for example, for .1 psi meas pressure, add .0314 to k1
k1ScaleFactor = zeros(1,numExperiments); 



%%%%%%%%%%%%%%%%%%
% 3-21-16
%%%%%%%%%%%%%%%%%%
% oocyte collection at 12 hrs post hCG
% mech measurement and fertilization at 16 hrs
% 0.1 psi pressure used
% 10 oocytes total in measured group (E1-10 were also measured post
% fertilization which I think harmed them; E11-20 were only measured before
% fertilization)
% controls were group cultured
% controls: 33/42 fertilized, 10/42 blast formation
experimentType(1) = 1;
dateList{1} = '3-21-16';
numOocytes(1) = 20; 
oocyteNums{1} = 1:20;
fertInfo{1} = [NaN*zeros(1,10) 0 1 1 1 0 1 1 1 1 0];
blastForm{1} = [NaN*zeros(1,10) 0 1 0 1 0 0 0 1 1 0];
hatchInfo{1} = [NaN*zeros(1,10) 0 0 0 1 0 0 0 0 1 0];
maturationEnv{1} = zeros(1,20);
k1ScaleFactor(1) = .0314;
morphologyInfo{1} = [1 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
measHour{1} = 14*ones(1,20);
fileNameApp{1} = cell(1,20); 
[fileNameApp{1}{:}] = deal(''); 



%%%%%%%%%%%%%%%%%%
% 3-28-16
%%%%%%%%%%%%%%%%%%
% oocyte collection at 12 hrs post hCG
% mech measurement and fertilization at 16 hrs
% 0.1 psi pressure used
% 24 oocytes total in measured group, only measured before fertilization
% controls were group cultured
% controls: 24/40 fertilized, 20/40 blast formation
experimentType(2) = 1;
dateList{2} = '3-28-16';
numOocytes(2) = 24; 
oocyteNums{2} = 1:24;
fertInfo{2} = [zeros(1,5) 1 0 0 0 1 0 0 1 1 0 0 1 1 0 0 0 1 1 1];
blastForm{2} = [zeros(1,12) 1 1 0 0 0 1 0 0 0 1 1 1];
hatchInfo{2} = [zeros(1,12) 1 1 0 0 0 1 0 0 0 1 1 1];
maturationEnv{2} = zeros(1,24);
k1ScaleFactor(2) = .0314;
morphologyInfo{2} = [2 2 2 2 2 2 2 1 2 2 1 2 2 1 2 2 2 2 1 2 2 2 2 2];
measHour{2} = 14*ones(1,24);
fileNameApp{2} = cell(1,24); 
[fileNameApp{2}{:}] = deal(''); 


%%%%%%%%%%%%%%%%%%
% 10-8-15 
%%%%%%%%%%%%%%%%%%
% In vitro maturation only
% 0.2 psi used
% mech group was measured at 0,8,16,24 hrs (so only take _24 expt), where 0
% is 5 hrs post hCG.
% exp group was measured only once at 24 hrs, but brought outside incubator
%   while mech group was being measured
% inc group was kept inside incubator and only measured once at 24 hrs
experimentType(3) = 0;
dateList{3} = '10-8-15';
numOocytes(3) = 125;
oocyteNums{3} = [repmat(1:20, [1 4]), 21:65];
fertInfo{3} = NaN*ones(1,125);
blastForm{3} = NaN*ones(1,125);
hatchInfo{3} = NaN*ones(1,125);
maturationEnv{3} = [zeros(1,20) ...
    ones(1,5) 2*ones(1,5) ones(1,5) 2*ones(1,5) ...
    ones(1,5) 2*ones(1,5) ones(1,5) 2*ones(1,5) ...
    ones(1,5) 2*ones(1,5) ones(1,5) 2*ones(1,5) ...
    ones(1,16) 2*ones(1,12) 2*ones(1,9) ones(1,8)];
k1ScaleFactor(3) = 0;
morphologyInfo{3} = [1 1 1 1 1 1 0 1 1 0 0 1 1 1 1 0 1 1 0 0 ...
    2 1 2 1 1 1 1 1 1 1 0 -1 1 1 1 1 2 1 1 1 ...
    2 1 -1 2 2 2 2 2 2 1 1 -1 2 2 1 1 2 1 1 2 ...
    2 1 -1 2 2 2 2 2 2 1 1 -1 2 2 1 1 2 1 1 2 ...
    1 1 0 0 1 1 0 1 1 2 0 1 1 2 2 2 2 2 2 1 2 -1 2 2 1 2 2 2 ...
    2 2 2 1 1 1 0 1 1 2 2 2 1 2 1 1 -1];
measHour{3} = [6*ones(1,20) 14*ones(1,20) 20*ones(1,20) 26*ones(1,20) ...
    26*ones(1,45)];
fileNameApp{3} = cell(1,125);
[fileNameApp{3}{1:20}] = deal('');
[fileNameApp{3}{21:40}] = deal('_8');
[fileNameApp{3}{41:60}] = deal('_16');
[fileNameApp{3}{61:80}] = deal('_24');
[fileNameApp{3}{81:125}] = deal('');



%%%%%%%%%%%%%%%%%%
% 10-21-15
%%%%%%%%%%%%%%%%%%
% In vitro maturation only
% 0.2 psi used
% mech group was measured at 0,8,16,24 hrs (so only take _24 expt), where 0
% is 5 hrs post hCG.
% exp group was measured only once at 24 hrs, but brought outside incubator
%   while mech group was being measured
% inc group was kept inside incubator and only measured once at 24 hrs
experimentType(4) = 0;
dateList{4} = '10-21-15';
numOocytes(4) = 161; 
oocyteNums{4} = [repmat(1:30, [1 4]), 31:71];
fertInfo{4} = NaN*ones(1,161);
blastForm{4} = NaN*ones(1,161);
hatchInfo{4} = NaN*ones(1,161);
maturationEnv{4} = [zeros(1,30) ones(1,10) 2*ones(1,20) ...
    ones(1,10) 2*ones(1,20) ones(1,10) 2*ones(1,20) ...
    ones(1,11) 2*ones(1,10) 2*ones(1,10) ones(1,10)];
k1ScaleFactor(4) = 0;
morphologyInfo{4} = [...
    0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 2 0 1 1 1 1 1 1 1 0 ...
    0 1 1 1 1 0 1 1 1 1 1 0 2 1 1 2 0 2 2 2 2 0 -1 1 2 2 1 2 1 2 ...
    0 2 -1 1 1 0 2 1 1 2 2 0 2 2 2 2 0 2 2 2 2 0 -1 2 2 2 2 2 2 2 ...
    0 2 -1 1 1 0 2 1 -1 2 2 0 2 2 1 2 0 2 2 2 2 0 -1 2 2 2 2 2 2 2 ...
    2 2 1 1 2 2 1 1 2 2 2 1 0 2 2 2 2 1 2 2 2 1 ...
    2 2 1 2 2 1 -1 2 1 1 1 1 1 2 -1 0 2 1 1];
measHour{4} = [6*ones(1,30) 14*ones(1,30) 20*ones(1,30) 26*ones(1,30) ...
    26*ones(1,41)];
fileNameApp{4} = cell(1,161);
[fileNameApp{4}{1:30}] = deal('');
[fileNameApp{4}{31:60}] = deal('_8');
[fileNameApp{4}{61:90}] = deal('_16');
[fileNameApp{4}{91:120}] = deal('_24');
[fileNameApp{4}{121:161}] = deal('');


% 12-2-15
% attempted IVF but failed
% 0.2 psi used
% oocytes matured in vivo and collected 16 hrs post hCG
% mech measurement done 16 hrs post hCG
experimentType(5) = 1;
dateList{5} = '12-2-15';
numOocytes(5) = 30; 
oocyteNums{5} = 1:30;
fertInfo{5} = zeros(1,30);
blastForm{5} = zeros(1,30);
hatchInfo{5} = zeros(1,30);
maturationEnv{5} = zeros(1,30);
k1ScaleFactor(5) = 0;
morphologyInfo{5} = [2 2 0 2 1 2 2 2 2 1 2 2 2 2 2 2 1 1 2 2 2 2 2 2 ...
    1 1 2 2 2 2];
measHour{5} = 20*ones(1,30);
fileNameApp{5} = cell(1,30); 
[fileNameApp{5}{:}] = deal(''); 



% 12-10-15
% attempted IVF but failed
% 0.2 psi used
% oocytes collected 3 hrs post hCG, matured in MM (1-20) or KSOM (21-40)
% mech measurement and fertilization done at 16 hrs post hCG
experimentType(6) = 1;
dateList{6} = '12-10-15';
numOocytes(6) = 80; 
oocyteNums{6} = [1:40 1:40];
fertInfo{6} = [NaN*zeros(1,40) zeros(1,40)];
blastForm{6} = [NaN*zeros(1,40) zeros(1,40)];
hatchInfo{6} = [NaN*zeros(1,40) zeros(1,40)];
maturationEnv{6} = [zeros(1,40) 2*ones(1,20) ones(1,20)];
k1ScaleFactor(6) = 0;
morphologyInfo{6} = [...
    1 1 0 0 0 0 0 1 1 0 1 0 1 1 1 1 1 0 1 0 ...
    1 1 0 0 0 0 1 0 0 1 1 0 1 0 0 1 0 0 1 0 ...
    2 2 1 0 1 2 1 1 2 1 2 1 1 2 2 2 2 1 2 0 ...
    1 1 0 1 1 0 1 1 1 2 1 0 1 0 1 1 1 1 1 1];
measHour{6} = [6*ones(1,40) 20*ones(1,40)];
fileNameApp{6} = cell(1,80);
[fileNameApp{6}{1:40}] = deal('');
[fileNameApp{6}{41:80}] = deal('_16');


% 1-13-16
% attempted IVF but failed
% 0.2 psi used
% oocytes collected 3 hrs post hCG, matured in MM (1-10) or KSOM (11-20)
% filenames ending in _16 correspond to mech measurement at 16 hrs post hCG
% mech measurement done at 5,16 hrs post hCG, fertilization done at 16 hrs post hCG
experimentType(7) = 1;
dateList{7} = '1-13-16';
numOocytes(7) = 40; 
oocyteNums{7} = [1:20 1:20];
fertInfo{7} = [NaN*zeros(1,20) zeros(1,20)];
blastForm{7} = [NaN*zeros(1,20) zeros(1,20)];
hatchInfo{7} = [NaN*zeros(1,20) zeros(1,20)];
maturationEnv{7} = [zeros(1,20) 2*ones(1,10) ones(1,10)]; 
k1ScaleFactor(7) = 0;
morphologyInfo{7} = [...
    0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ...
    0 1 1 2 1 2 1 0 2 1 1 2 0 0 2 1 2 2 1 2];
measHour{7} = [6*ones(1,20) 20*ones(1,20)];
fileNameApp{7} = cell(1,40);
[fileNameApp{7}{1:20}] = deal('');
[fileNameApp{7}{21:40}] = deal('_16');



% 1-20-16
% attempted IVF but failed
% 0.2 psi used
% oocytes collected 16 hrs post hCG
% mech measurement done at 16 hrs post hCG, fertilization done at 16 hrs post hCG
experimentType(8) = 1;
dateList{8} = '1-20-16';
numOocytes(8) = 20; 
oocyteNums{8} = 1:20;
fertInfo{8} = zeros(1,20);
blastForm{8} = zeros(1,20);
hatchInfo{8} = zeros(1,20);
maturationEnv{8} = zeros(1,20);
k1ScaleFactor(8) = 0;
morphologyInfo{8} = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 1 2];
measHour{8} = 20*ones(1,20);
fileNameApp{8} = cell(1,20); 
[fileNameApp{8}{:}] = deal(''); 




% 2-2-16
% attempted IVF, only worked on controls
% 0.2 psi used
% oocytes collected 12 hrs post hCG
% mech measurement done at 14 hrs post hCG, fertilization done at 16 hrs post hCG
experimentType(9) = 0;
dateList{9} = '2-2-16';
numOocytes(9) = 15; 
oocyteNums{9} = 1:15;
fertInfo{9} = NaN*ones(1,15);
blastForm{9} = NaN*ones(1,15);
hatchInfo{9} = NaN*ones(1,15);
maturationEnv{9} = zeros(1,15);
k1ScaleFactor(9) = 0;
morphologyInfo{9} = [1 1 2 2 2 2 2 2 2 2 1 1 2 2 2];
measHour{9} = 14*ones(1,15);
fileNameApp{9} = cell(1,15); 
[fileNameApp{9}{:}] = deal(''); 



% 5-6-16
% 3 hr post hCG collection
% 0.1 psi used
% morphology measured at 6 and 14 hr mark ('' and '_8')
% oocytes collected 1 hr post hCG, 30 measured + 12 controls
% mech measurement done every 4 hrs
experimentType(10) = 0;
dateList{10} = '5-6-16';
numOocytes(10) = 72; 
oocyteNums{10} = [1:30 1:30 101:112];
fertInfo{10} = NaN*ones(1,72);
blastForm{10} = NaN*ones(1,72);
hatchInfo{10} = NaN*ones(1,72);
maturationEnv{10} = [zeros(1,30) 2*ones(1,15) ones(1,27)];
k1ScaleFactor(10) = 0.0314;
morphologyInfo{10} = [...
    1 0 1 1 0 1 0 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 ...
    2 1 2 2 2 2 1 2 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 1 2 1 2 ...
    1 2 1 2 1 2 2 1 0 1 2 1];
measHour{10} = [6*ones(1,30) 14*ones(1,42)];
fileNameApp{10} = cell(1,72);
[fileNameApp{10}{1:30}] = deal('');
[fileNameApp{10}{31:60}] = deal('_12');
[fileNameApp{10}{61:72}] = deal('');



% 5-19-16
% Did IVF
% 0.1 psi used
% morphology measured at 10 hr mark
% oocytes collected 8 hrs post hCG, were mostly immature
% mech measurement done at 10 hrs post hCG, fertilization done at 11 hrs post hCG
experimentType(11) = 1;
dateList{11} = '5-19-16';
numOocytes(11) = 35; 
oocyteNums{11} = 1:35;
fertInfo{11} = [0 0 0 0 1 0 0 1 zeros(1,27)];
blastForm{11} = zeros(1,35);
hatchInfo{11} = zeros(1,35);
maturationEnv{11} = zeros(1,35);
k1ScaleFactor(11) = 0.0314;
morphologyInfo{11} = [0 0 1 1 2 1 1 1 1 0 0 1 1 1 0 1 1 0 1 0 ...
    1 1 0 0 0 1 0 0 0 1 0 0 1 1 0];
measHour{11} = 10*ones(1,35);
fileNameApp{11} = cell(1,35); 
[fileNameApp{11}{:}] = deal(''); 


% 6-6-16
% Did IVF
% 0.1 psi used
% morphology measured at 14 hr mark
% oocytes collected 12 hrs post hCG, were mostly immature
% mech measurement done at 14 hrs post hCG, fertilization done at 15 hrs post hCG
experimentType(12) = 1;
dateList{12} = '6-6-16';
numOocytes(12) = 40; 
oocyteNums{12} = 1:40;
fertInfo{12} = [0 0 1 1 1 0 1 1 0 1 0 0 0 1 0 1 1 1 0 0 1 1 1 1 1 1 ...
    0 1 0 0 0 0 1 1 1 1 1 1 0 1];
blastForm{12} = [0 0 1 1 1 0 0 1 0 1 0 0 0 1 0 1 1 0 0 0 1 0 1 0 1 1 ...
    0 1 0 0 0 0 1 1 1 0 1 0 0 1];
hatchInfo{12} = [0 0 1 0 0 0 0 1 0 1 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 1 ...
    0 1 0 0 0 0 0 0 1 0 0 0 0 0];
maturationEnv{12} = zeros(1,40);
k1ScaleFactor(12) = 0.0314;
morphologyInfo{12} = [2 2 2 2 2 2 2 2 1 2 2 1 2 2 2 2 2 1 1 1 2 2 2 2 2 2 ...
    1 1 1 2 1 2 2 2 2 2 2 1 1 2];
measHour{12} = 14*ones(1,40);
fileNameApp{12} = cell(1,40); 
[fileNameApp{12}{:}] = deal(''); 



outputStruct.experimentType = experimentType;
outputStruct.numExperiments = numExperiments;
outputStruct.dateList = dateList;
outputStruct.numOocytes = numOocytes;
outputStruct.oocyteNums = oocyteNums;
outputStruct.fertInfo = fertInfo;
outputStruct.blastForm = blastForm;
outputStruct.hatchInfo = hatchInfo;
outputStruct.maturationEnv = maturationEnv;
outputStruct.k1ScaleFactor = k1ScaleFactor;
outputStruct.morphologyInfo = morphologyInfo;
outputStruct.measHour = measHour;
outputStruct.fileNameApp = fileNameApp;

