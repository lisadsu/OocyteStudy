% load outcome data clinical study

function outputStruct = loadDataClinicalStudy()

% outcomeInfo: 10 = transfer, no hCG rise; 11 = transfer, hCG rise, no pregnancy; 
%              12 = transfer, confirmed pregnancy ultrasound
%              2 = unknown, use for prediction only
%              0 = no blast, 1 = blast

% still need to look up day 3 morphology for 11, 16, 22, 24


outputStruct = struct;

numParticipants = 33;
participantIDs = cell(1,numParticipants);
numEmbryos = zeros(1,numParticipants);
outcomeInfo = cell(1,numParticipants); % 
d3M = cell(1,numParticipants);
blastM = cell(1,numParticipants);
zygoteM = cell(1,numParticipants); % graded by me, 1 = worst, 2 = ok, 3 = best
ICSI = cell(1,numParticipants);
PGD = cell(1,numParticipants);
gender = cell(1,numParticipants); % 0 = male, 1 = female
patientAge = zeros(1,numParticipants);
patientBMI = zeros(1,numParticipants);


% EMECH-001
% 10-3-15
% E1 and E2 were non-ferts
participantIDs{1} = 'MECH001';
numEmbryos(1) = 7; 
outcomeInfo{1} = [NaN NaN 1 1 1 0 1 0 0];
zygoteM{1} = [2 2 3 3 2 3 3 3 3];
d3M{1} = {'4III', '5III', '7II', '8I', '7III', '7II', '8II', '8I', '7III'}; 
blastM{1} = {'', '', '6BB', '5AA', '5BB', '', '5BB', '', ''};
ICSI{1} = ones(1,9);
PGD{1} = [NaN NaN 1 1 0 NaN 0 NaN NaN];
gender{1} = [NaN NaN 0 0 0 NaN 0 NaN NaN];
patientAge(1) = 42;
patientBMI(1) = 43.52;


% EMECH-002
% exited from study due to low oocyte count
participantIDs{2} = 'MECH002';
numEmbryos(2) = 0; 
outcomeInfo{2} = [];
zygoteM{2} = [];
d3M{2} = {};
blastM{2} = {};
ICSI{2} = [];
PGD{2} = [];
gender{2} = [];
patientAge(2) = NaN;
patientBMI(2) = NaN;


% EMECH-003
% 10-23-15
% all blasts frozen for PGD and future transfer
participantIDs{3} = 'MECH003';
numEmbryos(3) = 14; 
outcomeInfo{3} = [1 1 1 1 1 1 1 1 1 0 1 1 0 0];
zygoteM{3} = [2 3 3 2 2 2 2 1 2 1 2 2 1 2];
d3M{3} = {'8I', '7I', '8I', '8II', '8II', '8I', '8I', '5III', '5II', ...
    '5III', '4III', '8II', '7III', '8I'};
blastM{3} = {'5BB', '5BA', '3BB', '3BB', '4BA', '5AA', '5BB', '5BB', ...
    '5BB', '', '2CC', '5AB', '5CB', ''};
ICSI{3} = ones(1,14);
PGD{3} = [1 1 1 1 1 1 0 0 1 NaN NaN 1 NaN NaN];
gender{3} = [0 1 0 0 1 1 1 1 0 NaN NaN 0 NaN NaN];
patientAge(3) = 34;
patientBMI(3) = 28.18;


% EMECH-004
% 10-23-15
% all discarded due to poor quality
participantIDs{4} = 'MECH004';
numEmbryos(4) = 13;
outcomeInfo{4} = [1 0 1 0 0 0 0 0 0 0 0 1 0];
zygoteM{4} = [3 2 3 3 3 1 3 2 3 3 2 2 2];
d3M{4} = {'8III', '6III', '8II', '8II', '8I', '4IV', '5III', '6III', ...
    '4IV', '8II', '5III', '8II', '7III'};
blastM{4} = {'2CC', '', '1CC', '1', '', '', '', '', '', '', '', '5CC', ''};
ICSI{4} = ones(1,13);
PGD{4} = NaN*ones(1,13);
gender{4} = NaN*ones(1,13);
patientAge(4) = 44;
patientBMI(4) = 31.57;


% EMECH-005
% 10-25-15
participantIDs{5} = 'MECH005';
numEmbryos(5) = 10;
outcomeInfo{5} = [1 0 1 1 0 1 0 1 0 1];
zygoteM{5} = [2 2 3 2 2 3 3 3 2 3];
d3M{5} = {'8II', '4I', '7I', '7I', '5V', '7II', '7II', '8II', '4V', '7II'};
blastM{5} = {'5BB', '', '6BB', '6AA', '', '5BB', '', '5BB', '', '5BB'};
ICSI{5} = [0 0 0 ones(1,7)];
PGD{5} = NaN*ones(1,10);
gender{5} = NaN*ones(1,10);
patientAge(5) = 30;
patientBMI(5) = 24.26;


% EMECH-006
% 11-9-15, last 2 embryos were 1PNs
% transferred E3, initially became pregnant but miscarried
% transferred E4, will have follow up ultrasound on 3/18
participantIDs{6} = 'MECH006';
numEmbryos(6) = 7;
outcomeInfo{6} = [0 0 11 12 0 1 0];
zygoteM{6} = [3 2 3 3 3 3 3];
d3M{6} = {'4I', '8III', '6II', '8II', '9III', '7II', '4IV'};
blastM{6} = {'', '', '6BB', '5BC', '', '3BC', ''};
ICSI{6} = ones(1,7);
PGD{6} = [NaN NaN 1 1 NaN 1 NaN];
gender{6} = [NaN NaN 1 1 NaN 1 NaN];
patientAge(6) = 34;
patientBMI(6) = 21.77;


% EMECH-007
% only 1 follicle, exited from study
participantIDs{7} = 'MECH007';
numEmbryos(7) = 0;
outcomeInfo{7} = [];
zygoteM{7} = [];
d3M{7} = {};
blastM{7} = {};
ICSI{7} = [];
PGD{7} = [];
gender{7} = [];
patientAge(7) = NaN;
patientBMI(7) = NaN;


% EMECH-008
% 11-11-15
% no normal PGS
participantIDs{8} = 'MECH008';
numEmbryos(8) = 11;
outcomeInfo{8} = [0 0 0 0 1 0 0 0 1 0 0];
zygoteM{8} = [3 2 3 3 3 3 3 3 3 3 2];
d3M{8} = {'3II', '3II', '8IV', '8II', '9II', '7II', '4V', '8II', ...
    '8II', '7IV', '8I'};
blastM{8} = {'', '', '', '', '6BB', '', '', '', '1', '', ''};
ICSI{8} = zeros(1,11);
PGD{8} = [NaN NaN NaN NaN 0 NaN NaN NaN 0 NaN NaN];
gender{8} = [NaN NaN NaN NaN 0 NaN NaN NaN 0 NaN NaN];
patientAge(8) = 36;
patientBMI(8) = 34.61;


% EMECH-009
% 12-2-15
participantIDs{9} = 'MECH009';
numEmbryos(9) = 13;
outcomeInfo{9} = [1 0 1 1 0 1 0 1 1 1 0 1 0];
zygoteM{9} = [3 3 3 3 3 3 3 3 3 3 3 2 3];
d3M{9} = {'7I', '4I', '8I', '7I', '6III', '8I', '8I', '7I', '8I', ...
    '7I', '8I', '6III', '4IV'};
blastM{9} = {'5BB', '', '5BC', '2', '', '5BB', '', '6BB', '5AB', ...
    '5BB', '', '2', ''};
ICSI{9} = zeros(1,13);
PGD{9} = [0 NaN 0 0 NaN 0 NaN 0 1 0 NaN 0 NaN];
gender{9} = [0 NaN 0 1 NaN 1 NaN 0 1 0 NaN 0 NaN];
patientAge(9) = 41;
patientBMI(9) = 29.97;


% EMECH-010
% 11-27-15
% E4 tried to make blast but no ICM
% E3 transferred at D3 due to poor embryo quality but no pregnancy
participantIDs{10} = 'MECH010';
numEmbryos(10) = 6;
outcomeInfo{10} = [0 0 NaN 1 0 0]; % E3 transferred at D3 but no pregnancy
zygoteM{10} = [1 1 1 1 1 1];
d3M{10} = {'2II', '7II', '7III', '6III', '5II', '6III'};
blastM{10} = {'', '', '', '3DB', '', ''};
ICSI{10} = ones(1,6);
PGD{10} = NaN*ones(1,6);
gender{10} = NaN*ones(1,6);
patientAge(10) = 32;
patientBMI(10) = 20.41;


% EMECH-011
% 12-4-15
% no blasts, all arrested at cleavage stage
participantIDs{11} = 'MECH011';
numEmbryos(11) = 9;
outcomeInfo{11} = [0 0 0 0 0 0 0 0 0];
zygoteM{11} = [2 2 3 1 1 1 3 2 1];
d3M{11} = {'4I', '7I', '4I', '4I', '7I', '4I', '5I', '6I', '4I'}; % fill this in
blastM{11} = {'', '', '', '', '', '', '', '', ''};
ICSI{11} = zeros(1,9);
PGD{11} = NaN*ones(1,9);
gender{11} = NaN*ones(1,9);
patientAge(11) = 42;
patientBMI(11) = 25.73;


% EMECH-012
% not enough oocytes
participantIDs{12} = 'MECH012';
numEmbryos(12) = 0;
outcomeInfo{12} = [];
zygoteM{12} = [];
d3M{12} = {};
blastM{12} = {};
ICSI{12} = [];
PGD{12} = [];
gender{12} = [];
patientAge(12) = NaN;
patientBMI(12) = NaN;


% EMECH-013
% 12-13-15
% Patient had E2 transferred, no pregnancy
% then had E20 transferred
% then had E3 and E16 transferred, resulted in 1 pregnancy so can't tell
% which is which as they are both female
participantIDs{13} = 'MECH013';
numEmbryos(13) = 21;
outcomeInfo{13} = [1 10 NaN 0 0 0 0 0 1 0 1 1 1 0 0 NaN 0 0 0 11 0];
zygoteM{13} = [3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 3 3 3 3 3 2];
d3M{13} = {'10I', '8III', '8II', '6IV', '7V', '6III', '8IV', '8IV', ...
    '7IV', '8IV', '7IV', '6V', '8III', '4V', '7III', '9II', '5III', ...
    '9III', '7II', '8II', '8II'};
blastM{13} = {'6BB', '5AA', '5BA', '', '', '', '', '', '5CB', '', '6CB', ...
    '5BC', '5AA', '', '', '5BB', '1', '1', '', '6AA', ''};
ICSI{13} = zeros(1,21);
PGD{13} = [1 1 1 NaN NaN NaN NaN NaN 0 NaN 0 1 0 NaN NaN 1 NaN NaN NaN 1 NaN];
gender{13} = [0 1 0 NaN NaN NaN NaN NaN 0 NaN 0 1 1 NaN NaN 0 NaN NaN NaN 1 NaN];
patientAge(13) = 35;
patientBMI(13) = 21.99;


% EMECH-014
% 2-4-16
% patient had #8 transferred
participantIDs{14} = 'MECH014';
numEmbryos(14) = 13;
outcomeInfo{14} = [0 0 0 1 1 1 1 12 1 0 0 1 1];
zygoteM{14} = [3 3 3 3 3 3 3 3 3 3 3 3 3];
d3M{14} = {'9I', '8I', '8I', '8I', '8I', '8II', '3IV', '2IV', '9II', ...
    '8I', '8II', '4II', '3III'}; 
blastM{14} = {'', '', '', '6AA', '5BB', '5AB', '5BB', '6AB', '2', ...
    '', '', '6AB', '5AB'};
ICSI{14} = [zeros(1,8) ones(1,5)];
PGD{14} = [NaN NaN NaN 1 1 0 0 1 0 NaN NaN 1 0];
gender{14} = [NaN NaN NaN 1 1 0 0 0 0 NaN NaN 1 1];
patientAge(14) = 35;
patientBMI(14) = 28.34;


% EMECH-015
% not started cycle yet
participantIDs{15} = 'MECH015';
numEmbryos(15) = 0;
outcomeInfo{15} = [];
zygoteM{15} = [];
d3M{15} = {};
blastM{15} = {};
ICSI{15} = [];
PGD{15} = [];
gender{15} = [];
patientAge(15) = NaN;
patientBMI(15) = NaN;


% EMECH-016
% 2-21-16
participantIDs{16} = 'MECH016';
numEmbryos(16) = 11;
outcomeInfo{16} = [0 1 0 0 0 0 0 1 0 1 0];
zygoteM{16} = NaN*ones(1,11);
d3M{16} = {'', '', '', '', '', '', '', '', '', '', ''};
blastM{16} = {'', '5CB', '', '', '', '', '', '5BB', '', '5AA', ''};
ICSI{16} = ones(1,11);
PGD{16} = [NaN 0 NaN NaN NaN NaN NaN 0 NaN 1 NaN];
gender{16} = [NaN 0 NaN NaN NaN NaN NaN 0 NaN 0 NaN];
patientAge(16) = 43;
patientBMI(16) = 23.38;


% EMECH-017
% 3-2-16
% had E5 transferred, positive US at 6wks but miscarriage at 9wks
participantIDs{17} = 'MECH017';
numEmbryos(17) = 7;
outcomeInfo{17} = [1 0 1 0 11 1 0];
zygoteM{17} = NaN*ones(1,7);
d3M{17} = {'9II', '9III', '8I', '8II', '8II', '8III', '7II'};
blastM{17} = {'5BB', '4CC', '5BB', '', '6AB', '6BB', ''};
ICSI{17} = zeros(1,10);
PGD{17} = [0 NaN 1 NaN 1 0 NaN];
gender{17} = [0 NaN 1 NaN 0 1 NaN];
patientAge(17) = 40;
patientBMI(17) = 33.87;


% EMECH-018
% 3-16-16 
% not enough eggs retrieved
participantIDs{18} = 'MECH018';
numEmbryos(18) = 0;
outcomeInfo{18} = [];
zygoteM{18} = [];
d3M{18} = {};
blastM{18} = {};
ICSI{18} = [];
PGD{18} = [];
gender{18} = [];
patientAge(18) = NaN;
patientBMI(18) = NaN;


% EMECH-019
% 3-26-16
% E4 transferred, positive US at 6wks, negative at 8wks
participantIDs{19} = 'MECH019';
numEmbryos(19) = 6;
outcomeInfo{19} = [1 0 0 11 0 0];
zygoteM{19} = NaN*ones(1,6);
d3M{19} = {'7I', '3I', '3II', '8II', '8II', '4III'};
blastM{19} = {'5AB', '', '', '5BB', '', ''};
ICSI{19} = ones(1,6);
PGD{19} = [0 NaN NaN 1 NaN NaN];
gender{19} = [0 NaN NaN 1 NaN NaN];
patientAge(19) = 40;
patientBMI(19) = 27.72;


% EMECH-020
% 3-29-16
% no normal ferts
participantIDs{20} = 'MECH020';
numEmbryos(20) = 0;
outcomeInfo{20} = [];
zygoteM{20} = [];
d3M{20} = {};
blastM{20} = {};
ICSI{20} = [];
PGD{20} = [];
gender{20} = [];
patientAge(20) = NaN;
patientBMI(20) = NaN;


% EMECH-021
% 4-2-16
participantIDs{21} = 'MECH021';
numEmbryos(21) = 6;
outcomeInfo{21} = [0 1 0 1 12 0]; 
zygoteM{21} = NaN*ones(1,6);
d3M{21} = {'3II', '8II', '3V', '6III', '5III', '5III'};
blastM{21} = {'', '2', '', '6AA', '5BB', ''};
ICSI{21} = ones(1,6);
PGD{21} = [NaN 1 NaN 0 1 NaN];
gender{21} = [NaN 0 NaN 0 1 NaN];
patientAge(21) = 35;
patientBMI(21) = 23.51;


% EMECH-022
% 4-4-16
participantIDs{22} = 'MECH022';
numEmbryos(22) = 9;
outcomeInfo{22} = [1 1 0 0 1 0 1 0 0];
zygoteM{22} = NaN*ones(1,9);
d3M{22} = {'', '', '', '', '', '', '', '', ''};
blastM{22} = {'4BA', '4BB', '', '', '5AB', '', '5BB', '', ''};
ICSI{22} = zeros(1,9);
PGD{22} = [0 0 NaN NaN 1 NaN 0 NaN NaN];
gender{22} = [0 0 NaN NaN 0 NaN 1 NaN NaN];
patientAge(22) = 35;
patientBMI(22) = 20.86;


% EMECH-023
% not measured yet
participantIDs{23} = 'MECH023';
numEmbryos(23) = 0;
outcomeInfo{23} = [];
zygoteM{23} = [];
d3M{23} = {};
blastM{23} = {};
ICSI{23} = [];
PGD{23} = [];
gender{23} = [];
patientAge(23) = NaN;
patientBMI(23) = NaN;


% EMECH-024
% 7-2-16
% last 3 are non-ferts
participantIDs{24} = 'MECH024';
numEmbryos(24) = 3;
outcomeInfo{24} = [1 1 0];% 0 0 0];
zygoteM{24} = NaN*ones(1,3);
d3M{24} = {'', '', ''};%, '', '', ''};
blastM{24} = {'5CC', '4BB', ''};%, '', '', ''};
ICSI{24} = ones(1,3);
PGD{24} = [0 0 NaN];% NaN NaN NaN];
gender{24} = [1 1 NaN];% NaN NaN NaN];
patientAge(24) = 43;
patientBMI(24) = 21.2;


% EMECH-025
% 7-7-16
% transferred E1 on 8/13/16
participantIDs{25} = 'MECH025';
numEmbryos(25) = 7;
outcomeInfo{25} = [10 0 1 1 1 1 0];
zygoteM{25} = NaN*ones(1,7);
d3M{25} = {'7III', '7I', '9I', '7I', '7IV', '9III', ''};
blastM{25} = {'5BB', '', '5BB', '5AB', '5BC', '5BB', ''};
ICSI{25} = ones(1,7);
PGD{25} = [1 NaN 1 0 0 0 NaN];
gender{25} = [1 NaN 0 0 0 1 NaN];
patientAge(25) = 36;
patientBMI(25) = 40.7;


% EMECH-026
% 7-9-16
% 6 thru 11 are non-ferts
participantIDs{26} = 'MECH026';
numEmbryos(26) = 5;
outcomeInfo{26} = [0 1 1 1 1];
zygoteM{26} = NaN*ones(1,5);
d3M{26} = {'5I', '7I', '7I', '8I', '8I'};
blastM{26} = {'1', '5AB', '5BB', '5BC', '6AA'};
ICSI{26} = zeros(1,5);
PGD{26} = NaN*zeros(1,5);
gender{26} = NaN*zeros(1,5);
patientAge(26) = 32;
patientBMI(26) = 17.8;


% EMECH-027
% 8-6-16
% transferred E6 in August 2016
participantIDs{27} = 'MECH027';
numEmbryos(27) = 10;
outcomeInfo{27} = [1 1 1 1 1 12 0 0 1 0];
zygoteM{27} = NaN*ones(1,10);
d3M{27} = {'9II', '8I', '8II', '6II', '9III', '8II', '8I', '9III', '7III', '8III'};
blastM{27} = {'5BB', '6AA', '6AA', '6BB', '5BB', '5BB', '', '', '6BB', ''};
ICSI{27} = ones(1,10);
PGD{27} = NaN*zeros(1,10);
gender{27} = NaN*zeros(1,10);
patientAge(27) = 29;
patientBMI(27) = 27.2;


% EMECH-028
% 8-19-16
% had E3 transferred, no hCG rise
participantIDs{28} = 'MECH028';
numEmbryos(28) = 6;
outcomeInfo{28} = [1 1 10 1 0 0];
zygoteM{28} = NaN*ones(1,6);
d3M{28} = {'8I', '9I', '7I', '8I', '', ''};
blastM{28} = {'5BB', '5AB', '5BB', '5AB', '', ''};
ICSI{28} = ones(1,6);
PGD{28} = [0 0 1 1 NaN NaN];
gender{28} = [0 1 0 0 NaN NaN];
patientAge(28) = 37;
patientBMI(28) = 37.3;


% EMECH-029
% 8-25-16
% banking
participantIDs{29} = 'MECH029';
numEmbryos(29) = 10;
outcomeInfo{29} = [1 0 1 NaN 0 0 1 0 0 0];
zygoteM{29} = NaN*ones(1,10);
d3M{29} = {'8II', '7II', '8II', '', '8I', '4II', '8III', '4III', '8II', '9II'};
blastM{29} = {'3BB', '1', '3AB', '', '', '', '3BB', '', '', ''};
ICSI{29} = zeros(1,10);
PGD{29} = NaN*zeros(1,10);
gender{29} = NaN*zeros(1,10);
patientAge(29) = 33;
patientBMI(29) = 22.8;

% EMECH-030
% not enough eggs
participantIDs{30} = 'MECH030';
numEmbryos(30) = 0;
outcomeInfo{30} = [];
zygoteM{30} = {};
d3M{30} = {};
blastM{30} = {};
ICSI{30} = [];
PGD{30} = [];
gender{30} = [];
patientAge(30) = NaN;
patientBMI(30) = NaN;

% EMECH-031
% 9-13-16
participantIDs{31} = 'MECH031';
numEmbryos(31) = 18;
outcomeInfo{31} = [1 1 0 1 1 1 1 0 0 1 0 1 0 0 0 1 0 1];
zygoteM{31} = NaN*ones(1,18);
d3M{31} = {'8II', '7II', '5II', '8I', '9I', '8I', '7I', '5III', '5II',...
           '6II', '6III', '8II', '9II', '8III', '9III', '9III', '7III', '4II'};
blastM{31} = {'6BA', '4BB', '', '5BB', '5BB', '5BB', '6BB', '', '', ...
              '6BA', '', '5BB', '', '', '', '5BC', '', '5BB'};
ICSI{31} = zeros(1,18);
PGD{31} = [0 0 NaN 1 1 0 1 NaN NaN 1 NaN 1 NaN NaN NaN 1 NaN 1];
gender{31} = [0 1 NaN 0 1 1 0 NaN NaN 0 NaN 1 NaN NaN NaN 0 NaN 0];
patientAge(31) = 32;
patientBMI(31) = 19.0;

% EMECH-032
% 9-15-16
participantIDs{32} = 'MECH032';
numEmbryos(32) = 5;
outcomeInfo{32} = [0 1 0 1 1];
zygoteM{32} = NaN*ones(1,5);
d3M{32} = {'', '', '', '', ''};
blastM{32} = {'', '5AA', '', '5BA', '4BC'};
ICSI{32} = ones(1,5);
PGD{32} = [NaN 1 NaN 1 0];
gender{32} = [NaN 1 0 NaN 1];
patientAge(32) = 38;
patientBMI(32) = 19.85;

% EMECH-033
% 9-18-16
participantIDs{33} = 'MECH033';
numEmbryos(33) = 8;
outcomeInfo{33} = [0 0 1 0 0 0 0 0];
zygoteM{33} = NaN*ones(1,8);
d3M{33} = {'6III', '4II', '9I', '8I', '9II', '3III', '8II', ''};
blastM{33} = {'', '', '6BB', '', '', '', '', ''};
ICSI{33} = zeros(1,8);
PGD{33} = NaN*ones(1,8);
gender{33} = NaN*ones(1,8);
patientAge(33) = 43;
patientBMI(33) = 38.01;


outputStruct.numParticipants = numParticipants;
outputStruct.participantIDs = participantIDs;
outputStruct.numEmbryos = numEmbryos;
outputStruct.outcomeInfo = outcomeInfo;
outputStruct.zygoteM = zygoteM;
outputStruct.d3M = d3M;
outputStruct.blastM = blastM;
outputStruct.ICSI = ICSI;
outputStruct.PGD = PGD;
outputStruct.gender = gender;
outputStruct.patientAge = patientAge;
outputStruct.patientBMI = patientBMI;






