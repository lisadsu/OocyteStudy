% add new experiments to appropriate section

function saveNewMorphology(type)

% ==================================
% ============ HUMAN ===============
% ==================================
if isequal(type, 'human')
    
%     clear all;
    % 1 = no blast
    % 2 = poor quality blast
    % 3 = average quality blast
    % 4 = good quality blast
    
    TimeLapseParams = xlsread('C:\Users\Livia\Desktop\IVF\Processed Data\Human\cell cycle parameters.xlsx');

    % 3-8-13 cleavage stage
    morphology3_20_13 = [1 2 4 4];
    timeLapse3_20_13 = [];
    
    % 4-2-13 cleavage stage
    morphology4_2_13 = [1 1 1 3];
    timeLapse4_2_13 = [];
    
    % 4-17-13 cleavage stage
    morphology4_17_13 = 1*ones(1,12);
    morphology4_17_13([7 9 11]) = 2; % poor blast
    morphology4_17_13([2 5 6]) = 3; % average blast
    timeLapse4_17_13 = [];
    % morphology4_17_13([1 9]) = 4; % good blast
    % morphology4_17_13(7) = 4; % poor quality blast
    % morphology4_17_13(11) = 4;
    
    % all of these died because of incubator, don't bother including
    morphology4_30_13 = 2*ones(1,12);
    morphology4_30_13([2 5 8]) = 4;
    morphology4_30_13([1 3 4 7 10 12]) = 3;
    timeLapse4_30_13 = [];
    
    % 2PN stage
    % #s 1-9
    % pipSize 241
    morphology5_13_13 = 1*ones(1,9);
    morphology5_13_13([6]) = 2; % poor quality
    morphology5_13_13([7 9]) = 3; % medium quality blast
    morphology5_13_13([1 2]) = 4; % good quality blast
    morphology5_13_13([3 8]) = NaN; % totally aspirated :(

    timeLapse5_13_13 = zeros(length(morphology5_13_13),3);
    timeLapse5_13_13(isnan(morphology5_13_13),:) = NaN;
    timeLapse5_13_13(~isnan(morphology5_13_13),:) = TimeLapseParams(1:7,2:4);
    
    % #s 10-21
    % pipSize 229
    morphology5_29_13 = 1*ones(1,12);
    morphology5_29_13([5]) = 2; % poor quality blast
    morphology5_29_13([1]) = 3; % medium quality blast
    morphology5_29_13([7 11]) = 4; % good blast
    morphology5_29_13(11) = NaN; % this one is an outlier, ignore for graphing
    
    timeLapse5_29_13 = zeros(length(morphology5_29_13),3);
    timeLapse5_29_13(isnan(morphology5_29_13),:) = NaN;
    timeLapse5_29_13(~isnan(morphology5_29_13),:) = TimeLapseParams([8:17 19],2:4);
    
    % #s 22-32
    % pipSize 234
    morphology6_28_13 = 1*ones(1,11);
    morphology6_28_13(1) = 2;
    morphology6_28_13([5 7 9]) = 3;
    morphology6_28_13(10) = 4;
    morphology6_28_13([2 4]) = NaN; % totally aspirated :(
    timeLapse6_28_13 = zeros(length(morphology6_28_13),3);
    timeLapse6_28_13(isnan(morphology6_28_13),:) = NaN;
    timeLapse6_28_13(~isnan(morphology6_28_13),:) = TimeLapseParams(20:28,2:4);
    
    % #s 33-45
    % pipSize 241
    morphology7_30_13 = 1*ones(1,13);
    morphology7_30_13([4 7]) = 2;
    morphology7_30_13([1 8]) = 3;
    morphology7_30_13([2 9]) = 4;
    
    timeLapse7_30_13 = zeros(length(morphology7_30_13),3);
    timeLapse7_30_13(isnan(morphology7_30_13),:) = NaN;
    timeLapse7_30_13(~isnan(morphology7_30_13),:) = TimeLapseParams(29:41,2:4);
    
    % pipSize 242
    morphology11_5_13 = 1*ones(1,8);
    morphology11_5_13([2 4 8]) = 4;
%     morphology11_5_13([2 9]) = 4;
    
    % no time lapse params measured, set all to NaNs
    timeLapse11_5_13 = NaN*zeros(length(morphology11_5_13),3);
    timeLapse11_5_13(isnan(morphology11_5_13),:) = NaN;
    timeLapse11_5_13(~isnan(morphology11_5_13),:) = NaN*TimeLapseParams(42:49,2:4);
    
    % pipSize 236
    morphology11_8_13 = 1*ones(1,21);
    morphology11_8_13([8 10]) = 4;
    morphology11_8_13([19]) = 3;
    morphology11_8_13([9]) = 2;
%     morphology11_8_13([1 7 13 14 16 17]) = 5; % didn't cleave
    timeLapse11_8_13 = zeros(length(morphology11_8_13),3);
    timeLapse11_8_13(isnan(morphology11_8_13),:) = NaN;
    timeLapse11_8_13(~isnan(morphology11_8_13),:) = TimeLapseParams(50:70,2:4);
    
    % used for RNA-seq
    % pipSize 259
    morphology1_14_14 = NaN*ones(1,16);
    morphology1_14_14([1 7 8 11 15 16]) = NaN; % did not amplify, left with 10
    timeLapse1_14_14 = NaN*zeros(length(morphology1_14_14),3);

    % used for RNA-seq
    % pipSize 281 for first three, 255 thereafter
    morphology2_13_14 = NaN*ones(1,13);
    morphology2_13_14(13) = NaN; % did not amplify, left with 12
    timeLapse2_13_14 = NaN*zeros(length(morphology2_13_14),3);

    % pipSize 275
    morphology3_19_14 = 1*ones(1,20);
    morphology3_19_14([12 14]) = 2;
    morphology3_19_14([11]) = 3;
    morphology3_19_14([18 19]) = 4;
    timeLapse3_19_14 = zeros(length(morphology3_19_14),3);
    timeLapse3_19_14(isnan(morphology3_19_14),:) = NaN;
    timeLapse3_19_14(~isnan(morphology3_19_14),:) = TimeLapseParams(71:90,2:4);
    
    % pipSize 254
    % for cortical granule staining
    morphology6_9_14 = 5*ones(1,9);
    timeLapse6_9_14 = NaN*zeros(length(morphology6_9_14),3);  
    
    % pipSize 242
    % for cortical granule staining
    morphology6_23_14 = 5*ones(1,9);
    timeLapse6_23_14 = NaN*zeros(length(morphology6_23_14),3);  
    
    save('morphologyHuman.mat');
    
    % ==================================
    % ======== MOUSE OOCYTE ============
    % ==================================
elseif isequal(type, 'mouse oocyte')
    
    % divided at least once = dark blue
    % made it to 8-cells = light blue
    % made it to blast = green
    morphology5_14_13 = 1*ones(1,24); % didn't cleave
    morphology5_14_13([3 4 13 14 17 19]) = 4; % made it to blast
    morphology5_14_13([1 8 10 11 12 15 16 22 23]) = 2; % divided at least once
    morphology5_14_13([1 8 11 12 16 23]) = 4; % made it to 8-cells
    timeLapse5_14_13 = NaN*zeros(length(morphology5_14_13),3);

    % divided at least once = dark blue
    % made it to 8-cells = light blue
    % made it to blast = green
    morphology5_30_13 = 1*ones(1,20);
    morphology5_30_13([4 7 9 11 14 16 17 18 20]) = 4; % blast
    morphology5_30_13([2 5 10 13 19]) = 2; % no blast
    timeLapse5_30_13 = NaN*zeros(length(morphology5_30_13),3);
    
    % good strain
    % pipSize 129
    morphology2_11_14 = 2*ones(1,35); % didn't cleave
%     morphology2_11_14([2:5 8 14:16 18 20:22 24:29 32:34]) = 4;
%     morphology2_11_14([1 11 12 31]) = 1;
    morphology2_11_14([1 7 11 13 17 23 30 31 35]) = 4; % 2-8 cells
    morphology2_11_14([9 10 19]) = 4; % > 8 cells but no blast
    morphology2_11_14([5 16 27]) = 4; % poor quality blast
    morphology2_11_14([2:4 8 14 15 18 20:22 24:26 28:29 32:34]) = 4; % good quality blast
    timeLapse2_11_14 = NaN*zeros(length(morphology2_11_14),3);

    % good strain
    % pipSize 131
%     morphology2_27_14 = 2*ones(1,41);
%     morphology2_27_14([27:41]) = NaN;
%     morphology2_27_14([2 4 6 8:13 17 18 21 23 30 32]) = 4;
%     morphology2_27_14([7 22 26]) = 1;
    morphology2_27_14 = 2*ones(1,41); % didn't cleave
    morphology2_27_14([27:41]) = NaN;
    morphology2_27_14([1 5 15 20 24 25]) = 4; % 2-8 cells
    morphology2_27_14([14 16 19]) = 4; % > 8 cells but no blast
    morphology2_27_14([9 18 38 39]) = 4; % poor quality blast
    morphology2_27_14([2 4 6 8 10:13 17 21 23 30 32]) = 4; % good quality blast
    timeLapse2_27_14 = NaN*zeros(length(morphology2_27_14),3);

    % GV
    % pipSize 134
    morphology4_8_14 = 1*ones(1,34); % GV
    morphology4_8_14(15:34) = 3; % M1
    timeLapse4_8_14 = [];
    
    % GV
    % pipSize 125
    morphology4_11_14 = 1*ones(1,42); % GV
    timeLapse4_11_14 = [];
    
    % M1
    % pipSize 130
    morphology4_24_14 = 3*ones(1,52); % M1
    timeLapse4_24_14 = [];

    % pipSize 126
    % for before and after IVF measurements
    % 39 was lost and should be NaN
    morphology7_27_14 = 5*ones(1,54); %5
    morphology7_27_14([8 9 10 11 14 15 16 17 19 21 24 25 32 35 37 41]) = 6; % 6: fragmented or didn't divide at all
%     morphology7_27_14([18 29 36 43]) = 7; % 7: blastocyst
    timeLapse7_27_14 = [];    
    
    % pipSize 125
    % for before and after IVF measurements, then CG imaging
    morphology11_26_14 = 5*ones(1,53);
    timeLapse11_26_14 = [];
    
    save('morphologyMouseOocyte.mat');
    
    
    % ==================================
    % ======== MOUSE EMBRYO ============
    % ==================================
elseif isequal(type, 'mouse embryo')
    
%     clear all;
    TimeLapseParams = xlsread('C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo\cell cycle parameters.xlsx');

    morphology621 = [3 4 4 2 4 2 2 4 4 4 2 4 4 2 3];
    timeLapse621 = [];
    morphology627 = [4 4 4 4 2 2 4 4 2 4 2 2 2];
    timeLapse627 = [];
    
    morphology712 = 4*ones(1,20);
    morphology712(1) = 2;
    morphology712(14) = 2;
    timeLapse712 = [];
    
    morphology726 = 4*ones(1,20);
    morphology726(7) = 2;
    morphology726(10:11) = 2;
    morphology726(15) = 2;
    timeLapse726 = [];
    
    morphology1015 = 4*ones(1,15);
    morphology1015(5:6) = 2; % 2
    morphology1015(9) = 2; % 3
    timeLapse1015 = [];
    morphology1024 = 4*ones(1,19);
    morphology1024(15) = 3;
    morphology1024(19) = 3;
    timeLapse1024 = [];
    
    morphology111 = 4*ones(1,15);
    morphology111(14) = 2; % 2
    morphology111(1:2) = 2; % 3
    morphology111(10) = 2; % 3
    timeLapse111 = [];
    
    morphology115 = 4*ones(1,15);
    morphology115(9) = 2;
    morphology115(14) = 2;
    morphology115(6:7) = 2;
    morphology115(13) = 2;
    timeLapse115 = [];
    
    morphology1119 = 4*ones(1,22);
    morphology1119(9) = 3;
    morphology1119([8 14 15 16 17 19 20 21 22]) = 2;
    
    cellEnter1119 = ones(1,22);
    cellEnter1119([5 11 13 19 21 22]) = 0;
    cellEnter1119([2 7 9 10 14 20]) = .5;
    
    timeLapse1119 = zeros(length(morphology1119),3);
    timeLapse1119(isnan(morphology1119),:) = NaN;
    order1119 = [11:20 1:10 21:22];
    timeLapse1119(~isnan(morphology1119),:) = TimeLapseParams(order1119,:);
    
    morphology1130 = 4*ones(1,37);
    morphology1130([2 6 11 14 16:20 28 30 31 34]) = 2;
    morphology1130([22]) = 3; % poor quality blast
    
    cellEnter1130 = ones(1,37);
    cellEnter1130([11 14 28]) = 0;
    cellEnter1130([1 3 6 12 16 19 20 24 34 37]) = .5;
    
    timeLapse1130 = zeros(length(morphology1130),3);
    timeLapse1130(isnan(morphology1130),:) = NaN;
    timeLapse1130(~isnan(morphology1130),:) = TimeLapseParams(23:59,:);
    
    morphology1213 = 4*ones(1,32);
    morphology1213([3 7 19 24 27]) = 3;
    morphology1213([6 9 14 18 22 31 32]) = 2;
    
    % .5 means there was a delay
    % 0 means cell didn't enter at all
    cellEnter1213 = ones(1,32);
    cellEnter1213([6 17 18 25]) = 0;
    cellEnter1213([14 22]) = .5;
    
    timeLapse1213 = zeros(length(morphology1213),3);
    timeLapse1213(isnan(morphology1213),:) = NaN;
    timeLapse1213(~isnan(morphology1213),:) = [TimeLapseParams(61:90,:) ; ...
        0 0 0 ; 0 0 0];
    
    morphology1220 = 4*ones(1,32);
    morphology1220([2 4 6 8 16 19 21 23 25 29 31 32]) = 2;
    
    cellEnter1220 = ones(1,32);
    cellEnter1220([6 31 32]) = .5;

    timeLapse1220 = zeros(length(morphology1220),3);
    timeLapse1220(isnan(morphology1220),:) = NaN;
    timeLapse1220(~isnan(morphology1220),:) = [TimeLapseParams(91:120,:) ; ...
        0 0 0 ; 0 0 0];
    
    morphology211 = 4*ones(1,15);
    morphology211(10) = 2;
    cellEnter211 = ones(1,15);
    
    timeLapse211 = zeros(length(morphology211),3);
    timeLapse211(isnan(morphology211),:) = NaN;
    timeLapse211(~isnan(morphology211),:) = TimeLapseParams(121:135,:);
    
    % some used for transfer
    morphology5_7_13 = NaN*ones(1,46);
    timeLapse5_7_13 = [];
    
%     morphology6_14_13 = 5*ones(1,55);
    % pipSize 132
    morphology6_14_13 = NaN*ones(1,55);
    morphology6_14_13([2 4 12 15 16 28 29 31 32 34 37 38 39 45]) = NaN; % transferred good
    morphology6_14_13([1 3 5 7 8 14 26 30 35 43 44 47 48 54]) = NaN; %transferred bad
    morphology6_14_13([9 10 13 17 18 20 21 22 23 24 25 33 36 41 46 49 50 51]) = 4;
    morphology6_14_13([11 19 27 40 42]) = 2;
    morphology6_14_13([6 52]) = 3;
    
    
    timeLapse6_14_13 = zeros(length(morphology6_14_13),3);
    timeLapse6_14_13(isnan(morphology6_14_13),:) = NaN;
    timeLapse6_14_13(~isnan(morphology6_14_13),:) = TimeLapseParams(136:160,:);
    
    morphology7_1_13 = 4*ones(1,36);
    morphology7_1_13([7 12 18 27 34 36]) = 2;
    timeLapse7_1_13 = [];
    
    % pipSize 127
    morphology7_12_13 = NaN*ones(1,44); % all 26 of the ones in culture survived
    morphology7_12_13([1 2 3 4 5 7 8 10 12 14 16 17 19 22 24 27 29 30]) = NaN; % used for transfer
%     morphology7_12_13([2 5 8 12 14 16 17 22 24]) = 5; % "bad" sample in RNA-seq
    timeLapse7_12_13 = [];
    
    % got blast formation info
    % pipSize = 134
    morphology8_9_13 = NaN*ones(1,59);
    morphology8_9_13([6 7 14 26 30 33 36 44 49 54 57 59]) = NaN; % transferred good
    morphology8_9_13([1 2 5 8 13 21 25 27 34 37 42 48]) = NaN; % transferred bad
    morphology8_9_13([4 11 12 15 16 17 18 20 22 23 28 29 31 35 38 39 43 45 47 53]) = 4; %4
    morphology8_9_13([32 41 46 52 56]) = 2; %2
    morphology8_9_13([12 28 53 22 45]) = 3; %3
    timeLapse8_9_13 = [];
    
    % used for RNA-seq
    % pipSize 128
    morphology9_23_13 = NaN*ones(1,48);
%     morphology9_23_13([4 9 10 15 16 22 23 28 44]) = 5; % "bad2" sample in RNA-seq
%     morphology9_23_13([3 5 6 18 20 25 26 30 32]) = 5; % "good1" sample in RNA-seq
%     morphology9_23_13([33 34 36 37 39 40 41 42 45]) = 5; % "good2" sample in RNA-seq
    timeLapse9_23_13 = [];
    
    % used for transfer experiment
    morphology2_19_14 = 4*ones(1,58);
    morphology2_19_14([14 21 48]) = 2;
    morphology2_19_14([4 6 7 10 16 17 19 23 24 27 52 53 56 57]) = NaN; %transferred good
    morphology2_19_14([20 22 25 28 29 30 31 33 35 42 43 45 46 47]) = NaN; %transferred bad
    timeLapse2_19_14 = [];
    
    morphology2_20_14 = NaN*ones(1,63); % transferred, predicted bad
    morphology2_20_14([4 5 16 19 25 29 33 34 41 42 48 49 50 54 56]) = NaN; % transferred predicted good
    morphology2_20_14([1 7 10 11 17 18 20 21 23 26 27 31 32 43 44 45 51 52 57 60 63]) = 4;
    morphology2_20_14([3 6 13 15 22 38 40 47 55 59 61 62]) = 2;
    timeLapse2_20_14 = [];
    
    % pipSize 133 
    % for CG stain
    morphology6_17_14 = NaN*ones(1,16); % 5
    timeLapse6_17_14 = [];
    
    % pipSize 133
    % for CG stain
    morphology6_27_14 = NaN*ones(1,79); % 5
    timeLapse6_27_14 = [];
    
    % pipSize 126
    % for before and after IVF measurements
    % CBA strain
    % 39 was lost and should be NaN    
    morphology7_27_14 = 5*ones(1,54);
%     morphology7_27_14([8 9 10 11 14 15 16 17 19 21 24 25 32 35 37 41]) = 6; % fragmented or didn't divide at all
    morphology7_27_14([18 29 36 43]) = 5; % blastocyst
    timeLapse7_27_14 = [];

    
    % pipSize 126
    % for after IVF measurements
    % CBA strain
    morphology8_3_14 = 6*ones(1,48); % ?? 6
    morphology8_3_14([11 12 19 22 26 27 30 37 38 40 44]) = 6; % 6: fragmented or didn't divide at all
    morphology8_3_14([2 3 6 15 16 24 25 31 33 34 36 46]) = 7; % 7: 2-cells
    timeLapse8_3_14 = [];    
    
    % pipSize 123
    % for after IVF measurements
    % CBA strain
    morphology8_14_14 = 6*ones(1,49); % 6: no blast
    morphology8_14_14([4 5 8 9 11 18 22 36 40 43]) = 7; % 7: blastocyst %38 also
    timeLapse8_14_14 = [];
    
    % pipSize 124
    % for after IVF measurements
    % CBA strain
    morphology9_5_14 = 6*ones(1,31); % 6: no blast
    morphology9_5_14([2 12]) = 7; % 7: blast
    timeLapse9_5_14 = []; 
    
    % pipSize 124 (1-29), 125 (31-84), 122 (85-104)
    % for after microinjection, IVF measurements
    % CBA strain
    morphology9_24_14 = NaN*ones(1,104);
    morphology9_24_14(1:34) = NaN;
    morphology9_24_14(35:65) = NaN;
    morphology9_24_14(66:84) = NaN;
    timeLapse9_24_14 = [];  
    
    % pipSize 124
    % for before and after IVF measurements, then CG imaging
    morphology11_26_14 = 5*ones(1,53);
    timeLapse11_26_14 = [];    
    
    save('morphologyMouseEmbryo.mat');
    
end

