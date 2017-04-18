% sort cell array in ascending order, with numbers to sort being between
% 'pattern1' and 'pattern2' in each element of cell array

function outputArray = FileNameSort(inputArray, pattern1, pattern2)


indexList = zeros(1,length(inputArray)); % holds index order of file names

% get index order of allFileNames
for i = 1:length(inputArray)
   
    currName = inputArray{i};
    firstIndex = strfind(currName, pattern1) + length(pattern1);
    lastIndex = strfind(currName, pattern2) - 1;
    indexList(i) = str2num(currName(firstIndex:lastIndex));
    
end

% sort allFileNames
[~, indexSort] = sort(indexList, 2, 'ascend');
outputArray = inputArray(indexSort);

















