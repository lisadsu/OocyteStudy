function embryoInfo = buildInfoStruct(fileToSave, numEmbryos)

if exist(fileToSave, 'file')
    load(fileToSave);
else
    
    embryoInfo = struct();
    
    for i = 1:numEmbryos
        fieldName = ['E' num2str(i)]
        embryoInfo.(fieldName) = struct();
    end
    
end




