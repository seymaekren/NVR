function changeLabilityInfoInModelData

addpath('/home/home4/apaydin/Mist/Matlab')

for fileIndex = 100:100:10000
  if (fileIndex < 1000)
    parsedStrideFilename = sprintf('ModelDataGeneration/StrideFiles/Parsed/Coords000%d.stride.parsed', fileIndex);
  elseif (fileIndex < 10000)
    parsedStrideFilename = sprintf('ModelDataGeneration/StrideFiles/Parsed/Coords00%d.stride.parsed', fileIndex);
  elseif (fileIndex == 10000)
    parsedStrideFilename = sprintf('ModelDataGeneration/StrideFiles/Parsed/Coords010000.stride.parsed');
  end
  
  modelDataFilename = sprintf('ModelDataGeneration/ModelDataFiles/ModelData%d', fileIndex);
  newModelDataFilename = sprintf('ModelDataGeneration/ModelDataFiles/NewModelData%d', fileIndex);
  replaceLabilityField(newModelDataFilename, modelDataFilename, parsedStrideFilename);
end

