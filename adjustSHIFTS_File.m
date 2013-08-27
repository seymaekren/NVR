function adjustSHIFTS_File(modelDataFilename)

%based on the template file modelDataFilename
%and the SHIFTS file, where
%SHIFTS file is defined below, generates new SHIFTS file
%which contains in the order corresponding to the model data file
%the same information as in the original SHIFTS file.

for modeIndex = 7:11
  for modelIndex = 1:11
    shiftsFilename       = sprintf('MySHIFTS.%d.model%d', modeIndex, modelIndex);
    adjustOneSHIFTS_File          (modelDataFilename, shiftsFilename);
  end
end
