function generateModelData  (templateModelFile)

%based on the template file modelDataFilename
%and the normal mode coordinates file  where
%normalModeCoordsFilename is defined below, generates new model file
%which contains updated H and NH bond info, the rest being obtained
%from the template file.
%The output filename is Modex.all.

addpath('/home/home4/apaydin/Mist/Matlab/Other');
addpath('/home/home4/apaydin/Mist/Matlab/General');

for modeIndex = 7:11

  normalModeCoordsFilename = sprintf('Mode%d.parsedPDB', modeIndex);

  modelDataFilename        = sprintf('Mode%d.all', modeIndex);

  generateOneModelData (templateModelFile, normalModeCoordsFilename, modelDataFilename);

end

  
