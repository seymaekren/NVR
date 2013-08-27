function autoAdjustAnd_SSE_InfoToSHIFTX_File(templateSHIFTX_File, shiftxFilename)

%based on the template SHIFTX file,
%and the SHIFTX file, where SHIFTX file is defined below, generates new SHIFTX file
%which contains in the order corresponding to the model data file
%the same information as in the original SHIFTX file, to which SSE
%info is added from the modeldata file.

addpath('/home/home4/apaydin/Mist/Matlab/General');

[baseModelFileResidueIDs baseModelFileResidueNames  sStruct unusedVar_Cs1 unusedVar_Cs2 unusedVar_Cs3 unusedVar_Cs4 unusedVar_Cs5 unusedVar_Cs6] =textread(templateSHIFTX_File,'%d %s %s %f %f %f %f %f %f');


outputShiftxFilename = sprintf('%s.adjusted', shiftxFilename);
[resIndex resName cs1 cs2 cs3 cs4 cs5 cs6] = textread(shiftxFilename, ...
						  '%d %s %f %f %f %f %f %f');

writeShiftxFile               (outputShiftxFilename, resIndex, ...
			       resName, sStruct, cs1,cs2,cs3,cs4,cs5,cs6,...
			       baseModelFileResidueIDs, baseModelFileResidueNames);
