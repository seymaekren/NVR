function adjustAndSSE_InfoToSHIFTX_File(templateSHIFTX_File, shiftxFilename)

%based on the template SHIFTX file,
%and the SHIFTX file, generates new SHIFTX file
%which contains in the order corresponding to the template data file
%the same information as in the original SHIFTX file, to which SSE
%info is added from the template data file.


[baseModelFileResidueIDs baseModelFileResidueNames  sStruct unusedVar_Cs1 unusedVar_Cs2 unusedVar_Cs3 unusedVar_Cs4 unusedVar_Cs5 unusedVar_Cs6] =textread(templateSHIFTX_File,'%d %s %s %f %f %f %f %f %f');
 

outputShiftxFilename = sprintf('%s.adjusted', shiftxFilename);
[resIndex resName cs1 cs2 cs3 cs4 cs5 cs6] = textread(shiftxFilename, ...
						  '%d %s %f %f %f %f %f %f');

%resIndex = resIndex - 11; %this is for 1CMZ/1DK8

writeShiftxFile               (outputShiftxFilename, resIndex, ...
			       resName, sStruct, cs1,cs2,cs3,cs4,cs5,cs6,...
			       baseModelFileResidueIDs, baseModelFileResidueNames);
