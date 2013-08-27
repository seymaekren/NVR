function ...
    adjustAndSSE_InfoToSHIFTX_FileWithoutTemplateShiftxFile(shiftxFilename, myinputFilename)

%based on the template SHIFTX file,
%and the SHIFTX file, generates new SHIFTX file
%which contains in the order corresponding to the template data file
%the same information as in the original SHIFTX file, to which SSE
%info is added from the template data file.


[baseModelFileResidueIDs baseModelFileResidueNames  unusedVar1 ...
 unusedVar2 unusedVar3 unusedVar4 unusedVar5 unusedVar6 unusedVar7 ...
 sStruct unusedVar8 unusedVar9 unusedVar10 unusedVar11 unusedVar12] ...
    =textread(myinputFilename,'%d %s %f %f %f %f %f %f %f %s %s %d %f %f %f');
 

outputShiftxFilename = sprintf('%s.adjusted', shiftxFilename);
[resIndex resName cs1 cs2 cs3 cs4 cs5 cs6] = textread(shiftxFilename, ...
						  '%d %s %f %f %f %f %f %f');

%resIndex = resIndex - 11; %this is for 1CMZ/1DK8

%resIndex = resIndex - 78;  %this is for 1CMZ

writeShiftxFile               (outputShiftxFilename, resIndex, ...
			       resName, sStruct, cs1,cs2,cs3,cs4,cs5,cs6,...
			       baseModelFileResidueIDs, baseModelFileResidueNames);
