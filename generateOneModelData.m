function generateOneModelData(templateModelFile, coordsFilename, outputFilename)

%based on the template file templateModelFile
%which contains updated H and NH bond info, the rest being obtained
%from the template file.

addpath('/home/home4/apaydin/Mist/Matlab/Other');
addpath('/home/home4/apaydin/Mist/Matlab/FileProcessing');
addpath('/home/home4/apaydin/Mist/NVR/Routines');


[vectors,types,baseModelFileResidueIDs,sStruct, H_Bond, allDists,I_AllDists,HSQCDATA] = loaddata(templateModelFile);

[N_Coords, H_Coords, parsedPDBFileResidueIDs] = readParsedPDBFile(coordsFilename); 

[NH_BondVectors]                              = computeNormalizedAmideBondVectors(N_Coords, H_Coords);

writeModelData(outputFilename, NH_BondVectors, H_Coords, parsedPDBFileResidueIDs, types, ...
	       baseModelFileResidueIDs, sStruct, H_Bond);
