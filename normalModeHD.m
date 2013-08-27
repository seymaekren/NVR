function normalModeHD (inputModeIndex, inputModelIndex)

%normalModeHD: This script runs HD on 1D3Z homology models. 

dbstop if error
dbstop if warning

%addpath('/home/home4/apaydin/Mist/NVR/hd')
%addpath('/home/home4/apaydin/Mist/NVR/nvr')
%addpath('/home/home4/apaydin/Mist/NVR/tenest')
%addpath('/home/home4/apaydin/Mist/NVR/Routines');
%addpath('/home/home4/apaydin/Mist/Matlab/ScoringFunctions');
%addpath('/home/home4/apaydin/Mist/Matlab/Other');

b_printAssignments = 1;

useOrigData        = 1;

[HSQCDATA, NOES] = loadUbqData(useOrigData);


modeIndex        = inputModeIndex;

modelIndex       = inputModelIndex;

outputFilename   = sprintf('Mode%d_Model%d_HD_vs_NMA.txt', ...
			   modeIndex, modelIndex);			
  
if (checkIfFileExists(outputFilename))
  error ('please erase the %s\n',outputFilename);
end

fout                    = fopen(outputFilename, 'w');

fprintf(fout, 'useOrigData = %d\n', useOrigData);

SHIFTX_Filename   = sprintf('SHIFTX/MySHIFTX.%d.model%d',modeIndex, modelIndex);

SHIFTS_Filename   = sprintf('SHIFTS/MySHIFTS.%d.model%d',modeIndex, modelIndex);

modelDataFilename       = sprintf('ModelDataGeneration/ModelDataFiles/Mode%d.coords%d',  modeIndex, modelIndex);

[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata(modelDataFilename);

[HD_SCORE, assignmentAccuracy, assignments] = ...%HD
    NVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6), ...
	HSQCDATA(:,1),NOES,VECTORS,TYPES,...
	RESNUMS,SSTRUCT, HBOND, ALLDISTS,...%IALLDISTS, ...
	SHIFTS_Filename, SHIFTX_Filename);

fprintf(fout, '%d %f %f\n', modelIndex, HD_SCORE, assignmentAccuracy);

if (b_printAssignments)
  
  assignmentsOutputFilename = ...
      sprintf('assignmentsForMode%d_Model%d.txt',modeIndex, ...
	      modelIndex);
  printAssignments(assignmentsOutputFilename, assignments);
  
end

fclose(fout);



