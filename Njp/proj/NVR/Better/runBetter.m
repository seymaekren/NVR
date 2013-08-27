function runBetter

%(useOrigData, useMBM_EM, SHIFTS_Filename, SHIFTX_Filename)



addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/Routines')
addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/Matlab');
addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/Matlab/FileProcessing');

dbstop if error
dbstop if warning

useOrigData     = 1;
SHIFTS_Filename = 'SHIFTS.txt';
SHIFTX_Filename = 'SHIFTX.txt';

[NOES, RDC1, RDC2, HDEXCHANGE, peakIDs, HSHIFTS, ...
 NSHIFTS, HN_SHIFTS, RDCs] = readNMR_Data(useOrigData);

if (useOrigData)

else
  %use distrib data
 
  [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = ...
      loadmodeldata('modeldata.txt');
 
end



betterHD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS, ...
 	 TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
 	 SHIFTS_Filename, SHIFTX_Filename, useMBM_EM);
 


% $$$ [assignmentAccuracy,NVR_SCORE,oldNVR_SCORE]= betterNVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:, 1),NOES,VECTORS, ...
% $$$ 						  TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,SHIFTS_Filename, SHIFTX_Filename);



