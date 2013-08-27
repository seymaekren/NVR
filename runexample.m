function runexample(type)

%runexample: This script runs one of 3 tests on the sample data provided in the distribution. 
%Input:  type: if type == 1, it will run an example of the tensor estimation algorithm
%					on rdc data from ubiquitin using the ubiquitin model 1UBI. It estimates
%					the tensors, shows them, shows the 'real' tensors, and then returns a
%					percentile of accuacy (Higher numbers are better). 
%
%					if type == 2, it will run an example of the NVR assignment algorithm. 
%					if type == 3, it will run an example of the HD homology detection algorithm. 


% $$$ if (nargin == 0)
% $$$     fprintf(1, 'usage: runexample(type, useOrigData)\n');
% $$$     fprintf(1, 'useOrigData == 1 by default\n');
% $$$     fprintf(1, 'type == 2 for NVR, type == 3 for HD\n');
% $$$     fprintf(1, 'useOrigData == 1 corresponds to using this program with data files that were not originally included with the NVR/HD\n');
% $$$     return;
% $$$ elseif (nargin == 1)    
% $$$     useOrigData = 0; %default value
% $$$ end 

useOrigData = 1;

[HSQCDATA, NOES] = loadUbqData(useOrigData);

%[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = ...
%      loadmodeldata('modeldata.m');



modeIndex = 7; modelIndex = 6;

SHIFTX_Filename   = sprintf('SHIFTX/MySHIFTX.%d.model%d',modeIndex, modelIndex);

SHIFTS_Filename   = sprintf('SHIFTS/MySHIFTS.%d.model%d',modeIndex, modelIndex);

%      modelDataFilename = sprintf('ModelDataGeneration/ModelDataFiles/Mode%d.coords%d',  modeIndex, modelIndex);

%      [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata(modelDataFilename);

[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ignoredHSQCDATA] = loaddata('myinput.m');

if(type==1)
   %estimate tensors
   [S1]=NVR_TENEST(HSQCDATA(:,4),VECTORS);
   [S2]=NVR_TENEST(HSQCDATA(:,5),VECTORS);
   
   T1 = load('REALTENM1');T1=T1.S;
   T2 = load('REALTENM2');T2=T2.S;
   
   ESTIMATED_TENSOR_1 = S1
   ACTUAL_TENSOR_1 = T1
   
   ESTIMATED_TENSOR_2 = S2
   ACTUAL_TENSOR_2 = T2
   
   fprintf('Percentile of accuracy medium 1: %f\n', NVR_COMP_TEN(S1, T1)*100); 
   fprintf('Percentile of accuracy medium 2: %f\n', NVR_COMP_TEN(S2, T2)*100); 
   
elseif(type==2)
  t0 = clock;
  
  fout = fopen('NVRresults.txt', 'w');  
 
  if (useOrigData)
    [assignmentAccuracy] = NVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT,...
			       HBOND,...
			       ALLDISTS,'PREDICTEDSHIFTS.m', ...
			       'SHIFTX.txt');
  else
    [assignmentAccuracy] = NVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT,...
						      HBOND, ALLDISTS,'SHIFTS.m','SHIFTX.m');
  end
  fprintf(' elapsed %d seconds\n',etime(clock,t0));
  

  dirName = pwd;
  fprintf(fout, '%s %f\n', dirName, assignmentAccuracy);
  fclose(fout);

elseif(type==3)

  fprintf('Running Example 1: A Ubiquitin Model (1UBI) on ubiquitin data, bb RMSD = 0.6 Angstroms\n\n');
   
%  cd('hd_1ubi');
  %model 1, an actual ubiqutin structure
  [HD_SCORE,assignmentAccuracy, assignments, weightedHD_Score] = ...
      HD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS,SHIFTS_Filename,SHIFTX_Filename);
   
%  printAssignments ('possiblyBuggyAssignments.txt', assignments);
  
  keyboard
   
   
% $$$   cd('../hd_1h8c');
% $$$   fprintf('Running Example 2: A Homolog of ubiquitin (1H8C) on ubiquitin data, bb RMSD = 1.8 Angstroms\n\n');
% $$$   [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] =  loadmodeldata('modeldata2.m');
% $$$   HD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS);
% $$$  
% $$$   cd('../hd_1esr');
% $$$   fprintf('Running Example 3: A Non-Homolog of ubiquitin (1ESR) on ubiquitin data, bb RMSD = 5.9 Angstroms\n\n');
% $$$   [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] =   loadmodeldata('modeldata3.m');
% $$$   HD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS);
% $$$   cd('..');
   
end





