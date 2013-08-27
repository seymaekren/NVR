numFiles2Read = 1;
   
filenames = cell(numFiles2Read,1);
%filenames{1} = 'allPeaksAssignmentEnvironment-For1UBI-forBayesianScoring-withNVR_ScoringMatrices-CH_RDCs-simpleInitialize.mat';
filenames{1} = 'allPeaksAssignmentEnvironment-For1G6J-forBayesianScoring-withNVR_ScoringMatrices-CH_RDCs-simpleInitialize.mat';

%filenames{1} = 'allPeaksAfterInitialAssignments--For1UBI-WithNVR_Matrices-CH_RDCs-simpleInitialize.mat';
%   filenames{1} = 'allPeaksAssignmentEnvironment-For1UBI-forBayesianScoring-withNVR_ScoringMatrices.mat';
%   filenames{2} = 'allPeaksAssignmentEnvironment-For1UBQ-forBayesianScoring-withNVR_ScoringMatrices.mat';
%   filenames{1} = 'allPeaksAssignmentEnvironment-For1G6J-forBayesianScoring-withNVR_ScoringMatrices.mat';
%   filenames{3} = 'allPeaksAssignmentEnvironment-For1UD7-forBayesianScoring-withNVR_ScoringMatrices.mat';
%    filenames{1} =
%    'allPeaksAssignmentEnvironment-For1UBI-forBayesianScoring-withNVR_ScoringMatrices-CH_RDCs.mat';

%numVoters                        = 5;


%for fileIndex = 1:numFiles2Read
fileIndex = 1;
  %     fprintf(1, 'loading allPeaksAssignmentEnvironment-For1UBI-forBayesianScoring-withNVR_ScoringMatrices.mat');
  %     load ('allPeaksAssignmentEnvironment-For1UBI-forBayesianScoring-withNVR_ScoringMatrices.mat');
  fprintf(1, 'loading %s\n', filenames{fileIndex});
  load(filenames{fileIndex});
  
  labelFilename  = '1G6J_Labels.txt';
  vectorsFilename = '1G6J_VectorsWithExtraColumnForMissingEntries.txt';
  
  printSVM_Information(MASTER,differenceMatrixH_SHIFTX, ...
		       differenceMatrixN_SHIFTX, differenceMatrixH_SHIFTS, ...
		       differenceMatrixN_SHIFTS,differenceMatrix_RDC1,differenceMatrix_RDC2, labelFilename, vectorsFilename); 
%end  
  
  
