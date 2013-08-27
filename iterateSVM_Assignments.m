numFiles2Read = 1;
   
filenames     = cell(numFiles2Read,1);
%filenames{1}  = 'allPeaksAssignmentEnvironment-ForHSRI-forBayesianScoring-withNVR_ScoringMatrices-CH_RDCs-simpleInitialize.mat';
%filenames{1}  = 'assignmentEnvironmentAfterSecondRoundOfSVM.mat';
filenames{1} = 'allPeaksAssignmentEnvironment-ForFF2-forBayesianScoring-withNVR_ScoringMatrices-CH_RDCs-simpleInitialize.mat';
load (filenames{1});
%load mat file that includes processedDifferenceMatrices and ...
%    ASSIGNTABLE
%initializeVoters

NEGATIVE_LARGE_VALUE = -1000;

if (0)
  numVoters              = 6;
  dummyVoter             = initialize6Voters(differenceMatrixH_SHIFTX(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)), differenceMatrixN_SHIFTX(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)), differenceMatrixH_SHIFTS(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)), differenceMatrixN_SHIFTS(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),differenceMatrix_RDC1(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),differenceMatrix_RDC2(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)));
%  hsriFirst10Assignments = load ('hsriFirst10AssignmentsBySVM.txt');
%  candidatePeakIndices   = hsriFirst10Assignments(:,2);

%[score label vector peakIndex residueIndex shiftx1 shiftx2 shifts1 ...
% shifts2 RDC1 RDC2 SSE1 SSE2 SSE3 usable] = textread('top1each.txt',' %f %d %d %d %d %f %f %f %f %f %f %d %d %d %d');

  [score label vector peakIndex residueIndex shiftx1 shiftx2 shifts1 ...
   shifts2 RDC1 RDC2 SSE1 SSE2 SSE3 usable] = textread('hsri4-m23.svmlight.top1each.predict',' %f %d %d %d %d %f %f %f %f %f %f %d %d %d %d');


  for i = 1:length(peakIndex)

    [maxScore, entryIndex] = max(score);
    assert (maxScore > NEGATIVE_LARGE_VALUE);
    score(entryIndex)      = NEGATIVE_LARGE_VALUE;
   
    relPeakIndex    = find(ROWIN == peakIndex(entryIndex));
    relResidueIndex = find(COLIN == residueIndex(entryIndex));
    
    if (isempty(relResidueIndex)) | (isempty(relPeakIndex))
      %this residue or peak has been assigned before
      continue;
    end
      
    
    
    
    assert (~isempty(relPeakIndex));
    
    newASSIGNTABLE         = updateAssignTable(ASSIGNTABLE, relPeakIndex, ...
					       relResidueIndex);
  
    fprintf(1, 'testing the assignment of %d to %d\n', peakIndex(entryIndex),residueIndex(entryIndex));
    
    scoreSoFar = 0;
  
    [impossibleToAssign, newMASTER, newASSIGNTABLE, newROWIN, newCOLIN, newScoreSoFar, voter, newNumVoters] = ...
	noePruneAndDoMoreAssignments(newASSIGNTABLE, MASTER, ...
				     ROWIN, ...
				     COLIN, ...
				     scoreSoFar,...
				     voter, numVoters,...
				     NOES, ...
				     IALLDISTS, ALLDISTS, NTH);%, NH_RDCS, ...
							       %				 CH_RDCS, VECTORS_NH,MB_ShiftScore, ...
							       %				 MB_RDC_Score, numCS, numRDC, ...
							       %				 HSQCDATA, csCoefficient, ...
							       %				 rdcCoefficient, noeCoefficient, ...
							       %
							       %				 numHN_NOES);
    if (~impossibleToAssign)
      MASTER      = newMASTER;
      ASSIGNTABLE = newASSIGNTABLE;
      ROWIN       = newROWIN;
      COLIN       = newCOLIN;
      fprintf(1, 'accepted this assignment.\n');
    else
      fprintf(1, 'did not accept this assignment.\n');
    end
								 
											   
  end
  
  save assignmentEnvironmentAfterThirdRoundOfSVM.mat
  
  
  numCorrect = 0; numAssignedPeaks = 0;
  for peakIndex = 1:size(MASTER,1)
    residueIndex = find(MASTER(peakIndex,:));
    assert (length(residueIndex) <= 1);
    if (peakIndex == residueIndex)
      numCorrect = numCorrect + 1;
    end
    if (~isempty(residueIndex))
      numAssignedPeaks = numAssignedPeaks + 1;
    end
  end
  
  fprintf(1, 'the aggregate assignment has %d correct out of %d positions.\n',numCorrect, numAssignedPeaks);
  
							     
  numAssignments = sum(sum(MASTER));
  [assignmentPeakIndices,assignmentResidueIndices] = find(MASTER);
  fprintf(1, 'the assignments made after the call to noePrune are:\n');
  for assignmentIndex = 1:numAssignments
    fprintf(1, '%f %f\n',assignmentPeakIndices(assignmentIndex), assignmentResidueIndices(assignmentIndex));
  end
end %if(1)
labelFilename   = 'ff2Labels_withSSEs_assigntableFiltered.txt';
vectorsFilename = 'ff2Vectors_withPeakAndResidueIndicesAndSSEs_assigntableFiltered.txt';
printSVM_Information(MASTER,ASSIGNTABLE, ROWIN, COLIN, differenceMatrixH_SHIFTX, ...
		     differenceMatrixN_SHIFTX, ...
		     differenceMatrixH_SHIFTS, ...
		     differenceMatrixN_SHIFTS, ...
		     differenceMatrix_RDC1,differenceMatrix_RDC2, alphaHelixMatrix,betaStrandMatrix,coilMatrix, labelFilename, vectorsFilename);
