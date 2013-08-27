function [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments]=divideAndConquer(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
						  residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
						  peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
						  MB_RDC_Score, ...
						  numCS, numRDC, ...
						  HSQCDATA, ...
						  csCoefficient, ...
						  rdcCoefficient, noeCoefficient,numHN_NOES);

NUM_PEAKS_TO_ASSIGN_WITHOUT_RECURSION=100;

[numPeaks, numResidues] = size(ASSIGNTABLE);

fprintf(1, 'wanna assign %d peaks and %d residues in divideAndConquer\n',numPeaks,numResidues);
%keyboard
if (numPeaks > NUM_PEAKS_TO_ASSIGN_WITHOUT_RECURSION)
  
  voterL = voter;
  voterR = voter;
  
  for i = 1:numVoters
    voterL{i} = voterL{i}(1:floor(numPeaks/2),:);
    voterR{i} = voterR{i}(floor(numPeaks/2)+1:numPeaks,:);
  end
  
  ASSIGNTABLE_L = ASSIGNTABLE(1:floor(numPeaks/2),:);
  ASSIGNTABLE_R = ASSIGNTABLE(floor(numPeaks/2)+1:numPeaks,:);
  
  
  assert (length(peakIndices) == numPeaks);
  
  peakIndicesL  = peakIndices(1:floor(numPeaks/2));
  peakIndicesR  = peakIndices(floor(numPeaks/2)+1:numPeaks);
  
% $$$   MASTER_L      = []; MASTER_R = [];
% $$$ 
% $$$   totalNumPeaks = size(MASTER,1);
% $$$   
% $$$   for peakIndex = 1:totalNumPeaks
% $$$     relPeakIndexL = find(peakIndices(1:numPeaks/2) == peakIndex);
% $$$     relPeakIndexR = find(peakIndices(numPeaks/2+1:numPeaks) == peakIndex);
% $$$     
% $$$     if (isempty(relPeakIndexL)) & (isempty(relPeakIndexR))
% $$$       %already assigned peak
% $$$       assert (~isempty(find(MASTER(peakIndex,:))))
% $$$       MASTER_L    = [MASTER_L; MASTER(peakIndex,:)];
% $$$       MASTER_R    = [MASTER_R; MASTER(peakIndex,:)];
% $$$     elseif (~isempty(relPeakIndexL))
% $$$       assert (isempty(relPeakIndexR));
% $$$       MASTER_L    = [MASTER_L; MASTER(peakIndex,:)];
% $$$     elseif (~isempty(relPeakIndexR))
% $$$       assert (isempty(relPeakIndexL));
% $$$       MASTER_R    = [MASTER_R; MASTER(peakIndex,:)];
% $$$     else
% $$$       fprintf(1,'cant happen. peakIndex cannot be found in both left and right peakIndices.\n');
% $$$       keyboard;
% $$$     end
% $$$   end
  
  [foundMASTERsL, finalScoresL, assignmentAccuraciesL, ...
   totalNumAssignmentsL] = divideAndConquer(MASTER, voterL, numVoters, ...
					    ASSIGNTABLE_L, peakIndicesL, ...
					    residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
					    peakIDs, RESNUMS, RDC1, ...
					    RDC2, VECTORS, MB_ShiftScore, ...
					    MB_RDC_Score, ...
					    numCS, numRDC, ...
					    HSQCDATA, ...
					    csCoefficient, rdcCoefficient, noeCoefficient);

  [foundMASTERsR, finalScoresR, assignmentAccuraciesR, ...
   totalNumAssignmentsR] = divideAndConquer(MASTER, voterR, numVoters, ...
					    ASSIGNTABLE_R, peakIndicesR, ...
					    residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
					    peakIDs, RESNUMS, RDC1, ...
					    RDC2, VECTORS, MB_ShiftScore, ...
					    MB_RDC_Score, ...
					    numCS, numRDC, ...
					    HSQCDATA, ...
					    csCoefficient, rdcCoefficient, noeCoefficient);
  
  if (totalNumAssignmentsL == 0) | (totalNumAssignmentsR == 0)
    fprintf(1, 'either the left or right subtree returned no assignments.Returning...\n');
    totalNumAssignments  = 0;
    foundMASTERs         = cell(0,1);
    finalScores          = [];
    assignmentAccuracies = [];
    return;
  end
  
  fprintf(1, 'merging assignments...\n');
  [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments] =...
      mergeAssignments(foundMASTERsL, foundMASTERsR, finalScoresL, finalScoresR, totalNumAssignmentsL, totalNumAssignmentsR, peakIndices, ...
		       residueIndices, voter,numVoters,ASSIGNTABLE,NOES,ALLDISTS,IALLDISTS,NTH,peakIDs,RESNUMS,RDC1,RDC2,VECTORS);%,totalNumPeaks);
  fprintf(1, 'merged assignments\n');
%  keyboard;

    
else
  [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments] ...
      = a_star_assign(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
		      residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
		      peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
		    MB_RDC_Score, numCS, numRDC, HSQCDATA, csCoefficient, ...
		      rdcCoefficient, noeCoefficient,numHN_NOES);
  
end







function resolveConflicts()

    
    fprintf(1,'found %d conflicting pairs of peaks (assigned to',numConflicts);
    fprintf(1, ' the same residue).\n');
    fprintf(1, 'conflicting peak indices are:\n');
    for myIndex = 1:length(conflictingAbsPeakIndices)
      fprintf(1, '%d ',conflictingAbsPeakIndices(myIndex));
    end
    fprintf(1, '\n');
    fprintf(1, 'unassigned residue indices are:\n');
        for myIndex = 1:length(unassignedAbsResidueIndices)
      fprintf(1, '%d ',unassignedAbsResidueIndices(myIndex));
    end
    fprintf(1, '\n');
%    keyboard
    
    voterConflicting             = voter;
    
    conflictingRelPeakIndices    = zeros(length(conflictingAbsPeakIndices),1);
    unassignedRelResidueIndices  = zeros(length(unassignedAbsResidueIndices),1);
    
    for myIndex = 1:length(conflictingRelPeakIndices)
      conflictingRelPeakIndices(myIndex)   = find(peakIndices == ...
					 conflictingAbsPeakIndices(myIndex));
    end
    
    [conflictingRelPeakIndices]  = sort(conflictingRelPeakIndices);
    [conflictingAbsPeakIndices]  = sort(conflictingAbsPeakIndices);
    fprintf(1,'conflicting abs peak indices are: ');
    for myIndex = 1:length(conflictingAbsPeakIndices)
      fprintf(1, '%d ', conflictingAbsPeakIndices(myIndex));
    end
    fprintf(1, '\n');
%    keyboard
    for myIndex = 1:length(unassignedRelResidueIndices)
      unassignedRelResidueIndices(myIndex) = find(residueIndices == unassignedAbsResidueIndices(myIndex));
    end
    
    for i = 1:numVoters
      voterConflicting{i}       = voterConflicting{i}(conflictingRelPeakIndices,unassignedRelResidueIndices);
    end
  
    ASSIGNTABLE_Conflicting     = ASSIGNTABLE        (conflictingRelPeakIndices,unassignedRelResidueIndices);
    
    [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments] ...
	= divideAndConquer(mergedMASTER, voterConflicting, numVoters, ASSIGNTABLE_Conflicting, conflictingAbsPeakIndices, ...
			   unassignedAbsResidueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
			   peakIDs, RESNUMS, RDC1, RDC2, VECTORS);
  