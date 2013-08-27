function  [relCandidatePeakIndices, relCandidateResidueIndices, ...
	   scores, numCandidates] =  voteAccordingToBayesianScore(voter, numVoters,...
						  ASSIGNTABLE, ...
						  unassignedPeakIndices,unassignedResidueIndices);


MAX_NUMCANDIDATES          = 3;
relCandidatePeakIndices    = zeros(MAX_NUMCANDIDATES,1); 
relCandidateResidueIndices = zeros(MAX_NUMCANDIDATES,1); 
scores                     = zeros(MAX_NUMCANDIDATES,1); 
numCandidates              = 0; 



assert (~isempty(ASSIGNTABLE));

V         = voter{1};




V         = V.*ASSIGNTABLE;

origV     = V; 



V                           = V.*ASSIGNTABLE;

% $$$ CS_RDC_COMPONENTS           = ASSIGNTABLE * 0 ;



stopGeneratingCandidates    = 0;

assert (length(unassignedPeakIndices) == size(ASSIGNTABLE,1));
assert (length(unassignedResidueIndices) == size(ASSIGNTABLE,2));

% $$$ reducedV                                          = cell(1,1);





while ((numCandidates < MAX_NUMCANDIDATES) & (~stopGeneratingCandidates))
%while (~stopGeneratingCandidates)

  maxScore  = max(max(V));
  
  if (maxScore == 0)
    break;
  end
  
  %if (maxVotes ~= maxVotesBeforeMaskingWithASSIGNTABLE)
  %  fprintf(1, 'ASSIGNTABLE masked a more likely assignment\n');
  %  keyboard
  %end
  
  [bestRelPeakIndex,bestRelResidueIndex] = find(V == maxScore);
  

  
  if (length(bestRelPeakIndex) > 1)
    
    bestRelPeakIndex    = bestRelPeakIndex(1);
    bestRelResidueIndex = bestRelResidueIndex(1);
    
  end

  if (ASSIGNTABLE(bestRelPeakIndex,bestRelResidueIndex)==0)
    stopGeneratingCandidates = 1;
    break;
  end


  
  
  
  numCandidates                              = numCandidates + 1;
  relCandidatePeakIndices(numCandidates)     = bestRelPeakIndex;
  relCandidateResidueIndices (numCandidates) = bestRelResidueIndex;
  scores(numCandidates)                      = -maxScore;

  fprintf(1, 'assignment of %d to %d has %f score.\n',unassignedPeakIndices(bestRelPeakIndex),unassignedResidueIndices(bestRelResidueIndex),scores(numCandidates));
% $$$   fprintf(1, 'the corresponding CS_RDC_COMPONENT is %f\n', CS_RDC_COMPONENTS(bestRelPeakIndex,bestRelResidueIndex));
  %  fprintf(1, 'MB_ShiftScore for that entry = %f, MB_RDC_Score for the',MB_ShiftScore(peakIndices(bestRelPeakIndex),residueIndices(bestRelResidueIndex)));
%  fprintf(1, ' entry is %f\n', MB_RdcScore(peakIndices(bestRelPeakIndex),residueIndices(bestRelResidueIndex)));
%  keyboard
  
  V(bestRelPeakIndex, bestRelResidueIndex) = 0;
end

[sortedScores, newIndices]       = sort(scores(1:numCandidates));
sortedRelCandidatePeakIndices    = zeros(numCandidates, 1);
sortedRelCandidateResidueIndices = zeros(numCandidates, 1);

for i = 1:numCandidates
  sortedRelCandidatePeakIndices(i) = ...
      relCandidatePeakIndices(newIndices(i));
  sortedRelCandidateResidueIndices(i) = relCandidateResidueIndices(newIndices(i)); 
end

relCandidatePeakIndices    = sortedRelCandidatePeakIndices;
relCandidateResidueIndices = sortedRelCandidateResidueIndices;
scores                     = sortedScores;
% $$$ if (numCandidates == 1)
% $$$   fprintf(1, 'numCandidates == 1\n');
% $$$   keyboard
% $$$ end

