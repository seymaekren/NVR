%load allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat
%origPeakIndices = peakIndices;
%load ('patch1-3-allPeaks-FullMBScore-moreCandidatesReturned-MB_HFunction.mat');
%load('patch4-6-allPeaks-MBShiftAndRDC_Score-moreCandidatesReturned-correctHFunction.mat');
%load
%('patch1-3-allPeaks-FullMBScore-moreCandidatesReturned-MB_HFunction-scoringFullyInRankingCandidates.mat');

%load ('patch4-6-allPeaks-FullMBScoreWithNewWeights-moreCandidatesReturned-MB_HFunction-scoringFullyInRankingCandidates.mat');
%load ('patch1-3-allPeaks-FullMBScoreWithNewWeights-moreCandidatesReturned-MB_HFunction-scoringFullyInRankingCandidates.mat');

%load ('patch1-6_AllResidues.mat');

%load ('patch1-5-weirdShiftScoresMatrix.mat');
%load ('patchAllPeaks-FullMBScoreWithNewWeights-MB_HFunction-scoringUsingFullG_Plus_H_InRankingCandidates.mat');
%peakIndices     = [ 1     2     3     4     6     7];%origPeakIndices(1:6)
%peakIndices     = origPeakIndices(1:5);
%load ('patchAllPeaks-BayesianScore-SimpleHFunction-scoringUsingOnlyG_MaxFragmentSize5.mat');
load ('patchAllPeaks-BayesianScore-MaxFragmentSize10.mat');

csCoefficient = 0.22; rdcCoefficient = 0.05; noeCoefficient =0.73;
%csCoefficient = 2/77; rdcCoefficient = 70/77; noeCoefficient =5/77;

fprintf(1, 'coeffs: cs: %f rdc: %f noe: %f\n', csCoefficient, ...
	rdcCoefficient, noeCoefficient);

foundIt = 0; maxNumCorrectPeaks = 0;MB_Scores = []; overallMASTER= foundMASTERs{1}*0;
numPeaks = size(overallMASTER,1);

for resultIndex = 1:totalNumAssignments
  numCorrectPeaks = 0;
  overallMASTER   = overallMASTER + foundMASTERs{resultIndex};
%  for relPeakIndex = 1:length(origPeakIndices)
  for peakIndex = 1:numPeaks
%    peakIndex    = origPeakIndices(relPeakIndex);
    residueIndex = find(foundMASTERs{resultIndex}(peakIndex, :));
    if (isempty(residueIndex))
      continue;
    end
    %   assert (~isempty(residueIndex));
    if (peakIndex == residueIndex)
      numCorrectPeaks = numCorrectPeaks + 1;
    end
  end
%  if (numCorrectPeaks == length(origPeakIndices) | (numCorrectPeaks == sum(sum(foundMASTERs{resultIndex}))))
  if (numCorrectPeaks == sum(sum(foundMASTERs{resultIndex})))
    fprintf(1, 'found the correct assignment as #%d\n',resultIndex);
    fprintf(1, 'its score is %f\n', finalScores(resultIndex));
    if (~foundIt)
      foundPosition = resultIndex;
      foundIt = 1;
    end
% $$$     fprintf(1,'write return to continue.\n');
% $$$     figure
% $$$     keyboard
% $$$     break;
  end
  
  if (numCorrectPeaks > maxNumCorrectPeaks)
    maxNumCorrectPeaks = numCorrectPeaks;
    foundPosition = resultIndex;
  end

% $$$   MB_Scores(resultIndex) = computeMB_Score(foundMASTERs{resultIndex}, MB_ShiftScore, MB_RDC_Score, ...
% $$$ 					   numCS, numRDC, HSQCDATA, ...
% $$$ 					   ALLDISTS, NTH, csCoefficient, ...
% $$$ 					   rdcCoefficient, noeCoefficient);
end

if (~foundIt)
  fprintf(1, 'did not find the correct assignment among %d results\n',totalNumAssignments);
  fprintf(1, 'found max %d correct peaks (out of %d) at position %d.\n', ...
	  maxNumCorrectPeaks, sum(sum(foundMASTERs{foundPosition})),foundPosition);
%  fprintf(1,'write return to continue.\n');
%  figure;
%  keyboard
end

numPeaks = size(overallMASTER,1);
for peakIndex=1:numPeaks
  [maxNumConcurringVotes, residueIndex] = max(overallMASTER(peakIndex,:));
  if (maxNumConcurringVotes > 0)
    fprintf(1, 'peak %d has %d votes to be assigned to residue %d\n',peakIndex,maxNumConcurringVotes,residueIndex);
    if (peakIndex ~= residueIndex)
      fprintf(1, 'the correct assignment had %d votes.\n', overallMASTER(peakIndex,peakIndex));
    end
    %fprintf(1, 'its MB_ShiftScore is %f\n', MB_ShiftScore(peakIndex,residueIndex));
  end
end




MB_Scores = [];BayesianScores = [];
for resultIndex = 1:totalNumAssignments
 MB_Scores(resultIndex) = computeMB_Score(foundMASTERs{resultIndex}, MB_ShiftScore,MB_RDC_Score, ...
 					   numCS,numRDC, HSQCDATA, ...
 					  ALLDISTS,NTH, csCoefficient, ...
					  rdcCoefficient, noeCoefficient);

 BayesianScores(resultIndex) = computeBayesianScore(foundMASTERs{resultIndex});

end

Scores = BayesianScores;

[sortedScores, sortedIndices] = sort(-Scores);
newRank = find(sortedIndices == foundPosition);
fprintf(1, 'the result at %dth position becomes at %d th position.\n',foundPosition,newRank);
%fprintf(1, 'the coeffs are: cs: %f rdc: %f noe: %f\n', csCoefficient, ...
%	rdcCoefficient, noeCoefficient);
fprintf(1, 'the score becomes %f\n', sortedScores(newRank));

numCorrect = 0; numAssignedPeaks = 0;
for peakIndex = 1:size(foundMASTERs{totalNumAssignments},1)
  residueIndex = find(foundMASTERs{sortedIndices(1)}(peakIndex,:));
  assert (length(residueIndex) <= 1);
  if (peakIndex == residueIndex)
    numCorrect = numCorrect + 1;
  end
  if (~isempty(residueIndex))
    numAssignedPeaks = numAssignedPeaks + 1;
  end
end
  
fprintf(1, 'the best assignment has %d correct out of %d positions.\n',numCorrect, numAssignedPeaks);

figure;
plot(sortedScores,'*')
