function doIndividualPieceManually(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
				   residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
				   peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
				   MB_RDC_Score, ...
				   numCS, numRDC, HSQCDATA, ...
				   csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES);

relPeakIndicesToBeAllowed = 1:size(ASSIGNTABLE,1);filename = 'patchAllPeaks-BayesianScore-SimpleHFunction-scoringUsingOnlyG_MaxFragmentSize20.mat';
%relPeakIndicesToBeAllowed = 1:5;filename = 'patch1-5-FullMBScoreWithNewWeightsAndIntensities-scoringUsingFullG_Plus_H_InRankingCandidates.mat';

for i = 1:numVoters
  voter{i} = voter{i}(relPeakIndicesToBeAllowed,:);
end
ASSIGNTABLE               = ASSIGNTABLE(relPeakIndicesToBeAllowed,:);
origPeakIndices = peakIndices;
peakIndices     = peakIndices(relPeakIndicesToBeAllowed);



[foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments]=divideAndConquer(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
						  residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
						  peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
						  MB_RDC_Score, numCS, numRDC, HSQCDATA,csCoefficient, ...
						  rdcCoefficient, ...
						  noeCoefficient,numHN_NOES);

save (filename, 'foundMASTERs', 'finalScores', 'assignmentAccuracies', 'totalNumAssignments','origPeakIndices');
fprintf('saved variables.\n');
figure; %empty figure to signal end of computation.
keyboard
for resultIndex = 1:totalNumAssignments
  numCorrectPeaks = 0;
  for relPeakIndex = 1:length(relPeakIndicesToBeAllowed)
    peakIndex    = peakIndices(relPeakIndicesToBeAllowed(relPeakIndex));
    residueIndex = find(peakIndex, :);
    assert (~isempty(residueIndex));
    if (peakIndex ~= residueIndex)
      break;
    else
      numCorrectPeaks = numCorrectPeaks + 1;
    end
  end
  if (numCorrectPeaks == length(relPeakIndicesToBeAllowed))
    fprintf(1, 'found the correct assignment as #%d\n',resultIndex);
    break;
  end
end