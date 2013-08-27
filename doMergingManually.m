function doMergingManually(peakIndices, ...
			   residueIndices, voter,numVoters,ASSIGNTABLE,NOES,ALLDISTS,IALLDISTS,NTH,peakIDs,RESNUMS,RDC1,RDC2,VECTORS, MB_ShiftScore, ...
			    MB_RDC_Score, numCS, numRDC, HSQCDATA, csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES)

%filename1= 'patch1-4-Total_20Peaks-MoreNOEs.mat'; filename2 = 'patch5-7-Total_20Peaks-MoreNOEs.mat'; 
%relPeakIndicesToBeAllowed = 1:7;outputFilename = 'patch1-7_Total20Peaks.mat';
filename2 = 'patch1-3-allPeaks-FullMBScoreWithNewWeights-moreCandidatesReturned-MB_HFunction-scoringFullyInRankingCandidates.mat';
filename1 = 'patch4-6-allPeaks-FullMBScoreWithNewWeights-moreCandidatesReturned-MB_HFunction-scoringFullyInRankingCandidates.mat';
relPeakIndicesToBeAllowed = 1:6;outputFilename = 'patch1-6_AllResidues.mat';

for i = 1:numVoters
  voter{i} = voter{i}(relPeakIndicesToBeAllowed,:);
end
ASSIGNTABLE = ASSIGNTABLE(relPeakIndicesToBeAllowed,:);
peakIndices = peakIndices(relPeakIndicesToBeAllowed);


load (filename1);
foundMASTERsL = foundMASTERs;finalScoresL  = finalScores;assignmentAccuraciesL=assignmentAccuracies;totalNumAssignmentsL=totalNumAssignments;
load (filename2);
foundMASTERsR = foundMASTERs;finalScoresR  = finalScores;assignmentAccuraciesR=assignmentAccuracies;totalNumAssignmentsR=totalNumAssignments;
if (totalNumAssignmentsL == 0) | (totalNumAssignmentsR == 0)
  fprintf(1, 'either the left or right subtree returned no assignments.Returning...\n');
  totalNumAssignments  = 0;
  foundMASTERs         = cell(0,1);
  finalScores          = [];
  assignmentAccuracies = [];
else
  
  fprintf(1, 'merging assignments...\n');
  fprintf(1, 'noeCoefficient = %f csCoefficient = %f\n',noeCoefficient,csCoefficient);
  [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments] =...
      mergeAssignments(foundMASTERsL, foundMASTERsR, finalScoresL, finalScoresR, totalNumAssignmentsL, totalNumAssignmentsR, peakIndices, ...
		       residueIndices, voter,numVoters,ASSIGNTABLE,NOES,ALLDISTS,IALLDISTS,NTH,peakIDs,RESNUMS,RDC1,RDC2,VECTORS,MB_ShiftScore, MB_RDC_Score, numCS, numRDC, HSQCDATA, csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES);%,totalNumPeaks);
  fprintf(1, 'merged assignments\n');
  
end
save (outputFilename, 'foundMASTERs', 'finalScores', 'assignmentAccuracies', 'totalNumAssignments');
fprintf(1,'saved merge result.\n');
figure
keyboard

