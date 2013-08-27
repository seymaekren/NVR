function MC_AroundCorrectAssignment (NOES, ALLDISTS, NTH, startingMASTER);

%load allPeaksNewAssignmentEnvironment-forBayesianScoring.mat
%load allPeaksAssignmentsEnvironment.mat
%load allPeaksAfterInitialAssignments-WithMarsScores.mat
%load allPeaksAfterInitialAssignments-WithMBShiftScore.mat
%load allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score.mat
%load ('allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat');
%load ('allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-CorrectRDC_Tensor.mat');
%MB_Score            = (2 * MB_ShiftScore + 70 * MB_RDC_Score)/(2*numCS ...
%						  + 70 * numRDC);
%MB_Score            = (2 * MB_ShiftScore +
%MB_RDC_Score);%/(2*numCS ...
%MB_Score            = MB_RDC_Score;%/(2*numCS ...
%+ 70 * numRDC);
%MB_Score             = MB_ShiftScore/numCS;


load combinedScoringMatrix.txt

[numPeaks, numResidues] = size(combinedScoringMatrix);

printAssignmentAccuracy = 0;
%numVoters = 7;
%voter                            = cell(numVoters,1);
%voter                            =
%initializeVoters(voter,SSCP,SXCP,CP,TP,HDE, RP1, RP2);
%numVoters = 1;
%voter     = marsScore;
%voter     = MB_Score; %MB_ShiftScore;
%myMASTER  = ASSIGNTABLE * 0;
%MASTER    = myMASTER;


%candidateScores = MB_ShiftScore(1:3,:)*csCoefficient + 






%csCoefficient = 0.22; rdcCoefficient = 0.05; noeCoefficient =
%0.73;
%csCoefficient = 0.35; rdcCoefficient = 0.06; noeCoefficient =
%0.59;
%csCoefficient = 0.19; rdcCoefficient = 0.03; noeCoefficient = 0.78;



rand('twister',sum(100*clock))

%numPeaks  = size(MASTER,1);


MASTER     = startingMASTER; %zeros(numPeaks, numResidues);

%MASTER    = MASTER * 0;
%for i = 1:numPeaks
 % MASTER(i,i) = 1;
%end

%score      = computeMB_Score(MASTER, voter);
%score      = computeBayesianScore(MASTER);
%score      = computeMarsScore(MASTER);
%score      = computeMB_Score(MASTER, MB_ShiftScore,MB_RDC_Score,numCS,numRDC,HSQCDATA,ALLDISTS,NTH,csCoefficient,rdcCoefficient,noeCoefficient);

[score, accuracy, precision] = computeScore2(combinedScoringMatrix, ...
					     MASTER);

fprintf(1, 'starting matrix a.a. is %f\n', accuracy);
fprintf(1, 'starting score is %f\n', score);


keyboard


prevScore                    = score;


%keyboard

scores                  = [];
assignmentAccuracies    = [];
scoreIndex              = 1;
scores(1)               = score;
assignmentAccuracies(1) = accuracy; %computeAssignmentAccuracy(peakIDs, RESNUMS,MASTER, ...
				  %	  printAssignmentAccuracy);

%percentiles1          = [];
%percentiles2          =  [];
assert (noeCheck(MASTER, NOES, ALLDISTS, NTH) == 1);

%S1_full                                       = updateTen(MASTER, ...
%						  RDC1, VECTORS);

%S2_full                                       = updateTen(MASTER, ...
%						  RDC2, VECTORS);

%percentiles1(1) = NVR_COMP_TEN(S1_full, S1_full);
%percentiles2(1) = NVR_COMP_TEN(S2_full, S2_full);


iter      = 1;
maxNumIters = 10000;

while (iter < maxNumIters)
  while (iter < maxNumIters)
    newMASTER = MC_Move(MASTER);
    iter      = iter + 1;
    if (noeCheck(newMASTER, NOES, ALLDISTS, NTH) == 1)
%      score = computeMarsScore(newMASTER, voter, numVoters);
%      score = computeMB_Score(newMASTER, voter,numCS,numRDC);

 %       score                    = computeBayesianScore(newMASTER);
       [score, accuracy, precision] = computeScore2(combinedScoringMatrix, ...
						    newMASTER);
       %       score = computeMarsScore(newMASTER);
       %      score = computeMB_Score(MASTER, MB_ShiftScore,MB_RDC_Score,numCS,numRDC,HSQCDATA,ALLDISTS,NTH,csCoefficient,rdcCoefficient,noeCoefficient);
       if ((score ~= inf) & MC_Accept(prevScore,score))
	 break;
       end
    end
  end

  if (iter >= maxNumIters)
    break;
  end
  
  scoreIndex                     = scoreIndex + 1;
  scores(scoreIndex)             = score;
  assignmentAccuracies(scoreIndex) = accuracy; %computeAssignmentAccuracy(peakIDs, ...
					     %	  RESNUMS,newMASTER, printAssignmentAccuracy);

%  S1_full_partiallyCorrect = updateTen(newMASTER,RDC1,VECTORS);
%  S2_full_partiallyCorrect = updateTen(newMASTER,RDC2,VECTORS);
  
%  percentile1  = NVR_COMP_TEN(S1_full_partiallyCorrect, S1_full);
%  percentile2  = NVR_COMP_TEN(S2_full_partiallyCorrect, S2_full);

%  percentiles1(scoreIndex)      = percentile1;\
%  percentiles2(scoreIndex)      = percentile2;
  
  
  MASTER                         = newMASTER;
  prevScore                      = score;
  fprintf(1, 'iter = %d\n', iter);
%  keyboard
end  

if (0)
figure;
plot(assignmentAccuracies,'*');
title('assignment accuracy');
figure;
plot(scores, 'o');
title('scores');
figure;
plot(assignmentAccuracies,scores,'*');
correlationOfScore = corrcoef(assignmentAccuracies,scores);
correlationOfScore = correlationOfScore(1,2);
fprintf(1, 'the correlation of the score is %f\n', ...
	correlationOfScore);

end

%lowerScores = find(scores < 1756.17);
%lowerScores = find(scores < 2047.04);
%lowerScores = find(scores < 1418.04);
%lowerScores = find(scores < 3512.07);
lowerScores = find(scores < 1345.18);
scores(lowerScores)
assignmentAccuracies(lowerScores)
fprintf(1, 'found total %d assignments\n', scoreIndex);
keyboard
%figure;
%plot(assignmentAccuracy, percentiles1, '*');
%figure;
%plot(assignmentAccuracy, percentiles2, '*');