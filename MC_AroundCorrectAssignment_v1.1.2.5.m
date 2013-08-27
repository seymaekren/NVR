%load allPeaksNewAssignmentEnvironment-forBayesianScoring.mat
%load allPeaksAssignmentsEnvironment.mat
load allPeaksAfterInitialAssignments-WithMarsScores.mat
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
csCoefficient = 0.19; rdcCoefficient = 0.03; noeCoefficient = 0.78;



 rand('twister',sum(100*clock))






numPeaks  = size(MASTER,1);

MASTER    = MASTER * 0;
for i = 1:numPeaks
  MASTER(i,i) = 1;
end

%score     = computeMB_Score(MASTER, voter);
score      = computeBayesianScore(MASTER);
%score      = computeMarsScore(MASTER);
%score      = computeMB_Score(MASTER, MB_ShiftScore,MB_RDC_Score,numCS,numRDC,HSQCDATA,ALLDISTS,NTH,csCoefficient,rdcCoefficient,noeCoefficient);
prevScore  = score;


%keyboard

scores    = [];
assignmentAccuracy = [];
scoreIndex = 1;
scores(1) = score;
assignmentAccuracy(1) = computeAssignmentAccuracy(peakIDs, RESNUMS,MASTER, ...
						  printAssignmentAccuracy);
assert (noeCheck(MASTER, NOES, ALLDISTS, NTH) == 1);
iter      = 1;
maxNumIters = 1000;

while (iter < maxNumIters)
  while (iter < maxNumIters)
    newMASTER = MC_Move(MASTER);
    iter      = iter + 1;
    if (noeCheck(newMASTER, NOES, ALLDISTS, NTH) == 1)
%      score = computeMarsScore(newMASTER, voter, numVoters);
%      score = computeMB_Score(newMASTER, voter,numCS,numRDC);

        score = computeBayesianScore(newMASTER);
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
  assignmentAccuracy(scoreIndex) = computeAssignmentAccuracy(peakIDs, ...
						  RESNUMS,newMASTER, printAssignmentAccuracy);
  MASTER                         = newMASTER;
  prevScore = score;
  fprintf(1, 'iter = %d\n', iter);
%  keyboard
end  


figure;
plot(assignmentAccuracy,'*');
title('assignment accuracy');
figure;
plot(scores, 'o');
title('scores');
figure;
plot(assignmentAccuracy,scores,'*');
correlationOfScore = corrcoef(assignmentAccuracy,scores);
correlationOfScore = correlationOfScore(1,2);
fprintf(1, 'the correlation of the score is %f\n',correlationOfScore);