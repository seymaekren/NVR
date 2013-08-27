load allPeaksAfterInitialAssignments-WithMarsScores.mat
%load allPeaksAfterInitialAssignments-WithMBShiftScore.mat

%load allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score.mat
%load ('allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-CorrectRDC_Tensor.mat')

TEST_MB_SCORE = 0;


if (TEST_MB_SCORE)
  
  %MB_Score            = (2 * MB_ShiftScore + 70 * MB_RDC_Score)/(2*numCS ...
  %						  + 70 * numRDC);
  debugRDC_Score
  %keyboard
  MB_Score = MB_RDC_Score; absoluteIndexed = 1;
  %MB_Score = RP2; absoluteIndexed = 0
  %MB_Score = RP1; absoluteIndexed = 0;
  %MB_Score = MB_ShiftScore; absoluteIndexed = 1;
  %MB_Score = SSCP; absoluteIndexed = 0;
  
  
  %MB_Score            = MB_RDC_Score;
% $$$ MB_Score            = MB_ShiftScore;
  numPeaks            = size(MASTER,1);
% $$$ voterGetsItRight    = zeros(numPeaks,1);
% $$$ rankOfCorrectVote   = [];
% $$$ for peakIndex = 1:numPeaks
% $$$   score                             = MB_Score(peakIndex,peakIndex);
% $$$   [sortedScores,sortedIndices]      = sort(-MB_Score(peakIndex,:));
% $$$   rankOfCorrectVote                 = [rankOfCorrectVote find(sortedIndices == peakIndex)];
% $$$   maxScore                          = max(MB_Score(peakIndex,:));
% $$$   if (score == maxScore)
% $$$     score2 = secondMax(MB_Score(peakIndex,:));
% $$$     if (score > score2)
% $$$       voterGetsItRight(peakIndex) = 1;
% $$$     end
% $$$   end
% $$$ end% $$$ 
% $$$ figure;
% $$$ plot(rankOfCorrectVote,'*');
% $$$ figure;
% $$$ plot(voterGetsItRight,'*');
% $$$ correctResiduesForVoter = find(voterGetsItRight)
% $$$ keyboard

  rankOfCorrectVote                             = [];
  voterGetsItRight                              = zeros(numPeaks, 1);
  %[assignedPeakIndices, assignedResidueIndices] = find(MASTER);
  %for i = 1:length(assignedPeakIndices)
  %  voterGetsItRight(assignedPeakIndices(i)) = 1;
  %end
  
  peakIndices                     = ROWIN;
  residueIndices                  = COLIN;
  %peakIndices                      = 1:numPeaks;
  %residueIndices                   = 1:size(MASTER,2);
  %for relPeakIndex = 1:length(peakIndices)
  for relPeakIndex = 1:size(ASSIGNTABLE,1)
    
    if (absoluteIndexed)
      
      score                         = MB_Score(peakIndices(relPeakIndex), ...
					       peakIndices(relPeakIndex));
      possibleResidueIndices        = residueIndices(find(ASSIGNTABLE(relPeakIndex,:)));
      [sortedScores, sortedIndices] = sort(-MB_Score(peakIndices(relPeakIndex),possibleResidueIndices));
      
    else
      
      score                         = MB_Score(relPeakIndex, ...
					       relPeakIndex);
      possibleRelResidueIndices        = (find(ASSIGNTABLE(relPeakIndex,:)));
      [sortedScores, sortedIndices] = sort(-MB_Score(relPeakIndex,possibleRelResidueIndices));
    end
    
    thisRank                      = find(sortedScores == -score);
    rankOfCorrectVote             = [rankOfCorrectVote thisRank(1)];
    
    if (absoluteIndexed)
      
      maxScore                      = max(MB_Score(peakIndices(relPeakIndex),possibleResidueIndices));
    else
      
      maxScore                      = max(MB_Score(relPeakIndex,possibleRelResidueIndices));
    end
    
    
    if (score == maxScore) & (score > 0)
      
      if (absoluteIndexed)
	
	score2 = secondMax(MB_Score(peakIndices(relPeakIndex), ...
				    possibleResidueIndices));
      else
	score2 = secondMax(MB_Score(relPeakIndex, ...
				    possibleRelResidueIndices));
	
      end
      
      %x    keyboard
      if (score > score2)
	voterGetsItRight(peakIndices(relPeakIndex)) = 1;
      end
    end
    %keyboard
  end
  
  
  figure;
  plot(rankOfCorrectVote,'*');
  figure;
  plot(voterGetsItRight,'*');
  correctResiduesForVoter = find(voterGetsItRight)
  keyboard
  
end



testMARS_Score = 0;

if (testMARS_Score)

  numPeaks = size(ASSIGNTABLE,1);
  voterGetsItRight = zeros(numPeaks,1);
  rankOfCorrectVote = [];
  for peakIndex = 1:numPeaks
    score    = marsScore(peakIndex,peakIndex);
    possibleResidueIndices = find(ASSIGNTABLE(peakIndex,:));
    sortedScores = sort(marsScore(peakIndex,possibleResidueIndices));
    rankOfCorrectVote = [rankOfCorrectVote find(sortedScores == score)];
    minScore = min(marsScore(peakIndex,possibleResidueIndices));
    if (score == minScore)
      score2 = secondMin(marsScore(peakIndex,possibleResidueIndices));
      if (score < score2)
	voterGetsItRight(peakIndex) = 1;
      end
    end
  end
  
  figure;
  plot(rankOfCorrectVote,'*');
  figure;
  plot(voterGetsItRight,'*');
  correctResiduesForVoter = find(voterGetsItRight)
  
  keyboard
end

testORIGINAL_VOTERS = 1;

if (testORIGINAL_VOTERS)
  
  %load allPeaksAssignmentsEnvironment.mat  
  load allPeaksAfterInitialAssignments--For1UBI.mat
  numVoters                        = 7;
%  voter                            = cell(numVoters,1);
  voter                            = initializeVoters(SSCP,SXCP,CP,TP,HDE, RP1, RP2);
  numPeaks                         = size(voter{1},1);

% $$$   for i = 1:numVoters
% $$$     diffProb = zeros(numPeaks,1);
% $$$     for peakIndex = 1:numPeaks
% $$$       diffProb(peakIndex) = voter{i}(peakIndex,peakIndex)-max(voter{i}(peakIndex,:));
% $$$       if (diffProb(peakIndex) == 0)
% $$$ 	diffProb(peakIndex) = voter{i}(peakIndex,peakIndex) - secondMax(voter{i}(peakIndex,:));
% $$$       end
% $$$     end
% $$$     figure; plot(diffProb, '*');
% $$$   end
  
% $$$   voterGetsItRight = zeros(numPeaks,numVoters);

  [numRelPeaks,numRelResidues] = size(voter{1});
  
  entropy          = zeros(numRelPeaks, numVoters);
  information      = zeros(numRelPeaks, numVoters);
  
%  for peakIndex = 1:numPeaks

   for relPeakIndex = 1:numRelPeaks
     for voterIndex = 1:numVoters
%      for residueIndex = 1:numResidues
      for relResidueIndex = 1:numRelResidues
	p = voter{voterIndex}(relPeakIndex,relResidueIndex);
	assert ((p>=0) & (p <=1 ));
	if (p > 0)
	  entropy(relPeakIndex,voterIndex) = entropy(relPeakIndex,voterIndex) ...
	      - p * log(p);
	end
      end
      information(relPeakIndex,voterIndex) = log(numRelResidues) - entropy(relPeakIndex,voterIndex);
      continue;
% $$$       p = voter{voterIndex}(peakIndex,peakIndex);
% $$$       maxProb = max(voter{voterIndex}(peakIndex,:));
% $$$       if (p == maxProb)
% $$$ 	p2 = secondMax(voter{voterIndex}(peakIndex,:));
% $$$ 	if (p > p2)
% $$$ 	  voterGetsItRight(peakIndex,voterIndex) = 1;
% $$$ 	end
     end
    end
  end
  
% $$$   correctResiduesForVoter = cell(1,1);
% $$$   for voterIndex = 1:numVoters
% $$$     figure;
% $$$     plot(voterGetsItRight(:,voterIndex),'*');
% $$$     correctResiduesForVoter{voterIndex} = find(voterGetsItRight(:,voterIndex));
% $$$   end
% $$$ end
close all
figure; plot(1:numPeaks, information(:,1),'*-',1:numPeaks, ...
	     information(:,2),'*-',1:numPeaks, information(:,3),'*-', ...
	     1:numPeaks, information(:,4),'*-',1:numPeaks,...
	     information(:,5),'*-',1:numPeaks,...
	     information(:,6),'*-',1:numPeaks, information(:,7),'*-');
legend('SSCP', 'SXCP', 'CP','TP','HDE','RP1','RP2');
% $$$ for voterIndex = 1:numVoters
% $$$   figure; plot(information(:,voterIndex),'*-');
% $$$   switch (voterIndex)
% $$$    case 1,title('info SSCP');
% $$$    case 2,title('info SXCP');
% $$$    case 3,title('info CP');
% $$$    case 4,title('info TP');
% $$$    case 5,title('HDE');
% $$$    case 6,title('RP1');
% $$$    case 7,title('RP2');
% $$$   end
% $$$ end
  keyboard
  
coefficients        = zeros(7,1);
coefficientsCell    = cell(1,1);
coefficientsCell{1} = coefficients;
coefficientsIndex   = 1;

EPSILON = 1E-7;
averageWeightedDiffProbs = [];

%for coefficient1 = 0:0.1:1
for coefficient1 = 0:5:10
  remainingSum2 = 1-coefficient1;
  coefficients(1) = coefficient1;
%  for coefficient2 = 0:0.1:remainingSum2
   for coefficient2 = 0:5:10
    remainingSum3 = remainingSum2 - coefficient2;
    coefficients(2) = coefficient2;
    for coefficient3 = 0:5:10
    %    for coefficient3 = 0:0.1:remainingSum3
      remainingSum4 = remainingSum3 - coefficient3;
      coefficients(3) = coefficient3;
%      for coefficient4 = 0:0.1:remainingSum4
       for coefficient4 = 0:5:10
        remainingSum5 = remainingSum4 - coefficient4;
	coefficients(4) = coefficient4;
%	for coefficient5 = 0:0.1:remainingSum5
        for coefficient5 = 0:5:10
          remainingSum6 = remainingSum5 - coefficient5;
	  coefficients(5) = coefficient5;
%	  for coefficient6 = 0:0.1:remainingSum6
          for coefficient6 = 0:5:10
%	    coefficient7 = remainingSum6-coefficient6;
%	    assert ((coefficient1 >= -EPSILON) & (coefficient1 <= 1));
%	    assert ((coefficient2 >= -EPSILON) & (coefficient2 <= 1));
%	    assert ((coefficient3 >= -EPSILON) & (coefficient3 <= 1));
%	    assert ((coefficient4 >= -EPSILON) & (coefficient4 <= 1));
%	    assert ((coefficient5 >= -EPSILON) & (coefficient5 <= 1));
%	    assert ((coefficient6 >= -EPSILON) & (coefficient6 <= 1));
%	    assert ((coefficient7 >= -EPSILON) & (coefficient7 <= 1));
%	    assert (abs(coefficient1 + coefficient2 + coefficient3 + ...
%		     coefficient4 + coefficient5 + coefficient6 + ...
%		     coefficient7 - 1) < EPSILON);
	    
	    coefficients(6) = coefficient6;
	    
	    for coefficient7 = 0:5:10
	    coefficients(7) = coefficient7;

	    sum_coefficients = 0;
	    for coefficientIndex = 1:7
	      sum_coefficients = sum_coefficients + coefficients(coefficientIndex);
	    end
	    
	    if (sum_coefficients == 0)
	      continue;
	    end
	    
	    sumWeightedDiffProb = 0;
	      
	    for peakIndex = 1:numPeaks
	      for voterIndex = 1:numVoters
		maxProb = max(voter{voterIndex}(peakIndex,:));
		if (voter{voterIndex}(peakIndex,peakIndex)  < maxProb)
		    diffProb = voter{voterIndex}(peakIndex,peakIndex) ...
			- maxProb;
		else
		  secondMaxProb = secondMax(voter{voterIndex}(peakIndex,:));
		  diffProb      = voter{voterIndex}(peakIndex,peakIndex) - secondMaxProb;
		end
		weightedDiffProb    = diffProb * coefficients(voterIndex)/sum_coefficients;
		sumWeightedDiffProb = sumWeightedDiffProb + weightedDiffProb;
	      end
	    end
	    coefficientsCell{coefficientsIndex}         = coefficients;
	    averageWeightedDiffProb                     = sumWeightedDiffProb / numPeaks;
	    averageWeightedDiffProbs(coefficientsIndex) = ...
		averageWeightedDiffProb;
	    coefficientsIndex                           = coefficientsIndex + 1;
	    end
	    end %coefficient6
	end
      end
    end
  end
end

figure; plot(averageWeightedDiffProbs,'*');
title('average weighted diff probabilities');
xlabel('coefficient index');
ylabel('weighted average probability difference of the correct assignment');
[maxV, maxI]=  max(averageWeightedDiffProbs);




fprintf(1, 'the max. probability difference that distinguishes the correct assignment has an average of %f weighted probability and has the following coefficients:\n',maxV);

for i = 1:numVoters
  fprintf(1, '%f ',coefficientsCell{maxI}(i));
end

fprintf(1, '\n');


