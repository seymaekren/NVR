function bayesianScore = computeBayesianScore(inputMASTER)
%load allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat

dbstop if error;
dbstop if warning;



persistent meanCorrect meanIncorrect standardDeviationCorrect 
persistent standardDeviationIncorrect voter numVoters
persistent EPSILON LARGE_VALUE numPeaks numResidues
persistent USE_CDF USE_HIST
persistent N_histCorrect X_histCorrect N_histIncorrect X_histIncorrect

if (isempty(meanCorrect))

  USE_CDF = 0; USE_HIST = 1;
  
  fprintf(1, 'loading statistics...\n');
  
%  load ('~/Workdir/allPeaksAssignmentsEnvironment.mat');
  load ('~/Workdir/allPeaksNewAssignmentEnvironment-forBayesianScoring.mat');
  EPSILON                          = 1E-3;
  LARGE_VALUE                      = 1E5;
  [numPeaks,numResidues]           = size(MASTER);
  
  voter                            = initializeVoters(SSCP,SXCP,CP,TP,HDE, RP1, RP2);
%  numVoters                        = 5;
  numVoters                        = 7;


  fprintf(1,'numVoters = %f\n',numVoters);
  
  
  correctAssignmentProbability   = cell(numVoters, 1);
  incorrectAssignmentProbability = cell(numVoters, 1);
  for voterIndex = 1:numVoters
    correctAssignmentProbability{voterIndex}   = [];
    incorrectAssignmentProbability{voterIndex} = [];
  end
  
  for peakIndex=1:numPeaks
    correctResidueIndex = peakIndex;
    
    for voterIndex = 1:numVoters
      correctAssignmentProbability{voterIndex} =...
	  [correctAssignmentProbability{voterIndex} voter{voterIndex}(peakIndex,correctResidueIndex)];
      for residueIndex = 1:numResidues
	if (residueIndex == correctResidueIndex)
	  continue;
	end
	incorrectAssignmentProbability{voterIndex} = [incorrectAssignmentProbability{voterIndex} voter{voterIndex}(peakIndex,residueIndex)];
      end
    end
  end
  
  meanCorrect                = zeros(numVoters,1);
  meanIncorrect              = zeros(numVoters,1);
  standardDeviationCorrect   = zeros(numVoters,1);
  standardDeviationIncorrect = zeros(numVoters,1);
  
  for voterIndex = 1:numVoters
    meanCorrect(voterIndex)                   = mean(correctAssignmentProbability{voterIndex});
    standardDeviationCorrect(voterIndex)      = std (correctAssignmentProbability{voterIndex});
    meanIncorrect(voterIndex)                 = mean(incorrectAssignmentProbability{voterIndex});
    standardDeviationIncorrect(voterIndex)    = std (incorrectAssignmentProbability{voterIndex});
  end

  
  N_histCorrect = cell(numVoters,1); X_histCorrect = cell(numVoters, ...
						  1);
  N_histIncorrect = cell(numVoters,1); X_histIncorrect = cell(numVoters, ...
						  1);
  for voterIndex = 1:numVoters
    [N_histCorrect{voterIndex},X_histCorrect{voterIndex}] = ...
	hist(correctAssignmentProbability{voterIndex});
    
    [N_histIncorrect{voterIndex}, X_histIncorrect{voterIndex}] = hist(incorrectAssignmentProbability{voterIndex});
    
% $$$   figure;
% $$$   hist(correctAssignmentProbability{voterIndex});
% $$$   title('correct assignment probabilities');
% $$$ 
% $$$   figure;
% $$$   hist(incorrectAssignmentProbability{voterIndex});
% $$$   title('incorrect assignment probabilities');
% $$$ 
  end

end

bayesianScore = 0;


% $$$ 
% $$$ 
% $$$ %keyboard
% $$$ 
% $$$ distCorrect = cell(numVoters,1);
% $$$ distIncorrect = cell(numVoters,1);
% $$$ 
% $$$ for voterIndex = 1:numVoters
% $$$   distCorrect{voterIndex}   = [];
% $$$   distIncorrect{voterIndex} = [];
% $$$ 
% $$$   for x = 0:0.01:1
% $$$     distCorrect{voterIndex}  =[distCorrect{voterIndex}   normpdf(x,meanCorrect(voterIndex),standardDeviationCorrect(voterIndex))];
% $$$     distIncorrect{voterIndex}=[distIncorrect{voterIndex} normpdf(x,meanIncorrect(voterIndex), standardDeviationIncorrect(voterIndex))];
% $$$   end
% $$$ 
% $$$ end
% $$$ 
% $$$ x = [0:0.01:1];
% $$$ for voterIndex = 1:numVoters
% $$$   figure; plot(x,distCorrect{voterIndex},'r-*',x,distIncorrect{voterIndex},'-o')
% $$$   legend('correct assignment','incorrect assignment');
% $$$ end



relativeProbabilityCorrectGivenVoters = ones(numPeaks,numResidues);
pCorrect                              = 1/numResidues; 





for peakIndex = 1:numPeaks
%  for residueIndex = 1:numResidues
  residueIndex = find(inputMASTER(peakIndex,:));
  %assert (~isempty(residueIndex));
  
  if (isempty(residueIndex))
    continue;
  end
  
  num = zeros(numVoters,1); denum = zeros(numVoters,1);
    for voterIndex = 1:numVoters

      if (USE_CDF)
	num(voterIndex) = normcdf(voter{voterIndex}(peakIndex,residueIndex)+EPSILON,meanCorrect(voterIndex),standardDeviationCorrect(voterIndex))-...
	    normcdf(voter{voterIndex}(peakIndex,residueIndex)-EPSILON,meanCorrect(voterIndex),standardDeviationCorrect(voterIndex));
      elseif (USE_HIST)
	bestBin = -1; minDifference = LARGE_VALUE;
	for histBinIndex = 1:length(X_histCorrect{voterIndex})
	  difference = abs(X_histCorrect{voterIndex}(histBinIndex)-voter{voterIndex}(peakIndex,residueIndex));
	  if (difference < minDifference)
	    bestBin = histBinIndex;
	    minDifference = difference;
	  end
	end
	
	assert (bestBin ~= -1);
	
	num(voterIndex)  = N_histCorrect{voterIndex}(bestBin)/sum(N_histCorrect{voterIndex});
      end
      
      relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex) = ...
	  relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex) * num(voterIndex);
      
      %      keyboard
      
      if (USE_HIST)
      
       bestBin = -1; minDifference = LARGE_VALUE;
        for histBinIndex = 1:length(X_histIncorrect{voterIndex})
 	 difference = abs(X_histIncorrect{voterIndex}(histBinIndex)-voter{voterIndex}(peakIndex,residueIndex));
 	 if (difference < minDifference)
 	   bestBin = histBinIndex;
 	   minDifference = difference;
 	 end
        end
       
        assert (bestBin ~= -1);
       
        
        %       denum(voterIndex) = num(voterIndex) * pCorrect + (1-pCorrect)*N_histIncorrect{voterIndex}(bestBin)/sum(N_histIncorrect{voterIndex});
        denum(voterIndex) = N_histIncorrect{voterIndex}(bestBin)/sum(N_histIncorrect{voterIndex});
       
      elseif (USE_CDF)
	
	denum(voterIndex) = (normcdf(voter{voterIndex}(peakIndex, ...
						       residueIndex)+EPSILON, ...
				     meanIncorrect(voterIndex), ...
				     standardDeviationIncorrect(voterIndex))-normcdf(voter{voterIndex}(peakIndex, ...
						  residueIndex)-EPSILON, ...
						  meanIncorrect(voterIndex), ...
						  standardDeviationIncorrect(voterIndex)));
	
      end
%      denum(voterIndex) = (normcdf(voter{voterIndex}(peakIndex, ...
% 						    residueIndex)+EPSILON, ...
%				  meanCorrect(voterIndex), ...
%				  standardDeviationCorrect(voterIndex))-normcdf(voter{voterIndex}(peakIndex, ...
%						    residueIndex)-EPSILON, ...
%						  meanCorrect(voterIndex), ...
%						  standardDeviationCorrect(voterIndex))) * ...
%	  pCorrect + (normcdf(voter{voterIndex}(peakIndex, ...
%					       residueIndex)+EPSILON, ...
%			     meanIncorrect(voterIndex), ...
%			     standardDeviationIncorrect(voterIndex))-normcdf(voter{voterIndex}(peakIndex, ...
%						  residueIndex)-EPSILON, ...
%						  meanIncorrect(voterIndex), ...
%						  standardDeviationIncorrect(voterIndex)))*(1-pCorrect);
      
     if (denum(voterIndex) < EPSILON)
       denum(voterIndex) = EPSILON;
     end

      assert ((denum(voterIndex) > 0) & (denum(voterIndex) <= 1));
      assert ((num(voterIndex) >= 0) & (num(voterIndex)<=1));
%      pCorrectGivenVoterValue =  (num(voterIndex) * pCorrect/ ...
%				  denum(voterIndex));
%      assert ((pCorrectGivenVoterValue >= 0) & (pCorrectGivenVoterValue<=1));

       relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex) = relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex)/denum(voterIndex);
%      keyboard
    end
    relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex) = ...
	relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex) * pCorrect/(1-pCorrect);
%  end
   bayesianScore = bayesianScore +  relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex);
end

% $$$ relativeProbabilityCorrectGivenVotersForIncorrectAssignments = [];
% $$$ 
% $$$ for peakIndex = 1:numPeaks
% $$$   for residueIndex = 1:numResidues
% $$$     if (peakIndex == residueIndex)
% $$$       continue;
% $$$     end
% $$$     relativeProbabilityCorrectGivenVotersForIncorrectAssignments = ...
% $$$ 	[relativeProbabilityCorrectGivenVotersForIncorrectAssignments relativeProbabilityCorrectGivenVoters(peakIndex,residueIndex)];
% $$$   end
% $$$ end
% $$$ 
% $$$ figure; plot(diag(relativeProbabilityCorrectGivenVoters),'*');
% $$$ title('relativeProbabilityCorrect for correct assignments');
% $$$ 
% $$$ figure; plot(relativeProbabilityCorrectGivenVotersForIncorrectAssignments,'*')
% $$$ title('relative probability correct for incorrect assignments');
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ fprintf(1, 'mean of relative prob for correct assignments=%f\n', ...
% $$$ 	mean(diag(relativeProbabilityCorrectGivenVoters)));
% $$$ fprintf(1, 'std of relative prob for correct assignments=%f\n', ...
% $$$ 	std(diag(relativeProbabilityCorrectGivenVoters)));
% $$$ 
% $$$ fprintf(1, 'mean of relative prob for incorrect assignments=%f\n', ...
% $$$ 	mean(relativeProbabilityCorrectGivenVotersForIncorrectAssignments));
% $$$ fprintf(1, 'std of relative prob for correct assignments=%f\n', ...
% $$$ 	std(relativeProbabilityCorrectGivenVotersForIncorrectAssignments));
% $$$ 
% $$$ keyboard