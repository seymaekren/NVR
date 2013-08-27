function  newH_Value = computeH_Value(newAssigntable, newVoter, ...
				  newNumVoters);

newH_Value = 0;

for peakIndex = 1:size(newAssigntable,1)

  residueIndices = find(newAssigntable(peakIndex,:));

  %the peak has at least one residue it can be assigned to.
  
  assert (~isempty(residueIndices));

  minScore = 0;
   
  for j = 1:length(residueIndices)

    residueIndex = residueIndices(j);
  
    score = 0;
    
    for voterIndex = 1:newNumVoters
%      score = score - log(newVoter{voterIndex}(peakIndex,
%      residueIndex));
      score = score - (newVoter{voterIndex}(peakIndex, residueIndex));
    end
  
    if (j == 1) | (score < minScore)
      minScore = score;
    end

%    fprintf(1, 'score of assigning peak #%d to residue#%d is %f\n',peakIndex,residueIndex,score);
    
  end
 
%  fprintf(1, 'score of assigning peak #%d is %f\n',peakIndex,score);
  
  newH_Value = newH_Value + minScore;
  
end