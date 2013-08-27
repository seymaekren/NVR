
function  newH_Value = computeBetterH_Value(newAssigntable, newVoter, ...
					    newNumVoters)%, unassignedPeakIndices, ...
%					    unassignedResidueIndices);

EPSILON    = 1E-5;

LARGE_VALUE = 1E9;

%assert (size(newAssigntable,1) == length(unassignedPeakIndices));
%assert (size(newAssigntable,2) == length(unassignedResidueIndices));

newH_Value = 0;

minusCombinedProbabilityMatrix = newAssigntable;

for voterIndex = 1:newNumVoters
  
  %minusLogCombinedProbabilityMatrix = minusLogCombinedProbabilityMatrix.*newVoter{voterIndex};
  minusCombinedProbabilityMatrix = minusCombinedProbabilityMatrix.*newVoter{voterIndex};
  
end

for peakIndex = 1:size(minusCombinedProbabilityMatrix,1)
  for residueIndex = 1:size(minusCombinedProbabilityMatrix,2)
%    if (minusCombinedProbabilityMatrix(peakIndex,residueIndex) ~= 0)
     if (newAssigntable(peakIndex,residueIndex) ~= 0)
      %      minusLogCombinedProbabilityMatrix(peakIndex,residueIndex) = - ...
      %
      %	  log(minusLogCombinedProbabilityMatrix(peakIndex,residueIndex));

      assert ( minusCombinedProbabilityMatrix(peakIndex,residueIndex)>=0);
      minusCombinedProbabilityMatrix(peakIndex,residueIndex) = - ...
      	  (minusCombinedProbabilityMatrix(peakIndex,residueIndex));
    else
      minusCombinedProbabilityMatrix(peakIndex,residueIndex) = LARGE_VALUE;
    end
  end
end
%[h,c] = hungarian(-combinedProbabilityMatrix');

squareMinusLogCombinedProbabilityMatrix = ...
    LARGE_VALUE*ones(max(size(minusCombinedProbabilityMatrix)));
squareMinusLogCombinedProbabilityMatrix(1:size(minusCombinedProbabilityMatrix, ...
					       1),1: ...
					size(minusCombinedProbabilityMatrix,2)) ...
    = minusCombinedProbabilityMatrix;


[h,testH_Value] = hungarian(squareMinusLogCombinedProbabilityMatrix');

newH_Value      = 0;

for peakIndex = 1:size(newAssigntable,1)
  newH_Value = newH_Value + minusCombinedProbabilityMatrix(peakIndex,h(peakIndex));
%  assert (minusCombinedProbabilityMatrix(peakIndex,h(peakIndex)) ~= ...
%	  LARGE_VALUE);
%  fprintf(1, 'assigning %d to %d with score %f in betterH_ValueComputation\n',unassignedPeakIndices(peakIndex),unassignedResidueIndices(h(peakIndex)),minusCombinedProbabilityMatrix(peakIndex,h(peakIndex)));
end


assert (abs(testH_Value - (newH_Value+LARGE_VALUE* ...
			   (size(squareMinusLogCombinedProbabilityMatrix,1) ...
			    - size(newAssigntable,1)))) < EPSILON);

fprintf(1, 'computed better H value. It is %f\n',newH_Value);