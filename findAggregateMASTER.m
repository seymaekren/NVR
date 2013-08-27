function aggregateMASTER = findAggregateMASTER(overallMASTER, ...
					       totalNumAssignments)



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


%confidenceThreshold              = 0;
[totalNumPeaks,totalNumResidues] = size(overallMASTER);

if (totalNumPeaks ~= totalNumResidues)
  squareOverallMaster = zeros(max(totalNumPeaks,totalNumResidues));
  squareOverallMaster(1:totalNumPeaks,1:totalNumResidues) = ...
      overallMASTER;
else
  squareOverallMaster = overallMASTER;
end

[h,C] = hungarian(-squareOverallMaster');
debugCost = 0;
for peakIndex = 1:length(h)
  debugCost = debugCost + squareOverallMaster(peakIndex,h(peakIndex));
end

assert (debugCost == -C);

aggregateMASTER                  = zeros(totalNumPeaks, totalNumResidues);
for peakIndex= 1:totalNumPeaks
  if (overallMASTER(peakIndex,h(peakIndex)) ~= 0)
    aggregateMASTER(peakIndex,h(peakIndex)) = 1;
  end
end


% $$$ for i = 1:totalNumPeaks
% $$$   [maxVotes, residueIndex] = max(overallMASTER(i,:));
% $$$   if (maxVotes > totalNumAssignments*confidenceThreshold)
% $$$     aggregateMASTER(i,residueIndex) = 1;
% $$$   end
% $$$ end