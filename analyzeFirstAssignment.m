numPeaks = size(foundMASTERs{1},1); foundIt = 0; maxNumCorrectPeaks ...
    = 0;

for resultIndex = 1:totalNumAssignments

  numCorrectPeaks = 0;

  for peakIndex = 1:numPeaks

    residueIndex = find(foundMASTERs{resultIndex}(peakIndex, :));

    if (isempty(residueIndex))
      continue;
    end

    if (peakIndex == residueIndex)
      numCorrectPeaks = numCorrectPeaks + 1;
    else
%      fprintf(1, 'peak %d is incorrectly assigned to residue %d\n',peakIndex,residueIndex);
    end
  
  end

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
  
  fprintf(1, 'numCorrectPeaks = %d\n',numCorrectPeaks);
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
%  fprintf(1, 'did not find the correct assignment among %d results\n',totalNumAssignments);
  fprintf(1, 'found max %d correct peaks (out of %d) at position %d.\n', ...
	  maxNumCorrectPeaks, sum(sum(foundMASTERs{foundPosition})),foundPosition);
%  fprintf(1,'write return to continue.\n');
%  figure;
%  keyboard
end