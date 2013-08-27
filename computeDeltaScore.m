function score = computeDeltaScore(origMASTER, newMASTER, voter, ...
				   numVoters, peakIndices, residueIndices);

[numPeaks,numResidues] = size(origMASTER);
score = 0;
for peakIndex = 1:numPeaks
  origResidue = find(origMASTER(peakIndex,:));
  newResidue  = find(newMASTER(peakIndex,:));

  if (isempty(origResidue) & ~isempty(newResidue))

    relPeakIndex = find(peakIndices == peakIndex);
    relResidueIndex = find(residueIndices == newResidue);
    
    assert ((~isempty(relPeakIndex)) & (~isempty(relResidueIndex)));
    
    for voterIndex = 1:numVoters
      score = score - log(voter{voterIndex}(relPeakIndex,relResidueIndex));
    end
  
  elseif (isempty(origResidue) & isempty(newResidue))
    %do nothing
  elseif (~isempty(origResidue) & ~isempty(newResidue))
    assert (origResidue == newResidue);
  else
    fprintf(1, 'we dont accept cases where the peak is assigned inthe orig matrix and is unassigned in the new matrix\n');
    assert(0);
  end
end
