function g = computeG(MASTER, voter)

[numPeaks,numResidues]  = size(MASTER);
assert (size(voter,1) == numPeaks);
assert (size(voter,2) == numResidues);

g = 0;

for peakIndex = 1:numPeaks
  residueIndex = find(MASTER(peakIndex,:));

  if (~isempty(residueIndex))
    assert (length(residueIndex) == 1);
    g = g - voter(peakIndex,residueIndex);
    fprintf(1, 'assignment of %d to %d. score contrib=  %f\n',peakIndex,residueIndex,voter(peakIndex,residueIndex));
  end

end