function score = computeScore(MASTER, voters, numVoters)

score = 0;
numPeaks = size(MASTER,1);

for i = 1:numVoters
  for peakIndex = 1:numPeaks
    residueIndex = find(MASTER(peakIndex,:));
    assert (length(residueIndex) == 1);
    score = score - log(voters{i}(peakIndex,residueIndex));
  end
end