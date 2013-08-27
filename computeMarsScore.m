function score = computeMarsScore(thisMASTER)

persistent marsScore

if (isempty(marsScore))

  fprintf(1, 'initializing mars score...\n');
  load ('~/Workdir/allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat');
  %load('~/Workdir/bestMarsScoreWithFullyCorrectRDCs.mat');

  %marsScore                        = 3.3 * marsShiftScore + marsRdcScore;
  
end
  
score    = 0;
numPeaks = size(thisMASTER,1);


for peakIndex = 1:numPeaks
  residueIndex = find(thisMASTER(peakIndex,:));
  assert (length(residueIndex) == 1);
  score = score + marsScore(peakIndex,residueIndex);
end



