function [score, accuracy, precision] = computeScore2(combinedScoringMatrix, MASTER);

%MASTER should be numPeaks x numResidues

persistent MAX_VALUE;

if (isempty(MAX_VALUE))
  MAX_VALUE = 1E+8;
  fprintf(1, 'set MAX_VALUE to 1E+8 in computeScore2\n');
end

count = 0;
[numRows, numColumns]   = size(MASTER);
score = 0; numAssigned  = 0;

for i = 1:numRows
  if (MASTER(i,i) == 1)
    count = count + 1;
  end
  residue = find(MASTER(i,:));
  if (~isempty(residue))
    numAssigned = numAssigned + 1;
    if (combinedScoringMatrix(i,residue) > MAX_VALUE)
      score     = inf;
      accuracy  = -1;
      precision = -1;
      return;
    end
    score       = score       + combinedScoringMatrix(i,residue);
  end
end
accuracy = count/numRows;
precision = count/numAssigned;

%fprintf(1, 'assignment accuracy = %f\n', count/numRows);
%fprintf(1, 'precision = %f\n', count/numAssigned);
%fprintf(1, 'score = %f\n', score);


