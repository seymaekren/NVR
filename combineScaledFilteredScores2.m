function combineScaledFilteredScores2

%this one uses scores, not probabilities. negates them.

fprintf(1, 'please cd to the directory where the files are.\n');
keyboard

[label_correct   s1c ] = textread('filteredCorrectScores.txt.scale.predict.polynomial.balanced.withScores','%d %f');
[label_incorrect s1i ] = textread('filteredIncorrectScores.txt.scale.predict.polynomial.balanced.withScores','%d %f');


fcsi                                  = load ('filteredCorrectScoreIndices.txt');
fisi                                  = load ('filteredIncorrectScoreIndices.txt');

load sizeOfScoringMatrix.txt

outputFilename                        = 'combinedScoringMatrix_scaled_polynomial_balanced_withScores.txt';

numPeaks                              = sizeOfScoringMatrix(1,1);
numResidues                           = sizeOfScoringMatrix(1,2);

fprintf(1, 'numPeaks = %d numResidues = %d\n', numPeaks, numResidues);

combinedScoringMatrix = zeros(numPeaks,numResidues);
unfilteredEntries     = zeros(numPeaks,numResidues);

VERY_LARGE_VALUE      = 1E+9;

for i = 1:size(fcsi,1)
  unfilteredEntries(fcsi(i,1),fcsi(i,2)) = 1;
end

for i = 1:size(fisi,1)
  unfilteredEntries(fisi(i,1),fisi(i,2)) = 1;
end

for i = 1:numPeaks
  for j = 1:numResidues
    if (unfilteredEntries(i,j) == 0)
      combinedScoringMatrix(i,j) = VERY_LARGE_VALUE;
    end
  end
end

for i = 1:length(fcsi)
  peakIndex     = fcsi(i,1); residueIndex = fcsi(i,2);
  combinedScore = computeMinusScore(i, s1c);
  combinedScoringMatrix(peakIndex,residueIndex) = combinedScore;
end

for i = 1:length(fisi)
  peakIndex     = fisi(i,1); residueIndex = fisi(i,2);
  combinedScore = computeMinusScore(i, s1i);
  combinedScoringMatrix(peakIndex,residueIndex) = combinedScore;
end

fid = fopen(outputFilename,'w');
for i = 1:numPeaks
  for j = 1:numResidues-1
    fprintf(fid, '%f\t', combinedScoringMatrix(i,j));
  end
  fprintf(fid, '%f\n', combinedScoringMatrix(i,numResidues));
end
fclose(fid);

fprintf(1, 'the computation is finished. Please check %s\n', outputFilename);
keyboard

function combinedScore =  computeMinusScore(indexInFile, s1)


  combinedScore = -(s1(indexInFile));



