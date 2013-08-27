function combineScaledFilteredScoresRBF

fprintf(1, 'please cd to the directory where the files are.\n');
keyboard

[label_correct   p1c p2c] = textread('filteredCorrectScores_1UBQ.txt.scale.predict','%d %f %f');
[label_incorrect p1i p2i] = textread('filteredIncorrectScores_1UBQ.txt.scale.predict','%d %f %f');


fcsi                                  = load ('filteredCorrectScoreIndices.txt');
fisi                                  = load ('filteredIncorrectScoreIndices.txt');

load sizeOfScoringMatrix.txt

outputFilename                        = 'combinedScoringMatrix_scaled_rbf_balanced.txt';

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
  combinedScore = computeMinusLogScore(i, p1c);
  combinedScoringMatrix(peakIndex,residueIndex) = combinedScore;
end

for i = 1:length(fisi)
  peakIndex     = fisi(i,1); residueIndex = fisi(i,2);
  combinedScore = computeMinusLogScore(i, p1i);
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

function combinedScore =  computeMinusLogScore(indexInFile, p1)


  combinedScore = -log(p1(indexInFile));



