function combineScaledFilteredScores

fprintf(1, 'please cd to the directory where the files are.\n');
keyboard

[alpha_incorrect xi1s xi2s xi3s xi4s] = textread('filteredIncorrectScores.txt.scale','%f %s %s %s %s');
[alpha_correct   xc1s xc2s xc3s xc4s] = textread('filteredCorrectScores.txt.scale','%f %s %s %s %s');

fcsi                                  = load ('filteredCorrectScoreIndices.txt');
fisi                                  = load ('filteredIncorrectScoreIndices.txt');

load sizeOfScoringMatrix.txt

outputFilename                        = 'combinedScoringMatrix_scaled.txt';

numPeaks                              = sizeOfScoringMatrix(1,1);
numResidues                           = sizeOfScoringMatrix(1,2);

fprintf(1, 'numPeaks = %d numResidues = %d\n', numPeaks, numResidues);

combinedScoringMatrix = zeros(numPeaks,numResidues);
unfilteredEntries     = zeros(numPeaks,numResidues);

weight                = -[-3.7656   -1.9749    0.1364   -6.5261];

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
  combinedScore = computeCombinedWeightedScore(i, xc1s,xc2s,xc3s,xc4s, ...
					       weight);
  combinedScoringMatrix(peakIndex,residueIndex) = combinedScore;
end

for i = 1:length(fisi)
  peakIndex     = fisi(i,1); residueIndex = fisi(i,2);
  combinedScore = computeCombinedWeightedScore(i, xi1s,xi2s,xi3s,xi4s, ...
					       weight);
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

function combinedScore =  computeCombinedWeightedScore(indexInFile, xi1s,xi2s,xi3s,xi4s, ...
						  weight);
  [T,R] = strtok(xi1s(indexInFile),':');
  xi1   = strtok(R,':');
  xi1   = str2num(xi1{1});
  
  [T,R] = strtok(xi2s(indexInFile),':');
  xi2   = strtok(R,':');
  xi2   = str2num(xi2{1});
  
  
  [T,R] = strtok(xi3s(indexInFile),':');
  xi3   = strtok(R,':');
  xi3   = str2num(xi3{1});
  
  
  [T,R] = strtok(xi4s(indexInFile),':');
  xi4   = strtok(R,':');
  xi4   = str2num(xi4{1});
  
  
  xi    = [xi1 xi2 xi3 xi4];
  combinedScore = dot(xi,weight);



