clear all;
%combinedScoringMatrix = load ('/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/100PercentCorrectAssignments/combinedScoringMatrix.txt');
%MASTER = load ('/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/hSRI_100PercentCorrect.txt');

%load combinedScoringMatrix.txt
combinedScoringMatrix = load ('combinedScoringMatrix_3RDC_Matrices.txt');
%MASTER = load ('3LYZ.txt');
[numPeaks,numResidues] = size(combinedScoringMatrix);
MASTER  = zeros(numPeaks,numResidues);
for i = 1:numPeaks
  MASTER(i,i) = 1;
end



count = 0;
score = 0; numAssigned  = 0;
for i = 1:numPeaks
  if (MASTER(i,i) == 1)
    count = count + 1;
  end
  residue = find(MASTER(i,:));
  if (~isempty(residue))
    numAssigned = numAssigned + 1;
    score = score + combinedScoringMatrix(i,residue);
  end
end
fprintf(1, 'assignment accuracy = %f\n', count/numPeaks);
fprintf(1, 'precision = %f\n', count/numAssigned);
fprintf(1, 'score = %f\n', score);


