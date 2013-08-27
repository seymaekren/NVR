
%combinedScoringMatrix = load ('/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/100PercentCorrectAssignments/combinedScoringMatrix.txt');
%MASTER = load ('/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/hSRI_100PercentCorrect.txt');

load combinedScoringMatrix.txt
MASTER = load ('193L.txt');

count = 0;
[numRows, numColumns] = size(MASTER);
MASTER = MASTER(1:numRows, 1:numColumns - 3);
fprintf(1, 'reduced numColumns of MASTER to %d\n', size(MASTER,2));
score = 0; numAssigned  = 0;
for i = 1:numRows
  if (MASTER(i,i) == 1)
    count = count + 1;
  end
  residue = find(MASTER(i,:));
  if (~isempty(residue))
    numAssigned = numAssigned + 1;
    score = score + combinedScoringMatrix(i,residue);
  end
end
fprintf(1, 'assignment accuracy = %f\n', count/numRows);
fprintf(1, 'precision = %f\n', count/numAssigned);
fprintf(1, 'score = %f\n', score);


