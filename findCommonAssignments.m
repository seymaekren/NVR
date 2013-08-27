hSRI_withoutRDC = load ('/home2/apaydin/speech/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/WithoutRDC/hSRI.txt');


hSRI_withRDC = load ('/home2/apaydin/speech/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/WithRDC_FirstRound/hSRI.txt');


[numPeaks,numColumns] = size(hSRI_withoutRDC);

withoutRDC  = hSRI_withoutRDC(:, numColumns-2:numColumns-1); 
withRDC     = hSRI_withRDC   (:, numColumns-2:numColumns-1); 

consensusMatrix = zeros(numPeaks, numColumns);


for i = 1:numPeaks
  residue1 = withoutRDC(i,2);
  residue2 = withRDC(i,2);
  if (residue1 == residue2)
    consensusMatrix(i, residue1)       = 1;
    consensusMatrix(i, numColumns - 1) = residue1;
    fprintf(1, 'i = %d residue1 = residue2 = %d\n', i,residue1);
  else
    consensusMatrix(i, numColumns - 1) = 0;
  end
  
  consensusMatrix(i, numColumns - 2)   = i;

  if ((i == residue1) & (residue1 == residue2))
    consensusMatrix(i, numColumns) = 1;
  else
    consensusMatrix(i, numColumns) = 0;
  end
end

fid = fopen('consensusMatrix.txt','w');
fprintf(1, 'check out consensusMatrix.txt \n');
for i = 1:numPeaks
  for j = 1:numColumns
    fprintf(fid, '%d ', consensusMatrix(i,j));
  end
  fprintf(fid, '\n');
end

fclose(fid);