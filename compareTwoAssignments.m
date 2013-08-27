hSRI_withoutRDC = load ('/home2/apaydin/speech/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/consensusMatrix.txt');


hSRI_withRDC = load ('/home2/apaydin/speech/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/WithRDC_FourthRound/hSRI.txt');


[numPeaks,numColumns] = size(hSRI_withoutRDC);

withoutRDC  = hSRI_withoutRDC(:, numColumns-2:numColumns-1); 
withRDC     = hSRI_withRDC   (:, numColumns-2:numColumns-1); 

consensusMatrix = zeros(numPeaks, numColumns);

same = 0; different = 0;
for i = 1:numPeaks
  residue1 = withoutRDC(i,2);
  residue2 = withRDC   (i,2);
  if ((residue1 ~= 0) & (residue2 ~= 0))
    if (residue1 == residue2)
      same = same + 1;
    else
      different  = different + 1;
    end
  end
end

fprintf(1, 'same assignments=   %d different assignments = %d\n',same,different);