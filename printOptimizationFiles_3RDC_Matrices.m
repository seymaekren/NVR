function printOptimizationFiles_3RDC_Matrices (CP, SXCP, SSCP, TP, HDE,  ...
					       RP1, RP2, RP3, NOES, ALLDISTS, COLIN, ...
					       NTH)

VERY_LARGE_NUMBER     = 1E9;
combinedScoringMatrix = zeros(size(CP));
fid                   = fopen('combinedScoringMatrix_3RDC_Matrices.txt','w');
fprintf(1, 'printing to combinedScoringMatrix_3RDC_Matrices.txt\n');

for rowIndex = 1:size(combinedScoringMatrix,1)

  checkEntry(CP(rowIndex,rowIndex),   rowIndex, 'CP');
  checkEntry(SXCP(rowIndex,rowIndex), rowIndex, 'SXCP');
  checkEntry(SSCP(rowIndex, rowIndex),rowIndex, 'SSCP');
  checkEntry(TP(rowIndex,rowIndex),   rowIndex, 'TP');
  checkEntry(HDE(rowIndex,rowIndex),  rowIndex, 'HDE');
  checkEntry(RP1(rowIndex,rowIndex),  rowIndex, 'RP1');
  checkEntry(RP2(rowIndex,rowIndex),  rowIndex, 'RP2');
  checkEntry(RP3(rowIndex,rowIndex),  rowIndex, 'RP3');

  for columnIndex = 1:size(combinedScoringMatrix,2)
    
    if ((CP(rowIndex,columnIndex)   == 0) | ...
	(SXCP(rowIndex,columnIndex) == 0) | ...
	(SSCP(rowIndex,columnIndex) == 0) | ...
	(TP(rowIndex,columnIndex)   == 0) | ...
	(HDE(rowIndex,columnIndex)  == 0) | ...
	(RP1(rowIndex,columnIndex)  == 0) | ...
	(RP2(rowIndex,columnIndex)  == 0) | ...
	(RP3(rowIndex,columnIndex)  == 0))
      
      combinedScoringMatrix(rowIndex,columnIndex) = ...
	  VERY_LARGE_NUMBER;
      
    else
      
%      fprintf(fid, 'before: %f ', combinedScoringMatrix(rowIndex, columnIndex));
      
      combinedScoringMatrix(rowIndex,columnIndex) = - log(CP(rowIndex,columnIndex)) ...
	  -log(SXCP(rowIndex, columnIndex)) ...
	  -log(SSCP(rowIndex,columnIndex)) ...
	  -log(TP(rowIndex, columnIndex))...
	  -log(HDE(rowIndex,columnIndex))...
	  -log(RP1(rowIndex,columnIndex))...
	  -log(RP2(rowIndex,columnIndex))...
	  -log(RP3(rowIndex,columnIndex));
      
      %      fprintf(fid, 'after once: %f ', combinedScoringMatrix(rowIndex, columnIndex));
      
      
    end
    
    fprintf(fid, '%f ', combinedScoringMatrix(rowIndex, columnIndex));
    
  end
  fprintf(fid, '\n');
end

fclose(fid);
  
printNOE_List(NOES);
printBinaryDistances(COLIN, NTH, ALLDISTS);
fprintf(1, 'printed all 3 files required for optimization.\n');



