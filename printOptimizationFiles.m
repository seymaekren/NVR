function printOptimizationFiles (CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
				 RP1, RP2, NOES, ALLDISTS, COLIN, ...
				 NTH, useBayesianMatrix)

VERY_LARGE_NUMBER     = 1E9;
combinedScoringMatrix = zeros(size(CP));
fid                   = fopen('combinedScoringMatrix.txt','w');
fprintf(1, 'printing to combinedScoringMatrix.txt\n');

for rowIndex = 1:size(combinedScoringMatrix,1)
%  fprintf(fid, '%d: ', rowIndex);
  
  
  
  for columnIndex = 1:size(combinedScoringMatrix,2)
  
    if (rowIndex == columnIndex)
    
      checkEntry(CP(rowIndex,rowIndex),   rowIndex, 'CP');
      checkEntry(SXCP(rowIndex,rowIndex), rowIndex, 'SXCP');
      checkEntry(SSCP(rowIndex, rowIndex),rowIndex, 'SSCP');
      checkEntry(TP(rowIndex,rowIndex),   rowIndex, 'TP');
      checkEntry(HDE(rowIndex,rowIndex),  rowIndex, 'HDE');
      if (refineWithRDCs)
	checkEntry(RP1(rowIndex,rowIndex),  rowIndex, 'RP1');
	checkEntry(RP2(rowIndex,rowIndex),  rowIndex, 'RP2');
      end
    
    end
    
    
    if ((CP(rowIndex,columnIndex)   == 0) | ...
	(SXCP(rowIndex,columnIndex) == 0) | ...
	(SSCP(rowIndex,columnIndex) == 0) | ...
	(TP(rowIndex,columnIndex)   == 0) | ...
	(HDE(rowIndex,columnIndex)  == 0))
      
      combinedScoringMatrix(rowIndex,columnIndex) = ...
	  VERY_LARGE_NUMBER;
      
    elseif (refineWithRDCs & ((RP1(rowIndex,columnIndex) == 0) | (RP2(rowIndex,columnIndex)  == 0)))
      
      combinedScoringMatrix(rowIndex,columnIndex) = VERY_LARGE_NUMBER;
    
    else
      
%      fprintf(fid, 'before: %f ', combinedScoringMatrix(rowIndex, columnIndex));
      
      combinedScoringMatrix(rowIndex,columnIndex) = - log(CP(rowIndex,columnIndex)) ...
	  -log(SXCP(rowIndex, columnIndex)) ...
	  -log(SSCP(rowIndex,columnIndex)) ...
	  -log(TP(rowIndex, columnIndex))...
	  -log(HDE(rowIndex,columnIndex));
      
      %      fprintf(fid, 'after once: %f ', combinedScoringMatrix(rowIndex, columnIndex));
      
      
      if (refineWithRDCs)
	combinedScoringMatrix(rowIndex,columnIndex) = ...
	    combinedScoringMatrix(rowIndex, columnIndex) ...
	    - log(RP1(rowIndex,columnIndex)) ...
	    - log(RP2(rowIndex,columnIndex));
	
%	fprintf(fid, 'after twice: %f ', combinedScoringMatrix(rowIndex, columnIndex));
      end
    end
    
    fprintf(fid, '%f ', combinedScoringMatrix(rowIndex, columnIndex));
    
    if (useBayesianMatrix)
      fprintf(fid, '%f ', bayesianMatrix(rowIndex, columnIndex));
    end
  end
  fprintf(fid, '\n');
end

fclose(fid);
    

printNOE_List(NOES);
printBinaryDistances      (COLIN, NTH, ALLDISTS);
printBinaryDistancesForQAP(COLIN, NTH, ALLDISTS);



