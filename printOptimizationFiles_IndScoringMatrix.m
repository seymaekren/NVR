function printOptimizationFiles_IndScoringMatrix (CP, NOES, ALLDISTS, COLIN, ...
						  NTH, ...
						  scoringMatrixFilename, ...
						  MATRIX2, MATRIX3)

printingRDC = 0; printingCS = 0;
if (nargin == 7)
  fprintf(1, 'printing RDC scoring matrices\n');
  printingRDC = 1;
%  keyboard
elseif (nargin == 8)
  fprintf(1, 'printing CS scoring matrices\n');
  printingCS = 1;
%  keyboard
end


VERY_LARGE_NUMBER     = 1E9;
combinedScoringMatrix = zeros(size(CP));
fid                   = fopen(scoringMatrixFilename,'w');
fprintf(1, 'printing to %s\n',scoringMatrixFilename);

for rowIndex = 1:size(combinedScoringMatrix,1)
%  fprintf(fid, '%d: ', rowIndex);

  checkEntry(CP(rowIndex,rowIndex),   rowIndex, 'CP');

  if (printingRDC == 1) 
    checkEntry(MATRIX2(rowIndex,rowIndex), rowIndex, 'MATRIX2');
  end
  
  if (printingCS == 1) 
    checkEntry (MATRIX3(rowIndex,rowIndex), rowIndex, 'MATRIX3');
  end
  
  for columnIndex = 1:size(combinedScoringMatrix,2)
    
    if ((CP(rowIndex,columnIndex)   == 0))
      combinedScoringMatrix(rowIndex,columnIndex) = ...
	  VERY_LARGE_NUMBER;
    elseif ((printingRDC == 1) & (MATRIX2(rowIndex,columnIndex) == 0) )
      combinedScoringMatrix(rowIndex,columnIndex) = ...
	  VERY_LARGE_NUMBER;
    elseif ((printingCS == 1) & (MATRIX3(rowIndex,columnIndex) == 0) )
      combinedScoringMatrix(rowIndex,columnIndex) = ...
	  VERY_LARGE_NUMBER;
    else
      combinedScoringMatrix(rowIndex,columnIndex) = - log(CP(rowIndex,columnIndex));
      if (printingRDC == 1)
	combinedScoringMatrix(rowIndex,columnIndex) = combinedScoringMatrix(rowIndex,columnIndex) - log(MATRIX2(rowIndex,columnIndex));
      elseif (printingCS == 1)
	combinedScoringMatrix(rowIndex,columnIndex) = ...
	    combinedScoringMatrix(rowIndex,columnIndex) - log(MATRIX2(rowIndex,columnIndex)) - log(MATRIX3(rowIndex,columnIndex));
      end
    end
    
    fprintf(fid, '%f ', combinedScoringMatrix(rowIndex, columnIndex));
    
  end
  fprintf(fid, '\n');
end

fclose(fid);
  

printNOE_List(NOES);
printBinaryDistances(COLIN, NTH, ALLDISTS);


fprintf(1, 'printed all 3 files required for optimization.\n');

return

