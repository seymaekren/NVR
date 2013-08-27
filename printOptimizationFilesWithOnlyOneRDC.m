function printOptimizationFilesWithOnlyOneRDC (CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
				 RP1, RP2, NOES, ALLDISTS, COLIN, ...
				 NTH, useBayesianMatrix)

VERY_LARGE_NUMBER     = 1E9;
combinedScoringMatrix = zeros(size(CP));
fid                   = fopen('combinedScoringMatrix.txt','w');
fprintf(1, 'printing to combinedScoringMatrix.txt\n');

for rowIndex = 1:size(combinedScoringMatrix,1)
%  fprintf(fid, '%d: ', rowIndex);
  if (CP(rowIndex,rowIndex) == 0) 
    fprintf(1, 'CP matrix entry empty for peak %d\n', ...
	    rowIndex);
  end
  
  if (SXCP(rowIndex,rowIndex) == 0) 
    fprintf(1, 'SXCP matrix empty for peak %d\n', rowIndex);
  end
  
  if (SSCP(rowIndex, rowIndex) == 0)
    fprintf(1, 'SSCP matrix empty for peak %d\n', rowIndex);
  end
  
  if (TP(rowIndex,rowIndex) == 0) 
    fprintf(1, 'TP matrix empty for peak %d\n', rowIndex);
  end
  
  if (HDE(rowIndex,rowIndex) == 0)
    fprintf(1, 'HDE matrix empty for peak %d\n', rowIndex);
  end
  
  if (refineWithRDCs)
    if (RP1(rowIndex,rowIndex) == 0) 
      fprintf(1, 'RP1 matrix empty for peak %d\n', rowIndex);
    end
    if (RP2(rowIndex,rowIndex) == 0)
      fprintf(1, 'RP2 matrix empty for peak %d\n', rowIndex);
    end
  end
  
  for columnIndex = 1:size(combinedScoringMatrix,2)
    
    if ((CP(rowIndex,columnIndex)   == 0) | ...
	(SXCP(rowIndex,columnIndex) == 0) | ...
	(SSCP(rowIndex,columnIndex) == 0) | ...
	(TP(rowIndex,columnIndex)   == 0) | ...
	(HDE(rowIndex,columnIndex)  == 0))
      
      combinedScoringMatrix(rowIndex,columnIndex) = ...
	  VERY_LARGE_NUMBER;
      
    elseif (refineWithRDCs & (RP2(rowIndex,columnIndex) == 0) )
      
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
  
fid = fopen('NOE_List.txt','w');
fprintf(1, 'printing to NOE_List.txt\n');

for peak1Index = 1:size(NOES,1)
  listOfNOEs = find(NOES(peak1Index,:));
  if (~isempty(listOfNOEs))
    fprintf(fid, '%d %d ', length(listOfNOEs), peak1Index);
    for peak2Index = 1:size(NOES,2)
      %	fprintf(fid, '%d ', NOES(ROWIN(rowIndex),
      %	COLIN(columnIndex)));
      if (NOES(peak1Index,peak2Index) == 1)
	fprintf(fid, ' %d ',peak2Index);
      end
    end
    fprintf(fid, '\n');
  end
end

fclose(fid);
    
fid = fopen('BinaryDistances.txt','w');
fprintf(1, 'printing to BinaryDistances.txt\n');
for residue1Index = 1:length(COLIN)
  for residue2Index =  1:length(COLIN)
    if (ALLDISTS(COLIN(residue1Index),COLIN(residue2Index)) > NTH)
      fprintf(fid, '0 ');
    else
      fprintf(fid, '1 ');
    end
  end	
  fprintf(fid, '\n');
end
fclose(fid);
fprintf(1, 'printed all 3 files required for optimization.\n');
return

