function printIndividualScores (CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
				 RP1, RP2)

VERY_LARGE_NUMBER     = 1E9;
fid_correct           = fopen('correctScores.txt','w');
fid_correct_indices   = fopen('correctScoreIndices.txt','w');
fid_incorrect         = fopen('incorrectScores.txt','w');
fid_incorrect_indices = fopen('incorrectScoreIndices.txt','w');
fid_sizeOfScoringMatrix = fopen('sizeOfScoringMatrix.txt','w');


fprintf(fid_sizeOfScoringMatrix, '%d %d\n', size(CP,1), size(CP,2));
fclose (fid_sizeOfScoringMatrix);


fprintf(1, 'printing to correctScores.txt and incorrectScores.txt\n');
fprintf(1, 'printing to correctScoreIndices.txt and incorrectScoreIndices.txt\n');

for rowIndex = 1:size(CP,1)

  checkDiagonalEntries(CP,SXCP,SSCP,TP,HDE,RP1,RP2,refineWithRDCs, rowIndex);
 
  for columnIndex = 1:size(CP,2)
  
    [minusLog_CP,minusLog_SXCP, minusLog_SSCP,minusLog_TP, ...
     minusLog_HDE,minusLog_RP1,minusLog_RP2] = computeMinusLogs(CP,SXCP,SSCP,TP,HDE,refineWithRDCs, RP1,RP2, rowIndex, columnIndex,VERY_LARGE_NUMBER);
    
    if (refineWithRDCs)
      if (rowIndex == columnIndex)
	fprintf(fid_correct,   '%f %f %f %f %f %f %f\n', minusLog_CP, ...
		minusLog_SXCP, minusLog_SSCP, minusLog_TP, minusLog_HDE, minusLog_RP1, minusLog_RP2);
	
	fprintf(fid_correct_indices, '%d %d\n', rowIndex, columnIndex);
      else
	fprintf(fid_incorrect, '%f %f %f %f %f %f %f\n', minusLog_CP, ...
		minusLog_SXCP, minusLog_SSCP, minusLog_TP, minusLog_HDE, minusLog_RP1, minusLog_RP2);
	fprintf(fid_incorrect_indices, '%d %d\n', rowIndex, columnIndex);
      end
    else
      if (rowIndex == columnIndex)
	fprintf(fid_correct,   '%f %f %f %f %f\n', minusLog_CP, minusLog_SXCP, minusLog_SSCP, minusLog_TP, minusLog_HDE);
      	fprintf(fid_correct_indices, '%d %d\n', rowIndex, columnIndex);
      else
	fprintf(fid_incorrect, '%f %f %f %f %f\n', minusLog_CP, minusLog_SXCP, minusLog_SSCP, minusLog_TP, minusLog_HDE);
	fprintf(fid_incorrect_indices, '%d %d\n', rowIndex, columnIndex);
      end
    end
  end
end

fclose(fid_correct);
fclose(fid_incorrect);
  

    


function minusLog = computeMinusLog(prob, VERY_LARGE_NUMBER);

assert ((prob>= 0) & (prob <= 1));
if (prob == 0)
  minusLog = VERY_LARGE_NUMBER;
else
  minusLog = -log(prob);
end

function checkDiagonalEntries(CP,SXCP,SSCP,TP,HDE,RP1,RP2,refineWithRDCs, rowIndex);

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

function [minusLog_CP,minusLog_SXCP,minusLog_SSCP,minusLog_TP, ...
	  minusLog_HDE,minusLog_RP1,minusLog_RP2] = ...
    computeMinusLogs(CP,SXCP,SSCP,TP,HDE,refineWithRDCs, RP1,RP2, rowIndex, columnIndex, ...
		     VERY_LARGE_NUMBER);
    
minusLog_CP   = computeMinusLog(CP(rowIndex, columnIndex), VERY_LARGE_NUMBER);
minusLog_SXCP = computeMinusLog(SXCP(rowIndex, columnIndex), ...
				 VERY_LARGE_NUMBER);
minusLog_SSCP = computeMinusLog(SSCP(rowIndex, columnIndex), ...
				 VERY_LARGE_NUMBER);
minusLog_TP   = computeMinusLog(TP(rowIndex, columnIndex), ...
				 VERY_LARGE_NUMBER);
minusLog_HDE  = computeMinusLog(HDE(rowIndex, columnIndex), ...
				 VERY_LARGE_NUMBER);

if (refineWithRDCs)
  minusLog_RP1  = computeMinusLog(RP1(rowIndex, columnIndex), ...
				   VERY_LARGE_NUMBER);
  minusLog_RP2  = computeMinusLog(RP2(rowIndex, columnIndex), ...
				   VERY_LARGE_NUMBER);
else
  minusLog_RP1  = [];
  minusLog_RP2  = [];
end