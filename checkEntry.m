function checkEntry(nonZeroValue, rowIndex, matrixName);

if (nonZeroValue == 0) 
    fprintf(1, '%s matrix entry empty for peak %d\n', ...
	    matrixName, rowIndex);
end