function computeAssignmentAccuracy2(MASTER)
count = 0;
[numRows, numColumns] = size(MASTER);
for i = 1:numRows
  if (MASTER(i,i) == 1)
    count = count + 1;
  end
end
fprintf(1, 'assignment accuracy = %f\n', count/numRows);