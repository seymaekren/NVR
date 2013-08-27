function retval = twoMatricesAreClose(matrix1, matrix2, EPSILON)

[sizeX1, sizeY1] = size(matrix1);
[sizeX2, sizeY2] = size(matrix2);

if (sizeX1 ~= sizeX2) | (sizeY1 ~= sizeY2)
  retval = 0;
  return ;
end

comparisonMatrix = find(abs(matrix1 - matrix2) > EPSILON);
if (isempty(comparisonMatrix))
  retval = 1;
  return;
else
  retval = 0;
  return;
end