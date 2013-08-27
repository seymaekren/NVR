load combinedScoringMatrix.txt
MASTER = load ('hSRI.txt');
count  = 0;
[numRows, numColumns] = size(MASTER);
MASTER = MASTER(1:numRows, 1:numColumns - 3);
fprintf(1, 'reduced numColumns of MASTER to %d\n', size(MASTER,2));
score = 0;
for i = 1:numRows
  if (MASTER(i,i) == 1)
    correct = 1;
    count = count + 1;
  else
    correct = 0;
  end
  residue = find(MASTER(i,:));
  score = combinedScoringMatrix(i,residue);
  if (correct == 1)
    plot(i,score, 'b+');
  else
    plot(i,score, 'r-');
  end
  hold on
end
fprintf(1, 'assignment accuracy = %f\n', count/numRows);



