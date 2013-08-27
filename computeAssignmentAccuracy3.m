clear all;
order = load ('order.m.CAM');
assignments = load ('AST_WITH_TOCSY.txt');
%assignments = load ('AST_SANS_TOCSY.txt');
[numRows, numColumns] = size(assignments)
numCorrect = 0; numIncorrect = 0;
for i = 1:numRows
  assignmentPeak = assignments(i, numColumns-2);
  assignmentResidue = assignments(i, numColumns-1);
  if (order(assignmentPeak,2) == -1)
    continue;
  end
  if (assignmentResidue == order(assignmentPeak,2))
    numCorrect = numCorrect + 1;
  else
    numIncorrect = numIncorrect + 1;
  end
end
fprintf(1, 'numCorrect = %d numIncorrect = %d\n',numCorrect,numIncorrect);