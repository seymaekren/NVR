load combinedScoringMatrix.txt
[numPeaks,numResidues] = size(combinedScoringMatrix);
%assignments = zeros(numPeaks,2);
%for i = 1:numPeaks
%  assignments(i,:) = i;
%end
assignments = load ('1G6J.txt');
assignments = assignments(:,size(assignments,2)-2:size(assignments, ...
						  2)-1);
score = 0;
for i = 1:size(assignments,1)
  score = score + combinedScoringMatrix(assignments(i,1),assignments(i,2));
  if (score > 1E8)
    fprintf(1, 'the assignment of %d to %d has a score of %f\n', assignments(i,1),assignments(i,2),score);
  end
end
fprintf(1, 'score = %f\n', score);
fprintf(1, 'score of correct assignment = %f\n', sum(diag(combinedScoringMatrix)));
fprintf('checkAssignments run completed.\n');