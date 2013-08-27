function compareAssignmentAccuracies

assignmentMatrixFilename1 = 'AST1.TXT';
assignmentMatrixFilename2 = 'AST3.TXT';

myMASTER_WithExtraColumns1 =     load (assignmentMatrixFilename1);
myMASTER_WithExtraColumns2 =     load (assignmentMatrixFilename2);

fprintf(1, 'reading %s...\n', assignmentMatrixFilename1);
fprintf(1, 'reading %s...\n', assignmentMatrixFilename2);
fprintf(1, 'the matrix read has %d rows and %d columns.\n', ...
	size(myMASTER_WithExtraColumns1,1), ...
	size(myMASTER_WithExtraColumns1,2));
fprintf(1, 'the matrix read has %d rows and %d columns.\n', ...
	size(myMASTER_WithExtraColumns2,1), size(myMASTER_WithExtraColumns2,2));

myAss1                  =     myMASTER_WithExtraColumns1(:, size(myMASTER_WithExtraColumns1,2)-1);
myAss2                  =     myMASTER_WithExtraColumns2(:, size(myMASTER_WithExtraColumns2,2)-1);


numCorrectAssignmentsInTheFileRead = 0;

for peakIndex = 1:length(myAss1)
  if (myAss1(peakIndex) == myAss2(peakIndex))
    fprintf(1, '%d %d\n',peakIndex,myAss1(peakIndex));
  end
end
%fprintf(1, 'assignment accuracy in the inputted file = %f\n',numCorrectAssignmentsInTheFileRead/size(myMASTER,1));
