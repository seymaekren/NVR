function analyzeOverallAssignments(overallMASTER, peakIDs, RESNUMS)

a = load('answerkey.m');
assignments=0;
correct=0;
incorr=0;
numResidues = size(overallMASTER,2);
numResiduesCheck = length(RESNUMS);
assert (numResidues  == numResiduesCheck);

for(i=1:size(overallMASTER,1))

  pk  = peakIDs(i);  %get the peak id
  foo = find(a(:,1)==pk);
  correctResidueIndex = a(foo,2);
  relResidueIndex = find(RESNUMS == correctResidueIndex);
  numAssignmentsThatAssignTheCorrectResidue = overallMASTER(i, ...
						  relResidueIndex);
  
  [highestNumberOfAssignments, bestAssignmentResidueIndex] = max(overallMASTER(i,:));  
  
  assert (sum(overallMASTER(i,:))~= 0);
  
  ratioOfAssignmentsThatAssignTheCorrectResidue = numAssignmentsThatAssignTheCorrectResidue/sum(overallMASTER(i,:));
  
  fprintf(1, 'peak %d, correct residue %d, %% correct in the set of assignments %f\n',pk,correctResidueIndex,100*ratioOfAssignmentsThatAssignTheCorrectResidue);
  if (relResidueIndex ~= bestAssignmentResidueIndex)
    fprintf(1, 'the highest number of assignments agreed on residue#%d\n',RESNUMS(bestAssignmentResidueIndex));
  end
end
