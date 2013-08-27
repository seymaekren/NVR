function resultMaster     = pickRandomAssignments(assignments, overallMASTER, ...
						  totalNumAssignments);


assert (length(assignments) == size(overallMASTER,1));
resultMaster   = overallMASTER * 0;

for peakIndex  = 1:length(assignments)
  residueIndex = assignments(peakIndex);
  numVotes     = overallMASTER(peakIndex,residueIndex);
  ratioOfVotes = numVotes/totalNumAssignments;
  assert ((ratioOfVotes >= 0) & (ratioOfVotes <= 1));
  uniformNumber= rand(1);
  if (uniformNumber < ratioOfVotes)
    resultMaster(peakIndex,residueIndex) = 1;
  end
end