assignments = load ('assignment_PR.out.parsed');
%assignments = load ('confidentlyAssignedResidues_MH.txt');

correct = 0; incorrect = 0;
for i = 1:size(assignments,1)
  residueId = assignments(i,1);
  peakId    = assignments(i,2);
%  if (peakId+709 == residueId)
  if (peakId == residueId)
    correct = correct + 1;
  else
    incorrect = incorrect + 1;
  end
  

end

  fprintf(1, 'a.a. = %f\n', correct/(correct + incorrect));