function [overallMASTER, totalNumAssignments, foundMASTERs, finalScores, assignmentAccuracies] = reportAndStoreAssignments(MASTER, score, overallMASTER, peakIDs, RESNUMS, foundMASTERs, ...
						  finalScores, ...
						  totalNumAssignments, assignmentAccuracies);

foundThisAssignmentBefore = 0;
for i=1:totalNumAssignments
  if (foundMASTERs{i} == MASTER)
    foundThisAssignmentBefore = 1;
    fprintf(1, 'found this assignment as assignment #%d\n', i);
    break;
  end
end

if (foundThisAssignmentBefore)
  return;
end

[assignmentAccuracy, assignments] = computeAssignmentAccuracy(peakIDs, RESNUMS, MASTER);

overallMASTER                     = overallMASTER + MASTER;

fprintf(1, 'score is %f\n', score);

totalNumAssignments = totalNumAssignments + 1;

foundMASTERs{totalNumAssignments} = MASTER;
finalScores(totalNumAssignments) = score;
assignmentAccuracies(totalNumAssignments) = assignmentAccuracy;
%keyboard
save foundMASTERsWithFullRDC_Score.mat foundMASTERs totalNumAssignments assignmentAccuracies finalScores;