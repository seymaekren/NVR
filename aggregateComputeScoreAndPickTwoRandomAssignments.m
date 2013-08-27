function [prevScore, aggregatedMatrix,resultMaster1,resultMaster2, ASSIGNTABLE1, unassignedPeakIndices1, unassignedResidueIndices1, ASSIGNTABLE2, unassignedPeakIndices2, unassignedResidueIndices2] = aggregateComputeScoreAndPickTwoRandomAssignments(overallMASTER,totalNumAssignments, ASSIGNTABLE, unassignedPeakIndices, unassignedResidueIndices, ...
						  MB_ShiftScore,MB_RDC_Score, ...
						  numCS, numRDC, HSQCDATA, ...
						  ALLDISTS, NTH, csCoefficient, ...
						  rdcCoefficient, ...
						  noeCoefficient);

LARGE_VALUE  = 1E15;
EPSILON      = 1E-5;

scoreMatrix                                 = -overallMASTER;%.*MB_ShiftScore./2;

[numPeaks,numResidues]                      = size(scoreMatrix);

squareScoreMatrix                           = ones(max(numPeaks,numResidues))*LARGE_VALUE;

squareScoreMatrix(1:numPeaks,1:numResidues) = scoreMatrix; 

[assignments,minusTotalNumberOfVotes]       = hungarian(squareScoreMatrix');

debugCost                                   = 0;

for peakIndex = 1:size(overallMASTER,1)
  debugCost  = debugCost + scoreMatrix(peakIndex,assignments(peakIndex));;
end

assert (abs(debugCost + LARGE_VALUE * (max(numPeaks,numResidues)- numPeaks) - minusTotalNumberOfVotes) < EPSILON);

aggregatedMatrix  = overallMASTER * 0;

for peakIndex = 1:size(overallMASTER,1)
  aggregatedMatrix(peakIndex,assignments(peakIndex)) = 1;
end


prevScore         = computeMB_Score(aggregatedMatrix, MB_ShiftScore,MB_RDC_Score, ...
				    numCS, numRDC, HSQCDATA, ...
				    ALLDISTS, NTH, csCoefficient, ...
				    rdcCoefficient, noeCoefficient);
assignments       = assignments    (1:size(overallMASTER,1));
resultMaster1     = pickRandomAssignments(assignments, overallMASTER, ...
					  totalNumAssignments);

voter = cell(1,1); numVoters = 0;

[ASSIGNTABLE1, voter, numVoters, unassignedPeakIndices1, unassignedResidueIndices1] = ...
        removeAssignedPeaksAndResidues(resultMaster1, ASSIGNTABLE, voter, numVoters, unassignedPeakIndices, ...
				   unassignedResidueIndices);

resultMaster2     = pickRandomAssignments(assignments, overallMASTER, ...
					  totalNumAssignments);

[ASSIGNTABLE2, voter, numVoters, unassignedPeakIndices2, unassignedResidueIndices2]= removeAssignedPeaksAndResidues(resultMaster2, ASSIGNTABLE, voter, numVoters, unassignedPeakIndices, unassignedResidueIndices);

