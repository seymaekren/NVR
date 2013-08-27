%function stochasticSearch

MAX_NUM_ITERATIONS = 100;

unassignedPeakIndices    = peakIndices;
unassignedResidueIndices = residueIndices;
%first making a descent


[foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments]=divideAndConquer(MASTER, voter, numVoters, ASSIGNTABLE, unassignedPeakIndices, ...
						  unassignedResidueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
						  peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
						  MB_RDC_Score, numCS, numRDC, HSQCDATA,csCoefficient, ...
						  rdcCoefficient, noeCoefficient,numHN_NOES);

%keyboard

assert (totalNumAssignments >= 1);
overallMASTER= foundMASTERs{1}*0;
for resultIndex = 1:totalNumAssignments
  overallMASTER = overallMASTER + foundMASTERs{resultIndex};
end

numIterations = 0;
while (1)

  [prevScore, aggregatedMatrix, resultMASTER1,resultMASTER2, ASSIGNTABLE1, unassignedPeakIndices1, unassignedResidueIndices1, ASSIGNTABLE2, unassignedPeakIndices2, unassignedResidueIndices2] = ...
      aggregateComputeScoreAndPickTwoRandomAssignments(overallMASTER, ...
						  totalNumAssignments, ...
						  ASSIGNTABLE, unassignedPeakIndices, ...
						  unassignedResidueIndices,...
						  MB_ShiftScore, MB_RDC_Score, ...
						  numCS, numRDC, HSQCDATA, ...
						  ALLDISTS, NTH, csCoefficient, ...
						  rdcCoefficient, noeCoefficient);

  numIterations = numIterations + 1;
  
  if (resultMASTER1 == MASTER) | (resultMASTER2 == MASTER) | ...
	(resultMASTER1 == resultMASTER2)
  else
    break;
  end
  
  if (numIterations > MAX_NUM_ITERATIONS)
      fprintf(1, 'did too many iterations but couldnt find ');
      fprintf(1, 'two random MASTERs out of the set of assignments.\n');
      keyboard
  end
  
end
  
keyboard
while (1)

  [foundMASTERs1, finalScores1, assignmentAccuracies1, totalNumAssignments1]=divideAndConquer(resultMASTER1, voter, ...
						  numVoters, ASSIGNTABLE1, unassignedPeakIndices1, ...
						  unassignedResidueIndices1, NOES, ALLDISTS, IALLDISTS, NTH, ...
						  peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
						  MB_RDC_Score, numCS, numRDC, HSQCDATA,csCoefficient, ...
						  rdcCoefficient, ...
						  noeCoefficient,numHN_NOES);

  [foundMASTERs2, finalScores2, assignmentAccuracies2, totalNumAssignments2]=divideAndConquer(resultMASTER2, voter, ...
						  numVoters, ASSIGNTABLE2, unassignedPeakIndices2, ...
						  unassignedResidueIndices2, NOES, ALLDISTS, IALLDISTS, NTH, ...
						  peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
						  MB_RDC_Score, numCS, numRDC, HSQCDATA,csCoefficient, ...
						  rdcCoefficient, ...
						  noeCoefficient,numHN_NOES);
  %  gather new overall matrix from both set of matrices;
  assert (totalNumAssignments2 >= 0);
  assert (totalNumAssignments1 >= 0);

  overallMASTER = foundMASTERs1{1}*0;
  for resultIndex = 1:totalNumAssignments1
    overallMASTER = overallMASTER + foundMASTERs1{resultIndex};
  end
  
  for resultIndex = 1:totalNumAssignments2
    overallMASTER = overallMASTER + foundMASTERs2{resultIndex};
  end

  numIterations = 0; 
  
  while (1)
  
    [newScore, aggregatedMatrix, resultMatrix1,resultMatrix2, ASSIGNTABLE1, unassignedPeakIndices1, unassignedResidueIndices1, ASSIGNTABLE2, unassignedPeakIndices2, unassignedResidueIndices2] = aggregateComputeScoreAndPickTwoRandomAssignments(overallMASTER,totalNumAssignments1+totalNumAssignments2, ...
						  ASSIGNTABLE, ...
						  unassignedPeakIndices, ...
						  unassignedResidueIndices, ...
						  MB_ShiftScore,MB_RDC_Score, ...
						  numCS, numRDC, HSQCDATA, ...
						  ALLDISTS, NTH, csCoefficient, ...
						  rdcCoefficient, ...
						  noeCoefficient);
  
    numIterations = numIterations + 1;
    
    if ((resultMatrix1 == resultMatrix2) | (resultMatrix1 == resultMASTER1) ...
	| (resultMatrix1 == resultMASTER2) | (resultMatrix2 == ...
					      resultMASTER1) | ...
	(resultMatrix2 == resultMASTER2))
    else
      resultMASTER1 = resultMatrix1; resultMASTER2 =resultMatrix2;
      break;
    end
  
    if (numIterations > MAX_NUM_ITERATIONS)
      fprintf(1, 'did too many iterations but couldnt find ');
      fprintf(1, 'two random MASTERs out of the set of assignments.\n');
      keyboard
    end
  end
  
  if (newScore <= prevScore)
    break;
  end

  prevScore = newScore;
end

computeAssignmentAccuracy(aggregatedMatrix, peakIDs, RESNUMS);