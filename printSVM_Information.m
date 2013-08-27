function printSVM_Information(MASTER,ASSIGNTABLE, ROWIN, COLIN, differenceMatrixH_SHIFTX, ...
			      differenceMatrixN_SHIFTX, ...
			      differenceMatrixH_SHIFTS, ...
			      differenceMatrixN_SHIFTS, ...
			      differenceMatrix_RDC1,differenceMatrix_RDC2, alphaHelixMatrix,betaStrandMatrix,coilMatrix, labelFilename, vectorsFilename)

[numPeaks,numResidues]         = size(MASTER);
  

numVoters                          = 9;
correctAssignmentProbability       = cell(numVoters, 1);
incorrectAssignmentProbability     = cell(numVoters, 1);
correctAssignmentsPeakIndices      = [];
incorrectAssignmentsPeakIndices    = [];
incorrectAssignmentsResidueIndices = [];

for voterIndex = 1:numVoters
  correctAssignmentProbability{voterIndex}   = [];
  incorrectAssignmentProbability{voterIndex} = [];
end

%  voter                            = initialize5Voters(SSCP,SXCP,CP, RP1, RP2);
%voter                   = initialize6Voters(differenceMatrixH_SHIFTX, differenceMatrixN_SHIFTX, differenceMatrixH_SHIFTS, differenceMatrixN_SHIFTS,differenceMatrix_RDC1,differenceMatrix_RDC2);
voter                   = initialize9Voters(differenceMatrixH_SHIFTX, differenceMatrixN_SHIFTX, differenceMatrixH_SHIFTS, differenceMatrixN_SHIFTS,differenceMatrix_RDC1,differenceMatrix_RDC2,alphaHelixMatrix,betaStrandMatrix,coilMatrix);
%  numVoters                        = 5;

for peakIndex=1:numPeaks

  if (~isempty(find(MASTER(peakIndex,:))))
    continue;
  end
  
  relPeakIndex = find(ROWIN == peakIndex);
  
  if (isempty(relPeakIndex))
    continue;
  end
  
  
  correctResidueIndex = peakIndex;
  
  correctRelResidueIndex = find(COLIN == correctResidueIndex);

  if (~isempty(correctRelResidueIndex)) & (ASSIGNTABLE(relPeakIndex,correctRelResidueIndex) == 1)
    correctAssignmentsPeakIndices = [correctAssignmentsPeakIndices peakIndex];
  end
  
  for voterIndex = 1:numVoters

    if (~isempty(correctRelResidueIndex)) & (ASSIGNTABLE(relPeakIndex,correctRelResidueIndex) == 1)
      correctAssignmentProbability{voterIndex} =...
	  [correctAssignmentProbability{voterIndex} ...
	   voter{voterIndex}(peakIndex,correctResidueIndex)];
    end
    
    for residueIndex = 1:numResidues
      if (residueIndex == correctResidueIndex)
	continue;
      end
      
      relResidueIndex = find(COLIN == residueIndex);
      
      if (isempty(relResidueIndex))
	continue;
      end
      
      if  (ASSIGNTABLE(relPeakIndex,relResidueIndex) == 1)
	incorrectAssignmentProbability{voterIndex} = ...
	    [incorrectAssignmentProbability{voterIndex} voter{voterIndex}(peakIndex,residueIndex)];
	if (voterIndex == 1)
	  incorrectAssignmentsPeakIndices = [incorrectAssignmentsPeakIndices ...
		    peakIndex];
	  incorrectAssignmentsResidueIndices = ...
	      [incorrectAssignmentsResidueIndices residueIndex];
	end
      end
    end
  end
end

f_labels  = fopen(labelFilename,'w');
f_vectors = fopen(vectorsFilename,'w');
fprintf(f_labels, 'item\tlabel\n');
fprintf(f_vectors,'item\tv1\tv2\tv3\tv4\tv5\tv6\tv7\tv8\tv9\tv10\tv11\tv12\n');

numCorrectVotes = length(correctAssignmentProbability{1});
assert (numCorrectVotes == length(correctAssignmentsPeakIndices));

for correctVoteIndex = 1:numCorrectVotes
  fprintf(f_vectors, '%d',correctVoteIndex);
  fprintf(f_labels, '%d' ,correctVoteIndex);
  fprintf(f_vectors, '\t%d\t%d',correctAssignmentsPeakIndices(correctVoteIndex), correctAssignmentsPeakIndices(correctVoteIndex));
  for voterIndex=1:numVoters
    fprintf(f_vectors, '\t%f',correctAssignmentProbability{voterIndex}(correctVoteIndex));
  end
  if ((correctAssignmentProbability{5}(correctVoteIndex) < 0) | ...
      (correctAssignmentProbability{6}(correctVoteIndex) < 0))
%I think this check is to make sure that the data is valid. Somehow
%if the RDCs are missing I think the 5th or 6th entries are negative.
    fprintf(f_vectors, '\t0');
  else
    fprintf(f_vectors, '\t1');
  end
  
  fprintf(f_vectors, '\n');
  fprintf(f_labels, '\t1\n');
end

numIncorrectVotes = length(incorrectAssignmentProbability{1});
for incorrectVoteIndex = 1:numIncorrectVotes
  fprintf(f_vectors, '%d',incorrectVoteIndex+numCorrectVotes);
  fprintf(f_labels, '%d',incorrectVoteIndex+numCorrectVotes);
  
  fprintf(f_vectors, '\t%d\t%d',incorrectAssignmentsPeakIndices(incorrectVoteIndex), incorrectAssignmentsResidueIndices(incorrectVoteIndex));
  
  for voterIndex=1:numVoters
    fprintf(f_vectors, '\t%f',incorrectAssignmentProbability{voterIndex}(incorrectVoteIndex));
  end
  
  if ((incorrectAssignmentProbability{5}(incorrectVoteIndex) < 0) | ...
      (incorrectAssignmentProbability{6}(incorrectVoteIndex) < 0))
    fprintf(f_vectors, '\t0');
  else
    fprintf(f_vectors, '\t1');
  end
  
  fprintf(f_vectors, '\n');
  fprintf(f_labels,  '\t-1\n');
end

fclose(f_vectors);
fclose(f_labels);