% $$$ MASTER = foundMASTERs{1};
% $$$ mypeakIndices                      = 1:size(MASTER,1);
% $$$ myresidueIndices                   = 1:size(MASTER,2);
% $$$ 
% $$$ [NP,impossibleToAssign]=NVR_NOE2PROB(MASTER,NOES,IALLDISTS,NTH,1: ...
% $$$ 				     size(MASTER,1), 1:size(MASTER, ...
% $$$ 						  2));
% $$$ if (~impossibleToAssign)
% $$$   differenceMatrix = MASTER - NP;
% $$$   [differentPeakIndices,differentResidueIndices] = find(differenceMatrix);
% $$$   for i = 1:length(differentPeakIndices)
% $$$     u = find(peakIndices == differentPeakIndices(i));
% $$$     v = find(residueIndices == differentResidueIndices(i));
% $$$     if (isempty(u))
% $$$       %the pruned peak is not unassigned, i.e. it is assigned already.
% $$$       assert(isempty(v));
% $$$       %the difference is not in one of the unassigned positions. It
% $$$       %is at one of the assigned positions which has been pruned.
% $$$       impossibleToAssign = 1;
% $$$       break;
% $$$     end
% $$$   end
% $$$ end
% $$$ 
% $$$ fprintf(1, 'impossibleToAssign is %d\n',impossibleToAssign);
% $$$ keyboard
% $$$ 
% $$$ [ASSIGNTABLE,voter,numVoters, peakIndices,residueIndices] = removeAssignedPeaksAndResidues(MASTER,ASSIGNTABLE,voter,numVoters,peakIndices,residueIndices);
% $$$ fprintf(1, 'removed assigned peaks and residues from relIndiced matrices\n');
% $$$ 
% $$$ score = 0;
% $$$ 
% $$$ [impossibleToAssign, MASTER, ASSIGNTABLE,  unassignedPeakIndices, unassignedResidueIndices, score, voter, numVoters] = ...
% $$$     noePruneAndDoMoreAssignments(ASSIGNTABLE, MASTER, ...
% $$$ 				 peakIndices, ...
% $$$ 				 residueIndices, ...
% $$$ 				 score, ...
% $$$ 				 voter, numVoters,...
% $$$ 				 NOES, IALLDISTS, ALLDISTS, NTH, ...
% $$$ 				 RDC1,RDC2,VECTORS);
% $$$ 
% $$$ fprintf(1, 'impossibleToAssign as in the merge code is %d\n',impossibleToAssign);
% $$$ keyboard

noeCompliantAssignments = [];
fprintf(1, 'structures that are stringently NOE compliant:\n');
for assignmentIndex = 1:totalNumAssignments
%for assignmentIndex = 108:108
  MASTER = foundMASTERs{assignmentIndex};
  noeCompliant = 1;
  for peakIndex = 1:size(MASTER,1)
    peaksWithNOE = find(NOES (peakIndex,:));
    if (isempty(peaksWithNOE))
      continue;
    end
    for peak2Index = 1:length(peaksWithNOE)
      if (peaksWithNOE(peak2Index) > size(MASTER,1))
	continue;
      end
      residueIndex1 = find(MASTER(peakIndex,:));
      residueIndex2 = find(MASTER(peaksWithNOE(peak2Index),:));
      if (ALLDISTS(residueIndex1,residueIndex2) > NTH)
	noeCompliant = 0;
	break;
      end
    end
    if (~noeCompliant)
      break;
    end
  end
  if (noeCompliant)
    fprintf(1, '%d ', assignmentIndex);
    noeCompliantAssignments = [noeCompliantAssignments assignmentIndex];
  end
end


figure; plot(finalScores(noeCompliantAssignments),'*-')          
hold on, plot(assignmentAccuracies(noeCompliantAssignments),'ro-')    
keyboard
EPSILON = 1E-6;
for assignmentIndex = 1:length(noeCompliantAssignments)
  MASTER = foundMASTERs{noeCompliantAssignments(assignmentIndex)};
  score  = 0;
  computedScore = finalScores(noeCompliantAssignments(assignmentIndex));

  for peakIndex = 1:size(MASTER,1)
    residueIndex = find(MASTER(peakIndex,:));
    relPeakIndex = find(peakIndices == peakIndex);
    if (isempty(relPeakIndex))
      %fprintf(1, 'peak not counted in score computation.\n');
      continue;
    end
    relResidueIndex = find(residueIndices == residueIndex);
    assert (~isempty(relResidueIndex));
    for voterIndex = 1:7
      score = score - log(voter{voterIndex}(relPeakIndex, relResidueIndex));
    end
  end
  assert (abs(score - computedScore) < EPSILON);
  fprintf(1, 'assert passed.\n');
end
