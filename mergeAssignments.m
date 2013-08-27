%how can we obtain multiple assignments from the merge SB
function [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments] =mergeAssignments(foundMASTERsL, ...
						foundMASTERsR, finalScoresL,finalScoresR,totalNumAssignmentsL,totalNumAssignmentsR,...
						  inputPeakIndices, inputResidueIndices, inputVoter, inputNumVoters,inputASSIGNTABLE,NOES,ALLDISTS,IALLDISTS,NTH,peakIDs,RESNUMS,RDC1,RDC2,VECTORS,MB_ShiftScore, MB_RDC_Score, numCS, numRDC, HSQCDATA, csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES)
%,totalNumPeaks)

%SCORE_DELTA_THRESHOLD = 100; %also is present in a_star_assign.m

assert (totalNumAssignmentsL >= 1);
assert (totalNumAssignmentsR >= 1);
assert (~isempty(foundMASTERsL{1}));
assert (~isempty(foundMASTERsR{1}));


[numPeaks,numResidues]      = size(foundMASTERsL{1});

printAssignmentAccuracy     = 0;

%thereIsConflict        = 0;
% $$$ peakIndexL             = 0;
% $$$ peakIndexR = 0;
%minScore                    = finalScoresL(1)+finalScoresR(1);%this
                                                              %is
                                                              %actually
                                                              %not
                                                              %the
                                                              %minScore.
                                                              %minScore
                                                              %varies
                                                              %(is
                                                              %not
                                                              %uniformly
                                                              %increasing
                                                              %due
                                                              %to conflicts)
                                                              %and
                                                              %could
                                                              %even
                                                              %be 0.SB
%score                       = minScore;
leftAssignmentIndex         = 1;
rightAssignmentIndex        = 1;
totalNumAssignments         = 0;
finalScores                 = zeros(0,1);
foundMASTERs                = cell(0,1);
assignmentAccuracies        = zeros(0,1);

assert (totalNumAssignmentsR >= 100);
assert (totalNumAssignmentsL >= 100);

%for (leftAssignmentIndex = 1:totalNumAssignmentsL)
for (leftAssignmentIndex = 1:100)
  fprintf(1, 'leftAssignmentIndex = %d\n', leftAssignmentIndex);
  for (rightAssignmentIndex = 1:100)
%for (leftAssignmentIndex = 44:45)
%  for (rightAssignmentIndex = 55:56)


    %    score = finalScoresL(leftAssignmentIndex) + ...
    %	    finalScoresR(rightAssignmentIndex);
    %   keyboard
    
    mergedMASTER                = zeros(numPeaks, numResidues);
    %unassignedPeakIndices    = [];
    
    
    numCommonAssignedPeaks      = findNumCommonAssignments(foundMASTERsL{leftAssignmentIndex}, foundMASTERsR{rightAssignmentIndex});
    
    mergedMASTER                = mergeMASTERs(foundMASTERsL{leftAssignmentIndex}, foundMASTERsR{rightAssignmentIndex});
    
    %    fprintf(1,' score is %f before removing conflicts.\n',score);
    
% $$$     [mergedMASTER,conflictingAbsPeakIndices,unassignedAbsResidueIndices,numConflicts,score] ...
% $$$ 	= removeConflicts(mergedMASTER,inputVoter,inputNumVoters, inputPeakIndices, ...
% $$$ 			  inputResidueIndices, score);
    
    thereAreConflicts            = areThereConflicts(mergedMASTER);
    
    if (thereAreConflicts)
  %    fprintf(1, 'found conflicts.\n');
      %     fprintf(1, 'found %d conflicts and score is now %f\n',numConflicts, ...
      %	    score);
      continue;
    else
   %   fprintf(1, 'found no conflicts.\n');
    end
    
    %    if (numConflicts > 0)
% $$$    if (thereAreConflicts)
% $$$      continue; %bypassing all conflicting assignments, targeting
% $$$ 	       %tulum.
% $$$    else
% $$$      fprintf(1, 'found an assignment without conflicts.\n');
% $$$    end

    numMergedAssignments = sum(sum(mergedMASTER));
    assert (numMergedAssignments >= numCommonAssignedPeaks);

    if (numMergedAssignments == numCommonAssignedPeaks)
      fprintf(1, 'cancelled all new assignments since they are conflicting.\n');
      fprintf(1, 'unassigned peak indices are:\n');
      for i = 1:length(inputPeakIndices)
	fprintf(1, '%d --> ',inputPeakIndices(i));
	j = find(mergedMASTER(inputPeakIndices(i),:));

	if (isempty(j))
	  fprintf(1, 'unassigned\n');
	else
	  fprintf(1, 'assigned. How can this be?\n');
	  keyboard
	end
      end
      fprintf(1, '\n');
      
%      keyboard;
      continue;
    end
    
    
%    if (score > SCORE_DELTA_THRESHOLD + minScore)
%	break;
%    end
    
    
    [ASSIGNTABLE,voter,numVoters, peakIndices,residueIndices] = removeAssignedPeaksAndResidues(mergedMASTER,inputASSIGNTABLE,inputVoter,inputNumVoters,inputPeakIndices,inputResidueIndices);
 %   fprintf(1, 'removed assigned peaks and residues from relIndiced matrices\n');

    
    if ((~isempty(find(peakIndices == 13))) & (~isempty(find(mergedMASTER(13,:)))))
      fprintf(1, 'peak 13 seems unassigned yet MASTER has en entry for that.\n');
      keyboard
    end
    
    if (isempty(ASSIGNTABLE))
  %    fprintf(1, 'ASSIGNTABLE is empty!\n');
%      keyboard
%impossibleToAssign = 0;
    end
%    else
      
      score = 0;
      [impossibleToAssign, mergedMASTER, ASSIGNTABLE,  unassignedPeakIndices, unassignedResidueIndices, score, voter, numVoters] = ...
	  noePruneAndDoMoreAssignments(ASSIGNTABLE, mergedMASTER, ...
				       peakIndices, ...
				       residueIndices, ...
				       score, ...
				       voter, numVoters,...
				       NOES, IALLDISTS, ALLDISTS, NTH, ...
				       RDC1,RDC2,VECTORS,MB_ShiftScore, ...
				       MB_RDC_Score, numCS, numRDC, ...
				       HSQCDATA, csCoefficient, ...
				       rdcCoefficient, noeCoefficient,numHN_NOES);
      
  %    fprintf(1, 'noe pruned. impossibleToAssign is %d\n',impossibleToAssign);
      if (impossibleToAssign == 0)
%	fprintf(1, 'found an NOE proof assignment\n');
	%keyboard
      end
%      if (impossibleToAssign)
%	fprintf(1, 'impossible to assign!\n');
%	keyboard
%      end
%    end
    
% $$$     if (numConflicts>0)
% $$$       %  ...
% $$$       %  should we not try to resolve the conflicts but return as is?SB
% $$$       %    resolveConflicts();
% $$$     end
    
    if (~impossibleToAssign)
%      fprintf(1, 'storing found assignment.\n');
%      numMergedAssignments =  sum(sum(mergedMASTER));
%      assert (numMergedAssignments ~= 0);
      totalNumAssignments                       = totalNumAssignments + 1;
      finalScores(totalNumAssignments)          = computeMB_Score(mergedMASTER, MB_ShiftScore,MB_RDC_Score, ...
						  numCS,numRDC, HSQCDATA,ALLDISTS,NTH, ...
						  csCoefficient, rdcCoefficient, noeCoefficient);
      assignmentAccuracies(totalNumAssignments) = computeAssignmentAccuracy(peakIDs, RESNUMS, mergedMASTER,printAssignmentAccuracy);
      foundMASTERs{totalNumAssignments}         = mergedMASTER;
      
    else
 %     fprintf(1, 'there are noe clashes in this one.\n');
%      keyboard
    end
    
 %   keyboard
    
  end
end




[sortedFinalScores,sortedIndices]   = sort(-finalScores(1:totalNumAssignments));
reorderedAssignmentAccuracies       = zeros(totalNumAssignments,1);
reorderedMASTERs                    = cell(totalNumAssignments,1);
for i = 1:totalNumAssignments
  reorderedAssignmentAccuracies(i)  = assignmentAccuracies(sortedIndices(i));
  reorderedMASTERs{i}               = foundMASTERs{sortedIndices(i)};
end

foundMASTERs                        = reorderedMASTERs;
assignmentAccuracies                = reorderedAssignmentAccuracies;
finalScores                         = -sortedFinalScores;

function numCommonAssignedPeaks      = findNumCommonAssignments(MASTER1, MASTER2);

numCommonAssignedPeaks = 0;
assert (size(MASTER1,1) == size(MASTER2,1));
for peakIndex = 1:size(MASTER1,1)
  residueIndex1  = find(MASTER1(peakIndex,:));
  residueIndex2  = find(MASTER2(peakIndex,:));
  if (isempty(residueIndex1) | isempty(residueIndex2))
    continue;
  end
  if (residueIndex1==residueIndex2)
    numCommonAssignedPeaks = numCommonAssignedPeaks + 1;
  end
end


function  mergedMASTER                = mergeMASTERs(MASTER1, MASTER2);

[numPeaks,numResidues] = size(MASTER1);
mergedMASTER           = zeros(numPeaks,numResidues);


for peakIndex = 1:numPeaks
    
% $$$   relPeakIndexL = find(peakIndices(1:numPeaks/2) == peakIndex);
% $$$   relPeakIndexR = find(peakIndices(numPeaks/2+1:numPeaks) == peakIndex);
% $$$   
% $$$   if (isempty(relPeakIndexL)) & (isempty(relPeakIndexR))
% $$$     peakIndexL = peakIndexL + 1;
% $$$     peakIndexR = peakIndexR + 1;
% $$$   elseif (isempty(peakIndexR))
% $$$       peakIndexL = peakIndexL + 1;
% $$$   elseif (isempty(peakIndexL))
% $$$     peakIndexR = peakIndexR + 1;
% $$$   end
% $$$   SB
    
  residueIndexL = find(MASTER1(peakIndex,:));
  residueIndexR = find(MASTER2(peakIndex,:));
  
  assert (length(residueIndexL)<=1);
  assert (length(residueIndexR)<=1);    
  
  if ((~isempty(residueIndexL)) & (~isempty(residueIndexR)))
    if (residueIndexL == residueIndexR)
      mergedMASTER(peakIndex,residueIndexL) = 1;
    else
      fprintf(1, 'this cant happen since the same peak is not assigned by two different calls to divideAndConquer\n');
      
      keyboard
      %     thereIsConflict = 1;
    end
  elseif (isempty(residueIndexR)) & (~isempty(residueIndexL)) 
    mergedMASTER(peakIndex, residueIndexL) = 1;
  elseif (isempty(residueIndexL)) & (~isempty(residueIndexR)) 
    mergedMASTER (peakIndex, residueIndexR) = 1;
  end
  
end

function   thereAreConflicts  = areThereConflicts(mergedMASTER)%, voter, numVoters, peakIndices, ...
				       %			residueIndices, score);
  
%conflictingAbsPeakIndices   = [];
%unassignedAbsResidueIndices = [];
thereAreConflicts            = 0;  
numResidues                  = size(mergedMASTER,2);
% $$$ deltaScore                  = 0;

for residueIndex = 1:numResidues
  
  assignedPeakIndices = find(mergedMASTER(:,residueIndex));
  
  if (length(assignedPeakIndices) > 1)

    assert (length(assignedPeakIndices) == 2);
    assert (assignedPeakIndices(1) < assignedPeakIndices(2));
    thereAreConflicts                      = 1;
    break;
  end
end

% $$$     numConflicts                           = numConflicts + 1;
% $$$     mergedMASTER(assignedPeakIndices,:)    = 0;
% $$$     conflictingAbsPeakIndices              = [conflictingAbsPeakIndices; assignedPeakIndices];
% $$$     
% $$$     relPeakIndex1 = find(peakIndices == assignedPeakIndices(1));
% $$$     relPeakIndex2 = find(peakIndices == assignedPeakIndices(2));
% $$$ 
% $$$     assert (~isempty(relPeakIndex1));
% $$$     assert (~isempty(relPeakIndex2));
% $$$     relResidueIndex = find(residueIndices == residueIndex);
% $$$     assert (~isempty(relResidueIndex));

    %keyboard
    
    
% $$$     for voterIndex=1:numVoters
% $$$       deltaScore = deltaScore + ...
% $$$ 	  log(voter{voterIndex}(relPeakIndex1,relResidueIndex));
% $$$       
% $$$       deltaScore = deltaScore + ...
% $$$ 	  log(voter{voterIndex}(relPeakIndex2,relResidueIndex));
% $$$       
% $$$     end
% $$$   end
% $$$   
% $$$   if ((length(assignedPeakIndices) > 1) | (isempty(assignedPeakIndices)))
% $$$ %    fprintf('residue %d is unassigned\n',residueIndex);
% $$$     %    keyboard
% $$$     unassignedAbsResidueIndices  = [unassignedAbsResidueIndices residueIndex];
% $$$   end
% $$$   
% $$$ 
% $$$ end
% $$$ score = score + deltaScore;