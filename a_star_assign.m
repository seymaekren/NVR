function [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments] ...
    = a_star_assign(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
		    residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
		    peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
		    MB_RDC_Score, numCS, numRDC, HSQCDATA,csCoefficient, ...
		    rdcCoefficient, noeCoefficient,numHN_NOES);
  

%would be good to make the queue a single structure, with all
%subcomponents. SB.

%could also be good to port this into C++. SB


MAXNUMNODES           = 1000;
MAXNUMSOLUTIONS       = 100;
MAX_FRAGMENT_SIZE     = 10; %number of amino acids such that if so
                           %many assignments are made, the
                           %procedure will stop.
numNodes              = 1;
nodeIndex             = 1;
nodeIndicesArray      = zeros(MAXNUMNODES, 1); 
next                  = zeros(MAXNUMNODES, 1);
foundAnAssignment     = 0;
minScore              = 0; %of the first found assignment
SCORE_DELTA_THRESHOLD = 1; %change in the assignment score that
                             %is allowed to accept a given
                             %assignment. ALSO DEFINED IN voteAccordingToMB_Score.m
totalNumAssignments   = 0;
foundMASTERs          = cell(MAXNUMSOLUTIONS,1);
%overallMASTER         = zeros(size(MASTER,1), size(MASTER,2));
finalScores           = zeros(MAXNUMSOLUTIONS, 1);
assignmentAccuracies  = zeros(MAXNUMSOLUTIONS, 1);

numPeaks              = size(MASTER,1);
overallMASTER         = MASTER * 0;

initiallyNumAssignedPeaks = sum(sum(MASTER));
[ASSIGNTABLES, MASTERS, voters, numVoterss, peakIndicess, residueIndicess, scoreSoFars, h_values] = make_queue(ASSIGNTABLE, MASTER, voter, numVoters, peakIndices, ...
						  residueIndices, MAXNUMNODES, MB_ShiftScore,MB_RDC_Score, ...
						  numCS,numRDC, HSQCDATA,ALLDISTS,NTH,csCoefficient, ...
						  rdcCoefficient, noeCoefficient,numHN_NOES);

while (1)
  
  %  if (isempty(nodes)) 
  if (nodeIndex == 0)
    break;
  end
  
  %node = remove_front(nodes);

  MASTER         = MASTERS{nodeIndex};
  ASSIGNTABLE    = ASSIGNTABLES{nodeIndex};
  voter          = voters{nodeIndex};
  numVoters      = numVoterss(nodeIndex);
  peakIndices    = peakIndicess{nodeIndex};
  residueIndices = residueIndicess{nodeIndex};
  scoreSoFar     = scoreSoFars(nodeIndex);
  h_value        = h_values(nodeIndex);
  score          = scoreSoFar + h_value;

  %could print the fetched node here. how many assignments it has,
  %what is its score, etc.  
  
  if (foundAnAssignment) & (score > minScore + SCORE_DELTA_THRESHOLD) %  if (stoppingCondition)
    fprintf(1, 'score = %f is too much higher than minScore = %f.\n',score,minScore);
    fprintf(1, 'stopping A* search.\n');
    break;
  end
  
  if (totalNumAssignments > MAXNUMSOLUTIONS)
    fprintf(1, 'found too many solutions. stopping A* search.\n');
    break;
  end
  
  numAssignedPeaks = sum(sum(MASTER));
  fprintf(1, 'current node is #%d has %d assigned peaks\n', nodeIndex,numAssignedPeaks);
  fprintf(1, 'its score is g = %f h = %f total = %f\n', scoreSoFar, ...
	  h_value, score);
  
  if ((numAssignedPeaks == numPeaks) | (isempty(ASSIGNTABLE)) | ...
				       (numAssignedPeaks - ...
					initiallyNumAssignedPeaks >= MAX_FRAGMENT_SIZE))

%	if (numAssignedPeaks - initiallyNumAssignedPeaks > ...
%	    MAX_FRAGMENT_SIZE)
%          fprintf(1, 'current solution has %d assignment\n',numAssignedPeaks);
%	  fprintf(1, 'the initial # assignments was %d\n',initiallyNumAssignedPeaks);
%	  keyboard
%	end
	  %  if  (goal_test(node))

%  if (numAssignedPeaks == MAX_FRAGMENT_SIZE) | (isempty(ASSIGNTABLE))%  if (goal_test(node))

    fprintf(1, 'found an assignment.\n');
    if (foundAnAssignment == 0)
      minScore = score;
    end
    foundAnAssignment = 1;
    
    [overallMASTER, totalNumAssignments, foundMASTERs, finalScores, assignmentAccuracies]       = reportAndStoreAssignments(MASTER, score, ...
						  overallMASTER, peakIDs, ...
						  RESNUMS, foundMASTERs, finalScores, ...
						  totalNumAssignments, ...
						  assignmentAccuracies);
    
%    keyboard
    
  else
    [ASSIGNTABLES,MASTERS, voters, numVoterss, peakIndicess, residueIndicess, scoreSoFars, h_values, next, numNodes]  = ...
	expand(ASSIGNTABLES, MASTERS, voters, numVoterss, ...
	       peakIndicess, residueIndicess, scoreSoFars, h_values, next, numNodes, nodeIndex, ...
	       ASSIGNTABLE, MASTER, voter, numVoters, peakIndices, ...
	       residueIndices, scoreSoFar, NOES, IALLDISTS, ALLDISTS, NTH, RDC1, RDC2, VECTORS,MB_ShiftScore, ...
	       MB_RDC_Score, numCS, numRDC, ...
	       HSQCDATA,csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES);
    
  end
  
  nodeIndex = next(nodeIndex);
  
end

% $$$ aggregateMASTER         = findAggregateMASTER(overallMASTER, totalNumAssignments);
% $$$ totalNumAssignments     = 1;
% $$$ foundMASTERs{1}         = aggregateMASTER;
% $$$ finalScores(1)          = computeDeltaScore(MASTERS{1}, ...
% $$$ 					aggregateMASTER, voters{1}, numVoterss(1),peakIndicess{1},residueIndicess{1});
% $$$ assignmentAccuracies(1) = computeAssignmentAccuracy(peakIDs, RESNUMS, aggregateMASTER);
% $$$ numAssignedPeaks        = sum(sum(aggregateMASTER));
% $$$ numInitialPeaks         = sum(sum(MASTERS{1}));
% $$$ 
% $$$ assert (numAssignedPeaks >= numInitialPeaks);
% $$$ if (numAssignedPeaks > numInitialPeaks)
% $$$   fprintf(1, 'assigned %d extra peaks\n',numAssignedPeaks - ...
% $$$ 	  numInitialPeaks); 
% $$$ else
% $$$   fprintf(1, 'num assignments remained the same.\n');
% $$$ end

fprintf(1, 'finished a_star_assign\n');
% $$$ if (totalNumAssignments > 0)
% $$$   figure; plot(assignmentAccuracies(1:totalNumAssignments), '*-'), ...
% $$$       hold on, plot(finalScores(1:totalNumAssignments),'ro-')
% $$$ else
% $$$   fprintf(1, 'could find no assignment.\n');
% $$$ end
%keyboard

%
function [ASSIGNTABLES, MASTERS, voters, numVoterss, peakIndicess, ...
	  residueIndicess, scoreSoFars, h_values] = make_queue(ASSIGNTABLE, ...
						  MASTER, voter, numVoters, peakIndices, residueIndices, ...
						  MAXNUMNODES, MB_ShiftScore,MB_RDC_Score, ...
						  numCS,numRDC, ...
						  HSQCDATA,ALLDISTS,NTH, ...
						  csCoefficient, ...
						  rdcCoefficient, ...
						  noeCoefficient, numHN_NOES)

ASSIGNTABLES       = cell(MAXNUMNODES,1);
MASTERS            = cell(MAXNUMNODES,1);
voters             = cell(MAXNUMNODES,1);
numVoterss         = zeros(MAXNUMNODES,1);
peakIndicess       = cell(MAXNUMNODES,1);
residueIndicess    = cell(MAXNUMNODES,1);
scoreSoFars        = zeros(MAXNUMNODES, 1);
h_values           = zeros(MAXNUMNODES, 1);



ASSIGNTABLES{1}    = ASSIGNTABLE;
MASTERS{1}         = MASTER;
voters{1}          = voter;
numVoterss(1)      = numVoters;
peakIndicess{1}    = peakIndices;
residueIndicess{1} = residueIndices;
h_value            = computeBetterH_Value(ASSIGNTABLE, voter, numVoters);
% $$$ h_value            = computeMB_H_Value(ASSIGNTABLE, MASTER, HSQCDATA,...
% $$$ 				       peakIndices, residueIndices,...
% $$$ 				       ALLDISTS, NTH, ...
% $$$ 				       MB_ShiftScore, ...
% $$$ 				       MB_RDC_Score, numCS, numRDC, ...
% $$$ 				       csCoefficient, rdcCoefficient, ...
% $$$ 				       noeCoefficient, numHN_NOES);

% $$$ score              = computeMB_Score(MASTER, MB_ShiftScore,MB_RDC_Score, ...
% $$$ 			numCS,numRDC, HSQCDATA,ALLDISTS,NTH,csCoefficient, ...
% $$$ 				      rdcCoefficient, noeCoefficient);

score               = 0;

%normalizingConstant = (csCoefficient*numCS +  rdcCoefficient*numRDC + noeCoefficient*numHN_NOES);
scoreSoFars(1)      = -score;%/normalizingConstant;
h_values(1)         = h_value;%/normalizingConstant;





function [ASSIGNTABLES,MASTERS, voters, numVoterss, peakIndicess, residueIndicess, scoreSoFars, h_values, next, numNodes]  =  expand (ASSIGNTABLES, MASTERS, voters, numVoterss, ...
						  peakIndicess, residueIndicess, scoreSoFars, h_values, ...
						  next, numNodes, ...
						  nodeIndex, ...
						  ASSIGNTABLE, MASTER, ...
						  voter, numVoters, ...
						  peakIndices, ...
						  residueIndices, ...
						  scoreSoFar, NOES, ...
						  IALLDISTS, ALLDISTS, NTH, RDC1, RDC2, VECTORS,MB_ShiftScore, ...
						  MB_RDC_Score, ...
						  numCS, numRDC, ...
						  HSQCDATA,...
                                                  csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES)




EPSILON = 1E-6;

%[relCandidatePeakIndices, relCandidateResidueIndices, scores, numCandidates] = ...
%    voteAccordingToMBM_EM(voter, numVoters, ASSIGNTABLE, peakIndices,residueIndices); %or
										      %some
										      %variant
										      %of
										      %mbm-em
%[relCandidatePeakIndices, relCandidateResidueIndices, scores, numCandidates] = ...
%    voteAccordingToLowestScore(voter, numVoters, ASSIGNTABLE, peakIndices,residueIndices); %or

[relCandidatePeakIndices, relCandidateResidueIndices, scores, numCandidates] = ...
    voteAccordingToBayesianScore(voter, numVoters, ASSIGNTABLE, peakIndices,residueIndices); %or

%[relCandidatePeakIndices, relCandidateResidueIndices, scores, numCandidates] = ...
%    voteAccordingToMarsScore(voter, numVoters, ASSIGNTABLE,peakIndices,residueIndices); %or

% $$$ [relCandidatePeakIndices, relCandidateResidueIndices, scores, numCandidates] = ...
% $$$     voteAccordingToMB_Score(MB_ShiftScore, ...
% $$$ 			    MB_RDC_Score, ...
% $$$ 			    csCoefficient, ...
% $$$ 			    rdcCoefficient, ...
% $$$ 			    noeCoefficient, ...
% $$$ 			    ASSIGNTABLE,peakIndices,residueIndices,MASTER,HSQCDATA,ALLDISTS,NTH,numCS,numRDC,numHN_NOES);




fprintf(1, 'expanding node#%d\n',nodeIndex);

							    
for i = 1:numCandidates
   fprintf(1, 'peakIndex : %d residueIndex:%d score:%f\n', peakIndices(relCandidatePeakIndices(i)),residueIndices(relCandidateResidueIndices(i)),scores(i));
end

%keyboard
  
for candidateIndex = 1:numCandidates

  
  fprintf(1, 'candidate #%d\n',candidateIndex);
  
  newAssigntable = updateAssignTable(ASSIGNTABLE, ...
				     relCandidatePeakIndices(candidateIndex), ...
				     relCandidateResidueIndices(candidateIndex)); 					 
  

  [impossibleToAssign, newMaster, newAssigntable, newUnassignedPeakIndices, newUnassignedResidueIndices, newScoreSoFar, newVoter, newNumVoters] = ...
      noePruneAndDoMoreAssignments(newAssigntable, MASTER, ...
				   peakIndices, ...
				   residueIndices, ...
				   scoreSoFar,...
				   voter, numVoters,...
				   NOES, ...
				   IALLDISTS, ALLDISTS, NTH, RDC1, ...
				   RDC2, VECTORS,MB_ShiftScore, ...
				   MB_RDC_Score, numCS, numRDC, ...
				   HSQCDATA, csCoefficient, ...
				   rdcCoefficient, noeCoefficient,numHN_NOES);

%  if ((sum(sum(newMaster))  - sum(sum(MASTER))) > 1)
%    fprintf(1, 'did %d assignments in noePruneAndDoMoreAssignments.\n',(sum(sum(newMaster))  - sum(sum(MASTER))));
%    fprintf(1, 'total num assignments = %d, before it was %d\n',sum(sum(newMaster)),sum(sum(MASTER)));
%    keyboard
%  end
  fprintf(1, 'peakIndex : %d residueIndex:%d score:%f\n', peakIndices(relCandidatePeakIndices(candidateIndex)),residueIndices(relCandidateResidueIndices(candidateIndex)),newScoreSoFar);

  if (impossibleToAssign)
    fprintf(1, 'impossible to assign. continuing.\n');
    continue;
  end
  
% newH_Value = computeH_Value(newAssigntable, newVoter,newNumVoters);


% newH_Value = computeMarsH_Value(newAssigntable, newVoter, newNumVoters);

  newH_Value = computeBetterH_Value(newAssigntable, newVoter, ...
				       newNumVoters);
  
%  fprintf(1, 'difference between betterH_Value and newH_Value is %f\n', betterH_Value-newH_Value);
 
%  assert (betterH_Value >= newH_Value-EPSILON);
  
%  keyboard

%  if ((peakIndices(relCandidatePeakIndices(candidateIndex)) == 1) & (residueIndices(relCandidateResidueIndices(candidateIndex))==44))
%    keyboard
%  end

% $$$   MB_H_Value = computeMB_H_Value(newAssigntable, newMaster, HSQCDATA,...
% $$$ 				 newUnassignedPeakIndices, newUnassignedResidueIndices,...
% $$$ 				 ALLDISTS, NTH, ...
% $$$ 				 MB_ShiftScore, ...
% $$$ 				 MB_RDC_Score, ...
% $$$ 				 numCS,numRDC,csCoefficient, rdcCoefficient, ...
% $$$ 				 noeCoefficient,numHN_NOES);
%  assert (MB_H_Value >= h_values(nodeIndex)); This is not
%  necessarily true
  
% $$$   S1                               = updateTen(MASTER,RDC1,VECTORS);
% $$$   S2                               = updateTen(MASTER,RDC2,VECTORS);
% $$$   
% $$$   RP1                              = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1,ROWIN,COLIN);
% $$$   RP2                              = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2,ROWIN,COLIN);
% $$$   
% $$$ 
% $$$   nvrVoter





  [numNodes, next, ASSIGNTABLES, MASTERS, voters, numVoterss, peakIndicess, residueIndicess, scoreSoFars, h_values] = insertIntoQueue(numNodes, nodeIndex, next, newMaster, newScoreSoFar, newH_Value, newVoter, newNumVoters,  newAssigntable, newUnassignedPeakIndices, newUnassignedResidueIndices, ...
						  ASSIGNTABLES, ...
						  MASTERS, voters, ...
						  numVoterss, ...
						  peakIndicess, ...
						  residueIndicess, ...
						  scoreSoFars, h_values);
end

%would be good to test this function by passing it a state that has
%been seen before.
function   [numNodes, next, ASSIGNTABLES, MASTERS, voters, numVoterss, peakIndicess, residueIndicess, scoreSoFars, h_values] = ...
    insertIntoQueue(numNodes, nodeIndex, next, newMaster, newScoreSoFar, newH_Value, newVoter, newNumVoters, ...
		    newAssigntable, ...
		    newUnassignedPeakIndices, newUnassignedResidueIndices, ...
		    ASSIGNTABLES, ...
		    MASTERS, voters, ...
		    numVoterss, ...
		    peakIndicess, ...
		    residueIndicess, ...
		    scoreSoFars, h_values)

EPSILON                                = 2E-5;
insertToTheBeginningOfEqualValuedNodes = 0;

newScore                               = newScoreSoFar + newH_Value;
foundThisNodeBefore                    = 0;


assert (nodeIndex > 0);
assert (nodeIndex <= numNodes);
%now find where to insert the new node, adjust next array.
assert (scoreSoFars(nodeIndex)+h_values(nodeIndex) <= newScore+EPSILON);
%the new node should have higher cost than the predecessor

while (nodeIndex > 0) & ((scoreSoFars(nodeIndex) + h_values(nodeIndex))<= ...
			 newScore + EPSILON) & (~insertToTheBeginningOfEqualValuedNodes)
 
  if (abs(scoreSoFars(nodeIndex) + h_values(nodeIndex) - newScore) < EPSILON)
    
    if (equal(ASSIGNTABLES{nodeIndex},newAssigntable) & equal(MASTERS{nodeIndex},newMaster) ...
	  & (equal(peakIndicess{nodeIndex} , newUnassignedPeakIndices)) & ...
	  (equal(residueIndicess{nodeIndex}, newUnassignedResidueIndices)) & (newNumVoters == numVoterss(nodeIndex)))
      %do we need to check the contents of the voters arrays as
      %well? SB
      %I think not since the contents of MASTER, ASSIGNTABLE,
      %unassignedPeakIndices, unassignedResidueIndices, numVoters
      %being the same should imply that the remaining entries in
      %the voter arrays should match as well.
      %can check whether voters{nodeIndex}{voterIndex} matches
      %newVoter{voterIndex}, voterIndex = 1..numVoters. But this
      %requires using an epsilon due to numerical inaccuracies.

      
      for voterIndex = 1:newNumVoters
	assert (twoMatricesAreClose(voters{nodeIndex}{voterIndex}, ...
				    newVoter{voterIndex}, EPSILON) ...
		== 1);
      end

      fprintf(1, 'found this node before\n');
      
      foundThisNodeBefore = 1;
      break;
    else
      newNumAssignedPeaks        = sum(sum(newMaster));
      numAssignedPeaks_nodeIndex = sum(sum(MASTERS{nodeIndex}));
      
      if (newNumAssignedPeaks > numAssignedPeaks_nodeIndex)
	fprintf(1, 'found that the score of node to be inserted');
	fprintf(1, ' equals the node in the queue but ');
	fprintf(1, ' the node to be inserted has more assigned peaks. ');
	fprintf(1, '%d vs. %d. Therefore inserting after node %d\n', newNumAssignedPeaks, numAssignedPeaks_nodeIndex,nodeIndex);
	insertToTheBeginningOfEqualValuedNodes = 1;
	%keyboard
      end
      
    end
  end
  
  currentNodeIndex = nodeIndex;
  nodeIndex        = next(nodeIndex);
end

if (foundThisNodeBefore)
  %insert nothing. exit.
  return;
end

assert (scoreSoFars(currentNodeIndex) + h_values(currentNodeIndex) ...
	<= newScore + EPSILON);

fprintf(1, 'inserting node#%d\n',numNodes+1);
fprintf(1, 'its g value is %f, h value is %f ',newScoreSoFar, ...
	newH_Value);
fprintf(1, 'its total score is %f\n', newScore);

fprintf(1, 'the previous node has g=%f h=%f ',scoreSoFars(currentNodeIndex),h_values(currentNodeIndex));
fprintf(1, 'total score = %f numAssignedPeaks = %d\n', scoreSoFars(currentNodeIndex)+ ...
	h_values(currentNodeIndex), sum(sum(MASTERS{currentNodeIndex})));

%if (insertToTheBeginningOfEqualValuedNodes)
%  keyboard
%end

next(currentNodeIndex)      = numNodes+1;
next(numNodes+1)            = nodeIndex;


MASTERS{numNodes+1}         = newMaster;
ASSIGNTABLES{numNodes+1}    = newAssigntable;

if (newNumVoters == 0)
  voters{numNodes+1}        = cell(0,1);
else
  for i = 1:newNumVoters
    voters{numNodes+1}{i}     = newVoter{i};
  end
end
numVoterss(numNodes+1)      = newNumVoters;
peakIndicess{numNodes+1}    = newUnassignedPeakIndices;
residueIndicess{numNodes+1} = newUnassignedResidueIndices;
scoreSoFars(numNodes+1)     = newScoreSoFar;
h_values(numNodes+1)        = newH_Value;

numNodes                    = numNodes + 1;










%keyboard

function retval = equal (matrix1, matrix2)

[sizeX1, sizeY1] = size(matrix1);
[sizeX2, sizeY2] = size(matrix2);

if (sizeX1 ~= sizeX2) | (sizeY1 ~= sizeY2)
  retval = 0;
  return ;
end

comparisonMatrix = (matrix1 ~= matrix2);

if isempty(find(comparisonMatrix))
  retval = 1;
  return;
else
  retval = 0;
  return;
end
  



function [indexOfNoesThatMatchHN, indexOfNoesThatMatchH2]= ...
    computeNoeHSQC_Correspondence(HSQCDATA); 

%make these static variables SB or compute once and pass it on.

[h_ppm n_ppm h_ppm2] = textread('dnns.txt','%f %f %f');
indexOfNoesThatMatchHN = cell(size(HSQCDATA,1),1);
indexOfNoesThatMatchH2 = cell(size(HSQCDATA,1),1);
for i = 1:size(HSQCDATA,1)
  indexOfNoesThatMatchHN{i} = [];
  indexOfNoesThatMatchH2{i} = [];
end

for noeIndex = 1:length(h_ppm)
  hn_closePeaks                               = findCloseToHN(h_ppm(noeIndex),n_ppm(noeIndex),  HSQCDATA);
  secondProtonClosePeaks                      = findCloseToH(h_ppm2(noeIndex), ...
						  HSQCDATA);
  for hnClosePeakIndex =1:length(hn_closePeaks)
    indexOfNoesThatMatchHN{hn_closePeaks(hnClosePeakIndex)} = [ ...
	indexOfNoesThatMatchHN{hn_closePeaks(hnClosePeakIndex)} noeIndex];
  end
  
  for secondProtonClosePeakIndex = 1:length(secondProtonClosePeaks)
    indexOfNoesThatMatchH2{secondProtonClosePeaks(secondProtonClosePeakIndex)} = [indexOfNoesThatMatchH2{secondProtonClosePeaks(secondProtonClosePeakIndex)} noeIndex];
  end
end

function numNOEsBetweenPeaks  = ...
    findNOEsBetweenThesePeaks(indexOfNoesThatMatchHN,indexOfNoesThatMatchH2,peakIndex1,peakIndex2)

numNOEsBetweenPeaks = 0;
for noeIndex = 1:length(indexOfNoesThatMatchHN{peakIndex1})
  noeIndexMatchingHN_OfPeakIndex1 = ...
      indexOfNoesThatMatchHN{peakIndex1}(noeIndex);
  for noe2Index = 1:length(indexOfNoesThatMatchH2{peakIndex2})
    noeIndexMatchingH2_OfPeak2 = ...
	indexOfNoesThatMatchH2{peakIndex2}(noe2Index);
    if (noeIndexMatchingHN_OfPeakIndex1 == ...
	noeIndexMatchingH2_OfPeak2)
      numNOEsBetweenPeaks = numNOEsBetweenPeaks + 1;
      break;
    end
  end
end

for noeIndex = 1:length(indexOfNoesThatMatchHN{peakIndex2})
  noeIndexMatchingHN_OfPeakIndex2 = ...
      indexOfNoesThatMatchHN{peakIndex2}(noeIndex);
  for noe2Index = 1:length(indexOfNoesThatMatchH2{peakIndex1})
    noeIndexMatchingH2_OfPeak1 = ...
	indexOfNoesThatMatchH2{peakIndex1}(noe2Index);
    if (noeIndexMatchingHN_OfPeakIndex2 == ...
	noeIndexMatchingH2_OfPeak1)
      numNOEsBetweenPeaks = numNOEsBetweenPeaks + 1;
      break;
    end
  end
end