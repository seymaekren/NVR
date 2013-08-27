function [overallMASTER,foundAnAssignment, ...
	  stopGeneratingAssignments, minScore, totalNumAssignments, foundMASTERs, finalScores, assignmentAccuracies] = assign (totalNumAssignments, finalScores, foundMASTERs, foundAnAssignment, minScore, MASTER, ...
			       scoreSoFar, voter, numVoters, ASSIGNTABLE, ...
			      unassignedPeakIndices, ...
						  unassignedResidueIndices,NOES, ALLDISTS, IALLDISTS, NTH, overallMASTER, peakIDs, RESNUMS, RDC1, RDC2, VECTORS, assignmentAccuracies)

THRESHOLD = 200; %corresponds to a 2-order-of-magnitude change in the
               %probability of assignments

%MAX_NUM_RETURNED_ASSIGNMENTS = 100;

numAssignmentsSoFar = sum(sum(MASTER));
numPeaks            = size(MASTER,1);

stopGeneratingAssignments = 0;


if (numAssignmentsSoFar == numPeaks)
  fprintf(1, 'found an assignment\n');
%if assigned all the peaks
   if (foundAnAssignment == 0)
     
      minScore            = scoreSoFar;
    
      [overallMASTER, totalNumAssignments, foundMASTERs, finalScores, assignmentAccuracies]       = reportAndStoreAssignments(MASTER, scoreSoFar, ...
						  overallMASTER, peakIDs, ...
						  RESNUMS, foundMASTERs, finalScores, ...
						  totalNumAssignments, ...
						  assignmentAccuracies);
    else
      
      %found an assignment before.
      if (scoreSoFar - minScore < THRESHOLD)
	[overallMASTER, totalNumAssignments, foundMASTERs, finalScores, assignmentAccuracies] = reportAndStoreAssignments(MASTER, scoreSoFar, ...
						  overallMASTER, ...
						  peakIDs, RESNUMS, ...
						  foundMASTERs, ...
						  finalScores, ...
						  totalNumAssignments, assignmentAccuracies);
      else
	fprintf(1, 'score too high. it is : %f\n', scoreSoFar);
	stopGeneratingAssignments = 1; %first time an assignment has
                                     %lower probability, we stop
                                     %generating new assignments.
      end
    
      if (scoreSoFar < minScore)
	%interesting. found an assignment that actually has higher
	%probability than the earliest assignment we returned.
	%this means that we do not get the assignments in the sorted
	%order of decreasing overall probabilities. This is rather
	%unexpected since in theory MBM_EM returns the candidates in
	%the order of decreasing probabilities. It may be that the
	%voting scheme MBM_EM uses and the probability calculation do
	%not match each other. 
%	keyboard
	minScore = scoreSoFar;
      end
    
    end
  
  foundAnAssignment = 1;      
      
      
else
  
  fprintf(1, 'assigned %d peaks so far\n', numAssignmentsSoFar);
  
  
  [relCandidatePeakIndices, relCandidateResidueIndices, scores, numCandidates] = ...
      voteAccordingToMBM_EM(voter, numVoters, ASSIGNTABLE);%, NOES, ...

 
  
  
  for i = 1:numCandidates
    fprintf(1, 'peakIndex : %d residueIndex:%d score:%f\n', unassignedPeakIndices(relCandidatePeakIndices(i)),unassignedResidueIndices(relCandidateResidueIndices(i)),scores(i));
  end
  
  
  
  for peakIndex = 1:numCandidates
    
    newAssigntable = updateAssignTable(ASSIGNTABLE, ...
				       relCandidatePeakIndices(peakIndex), ...
				       relCandidateResidueIndices(peakIndex)); 

%    fprintf(1, 'updated assigntable.\n');
    
    
    
    [impossibleToAssign, newMaster, newAssigntable, newUnassignedPeakIndices, newUnassignedResidueIndices, newScore, newVoter, newNumVoters] = ...
	noePruneAndDoMoreAssignments(newAssigntable, MASTER, ...
				     unassignedPeakIndices, ...
				     unassignedResidueIndices, ...
				     scoreSoFar, ...
				     voter, numVoters,...
				     NOES, ...
				     IALLDISTS, ALLDISTS, NTH, RDC1, ...
				     RDC2, VECTORS);
    
%    fprintf(1, 'called noeprune and do more assignments\n');
    
    
    if (impossibleToAssign)
      fprintf(1, 'impossible to assign. continuing.\n');
      continue;
    end

    %updates numVoters.
							
    %      newAssignmensMadeSoFar =  assignmentsMadeSoFar + ...
    %	  (peakIndices(relPeakIndex), residueIndices(relPeakIndex)); ...
    %prune assigntable
							
    %newProbabilityOfAssignmentsSoFar = probabilityOfAssignmentsSoFar * \pi_voterIndex=1..numVoters ...
    %	  voter{voterIndex}(peakIndices(relPeakIndex), ...
    %			    residueIndices(relPeakIndex));
							
    %      newAssignTable = update assigntable; % remove peakIndices(relPeakIndex), ...
    % residueIndices(relPeakIndex)
    %			  newVoters= update voters;
    %      newNumVoters   = update numVoters;
							
    %      update voters;
							
%    fprintf(1, 'calling assign recursively.\n');
[overallMASTER, foundAnAssignment, stopGeneratingAssignments, minScore, totalNumAssignments, foundMASTERs, finalScores, assignmentAccuracies] = ...
	assign(totalNumAssignments, finalScores, foundMASTERs, foundAnAssignment, minScore, newMaster, newScore, newVoter, newNumVoters, ...
	       newAssigntable, ...
	       newUnassignedPeakIndices, newUnassignedResidueIndices, ...
	       NOES, ALLDISTS,IALLDISTS, NTH, overallMASTER, peakIDs, ...
	       RESNUMS, RDC1, RDC2, VECTORS, assignmentAccuracies);
    
%    fprintf(1, 'returned from recursive assign call.\n');
    
    if (stopGeneratingAssignments)
      fprintf(1, 'stop generating assignments is true.\n');
      break;
    end
  end
end


% $$$ while (numUnassignedPeaks > 0)
% $$$ 
% $$$   %first updating voters
% $$$ 
% $$$   if (numAssignedPeaks >= 5)
% $$$     [voter,numVoters,S1,S2] = ...
% $$$ 	updateVotersAndAlignmentTensors(voter, ASSIGNTABLE,...
% $$$ 					MASTER, S1, S2, RDC1,RDC2, VECTORS);
% $$$   end
% $$$ 
% $$$    if (useMBM_EM)
% $$$      [bestRelPeakIndex,bestRelResidueIndex] = ...
% $$$  	voteAccordingToMBM_EM    (voter,     numVoters,ASSIGNTABLE,NOES,ALLDISTS,NTH,peakIndices,residueIndices);
% $$$    else
% $$$      [bestRelPeakIndex,bestRelResidueIndex] = ...
% $$$ 	 voteAccordingToSimpleSum (numVoters, voter);
% $$$    end
% $$$ 
% $$$ %  [bestRelPeakIndex,bestRelResidueIndex] = ...
% $$$ %      voteAccordingToKemeny    (voter,     numVoters,ASSIGNTABLE,NOES,ALLDISTS,NTH,peakIndices,residueIndices);
% $$$       
% $$$   bestPeakIndex    = peakIndices(bestRelPeakIndex);
% $$$   bestResidueIndex = residueIndices(bestRelResidueIndex);
% $$$ 
% $$$   fprintf(1,'votes decided to assign %d to %d\n', bestPeakIndex, ...
% $$$ 	  bestResidueIndex);
% $$$   [numUnassignedPeaks, numUnassignedResidues] = size(voter{1});  
% $$$    votes = zeros(numUnassignedResidues, 1);
% $$$ %  scores = votes;
% $$$   relPeakIndex = bestRelPeakIndex;
% $$$   
% $$$   for voterIndex = 1:numVoters
% $$$     weightedRankings                            = voter{voterIndex}(relPeakIndex,:);
% $$$     votes = votes + weightedRankings';
% $$$   end	
% $$$   figure; plot(votes', '*'); length(nonzeros(votes))
% $$$   [i,j] = max(votes)
% $$$   bestRelResidueIndex
% $$$   keyboard
% $$$   
% $$$   if (debugIncorrectAssignment)
% $$$ 
% $$$     debugIncorrectAssignmentFunction(voter, numVoters, ...
% $$$ 				     bestRelPeakIndex, ...
% $$$ 				     bestRelResidueIndex, bestPeakIndex,bestResidueIndex,residueIndices);
% $$$   end
% $$$ 
% $$$   
% $$$   ASSIGNTABLE             = updateAssignTable (ASSIGNTABLE,bestRelPeakIndex,bestRelResidueIndex);
% $$$   
% $$$   numPotentialAssignments = sum(sum(ASSIGNTABLE));
% $$$ 
% $$$   for i = 1:100
% $$$     %uncomment the following for NVR type NOE pruning
% $$$     %     [ASSIGNTABLE]        = noePrune(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,peakIndices,residueIndices);
% $$$     
% $$$     [ASSIGNTABLE]        = noePrune(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, peakIndices,residueIndices);
% $$$     [ASSIGNTABLE]        = noePrune(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 40, peakIndices,residueIndices);
% $$$ 
% $$$     [MASTER,ASSIGNTABLE] = doAssignments(MASTER,ASSIGNTABLE,peakIndices, ...
% $$$ 					 residueIndices);
% $$$     [ASSIGNTABLE,voter, peakIndices,residueIndices] ...
% $$$ 	= all_removeAssignedResiduesAndUpdateVoteBPG(ASSIGNTABLE,voter,numVoters,peakIndices,residueIndices);
% $$$   
% $$$ 
% $$$     newNumPotentialAssignments = sum(sum(ASSIGNTABLE));
% $$$     
% $$$     if(newNumPotentialAssignments == numPotentialAssignments)
% $$$       break;
% $$$     end
% $$$ 
% $$$     numPotentialAssignments = newNumPotentialAssignments;
% $$$     
% $$$   end
% $$$ 
% $$$   numAssignedPeaks      = sum(sum(MASTER));
% $$$   numUnassignedPeaks    = totalNumPeaks    - numAssignedPeaks;
% $$$ 
% $$$   
% $$$ end
% $$$ 
% $$$ 
% $$$ computeCorrectness(MASTER,RESNUMS, peakIDs);















% $$$ %finds unassigned peaks and residues, stores their indices in peakIndices
% $$$ %and residueIndices. However does not resize peakIndices and residueIndices, just removes
% $$$ %gaps. SO AT THE End of this function only the first consecutive x
% $$$ %entries in peakIndices and residueIndices are unassigned, but this function does
% $$$ %not explicitly return x.
% $$$ function [peakIndices,residueIndices]=removeGapsInUnassignedPeakResidueLists(MASTER,peakIndices,residueIndices)
% $$$ sct=1;
% $$$ for(i=1:size(MASTER,1))
% $$$    x = find(MASTER(i,:));
% $$$    if(length(x)==0)
% $$$       peakIndices(sct)=i;
% $$$       sct=sct+1;
% $$$    end
% $$$ end
% $$$ sct=1;
% $$$ for(i=1:size(MASTER,2))
% $$$    x = find(MASTER(:,i));
% $$$    if(length(x)==0)
% $$$       residueIndices(sct)=i;
% $$$       sct=sct+1;
% $$$    end
%end







%pseudocode:
%  for each peak :
%    for each voter :
%      get votes, sort them.
%      select the one having highest prob. increase that residue's vote count.
%    find the residue that had the highest votes, and its associated
%    probability (multiplication of all probabilities).

%  select the peak-residue assignment with most votes and highest
%  probability.
%  do that assignment.
%  remove that peak and residue from the resp. lists.
%  normalize the probability of each voter for each peak.
function [bestRelPeakIndex,bestRelResidueIndex] = voteAccordingToSimpleSum(numVoters,voter)


[numUnassignedPeaks, numUnassignedResidues] = size(voter{1});

tentativeAssignResidue = zeros(numUnassignedPeaks,1);
% $$$   tentAssignProb = zeros(numUnassignedPeaks,1);
tentativeAssignScore = zeros(numUnassignedPeaks,1);

%assert (voteAccordingToSimpleSumOfVotersVotes);

for relPeakIndex = 1:numUnassignedPeaks
  votes = zeros(numUnassignedResidues, 1);
%  scores = votes;
  
  for voterIndex = 1:numVoters
    weightedRankings                            = voter{voterIndex}(relPeakIndex,:);
    votes = votes + weightedRankings';
  end	
  
  [maxScore, winnerRelResidueIndex] = max(votes);
  
  winnerRelResidueIndices = find(votes == maxScore);
  
  tentativeAssignScore  (relPeakIndex)        = maxScore;
  
  assert (length(winnerRelResidueIndices) == 1);
  assert (winnerRelResidueIndices == winnerRelResidueIndex);
  
  tentativeAssignResidue(relPeakIndex)        = winnerRelResidueIndex;
% $$$     tentAssignProb(relPeakIndex)    = probabilityOfWinner;
  
end	

[highestScore, bestRelPeakIndex]      = max(tentativeAssignScore);

bestRelPeakIndices = find(tentativeAssignScore == highestScore);

assert (length(bestRelPeakIndices) == 1);
assert (bestRelPeakIndices == bestRelPeakIndex);

%  if (length(bestRelPeakIndices) == 1)
%    bestRelPeakIndex = bestRelPeakIndices;
%  else
%    [highestProbability, bestRelRelPeakIndex] = max(tentAssignProb(bestRelPeakIndices));

%    bestRelPeakIndex                          = bestRelPeakIndices(bestRelRelPeakIndex);
%  end

bestRelResidueIndex = tentativeAssignResidue(bestRelPeakIndex);
%keyboard
  

  



% $$$ for(i=1:size(V,1))
% $$$   x = find(V(i,:));
% $$$   if(length(x)==1)
% $$$     bestRelPeakIndex = i;
% $$$     bestRelResidueIndex = x;
% $$$     break;
% $$$     
% $$$ %actually here may not return early so as to do all assignments but
% $$$ %it may be safer to go slowly.    
% $$$ %    ASSIGNTABLE(i,:)=0;
% $$$ %    ASSIGNTABLE(i,x)=1;
% $$$   end
% $$$ end

function [bestRelPeakIndex,bestRelResidueIndex] =  voteAccordingToKemeny    (voter,     numVoters)

% $$$ for each peak;
% $$$   for each voter;
% $$$     
% $$$     find the best residue and second best residue;
% $$$     compute the probability difference between them;
% $$$     
% $$$   end
% $$$   
% $$$   if the number of times the best residue is selected the same for ...
% $$$ 	  all voters; record this case;
% $$$   end
% $$$   
% $$$ end
% $$$ 
% $$$ return the peaik that has the highest number of time s it is assigned ...
% $$$     to the same residue;


for peakIndex = 1:numPeaks
  for voterIndex = 1:numVoters
    
    [sortedVotes, residueIndices] = sort(-votes{voterIndex});
    
    margin(voterIndex)         = sortedVotes(1) - sortedVotes(2);
    residueIndex(voterIndex)   = residueIndices(1);
    
  end
  
  numConcurringVoters(peakIndex) = length(find(residueIndex == ...
					       mode(residueIndex)));
  
  sumOfMargin(peakIndex) = sum(margin(find(residueIndex == mode(residueIndex))));
end


bestPeakNumConcurringVoters = max  (numConcurringVoters);

bestPeakIndices             = find (numConcurringVoters == ...
				    bestPeakNumConcurringVoters);

[bestPeakIndexMargin,relBestPeakIndex] = ...
    max(sumOfMargin(bestPeakIndices));
bestPeakIndex = bestPeakIndices(relBestPeakIndex);


fprintf(1, 'vote according to Kemeny is currently unfinished.');
%bestResidueIndex = 





function voterPreferences = computeMargin (voterPreferences,M,peakIndex,voterIndex,assignedResidueIndex);

assignmentProb                                   = M(peakIndex, assignedResidueIndex);
[bestAssignmentProb,bestAssignmentResidueIndex]  = max(M(peakIndex,:));

if (bestAssignmentResidueIndex == assignedResidueIndex)
  %in case the resonance assignment corresponds to the highest
  %probability, which is hopefully the case, we would like to
  %compute the margin using the difference between the best
  %and second best assignment.
  
  M(peakIndex,assignedResidueIndex)              = 0;
  SecondBestAssignmentProb                       = max(M(peakIndex,:));
  M(peakIndex,assignedResidueIndex)              = assignmentProb;
  
  %voter conforms to MASTER for this peak.
  
  voterPreferences.margins(peakIndex,voterIndex) = assignmentProb - SecondBestAssignmentProb; 
  
  assert (voterPreferences.margins(peakIndex,voterIndex) >= 0);
else
  
  %voter does not conform to MASTER for this peak.
  
  voterPreferences.margins(peakIndex,voterIndex) = assignmentProb - bestAssignmentProb;

  assert (voterPreferences.margins(peakIndex,voterIndex) <= 0);

end %if


%-----
function debugIncorrectAssignmentFunction(voter, numVoters, ...
					  bestRelPeakIndex, ...
					  bestRelResidueIndex, bestPeakIndex,bestResidueIndex,residueIndices)

[numUnassignedPeaks, numAssignedResidues] = size(voter{1});
 
if (bestPeakIndex ~= bestResidueIndex)
  
  desiredResidueIndex    = bestPeakIndex;
  desiredRelResidueIndex = find(residueIndices == desiredResidueIndex);
  
  %for relPeakIndex = 1:numUnassignedPeaks
  relPeakIndex = bestRelPeakIndex;
  votes = zeros(numUnassignedResidues, 1);
  
  for voterIndex = 1:numVoters
    weightedRankings                            = voter{voterIndex}(relPeakIndex,:);
    
    [sortedWeightedRankings, relResidueIndices] = sort(-weightedRankings);
    
    rankOfCorrectAssignmentAccordingToVoter = find(relResidueIndices == desiredResidueIndex)
    rankOfIncorrectAssignmentAccordingToVoter = find(relResidueIndices == bestResidueIndex)
    
    differenceInProbabilityAccordingToThisVoter = ...
	sortedWeightedRankings(rankOfIncorrectAssignmentAccordingToVoter) - ...
	sortedWeightedRankings(rankOfCorrectAssignmentAccordingToVoter)
    
    %what if two peaks have the same value?
    %%%winnerRelResidueIndex            = relResidueIndices(1);
    %%%votes(winnerRelResidueIndex)     =
    %votes(winnerRelResidueIndex) + 1;
    
    votes = votes + weightedRankings';
  end	
  
  
  %end
  
  fprintf(1, 'the desired residue had a total sum of %f\n', ...
	  votes(desiredRelResidueIndex));
  
  sortedVotes = sort(-votes);
  sortedVotes(find(sortedVotes))
  fprintf(1, 'the sorted votes are as above\n');
  
  keyboard
end



% $$$ function [voter,numVoters,S1,S2,RP1,RP2] = ...
% $$$     updateVotersAndRDC_BipartiteGraphsAndAlignmentTensors( voter, ASSIGNTABLE,SSCP, SXCP, CP, TP, HDE, ...
% $$$ 						  MASTER, S1, S2, RDC1,RDC2,RP1,RP2, ...
% $$$ 						  VECTORS)
% $$$   numVoters = 5;
% $$$ 
% $$$ 
% $$$   voter{1}  = SSCP;
% $$$   voter{2}  = SXCP;
% $$$   voter{3}  = CP;
% $$$   voter{4}  = TP;
% $$$   voter{5}  = HDE;
% $$$ 
% $$$   numAssignedPeaks    = sum(sum(MASTER));
% $$$   
% $$$ %  if (useRDCs)
% $$$     if (numAssignedPeaks >= 5)
% $$$       S1  = updateTen(MASTER,S1,RDC1,VECTORS);
% $$$       S2  = updateTen(MASTER,S2,RDC2,VECTORS);
% $$$       RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
% $$$       RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
% $$$       
% $$$       numVoters = 7;
% $$$       voter{6} = RP1;
% $$$       voter{7} = RP2;
% $$$     end
% $$$ %  end



function computeCorrectness(MASTER,RESNUMS, peakIDs);

%now, compute correctness

a = load('answerkey.m');
assignments=0;
correct=0;
incorr=0;

for(i=1:size(MASTER,1))
   pk=peakIDs(i);%get the peak id
   x = find(MASTER(i,:));
    
   rn=RESNUMS(x);%get the residue it was assigned to 
   
   assignments(i,1)=pk;
   assignments(i,2)=rn;
  
   foo = find(a(:,1)==pk);
   if(rn==a(foo,2))
      correct=correct+1;
      incorr(i)=0;
    else
      incorr(i)=1;
   end

end



fprintf('\n');
fprintf('\n');
fprintf('Assignment Accuracy = %f percent \n',correct/size(MASTER,1)*100);
assignmentAccuracy = correct/size(MASTER,1)*100;
fprintf('\n');
fprintf('\n');
% $$$ fprintf('Assignments\n');
% $$$ fprintf('Peak ID -> Residue Number\n');
% $$$ for(i=1:size(assignments,1))
% $$$    if(incorr(i)==0)
% $$$       fprintf('%d	%d\n',assignments(i,1),assignments(i,2));
% $$$    else
% $$$       fprintf('*%d	%d  (correct=%d->%d)\n',assignments(i,1),assignments(i,2),a(i,1),a(i,2));
% $$$    end
% $$$ end
% $$$ fprintf('\n');
% $$$ fprintf('\n');
%keyboard
if(size(MASTER,1)~=size(MASTER,2))
   %find which residues aren't assigned
   ct=1;
   missing=[];
   
   for(i=1:length(RESNUMS))
      rn = RESNUMS(i);
      if(length(find(assignments(:,2)==rn))==0)
         missing(ct)=rn;
         ct=ct+1;
      end
   end
   
   fprintf('The following residues are missing peaks\n');
   for(i=1:length(missing))
     fprintf('%d\n',missing(i));
   end
end















function analyzeVoters(voter, numVoters, peakIndices, ...
		       residueIndices);
for i = 1:numVoters
  for peakIndex = 1:length(peakIndices)
    votes = voter{i}(peakIndex,:);
    [sortedVotes, sortedOrder] = sort(-votes)
    keyboard
  end
end
