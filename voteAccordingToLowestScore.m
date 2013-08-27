function  [relCandidatePeakIndices, relCandidateResidueIndices, scores, numCandidates] =  voteAccordingToLowestScore(voter, numVoters, ASSIGNTABLE,peakIndices,residueIndices);

%bestRelPeakIndex    = -1;
%bestRelResidueIndex = -1;

MAX_NUMCANDIDATES          = 100;
relCandidatePeakIndices    = zeros(MAX_NUMCANDIDATES,1); 
relCandidateResidueIndices = relCandidatePeakIndices;
scores                     = relCandidatePeakIndices;
numCandidates              = 0; 



assert (~isempty(ASSIGNTABLE));

if (numVoters == 5)
  V       = voter{1}.*voter{2}.*voter{3}.*voter{4}.*voter{5};%CP,SXCP,SSCP,TP,HDE,
  %ASSIGNTABLE);
elseif (numVoters == 7)
  V       = voter{1}.*voter{2}.*voter{3}.*voter{4}.*voter{5}.*voter{6}.*voter{7};%,...%CP,SXCP,SSCP,TP,RP1,RP2,HDE,
%	   ASSIGNTABLE);
elseif (numVoters == 1)
  V       = voter{1};
end

V         = V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));

%maxVotesBeforeMaskingWithASSIGNTABLE = max(max(V));

V         = V.*ASSIGNTABLE;

stopGeneratingCandidates = 0;



% $$$ combinedVoter = voter{1}.*voter{2}.*voter{3}.*voter{4} .*voter{5}.*voter{6} .*voter{7};
% $$$ for candidateIndex=1:100
% $$$   maxCombinedProbability         = max(max(combinedVoter));
% $$$   if (maxCombinedProbability == 0)
% $$$     break;
% $$$   end
% $$$   [relPeakIndex,relResidueIndex] = find(combinedVoter == ...
% $$$ 					maxCombinedProbability);
% $$$   fprintf(1, 'assignment of %d to %d has score %f\n',peakIndices(relPeakIndex),residueIndices(relResidueIndex),-log(maxCombinedProbability));
% $$$   combinedVoter(relPeakIndex,relResidueIndex) = 0.0;
% $$$ end

% $$$ for relPeakIndex=1:length(peakIndices)
% $$$   peakIndex        = peakIndices(relPeakIndex);
% $$$   relResidueIndex  = find(residueIndices == peakIndex);
% $$$   fprintf(1, 'vote of peak %d-residue%d assignment is %d\n',peakIndex,peakIndex,V(relPeakIndex,relResidueIndex));
% $$$   score = 0;
% $$$   assert (numVoters == 7);
% $$$   for voterIndex = 1:7
% $$$     assert (voter{voterIndex}(relPeakIndex,relResidueIndex)~=0);
% $$$     score = score - log(voter{voterIndex}(relPeakIndex,relResidueIndex));
% $$$   end
% $$$   fprintf(1, 'score of this assignment is %f\n', score);
% $$$ end
% $$$ 
% $$$ keyboard



while ((numCandidates < MAX_NUMCANDIDATES) & (~stopGeneratingCandidates))
%while (~stopGeneratingCandidates)

  maxCombinedProbability  = max(max(V));
  
  if (maxCombinedProbability == 0)
    break;
  end
  
  %if (maxVotes ~= maxVotesBeforeMaskingWithASSIGNTABLE)
  %  fprintf(1, 'ASSIGNTABLE masked a more likely assignment\n');
  %  keyboard
  %end
  
  [bestRelPeakIndex,bestRelResidueIndex] = find(V == maxCombinedProbability);
  
  if (length(bestRelPeakIndex) > 1)
    
    bestRelPeakIndex    = bestRelPeakIndex(1);
    bestRelResidueIndex = bestRelResidueIndex(1);
    
  
    
  end
  
  numCandidates                              = numCandidates + 1;
  relCandidatePeakIndices(numCandidates)     = bestRelPeakIndex;
  relCandidateResidueIndices (numCandidates) = bestRelResidueIndex;
  scores(numCandidates)                      = 0;

  for (i = 1:numVoters)
    voter_i_vote = voter{i}(bestRelPeakIndex, bestRelResidueIndex);
    if (voter_i_vote == 0)
      numCandidates = numCandidates - 1;
      stopGeneratingCandidates = 1;
      break;
    end
    scores(numCandidates) = scores(numCandidates) - ...
	log(voter_i_vote);
  end
  
% $$$   if (numVoters == 7)
% $$$     for i = 6:7
% $$$       voter_i_vote = voter{i}(bestRelPeakIndex, bestRelResidueIndex);
% $$$       if (voter_i_vote == 0)
% $$$ 	numCandidates = numCandidates - 1;
% $$$ 	stopGeneratingCandidates = 1;
% $$$ 	break;
% $$$       end
% $$$       scores(numCandidates) = scores(numCandidates) - 	log(voter_i_vote); 
% $$$     end
% $$$   end
  
  fprintf(1, 'assignment of %d to %d has %f score.\n',peakIndices(bestRelPeakIndex),residueIndices(bestRelResidueIndex),scores(numCandidates));
  
  V(bestRelPeakIndex, bestRelResidueIndex) = 0;
end

[sortedScores, newIndices]       = sort(scores(1:numCandidates));
sortedRelCandidatePeakIndices    = zeros(numCandidates, 1);
sortedRelCandidateResidueIndices = zeros(numCandidates, 1);

for i = 1:numCandidates
  sortedRelCandidatePeakIndices(i) = ...
      relCandidatePeakIndices(newIndices(i));
  sortedRelCandidateResidueIndices(i) = relCandidateResidueIndices(newIndices(i)); 
end

relCandidatePeakIndices    = sortedRelCandidatePeakIndices;
relCandidateResidueIndices = sortedRelCandidateResidueIndices;
scores                     = sortedScores;
% $$$ if (numCandidates == 1)
% $$$   fprintf(1, 'numCandidates == 1\n');
% $$$   keyboard
% $$$ end

function totalVotes = vote5(CP,SXCP,SSCP,TP,HDE)
totalVotes =CP*0;
for(i=1:5)
   m = nchoosek(1:5,i);
   for(j = 1:size(m,1))
      bipartiteGraphToVoteOn = CP*0+1;
      c = m(j,:);
      for(k=1:length(c))
         if(c(k)==1)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*CP;
         elseif(c(k)==2)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*SXCP;
         elseif(c(k)==3)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*SSCP;
         elseif(c(k)==4)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*TP;
         elseif(c(k)==5)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*HDE;
         end
      end
      votes = runHungarianAlgorithm(bipartiteGraphToVoteOn);
      votes = votes(1:size(totalVotes,1),1:size(totalVotes,2));
      totalVotes = totalVotes+votes;
   end
end

function totalVotes = vote(CP,SXCP,SSCP,TP,RP1,RP2,HDE)

totalVotes =CP*0; 
for(i=1:7)
   m = nchoosek(1:7,i);
   for(j = 1:size(m,1))
      bipartiteGraphToVoteOn = CP*0+1; 
      c = m(j,:);
      for(k=1:length(c))
         if(c(k)==1)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*CP;
         elseif(c(k)==2)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*SXCP;
         elseif(c(k)==3)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*SSCP;
         elseif(c(k)==4)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*TP;
         elseif(c(k)==5)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*RP1;
         elseif(c(k)==6)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*RP2;
         elseif(c(k)==7)
            bipartiteGraphToVoteOn = bipartiteGraphToVoteOn.*HDE;
         end
      end
      votes = runHungarianAlgorithm(bipartiteGraphToVoteOn);
      if (~isempty(votes))
	votes          = votes(1:size(totalVotes,1),1:size(totalVotes,2));
	totalVotes = totalVotes + votes;
      end
   end
end

% $$$ function [ASSIGNTABLE] = makeASSIGNTABLE_SquareAndNormalized(ASSIGNTABLE)
% $$$ 
% $$$ F=ones(max(size(ASSIGNTABLE)));
% $$$ F(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=ASSIGNTABLE;
% $$$ ASSIGNTABLE=F;   
% $$$ for(i=1:size(ASSIGNTABLE,1))
% $$$   ASSIGNTABLE(i,:)=ASSIGNTABLE(i,:)/sum(ASSIGNTABLE(i,:));
% $$$ end