function  [relCandidatePeakIndices, relCandidateResidueIndices, ...
	   scores, numCandidates] =  voteAccordingToMB_Score(MB_ShiftScore, ...
						  MB_RdcScore, ...
						  csCoefficient, ...
						  rdcCoefficient, ...
						  noeCoefficient, ...
						  ASSIGNTABLE, ...
						  peakIndices,residueIndices,MASTER,HSQCDATA,ALLDISTS,NTH,numCS,numRDC,numHN_NOES);

%bestRelPeakIndex    = -1;
%bestRelResidueIndex = -1;

MAX_NUMCANDIDATES          = 10;
SCORE_DELTA_THRESHOLD      = 0.01; %also defined in a_star_assign.m
relCandidatePeakIndices    = zeros(MAX_NUMCANDIDATES,1); 
relCandidateResidueIndices = zeros(MAX_NUMCANDIDATES,1); 
scores                     = zeros(MAX_NUMCANDIDATES,1); 
numCandidates              = 0; 



assert (~isempty(ASSIGNTABLE));

V         = csCoefficient * MB_ShiftScore + rdcCoefficient * MB_RdcScore;



%V         = V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
V         = V(peakIndices, residueIndices);

origV     = V; 

assert (length(peakIndices)    == size(ASSIGNTABLE,1));
assert (length(residueIndices) == size(ASSIGNTABLE,2));

%maxVotesBeforeMaskingWithASSIGNTABLE = max(max(V));

V                           = V.*ASSIGNTABLE;

% $$$ CS_RDC_COMPONENTS           = ASSIGNTABLE * 0 ;

noeScores                   = ASSIGNTABLE * 0 ;

stopGeneratingCandidates    = 0;

assert (length(peakIndices) == size(ASSIGNTABLE,1));
assert (length(residueIndices) == size(ASSIGNTABLE,2));

% $$$ reducedV                                          = cell(1,1);

% $$$ for relPeakIndex = 1:size(ASSIGNTABLE,1)
% $$$ 
% $$$   fprintf(1, 'relPeakIndex = %f\n', relPeakIndex);
% $$$   remainingRelPeakIndices        = [1:relPeakIndex-1,relPeakIndex+1:size(ASSIGNTABLE,1)];
% $$$   newRowReducedAssigntable       = ASSIGNTABLE(remainingRelPeakIndices,:);
% $$$   newRowReducedV                 = origV      (remainingRelPeakIndices,:);
% $$$   
% $$$   
% $$$   for relResidueIndex = 1:size(ASSIGNTABLE,2)
% $$$     if (ASSIGNTABLE(relPeakIndex,relResidueIndex) == 0) 
% $$$       continue;
% $$$     end
% $$$ 
% $$$ %    reducedV{1}                                       = origV;
% $$$ %    newV                                          = V;
% $$$ % $$$     newMASTER                                         = MASTER;
% $$$ % $$$     newMASTER(peakIndices(relPeakIndex),residueIndices(relResidueIndex)) ...
% $$$ % $$$ 	= 1;
% $$$ 
% $$$     remainingRelResidueIndices     = [1:relResidueIndex-1,relResidueIndex+1:size(ASSIGNTABLE,2)];
% $$$     
% $$$ %    newAssigntable                 = ASSIGNTABLE(remainingRelPeakIndices,remainingRelResidueIndices);
% $$$ %    reducedV{1}             	   = reducedV{1}
% $$$ %    (remainingRelPeakIndices,remainingRelResidueIndices);
% $$$ 
% $$$     newAssigntable                 = newRowReducedAssigntable(:,remainingRelResidueIndices);
% $$$     reducedV{1}             	   = newRowReducedV          (:,remainingRelResidueIndices);
% $$$ 
% $$$ 
% $$$ % $$$     newUnassignedPeakIndices       = ...
% $$$ % $$$ 	peakIndices(remainingRelPeakIndices);
% $$$ % $$$     newUnassignedResidueIndices    = residueIndices(remainingRelResidueIndices);
% $$$ % $$$     [newH_Value, CS_RDC_Component_Of_H_Value] = computeMB_H_Value(newAssigntable, newMASTER, HSQCDATA,...
% $$$ % $$$ 					 newUnassignedPeakIndices, newUnassignedResidueIndices,...
% $$$ % $$$ 					 ALLDISTS, NTH, ...
% $$$ % $$$ 					 MB_ShiftScore, ...
% $$$ % $$$ 					 MB_RdcScore, numCS, numRDC, ...
% $$$ % $$$ 					 csCoefficient, rdcCoefficient, ...
% $$$ % $$$ 					 noeCoefficient, numHN_NOES);
% $$$     betterH_Value                  = computeBetterH_Value(newAssigntable, ...
% $$$ 						  reducedV,1);
% $$$ % $$$     
% $$$ % $$$     debugH_Value                   = computeH_Value(newAssigntable, ...
% $$$ % $$$ 						    reducedV,1);
% $$$ % $$$     
% $$$ % $$$     assert (debugH_Value  == CS_RDC_Component_Of_H_Value);
% $$$ % $$$     
% $$$ % $$$     assert (betterH_Value >= debugH_Value);
% $$$     
% $$$     
% $$$ % $$$     CS_RDC_Component_Of_H_Value                     = CS_RDC_Component_Of_H_Value * -1;
% $$$ % $$$     CS_RDC_COMPONENTS(relPeakIndex,relResidueIndex) = CS_RDC_Component_Of_H_Value;
% $$$     %fprintf(1, 'the assignment of %d to %d has %f score and %f h_value.\n',peakIndices(relPeakIndex),residueIndices(relResidueIndex),V(relPeakIndex,relResidueIndex),CS_RDC_Component_Of_H_Value);
% $$$     V(relPeakIndex,relResidueIndex)                 = V(relPeakIndex,relResidueIndex)-betterH_Value;
% $$$   
% $$$   end
% $$$ end

if (noeCoefficient ~= 0)

  for relPeakIndex = 1:size(ASSIGNTABLE,1)
    for relResidueIndex = 1:size(ASSIGNTABLE,2)
      if (ASSIGNTABLE(relPeakIndex,relResidueIndex) == 0) 
	continue;
      end
      
      newMASTER                                     = MASTER;
      newMASTER (peakIndices,residueIndices)        = ASSIGNTABLE;
      newMASTER (peakIndices(relPeakIndex),:)       = 0;
      newMASTER (:,residueIndices(relResidueIndex)) = 0;
      newMASTER(peakIndices(relPeakIndex),residueIndices(relResidueIndex)) ...
	  = 1;
      
      [noeScores(relPeakIndex, relResidueIndex),numHN_NOES] = computeMB_NoeScore(HSQCDATA, ...
						  ALLDISTS, NTH, newMASTER);
      fprintf(1, 'assigning %d to %d has an NOE score of %f.\n',peakIndices(relPeakIndex),residueIndices(relResidueIndex),noeScores(relPeakIndex,relResidueIndex));
      fprintf(1, 'the corresponding shift and rdc score is %f\n', V(relPeakIndex,relResidueIndex));
							   %    keyboard
    end
  end
end
V = (V + noeCoefficient*noeScores)/(csCoefficient*numCS +  rdcCoefficient*numRDC + noeCoefficient*numHN_NOES);
%keyboard

V = V .* ASSIGNTABLE;

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

firstScore = max(max(V));

while ((numCandidates < MAX_NUMCANDIDATES) & (~stopGeneratingCandidates))
%while (~stopGeneratingCandidates)

  maxScore  = max(max(V));
  
  if (maxScore == 0)
    break;
  end
  
  %if (maxVotes ~= maxVotesBeforeMaskingWithASSIGNTABLE)
  %  fprintf(1, 'ASSIGNTABLE masked a more likely assignment\n');
  %  keyboard
  %end
  
  [bestRelPeakIndex,bestRelResidueIndex] = find(V == maxScore);
  

  
  if (length(bestRelPeakIndex) > 1)
    
    bestRelPeakIndex    = bestRelPeakIndex(1);
    bestRelResidueIndex = bestRelResidueIndex(1);
    
  end

  if (ASSIGNTABLE(bestRelPeakIndex,bestRelResidueIndex)==0)
    stopGeneratingCandidates = 1;
    break;
  end

  if (abs(firstScore - maxScore) > SCORE_DELTA_THRESHOLD)
    stopGeneratingCandidates = 1;
    break;
  end
  
  
  
  numCandidates                              = numCandidates + 1;
  relCandidatePeakIndices(numCandidates)     = bestRelPeakIndex;
  relCandidateResidueIndices (numCandidates) = bestRelResidueIndex;
  scores(numCandidates)                      = -maxScore;

  fprintf(1, 'assignment of %d to %d has %f MB score.\n',peakIndices(bestRelPeakIndex),residueIndices(bestRelResidueIndex),scores(numCandidates));
% $$$   fprintf(1, 'the corresponding CS_RDC_COMPONENT is %f\n', CS_RDC_COMPONENTS(bestRelPeakIndex,bestRelResidueIndex));
  %  fprintf(1, 'MB_ShiftScore for that entry = %f, MB_RDC_Score for the',MB_ShiftScore(peakIndices(bestRelPeakIndex),residueIndices(bestRelResidueIndex)));
%  fprintf(1, ' entry is %f\n', MB_RdcScore(peakIndices(bestRelPeakIndex),residueIndices(bestRelResidueIndex)));
  %  keyboard
  
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