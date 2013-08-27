function [HD_SCORE,assignmentAccuracy] = betterHD(peaks,rdcs,HDEXCHANGE, ...
						  peakIDs, NOES, ...
						  VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, SHIFTS_Filename,SHIFTX_Filename, useMBM_EM)

%HD: This program computes how well a given model fits  
%      Input:  peaks = a Nx2 matrix containing N hsqc peaks. column 1 is the H shift, column 2 is the N shift
%					rdcs = a Nx2 matrix containing N rdcs in 2 media medium, in an arbitrary order.
%					HDEXCHANGE= a Nx1 matrix containing N boolean values indicating whether that peak is a slow-exchanging peak. 1 = slow, 0 = fast
%					peakIDs= a Nx1 matrix containing N ids for the peak. eg peakIDs(1) = the id for peaks(1)
%								NOTE: the rows of peaks, rdcs, HDEXCHANGE, and NOES must have the same order as the
%								rows of peakIDs. 
%					NOES = a NxN matrix. If element NOES(i,j) = 1, then peakIDs(i) is involved in a Dnn with peakIDs(j)
%					VECTORS = A Mx3 matrix containing the normalized backbone amide bond vectors
%					TYPES = A Mx1 matrix containing the 3 letter amino acid code for the backbone amide bond vectors
%					RESNUMS = A Mx1 matrix containing the resiude. Note the rows in the parameters named 
%								 VECTORS, TYPES, SSTRUCT, HBOND, and ALLDISTS should all have the same order
% 								 as listed in the rows of RESNUMS
%					SSTRUCT = A Mx1 matrix containing secondary structure type each residue in the model. C=coil, B=beta, H=helix
%					HBOND = A Mx1 matrix containing the distances between all pairs of backbone amide protons, Y=yes, its involved in an H bond (or is not solvent accesible), N=means that it is essntially labile. 
%					ALLDISTS = A MxM matrix containing the distances between all pairs of backbone amide protons
%					IALLDISTS = A MxM matrix containing an altered version of the distances, that effectively ignores residues
%									involed in random coil
%		 Output:   


dbstop if error
dbstop if warning

[ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, RP1,...
 RP2, S1, S2, RDC1, RDC2, NTH] = initialize(peaks,rdcs,HDEXCHANGE, ...
					    peakIDs, NOES, ...
					    VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, SHIFTS_Filename, SHIFTX_Filename);


[totalNumPeaks,totalNumResidues] = size(MASTER);
numAssignedPeaks                 = 0;
numUnassignedPeaks               = totalNumPeaks;

numVoters           = 5;

%change in variable names:
peakIndices    = ROWIN;
residueIndices = COLIN;
voter          = cell(7,1);
voter          = initializeVoters(voter,SSCP,SXCP,CP,TP,HDE,RP1,RP2);

analyzeVoters(voter, numVoters, peakIndices, residueIndices);

%parameters of the program...

%useMBM_EM = 1
debugIncorrectAssignment = 0

while (numUnassignedPeaks > 0)

  %first updating voters

  if (numAssignedPeaks >= 5)
    [voter,numVoters,S1,S2] = ...
	updateVotersAndAlignmentTensors(voter, ASSIGNTABLE,...
					MASTER, S1, S2, RDC1,RDC2, VECTORS);
  end

  if (useMBM_EM)
    [bestRelPeakIndex,bestRelResidueIndex] = ...
	voteAccordingToMBM_EM    (voter,     numVoters,ASSIGNTABLE,NOES,ALLDISTS,NTH,peakIndices,residueIndices);
  else
    [bestRelPeakIndex,bestRelResidueIndex] = ...
	voteAccordingToSimpleSum (numVoters, voter);
  end
  
  bestPeakIndex    = peakIndices(bestRelPeakIndex);
  bestResidueIndex = residueIndices(bestRelResidueIndex);

  fprintf(1,'votes decided to assign %d to %d\n', bestPeakIndex, ...
	  bestResidueIndex);
[numUnassignedPeaks, numUnassignedResidues] = size(voter{1});  
   votes = zeros(numUnassignedResidues, 1);
%  scores = votes;
  relPeakIndex = bestRelPeakIndex;
  
  for voterIndex = 1:numVoters
    weightedRankings                            = voter{voterIndex}(relPeakIndex,:);
    votes = votes + weightedRankings';
  end	
  figure; plot(votes', '*'); length(nonzeros(votes))
  [i,j] = max(votes)
  bestRelResidueIndex
  keyboard
  
  if (debugIncorrectAssignment)

    debugIncorrectAssignmentFunction(voter, numVoters, ...
				     bestRelPeakIndex, ...
				     bestRelResidueIndex, bestPeakIndex,bestResidueIndex,residueIndices);
  end

  
  ASSIGNTABLE             = updateAssignTable (ASSIGNTABLE,bestRelPeakIndex,bestRelResidueIndex);
  
  numPotentialAssignments = sum(sum(ASSIGNTABLE));

  for i = 1:100
    %uncomment the following for NVR type NOE pruning
    %     [ASSIGNTABLE]        = noePrune(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,peakIndices,residueIndices);
    
    [ASSIGNTABLE]        = noePrune(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, peakIndices,residueIndices);
    [ASSIGNTABLE]        = noePrune(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 40, peakIndices,residueIndices);

    [MASTER,ASSIGNTABLE] = doAssignments(MASTER,ASSIGNTABLE,peakIndices, ...
					 residueIndices);
    [ASSIGNTABLE,voter, peakIndices,residueIndices] ...
	= all_removeAssignedResiduesAndUpdateVoteBPG(ASSIGNTABLE,voter,numVoters,peakIndices,residueIndices);
  

    newNumPotentialAssignments = sum(sum(ASSIGNTABLE));
    
    if(newNumPotentialAssignments == numPotentialAssignments)
      break;
    end

    numPotentialAssignments = newNumPotentialAssignments;
    
  end

  numAssignedPeaks      = sum(sum(MASTER));
  numUnassignedPeaks    = totalNumPeaks    - numAssignedPeaks;

  
end


computeCorrectness(MASTER,RESNUMS, peakIDs);




%keyboard


function   ASSIGNTABLE = updateAssignTable(ASSIGNTABLE,bestRelPeakIndex,bestRelResidueIndex)

  
ASSIGNTABLE  (bestRelPeakIndex,                   :) = 0;
ASSIGNTABLE  (:               , bestRelResidueIndex) = 0;
ASSIGNTABLE  (bestRelPeakIndex, bestRelResidueIndex) = 1;

%what to do if there is a peak left without a potential
%assignment in this case? restore its assignments to all residues
%uniformly could be an option.

ASSIGNTABLE                           = reinitializeAndNormalizeBPG_Entries_If_Necessary(ASSIGNTABLE);

ASSIGNTABLE                           = and (ASSIGNTABLE,ASSIGNTABLE);


%updates ASSIGNTABLE based on NOEs and currently assigned
%peak-residue pairs.
function [ASSIGNTABLE] = noePrune(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,peakIndices,residueIndices)
%prune via NOEs
%NP = NVR_NOE2PROB(ASSIGNTABLE,NOES,ALLDISTS,NTH,peakIndices,residueIndices);
mypeakIndices                      = 1:size(MASTER,1);
myresidueIndices                   = 1:size(MASTER,2);

assert (sum(sum(MASTER(peakIndices, residueIndices))) == 0);
MASTER(peakIndices,residueIndices) = ASSIGNTABLE;
NP                                 = NVR_NOE2PROB (MASTER,NOES,ALLDISTS,NTH, mypeakIndices,myresidueIndices);
NP                                 = NP           (peakIndices,residueIndices);
ASSIGNTABLE                        = and          (NP,NP);


%this function assigns those peaks which have only one residue
%assignment. It also then removes the potential assignability of
%that residue to other peaks. 

function [MASTER,ASSIGNTABLE]    = doAssignments(MASTER,ASSIGNTABLE,peakIndices,residueIndices);
%first, propagate contraints

madeAnAssignment                           = 1;

[numUnassignedPeaks,numUnassignedResidues] = size(ASSIGNTABLE);

peakIsUnassignableAtPresent                = ones(1,numUnassignedPeaks); %BECOMES 0 if the i'th peak has a
									 %single possible residue. All
									 %the other possible peak
									 %assignments to that residue are removed.
residueIsUnassigned                        = ones(1,numUnassignedResidues);
ASSIGNTABLE_IN                             = ASSIGNTABLE;

while (madeAnAssignment == 1)

  madeAnAssignment = 0;
  
  for (relPeakIndex = 1:numUnassignedPeaks)
  
    relResidueIndexOfPossibleResiduesForPeak_relPeakIndex = find(ASSIGNTABLE(relPeakIndex,:));
    
    if  (length(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)==1) 

      ASSIGNTABLE  (:            , relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)=0;
      
      %will have to change the above line, not doing the
      %assignment for these cases.SB

      ASSIGNTABLE  (relPeakIndex , relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)=1;

      if (peakIsUnassignableAtPresent(relPeakIndex) == 1)
	
	peakIsUnassignableAtPresent(relPeakIndex)  = 0;
	
	residueIsUnassigned(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex) = 0;
	
	madeAnAssignment                           = 1;
      
      end
      
    end
  
  end

  for(relPeakIndex=1:numUnassignedResidues)

    x=find(ASSIGNTABLE(:,relPeakIndex));

    if ((length(x)==1) & (numUnassignedPeaks==numUnassignedResidues))

      fprintf(1, 'I dont think the code reaches this point\n');
      ASSIGNTABLE(x,:)=0;
      ASSIGNTABLE(x,relPeakIndex)=1;
      if(residueIsUnassigned(relPeakIndex)==1)
	residueIsUnassigned(relPeakIndex)=0;
	madeAnAssignment=1;
      end
    
    end
  
  end

end


%the loop below restores the values of zero-rows in ASSIGNTABLE.
for (relPeakIndex=1:numUnassignedPeaks)

  if  (sum(ASSIGNTABLE(relPeakIndex,:))==0)

    ASSIGNTABLE(relPeakIndex,:)                           = ASSIGNTABLE_IN(relPeakIndex,:);

    relResidueIndexOfPossibleResiduesForPeak_relPeakIndex = find(ASSIGNTABLE(relPeakIndex,:));

    assert (~isempty(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)); %I hope this restoration helped to
									      %restore the values of ASSIGNTABLE.
% $$$     if (relResidueIndexOfPossibleResiduesForPeak_relPeakIndex <= numUnassignedPeaks) %I don't see why this is necessary.
% $$$       fprintf(1, 'code that I dont understand is reached\n');
% $$$       keyboard
% $$$       ASSIGNTABLE(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex,:)=ASSIGNTABLE_IN(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex,:);
% $$$     end
    
  end
  
end


%next, make assignments
for (relPeakIndex=1:numUnassignedPeaks)

  relResidueIndexOfPossibleResiduesForPeak_relPeakIndex = find(ASSIGNTABLE(relPeakIndex,:)); %what if there is another peak whose
											     %sole assignment is x?
   if (length(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)==1)

     MASTER(peakIndices(relPeakIndex),residueIndices(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex))=1;

%     fprintf(1, 'doAssignments: assigning %d to %d\n', peakIndices(relPeakIndex), residueIndices(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex));
   
   end

end


%updates peakIndices, residueIndices and d all 7 bpgs according to
%which rows of ASSIGNTABLE contain only one non-zero entry. (These
%correspond to the assignments made in doAssignments). Then also
%reduces the rows and columns of ASSIGNTABLE. Binarizes back
%ASSIGNTABLE. 

function [ASSIGNTABLE,voter, peakIndices, residueIndices]=     all_removeAssignedResiduesAndUpdateVoteBPG(ASSIGNTABLE, voter, ...
						  numVoters, peakIndices,residueIndices);


[unassignedRelPeakIndices, unassignedRelResidueIndices] = determineUnassignedPeaksAndResidues(ASSIGNTABLE);

for i = 1:numVoters
  voter{i}   = individualBPG_removeAssignedResiduesAndUpdate ...
      (ASSIGNTABLE,voter{i},unassignedRelPeakIndices, unassignedRelResidueIndices);
end


[ASSIGNTABLE] = individualBPG_removeAssignedResiduesAndUpdate ...
    (ASSIGNTABLE,ASSIGNTABLE,unassignedRelPeakIndices, unassignedRelResidueIndices);

ASSIGNTABLE               = and    (ASSIGNTABLE,ASSIGNTABLE);


peakIndices    =  peakIndices    (unassignedRelPeakIndices);
residueIndices =  residueIndices (unassignedRelResidueIndices);


%removes from consideration those peaks and residues (in
%BPG_TABLE,peakIndices and residueIndices) which are
%unambiguously assigned by ASSIGNTABLE. multiplies BPG_TO_UPDATE
%with ASSIGNTABLE. Keeps the portion of BPG_TO_UPDATE which
%contains more than one non-zero candidate for each peak. Also
%updates the relativePeakIndices and relativeResidueIndices as
%stored in peakIndices and residueIndices accordingly. Normalizes the rows of BPG_TO_UPDATE.
function [BPG_TO_UPDATE]=individualBPG_removeAssignedResiduesAndUpdate(ASSIGNTABLE,BPG_TO_UPDATE, ...
						  unassignedRelPeakIndices, unassignedRelResidueIndices)

if (isempty(BPG_TO_UPDATE))
  return
end

%ASSIGNTABLE_IN = ASSIGNTABLE;
%BPG_TO_UPDATE_IN   = BPG_TO_UPDATE;
BPG_TO_UPDATE      = BPG_TO_UPDATE.*ASSIGNTABLE;

%BPG_TO_UPDATE = reinitializeAndNormalizeBPG_Entries_If_Necessary(BPG_TO_UPDATE);%, BPG_TO_UPDATE_IN);

BPG_TO_UPDATE      = BPG_TO_UPDATE  (unassignedRelPeakIndices,unassignedRelResidueIndices);

BPG_TO_UPDATE      = reinitializeAndNormalizeBPG_Entries_If_Necessary(BPG_TO_UPDATE);


function [BPG_TO_UPDATE] = reinitializeAndNormalizeBPG_Entries_If_Necessary(BPG_TO_UPDATE)

for(i=1:size(BPG_TO_UPDATE,1))
% $$$    if(sum(BPG_TO_UPDATE(i,:))==0)
% $$$      try_erroneous_code = 0;
  
% $$$      if (try_erroneous_code)
% $$$      
% $$$        if(i<=size(BPG_TO_UPDATE_IN,1) & size(BPG_TO_UPDATE,2)==size(BPG_TO_UPDATE_IN,2))
% $$$ 	 numPeaks_in = length(BPG_TO_UPDATE_IN(:,1));
% $$$ 	 numPeaksRemainingToAssign = length(BPG_TO_UPDATE(:,1));
% $$$ 	 
% $$$ 	 if (numPeaks_in ~= numPeaksRemainingToAssign)
% $$$ 	   if (unassignedRelPeakIndices(i) ~= i)
% $$$ 	     fprintf(1, 'error. we want to update relative peak %d',i);
% $$$ 	     fprintf(1, 'using the information of peak %d',unassignedRelPeakIndices(i));
% $$$ 	     fprintf(1, 'in the input matrix\n');
% $$$ 	     keyboard
% $$$ 	   end
% $$$ 	 end
% $$$        end
% $$$ 	
% $$$ 	BPG_TO_UPDATE(i,:)=BPG_TO_UPDATE_IN(i,:); %is this
% $$$         %correct? the            %rows of
% $$$                                                    %BPG_TO_UPDATE
% $$$                                                    %no longer
% $$$                                                    %correspond to
% $$$                                                    %the rows of
% $$$                                                    %BPG_TO_UPDATE_IN IMO.
% $$$      end
     
  if (sum(BPG_TO_UPDATE(i,:))==0)

    fprintf(1, 'restoring the entries for bpg for peak relIndex %d\n',i);
%   keyboard
    BPG_TO_UPDATE(i,:)=1;
  
  end
% $$$    end

%  BPG_TO_UPDATE(i,:)=BPG_TO_UPDATE(i,:) / sum(BPG_TO_UPDATE(i,:));

end


% $$$ function [BPG_TO_UPDATE] = reinitializeAndNormalizeBPG_Entries_If_Necessary(BPG_TO_UPDATE, BPG_TO_UPDATE_IN)
% $$$ 
% $$$ for(i=1:size(BPG_TO_UPDATE,1))
% $$$   if(sum(BPG_TO_UPDATE(i,:))==0)
% $$$     %    if(i<=size(BPG_TO_UPDATE_IN,1))
% $$$     BPG_TO_UPDATE(i,:)=BPG_TO_UPDATE_IN(i,:);
% $$$     %    end
% $$$     if(sum(BPG_TO_UPDATE(i,:))==0)
% $$$       BPG_TO_UPDATE(i,:)=1;
% $$$     end
% $$$   end
% $$$   BPG_TO_UPDATE(i,:)=BPG_TO_UPDATE(i,:)/sum(BPG_TO_UPDATE(i,:));
% $$$ end


function [unassignedRelPeakIndices, unassignedRelResidueIndices] = determineUnassignedPeaksAndResidues(ASSIGNTABLE)

unassignedRelPeakIndices     = [1:size(ASSIGNTABLE,1)];
unassignedRelResidueIndices  = [1:size(ASSIGNTABLE,2)];


if (isempty(unassignedRelPeakIndices))
  return
end


for (relPeakIndex=1:size(ASSIGNTABLE,1))
  
   relResidueIndexOfPossibleResiduesForPeak_relPeakIndex  = find(ASSIGNTABLE(relPeakIndex,:));
   
   if (length(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)==1)
     
      unassignedRelPeakIndices(relPeakIndex)                    =0;
      unassignedRelResidueIndices(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)=0;   
      
   end
   
end

unassignedRelPeakIndices   =find(unassignedRelPeakIndices);
unassignedRelResidueIndices=find(unassignedRelResidueIndices);




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


%##################################################################################
%##  updateTen
function [SAUPE, flag]= updateTen(TABLE,ET,RDC,VECS)

SAUPE  = ET;

M = zeros(1,5);
vec = 0;
ct=1;
for(i=1:size(TABLE,1))
   x = find(TABLE(i,:));
   if(length(x)==1 & RDC(i)>-999)
      M(ct,1) = VECS(x,1).^2-VECS(x,3).^2;
      M(ct,2) = 2*(VECS(x,1)*VECS(x,2));
      M(ct,3) = 2*(VECS(x,1)*VECS(x,3));
      M(ct,4) = VECS(x,2).^2-VECS(x,3).^2;
      M(ct,5) = 2*(VECS(x,2)*VECS(x,3));
      vec(ct) =  RDC(i); 
      ct=ct+1;
   end
end
flag = 0;
if(size(M,1)<5)
   size(M,1);
   flag=1;
   return
end
%invert the matrix
if(size(M,1)==5 & size(M,2)==5)
   M = inv(M);   
else
   M = pinv(M);
end
%do the least squares fitting
s = M*vec';
SAUPE = zeros(3,3);
SAUPE(1,1) = s(1);
SAUPE(2,2) = s(4);
SAUPE(3,3) = -1*(s(1)+s(4));
SAUPE(1,2) = s(2);
SAUPE(2,1) = s(2);
SAUPE(1,3) = s(3);
SAUPE(3,1) = s(3);
SAUPE(2,3) = s(5);
SAUPE(3,2) = s(5);




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
  

  
function [bestRelPeakIndex,bestRelResidueIndex] = ...
    voteAccordingToMBM_EM(voter, numVoters, ASSIGNTABLE,NOES,ALLDISTS,NTH,peakIndices,residueIndices)

bestRelPeakIndex    = -1;
bestRelResidueIndex = -1;

assert (~isempty(ASSIGNTABLE));

if (numVoters == 5)
  V       = vote5(voter{1},voter{2},voter{3},voter{4},voter{5});%CP,SXCP,SSCP,TP,HDE,
  %ASSIGNTABLE);
else
  V       = vote (voter{1},voter{2},voter{3},voter{4},voter{5},voter{6},voter{7});%,...%CP,SXCP,SSCP,TP,RP1,RP2,HDE,
%	   ASSIGNTABLE);
end

V         = V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));

%maxVotesBeforeMaskingWithASSIGNTABLE = max(max(V));

V         = V.*ASSIGNTABLE;

maxVotes  = max(max(V));

%if (maxVotes ~= maxVotesBeforeMaskingWithASSIGNTABLE)
%  fprintf(1, 'ASSIGNTABLE masked a more likely assignment\n');
%  keyboard
%end

[bestRelPeakIndex,bestRelResidueIndex] = find(V == maxVotes);

if (length(bestRelPeakIndex) > 1)
  
  bestRelPeakIndex    = bestRelPeakIndex(1);
  bestRelResidueIndex = bestRelResidueIndex(1);

end

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


function voter = initializeVoters(voter,SSCP,SXCP,CP,TP,HDE,RP1,RP2)

voter        = cell(7,1);
voter{1}  = SSCP;
voter{2}  = SXCP;
voter{3}  = CP;
voter{4}  = TP;
voter{5}  = HDE;
voter{6}  = RP1;
voter{7}  = RP2;

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


function [voter,numVoters,S1,S2] = ...
    updateVotersAndAlignmentTensors(voter, ASSIGNTABLE,...
				    MASTER, S1, S2, RDC1,RDC2,...
				    VECTORS)
numVoters = 7;

S1  = updateTen(MASTER,S1,RDC1,VECTORS);
S2  = updateTen(MASTER,S2,RDC2,VECTORS);


voter{6} = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
voter{7} = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);


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

%assignmentType is 1 for making assignments. it is 2 to retain the
%top top candidates for each each row.
function assignmentsBipartiteGraph = runHungarianAlgorithm(bipartiteGraphToVoteOn)
% $$$ if(size(ASSIGNTABLE,1)-size(ASSIGNTABLE,2)~=0)
% $$$ 
% $$$   ASSIGNTABLE = makeASSIGNTABLE_SquareAndNormalized(ASSIGNTABLE);
% $$$ 
% $$$ end
% $$$ 
% $$$ ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);

if(size(bipartiteGraphToVoteOn,1)-size(bipartiteGraphToVoteOn,2)~=0)

  squareBipartiteGraphToVoteOn =   makeBipartiteGraphToVoteOnSquare(bipartiteGraphToVoteOn);
  squareBipartiteGraphToVoteOn =   reinitializeAndNormalizeBPG_Entries_If_Necessary(squareBipartiteGraphToVoteOn);

else
  squareBipartiteGraphToVoteOn =   bipartiteGraphToVoteOn;
end

assignmentsBipartiteGraph      =   squareBipartiteGraphToVoteOn*0;


assignmentVector = hungarian(squareBipartiteGraphToVoteOn*-1);

% $$$ elseif(assignmentType==2)
% $$$ 
% $$$   fprintf(1, 'currently, i dont think code is ever executed\n');
% $$$   keyboard
% $$$   
% $$$   sz=top;
% $$$    for(i=1:size(HM,1))
% $$$       [x ind]= sort(HM(i,:));
% $$$       in2=find(x);
% $$$       ind=ind(in2);x=x(in2);
% $$$       if(range(x)>0)
% $$$          if(length(x)<sz+1)
% $$$             assignments(i,ind)=1;
% $$$          else
% $$$             assignments(i,ind(length(ind)-sz:length(ind)))=1;
% $$$          end
% $$$       else
% $$$          assignments(i,ind)=1;
% $$$       end
% $$$       assignmentVector=[];   
% $$$    end
% $$$ else
% $$$   fprintf(1, 'I dont think this part of the code is reached\n');
% $$$   keyboard;
% $$$ %  assignmentVector=simp(HM);
% $$$ end
for(i=1:length(assignmentVector))
   assignmentsBipartiteGraph(i,assignmentVector(i)) = 1;
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


function [squareBipartiteGraphToVoteOn] = makeBipartiteGraphToVoteOnSquare (bipartiteGraphToVoteOn)

F=ones(max(size(bipartiteGraphToVoteOn)));
F(1:size(bipartiteGraphToVoteOn,1),1:size(bipartiteGraphToVoteOn,2))=bipartiteGraphToVoteOn;
squareBipartiteGraphToVoteOn=F;   


% $$$ 
% $$$ for(i=1:size(squareBipartiteGraphToVoteOn,1))
% $$$   if(sum(squareBipartiteGraphToVoteOn(i,:))==0)
% $$$     %i don't think the following line is reached but is ok to keep it.
% $$$     squareBipartiteGraphToVoteOn(i,:)=1;
% $$$   end
% $$$   squareBipartiteGraphToVoteOn(i,:)=squareBipartiteGraphToVoteOn(i,:)/sum(squareBipartiteGraphToVoteOn(i,:));
% $$$ end



function [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, RP1,...
	  RP2, S1, S2, RDC1, RDC2, NTH] = initialize(peaks,rdcs,HDEXCHANGE, ...
						  peakIDs, NOES, ...
						  VECTORS,TYPES, ...
						  RESNUMS,SSTRUCT, ...
						  HBOND, ALLDISTS,IALLDISTS, SHIFTS_Filename, SHIFTX_Filename);

HSHIFTS = peaks(:,1);
NSHIFTS = peaks(:,2);
RDC1    = rdcs(:,1);
RDC2    = rdcs(:,2);

%compute a tollerance for NOE distances
NTH=4.8;
mu=mean(mean(ALLDISTS));
if(mu-12.9>0)
   NTH=NTH+(mu-12.9); NTH=min(NTH,8);
end

ASSIGNTABLE = ones(length(HSHIFTS),size(VECTORS,1))/size(VECTORS,1);
OASSIGNTABLE=ASSIGNTABLE;
%these keep track of which peaks and residues are represented in the current
%matricies
ROWIN=1:size(ASSIGNTABLE,1);
COLIN=1:size(ASSIGNTABLE,2);
%This is the master assignment table
MASTER=ASSIGNTABLE*0;

fprintf('computing assignments...\n');
HDE = HD_HD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND);
ASSIGNTABLE = and(ASSIGNTABLE,HDE);

TP =   ASSIGNTABLE;
TP = NVR_TOCSY2PROB(peakIDs,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,NTH,ROWIN,COLIN);
ASSIGNTABLE = and(ASSIGNTABLE,TP );

[CP] = HD_CS2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,NTH,ROWIN,COLIN);
ASSIGNTABLE = and(ASSIGNTABLE,CP);

%prune via the program shiftx
[SXCP] = HD_SHIFTX2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			NOES,IALLDISTS,NTH,ROWIN,COLIN, SHIFTX_Filename);
ASSIGNTABLE = and(ASSIGNTABLE,SXCP);   

%prune via the program shifts
[SSCP] = HD_SHIFTS2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			NOES,IALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename, ...
			SHIFTX_Filename);
ASSIGNTABLE = and(ASSIGNTABLE,SSCP);

SSCP=SSCP.*ASSIGNTABLE;SXCP=SXCP.*ASSIGNTABLE;CP=CP.*ASSIGNTABLE;TP=TP.*ASSIGNTABLE;HDE=HDE.*ASSIGNTABLE;

%initialize the RDC variables
S1 = ones(3,3);S2 = S1;
S1 = updateTen(MASTER,S1,RDC1,VECTORS);
S2 = updateTen(MASTER,S2,RDC2,VECTORS);
RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);

function analyzeVoters(voter, numVoters, peakIndices, ...
		       residueIndices);
for i = 1:numVoters
  for peakIndex = 1:length(peakIndices)
    votes = voter{i}(peakIndex,:);
    [sortedVotes, sortedOrder] = sort(-votes)
    keyboard
  end
end
