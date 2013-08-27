function  [impossibleToAssign, MASTER, ASSIGNTABLE,  unassignedPeakIndices, unassignedResidueIndices, score, voter, numVoters] = ...
	noePruneAndDoMoreAssignments(ASSIGNTABLE, MASTER, ...
				     unassignedPeakIndices, ...
				     unassignedResidueIndices, ...
				     score, ...
				     voter, numVoters,...
				     NOES, IALLDISTS, ALLDISTS, NTH)%, ...
%				     RDC1,RDC2,VECTORS,MB_ShiftScore, ...
%				     MB_RDC_Score, numCS, numRDC, ...
%				     HSQCDATA, csCoefficient, ...
%				     rdcCoefficient, noeCoefficient,numHN_NOES);

%also
%updates numVoters.

%keyboard

numPotentialAssignments = sum(sum(ASSIGNTABLE));

for i = 1:100
  %uncomment the following for NVR type NOE pruning
  %     [ASSIGNTABLE]        = noePrune(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,peakIndices,residueIndices);
    
%  [ASSIGNTABLE, impossibleToAssign]        = noePrune(MASTER, ...
%						  ASSIGNTABLE,NOES,IALLDISTS, NTH, unassignedPeakIndices,unassignedResidueIndices);
  
%  if (impossibleToAssign)
%    fprintf(1, 'found impossibletoassign in noePrune\n');
%    return;
%  end
  
%  [ASSIGNTABLE, impossibleToAssign]        = noePrune(MASTER, ...
%
%						  ASSIGNTABLE,NOES,ALLDISTS, 40, unassignedPeakIndices,unassignedResidueIndices);

%keyboard

   [ASSIGNTABLE, impossibleToAssign]        = noePrune(MASTER, ASSIGNTABLE,NOES,ALLDISTS, NTH, unassignedPeakIndices,unassignedResidueIndices);
  
  if (impossibleToAssign)
    fprintf(1, 'found impossibletoassign in noePrune\n');
    return;
  end
  
  
%  fprintf(1, 'did noe pruning.\n');
  
  
  [impossibleToAssign, MASTER,ASSIGNTABLE, score] = doAssignments(MASTER,ASSIGNTABLE,unassignedPeakIndices, ...
						  unassignedResidueIndices, ...
						  voter, numVoters, ...
						  score);
%  update the score here and there? SB  
%  fprintf(1, 'called doAssignments\n');

  for i = 1:size(MASTER,2)
    if (sum(MASTER(:,i)) > 1)
      fprintf(1, 'warning. two or more peaks are assigned_to_the_same_residue\n');
      keyboard
    end
  end

  
  if (impossibleToAssign)
    fprintf(1, 'found impossibletoassign in doassignments\n');
    return;
  end

  numAssignedPeaks       = sum(sum(MASTER));

%  if (numAssignedPeaks == size(MASTER,1))
%    break;
%  end
  
%  fprintf(1,'keyboard command in noePruneAndDoMoreAssignments.m\n');
  
%  keyboard
  
  [impossibleToAssign, ASSIGNTABLE,voter, unassignedPeakIndices,unassignedResidueIndices] ...
      = all_removeAssignedResiduesAndUpdateVoteBPG(ASSIGNTABLE,voter,numVoters,unassignedPeakIndices,unassignedResidueIndices);
  
%  fprintf(1, 'called all_removeAssignedResiduesAndUpdateVoteBPG\n');

  if (isempty(ASSIGNTABLE))
%    fprintf(1, 'assigned all the available peaks. returning.\n');
%    keyboard
    break;
  end


  
  if (impossibleToAssign)
    fprintf(1, 'found impossibletoassign in all_removeassigned residues\n');
    return;
  end
  

  
  
% $$$   
% $$$   if (numAssignedPeaks >= 5)
% $$$     [voter,numVoters] = ...
% $$$       	updateVoters(voter, ASSIGNTABLE,...
% $$$ 		     MASTER, RDC1,RDC2, VECTORS);
% $$$   end
  
  newNumPotentialAssignments = sum(sum(ASSIGNTABLE));
  
  if(newNumPotentialAssignments == numPotentialAssignments)
    break;
  end
  
  numPotentialAssignments = newNumPotentialAssignments;
  
end

%keyboard
%score = computeMB_Score(MASTER, MB_ShiftScore,MB_RDC_Score, ...
%			numCS,numRDC, HSQCDATA,ALLDISTS,NTH, ...
%			csCoefficient, rdcCoefficient, noeCoefficient);

% $$$ numPeaks = size(MASTER,1);
% $$$ for peakIndex = 1:numPeaks
% $$$   residueIndex = find(MASTER(peakIndex,:));
% $$$   assert (length(residueIndex) <= 1);
% $$$   if (~isempty(residueIndex))
% $$$     score = score + csCoefficient*MB_ShiftScore(peakIndex,residueIndex) ...
% $$$ 	    + rdcCoefficient * MB_RDC_Score(peakIndex, residueIndex);
% $$$   end
% $$$ end
% $$$ score = score/(csCoefficient*numCS +  rdcCoefficient*numRDC + noeCoefficient*numHN_NOES);
% $$$ score = score * -1;






%this function assigns those peaks which have only one residue
%assignment. It also then removes the potential assignability of
%that residue to other peaks. 

function [impossibleToAssign, MASTER,ASSIGNTABLE, score]    = doAssignments(MASTER,ASSIGNTABLE,peakIndices,residueIndices, ...
		  voter, numVoters, score);


impossibleToAssign = 0;

madeAnAssignment   = 1;%initializing

[numUnassignedPeaks,numUnassignedResidues] = size(ASSIGNTABLE);


%first update ASSIGNTABLE
while (madeAnAssignment == 1)

  madeAnAssignment = 0;
  
  for (relPeakIndex = 1:numUnassignedPeaks)
  
    relResidueIndexOfPossibleResiduesForPeak_relPeakIndex = find(ASSIGNTABLE(relPeakIndex,:));
    
    if (isempty(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex))
      impossibleToAssign = 1;
      return;
    end
    
    
    if  (length(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)==1) 

      numPeakCandidates = length(find(ASSIGNTABLE(:, relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)));
      
      assert (numPeakCandidates >= 1);

      if (numPeakCandidates > 1)
      
	ASSIGNTABLE  (:            , relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)=0;
      
	ASSIGNTABLE  (relPeakIndex , relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)=1;
      
	madeAnAssignment                           = 1;
      end
      
    end
  
  end

  
%  fprintf(1, 'went through the peaks in doAssignments\n');
  
  for(relResidueIndex=1:numUnassignedResidues)

    x=find(ASSIGNTABLE(:,relResidueIndex));

    if ((length(x)==0) & (numUnassignedPeaks== ...
			  numUnassignedResidues))
      impossibleToAssign = 1;
      return;
    end
    
    if ((length(x)==1) & (numUnassignedPeaks==numUnassignedResidues))

      numResidueCandidates = length(find(ASSIGNTABLE(x,:)));
      
      assert (numResidueCandidates >= 1);
      
      if (numResidueCandidates > 1)
	ASSIGNTABLE(x,:)=0;
	ASSIGNTABLE(x,relResidueIndex)=1;
	
	madeAnAssignment=1;
      end
      
    end
  
  end
  
%  fprintf(1, 'went through the residues in doAssignments\n');

end

%fprintf(1, 'now updating master\n');

%next, make assignments (i.e. update MASTER)
for (relPeakIndex=1:numUnassignedPeaks)
%no need to go over the columns (residues) as well since if there
%is a unique assignment for a residue, the corresponding peak also
%has a unique assignment for that residue.
  relResidueIndexOfPossibleResiduesForPeak_relPeakIndex = find(ASSIGNTABLE(relPeakIndex,:)); %what if there is another peak whose
											     %sole assignment is x?
  if (length(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)==1)

    
%    ASSIGNTABLE  (:            , ...
%		  relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)=0; 
%    making sure that the residue is only assigned to one peak.
    
%    ASSIGNTABLE(relPeakIndex, ...
%		relResidueIndexOfPossibleResiduesForPeak_relPeakIndex)= 1;
    
    MASTER(peakIndices(relPeakIndex),residueIndices(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex))=1;

    for voterIndex = 1:numVoters
%      score = score - log(voter{voterIndex}(relPeakIndex, relResidueIndexOfPossibleResiduesForPeak_relPeakIndex));
       score = score - voter{voterIndex}(relPeakIndex, relResidueIndexOfPossibleResiduesForPeak_relPeakIndex);
    end
    
    fprintf(1, 'doAssignments: assigning %d to %d\n', peakIndices(relPeakIndex), residueIndices(relResidueIndexOfPossibleResiduesForPeak_relPeakIndex));
   
  end

end




%updates peakIndices, residueIndices and d all 7 bpgs according to
%which rows of ASSIGNTABLE contain only one non-zero entry. (These
%correspond to the assignments made in doAssignments). Then also
%reduces the rows and columns of ASSIGNTABLE. Binarizes back
%ASSIGNTABLE. 

function [impossibleToAssign, ASSIGNTABLE, voter, peakIndices, residueIndices]=     all_removeAssignedResiduesAndUpdateVoteBPG(ASSIGNTABLE, voter, ...
						  numVoters, peakIndices,residueIndices);


[unassignedRelPeakIndices, unassignedRelResidueIndices] = determineUnassignedPeaksAndResidues(ASSIGNTABLE);


%fprintf(1, 'keyboard in allremoveAssignedResiduesAndUpdateVoteBPG\n');
%keyboard

for i = 1:numVoters
  [voter{i}, impossibleToAssign]   = individualBPG_removeAssignedResiduesAndUpdate ...
      (ASSIGNTABLE,voter{i},unassignedRelPeakIndices, ...
       unassignedRelResidueIndices);
  if (impossibleToAssign)
    keyboard
    return;
  end
end


[ASSIGNTABLE, impossibleToAssign] = individualBPG_removeAssignedResiduesAndUpdate ...
    (ASSIGNTABLE,ASSIGNTABLE,unassignedRelPeakIndices, unassignedRelResidueIndices);

if (impossibleToAssign)
%  keyboard
  return;
end

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
function [BPG_TO_UPDATE, impossibleToAssign]=individualBPG_removeAssignedResiduesAndUpdate(ASSIGNTABLE,BPG_TO_UPDATE, ...
						  unassignedRelPeakIndices, unassignedRelResidueIndices)

%fprintf(1, 'keyboard in individualBPG_removeAssignedResiduesAndUpdate\n');
%keyboard

impossibleToAssign = 0;
if (isempty(BPG_TO_UPDATE))
%  fprintf(1 , 'input bpg is empty. impossible to assign\n');
%  impossibleToAssign = 1;
  return
end

%ASSIGNTABLE_IN = ASSIGNTABLE;
%BPG_TO_UPDATE_IN   = BPG_TO_UPDATE;
BPG_TO_UPDATE      = BPG_TO_UPDATE.*ASSIGNTABLE;

%BPG_TO_UPDATE = reinitializeAndNormalizeBPG_Entries_If_Necessary(BPG_TO_UPDATE);%, BPG_TO_UPDATE_IN);

BPG_TO_UPDATE      = BPG_TO_UPDATE  (unassignedRelPeakIndices,unassignedRelResidueIndices);

BPG_TO_UPDATE_IN   = BPG_TO_UPDATE;

[BPG_TO_UPDATE, impossibleToAssign]      = myRenormalize(BPG_TO_UPDATE);

BPG_TO_UPDATE = BPG_TO_UPDATE_IN; %undoing renormalization, keeping
                                  %probabilities intact.


%BPG_TO_UPDATE      = reinitializeAndNormalizeBPG_Entries_If_Necessary(BPG_TO_UPDATE);

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





% $$$ %##################################################################################
% $$$ %##  updateTen
% $$$ function [SAUPE, flag]= updateTen(TABLE,RDC,VECS)
% $$$ M = zeros(1,5);
% $$$ vec = 0;
% $$$ ct=1;
% $$$ for(i=1:size(TABLE,1))
% $$$    x = find(TABLE(i,:));
% $$$    if(length(x)==1 & RDC(i)>-999)
% $$$       M(ct,1) = VECS(x,1).^2-VECS(x,3).^2;
% $$$       M(ct,2) = 2*(VECS(x,1)*VECS(x,2));
% $$$       M(ct,3) = 2*(VECS(x,1)*VECS(x,3));
% $$$       M(ct,4) = VECS(x,2).^2-VECS(x,3).^2;
% $$$       M(ct,5) = 2*(VECS(x,2)*VECS(x,3));
% $$$       vec(ct) =  RDC(i); 
% $$$       ct=ct+1;
% $$$    end
% $$$ end
% $$$ flag = 0;
% $$$ if(size(M,1)<5)
% $$$    size(M,1);
% $$$    flag=1;
% $$$    return
% $$$ end
% $$$ %invert the matrix
% $$$ if(size(M,1)==5 & size(M,2)==5)
% $$$    M = inv(M);   
% $$$ else
% $$$    M = pinv(M);
% $$$ end
% $$$ %do the least squares fitting
% $$$ s = M*vec';
% $$$ SAUPE = zeros(3,3);
% $$$ SAUPE(1,1) = s(1);
% $$$ SAUPE(2,2) = s(4);
% $$$ SAUPE(3,3) = -1*(s(1)+s(4));
% $$$ SAUPE(1,2) = s(2);
% $$$ SAUPE(2,1) = s(2);
% $$$ SAUPE(1,3) = s(3);
% $$$ SAUPE(3,1) = s(3);
% $$$ SAUPE(2,3) = s(5);
% $$$ SAUPE(3,2) = s(5);