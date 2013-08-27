function HD_SCORE = computeKemenyScore ( RDC1, RDC2, VECTORS, ...
			  HDEXCHANGE, HBOND, peakIDs, HSHIFTS, ...
			  NSHIFTS, TYPES, SSTRUCT, NOES, ALLDISTS, IALLDISTS,...
			   SHIFTS_Filename, SHIFTX_Filename, ...
				     MASTER)

addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/Routines');

S1 = updateTen(MASTER,RDC1,VECTORS);
S2 = updateTen(MASTER,RDC2,VECTORS);

[RP1_A,RP2_A,HDE_A,TP_A,CP_A,SXCP_A,SSCP_A] = prepareBPGs (RDC1,VECTORS,S1,RDC2,S2,HDEXCHANGE,HBOND,peakIDs,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,ALLDISTS, SHIFTS_Filename,SHIFTX_Filename,MASTER);

MASTER = fixDoubleAssignments(MASTER,S1,RDC1,VECTORS,S2,RDC2,TP_A,CP_A,SXCP_A,SSCP_A,HDE_A,RP1_A,RP2_A);

%JOINT       = ren(RP1_A.*RP2_A.*HDE_A.*TP_A.*CP_A.*SXCP_A.*SSCP_A).*MASTER;
%JOINT       = [nonzeros(JOINT)'];


HD_SCORE    = computeKemenyScoreForOneBPG(RP1_A,MASTER) + computeKemenyScoreForOneBPG(RP2_A,MASTER);% + computeKemenyScore(HDE_A,MASTER) + computeKemenyScore(RP1_A,MASTER) + computeKemenyScore(RP1_A,MASTER) + computeKemenyScore(RP1_A,MASTER) + computeKemenyScore(RP1_A,MASTER);



% $$$ keyboard
% $$$ 5
% $$$ 
% $$$ pseudocode:
% $$$ say for rdc1:
% $$$   given assignments.
% $$$   for each peak and its assignment:
% $$$     find the kemeny score for that pairing. that is, find the pairwise ...
% $$$ 	  score for that residue such that that residue wins in all ...
% $$$ 	  cases.would be the kement score for that peak, and so on ...
% $$$ 	  for other peaks.
% $$$       similarly for the other bpg's.
% $$$ 	sum them up.
% $$$ 	that is your kemeny score.

actualScore = function computeKemenyScoreForOneBPG(BPG, MASTER)

actualScore = 0;
	
for peakIndex = 1:70
  
  actualAssignedResidueIndex = MASTER(peakIndex,:);
  
  for residueIndex = 1:72
    if (residueIndex == actualAssignedResidueIndex)
      continue;
    end
    
    actualScore = actualScore + BPG(peakIndex,actualAssignedResidueIndex) - BPG(peakIndex,residueIndex);
        
  end
end



