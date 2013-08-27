%load allPeaksAssignmentsEnvironment.mat  
%load allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat
load allPeaksNewAssignmentEnvironment.mat
%load debugRDC_Score.mat
% $$$ 
% $$$ if (isempty(MB_RDC_Score))
% $$$   [MASTER, newASSIGNTABLE, CP, SXCP, SSCP, TP, HDE, ROWIN, COLIN,marsRdcScore,MB_RDC_Score,numRDC]= ...
% $$$       AssignFirstFiveResidues(MASTER, ASSIGNTABLE, NOES, IALLDISTS, ...
% $$$ 			      NTH, ROWIN, COLIN, ALLDISTS, CP, SXCP, ...
% $$$ 			      SSCP, RP1, RP2, TP, HDE, S1, S2, RDC1, ...
% $$$ 			      RDC2, VECTORS,marsShiftScore);
% $$$   save debugRDC_Score.mat
% $$$ end
% $$$ 



% [MASTER, ASSIGNTABLE, CP, SXCP, SSCP, TP, HDE, ROWIN, COLIN,marsRdcScore,MB_RDC_Score,numRDC]= ...
%     AssignFirstFiveResidues(MASTER, ASSIGNTABLE, NOES, IALLDISTS, ...
% 			    NTH, ROWIN, COLIN, ALLDISTS, CP, SXCP, ...
% 			    SSCP, RP1, RP2, TP, HDE, S1, S2, RDC1, ...
% $$$ 			    RDC2, VECTORS);
% $$$ 
% $$$ fprintf(1, 'assigned first %d residues. These are:\n', ...
% $$$  	sum(sum(MASTER)));
% $$$ for i = 1:size(MASTER,1)
% $$$   v = find(MASTER(i,:));
% $$$   if isempty(v)
% $$$     continue;
% $$$   end
% $$$   for j = 1:length(v)
% $$$     fprintf(1, '%d %d\n', i, v(j));
% $$$   end
% $$$ end
for peakIndex = 1:size(MASTER,1)
   MASTER(peakIndex,peakIndex) = 1;
end
S1 = updateTen(MASTER,RDC1,VECTORS);
S2 = updateTen(MASTER,RDC2,VECTORS);
RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1,ROWIN, COLIN);
RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2,ROWIN, COLIN);


[MB_RDC_Score,numRDC]                             = computeMB_RDC_Score (ASSIGNTABLE,RDC1, RDC2, VECTORS, ...
 						  S1, S2);




%[MB_RDC_Score,numRDC]                             = computeMB_RDC_Score (ASSIGNTABLE,RDC1, RDC2, VECTORS, ...
%						  S1, S2);

%[S1,saupeWrongFlag]                               = updateTen(MASTER,RDC1,VECTORS);
%assert (saupeWrongFlag == 0);
%[S2,saupeWrongFlag]                               = updateTen(MASTER,RDC2,VECTORS);
%assert (saupeWrongFlag == 0);

%[MB_RDC_Score,numRDC]                             = computeMB_RDC_Score(ASSIGNTABLE,RDC1,RDC2, VECTORS,S1, S2);b