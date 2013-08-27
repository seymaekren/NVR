clear all;
automatedRun                      = 1;
v_refineWithRDCs                  = [0 0 0 0 0 1 1];
v_b_printIndividualScoringMatrix  = [1 1 1 1 1 1 0];
v_b_printCP                       = [0 1 0 0 0 0 0];
v_b_printSSCP                     = [1 0 0 0 0 0 0];
v_b_printSXCP                     = [0 0 1 0 0 0 0];
v_b_printTP                       = [0 0 0 1 0 0 0];
v_b_printCS                       = [0 0 0 0 1 0 0];
v_b_printRDC                      = [0 0 0 0 0 1 0];


% $$$ v_refineWithRDCs                  = [            1];
% $$$ v_b_printIndividualScoringMatrix  = [            0];
% $$$ v_b_printCP                       = [            0];
% $$$ v_b_printSSCP                     = [            0];
% $$$ v_b_printSXCP                     = [            0];
% $$$ v_b_printTP                       = [            0];
% $$$ v_b_printCS                       = [            0];
% $$$ v_b_printRDC                      = [            0];

for runIndex = 1:length(v_refineWithRDCs)

  refineWithRDCs = v_refineWithRDCs(runIndex);
  b_printIndividualScoringMatrix = ...
      v_b_printIndividualScoringMatrix(runIndex);
  b_printCP                      = v_b_printCP(runIndex);
  b_printSSCP                    = v_b_printSSCP(runIndex);
  b_printSXCP                    = v_b_printSXCP(runIndex);
  b_printTP                      = v_b_printTP(runIndex);
  b_printCS                      = v_b_printCS(runIndex);
  b_printRDC                     = v_b_printRDC(runIndex);
  runBetterAutomated(automatedRun, refineWithRDCs, ...
		     b_printIndividualScoringMatrix, b_printCP, ...
		     b_printSSCP, b_printSXCP, b_printTP, ...
		     b_printCS, b_printRDC);
end
