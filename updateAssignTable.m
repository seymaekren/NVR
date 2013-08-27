function   ASSIGNTABLE = updateAssignTable(ASSIGNTABLE,bestRelPeakIndex,bestRelResidueIndex)



assert ((bestRelPeakIndex>=1) & (bestRelPeakIndex <= ...
				 size(ASSIGNTABLE,1)));
assert ((bestRelResidueIndex >= 1) & (bestRelResidueIndex <= size(ASSIGNTABLE,2)));
ASSIGNTABLE  (bestRelPeakIndex,                   :) = 0;
ASSIGNTABLE  (:               , bestRelResidueIndex) = 0;
ASSIGNTABLE  (bestRelPeakIndex, bestRelResidueIndex) = 1;

%what to do if there is a peak left without a potential
%assignment in this case? restore its assignments to all residues
%uniformly could be an option.

%ASSIGNTABLE                           = reinitializeAndNormalizeBPG_Entries_If_Necessary(ASSIGNTABLE);

%ASSIGNTABLE                           = and (ASSIGNTABLE,ASSIGNTABLE);