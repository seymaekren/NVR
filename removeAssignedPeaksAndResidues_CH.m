function [ASSIGNTABLE,voter,numVoters, unassignedPeakIndicesCorrespondingToAssigntable,unassignedResidueIndicesCorrespondingToAssigntable,SSCP,SXCP,CP] ...
    = removeAssignedPeaksAndResidues_CH(MASTER,ASSIGNTABLE,voter,numVoters,unassignedPeakIndicesCorrespondingToAssigntable,unassignedResidueIndicesCorrespondingToAssigntable,SSCP,SXCP,CP);

%This function finds the overall unassigned peak and residue
%indices. It then finds which of these peaks and residues
%correspond to the rows and columns of ASSIGNTABLE, and retains
%those rows and columns in ASSIGNTABLE (and of the voter arrays in
%parallel). So for instance if there is a peak which is unassigned
% but which does not exist amongst the
%set of unassigned peaks, we can't do anything in terms of the
%reduction of ASSIGNTABLE. 


%keyboard
overallUnassignedPeakIndices = []; overallUnassignedResidueIndices = [];

[totalNumPeaks,totalNumResidues] = size(MASTER);
for i = 1:totalNumPeaks
  if isempty(find(MASTER(i,:)))
    overallUnassignedPeakIndices = [overallUnassignedPeakIndices i];
  end
end

for j = 1:totalNumResidues
  if isempty(find(MASTER(:,j)))
    overallUnassignedResidueIndices = [overallUnassignedResidueIndices j];
  end
end

relPeakIndices = []; relResidueIndices = [];

for peakIndex = 1:length(unassignedPeakIndicesCorrespondingToAssigntable)
  %assert (~isempty(find(peakIndices == ROWIN(i))));
  if (~isempty(find(overallUnassignedPeakIndices == unassignedPeakIndicesCorrespondingToAssigntable(peakIndex))))
    relPeakIndices = [relPeakIndices peakIndex];
  end
end

for residueIndex = 1:length(unassignedResidueIndicesCorrespondingToAssigntable)
%  assert (~isempty(find(residueIndices == COLIN(j))));
  if (~isempty(find(overallUnassignedResidueIndices == unassignedResidueIndicesCorrespondingToAssigntable(residueIndex))))
    relResidueIndices = [relResidueIndices residueIndex];
  end
end

%  [SSCP,SXCP,CP,TP,HDE, RP1, RP2, ASSIGNTABLE] = removeAssignedRowsAndColumns(SSCP,SXCP,CP,TP,HDE, RP1, RP2, ...
%						  ASSIGNTABLE,relPeakIndices, relResidueIndices);

%assert (numVoters == 7);
for i = 1:numVoters
  voter{i} =  voter{i}(relPeakIndices, relResidueIndices);
end
ASSIGNTABLE    = ASSIGNTABLE (relPeakIndices, relResidueIndices);

SSCP = SSCP(relPeakIndices, relResidueIndices);
SXCP = SXCP(relPeakIndices, relResidueIndices);
CP   = CP  (relPeakIndices, relResidueIndices);


unassignedPeakIndicesCorrespondingToAssigntable    = unassignedPeakIndicesCorrespondingToAssigntable(relPeakIndices);
unassignedResidueIndicesCorrespondingToAssigntable = unassignedResidueIndicesCorrespondingToAssigntable(relResidueIndices);
