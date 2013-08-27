function noesAreOK = noeCheck(MASTER, NOES, ALLDISTS, NTH)
mypeakIndices                      = 1:size(MASTER,1);
myresidueIndices                   = 1:size(MASTER,2);

[NP,impossibleToAssign]            = NVR_NOE2PROB (MASTER,NOES, ...
						  ALLDISTS,NTH, mypeakIndices,myresidueIndices);
if (impossibleToAssign)
  noesAreOK = 0;
  return;
end

noesAreOK = 1;

differenceMatrix = MASTER - NP;
[differentPeakIndices,differentResidueIndices] = ...
    find(differenceMatrix);
if (~isempty(differentPeakIndices))
  noesAreOK = 0;
  return;
end