%updates ASSIGNTABLE based on NOEs and currently assigned
%peak-residue pairs.
function [ASSIGNTABLE, impossibleToAssign] = noePrune(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,peakIndices,residueIndices)
%prune via NOEs
%NP = NVR_NOE2PROB(ASSIGNTABLE,NOES,ALLDISTS,NTH,peakIndices,residueIndices);

%keyboard
mypeakIndices                      = 1:size(MASTER,1);
myresidueIndices                   = 1:size(MASTER,2);

assert (sum(sum(MASTER(peakIndices, residueIndices))) == 0);
                                                             %%SB.
                                                             %ACTUALLY
                                                             %PEAKINDICES,:
                                                             %SHOULD
                                                             %BE EMPTY.
MASTER(peakIndices,residueIndices) = ASSIGNTABLE;
[NP,impossibleToAssign]            = NVR_NOE2PROB (MASTER,NOES,ALLDISTS,NTH, mypeakIndices,myresidueIndices);

if (~impossibleToAssign)
  differenceMatrix = MASTER - NP;
  [differentPeakIndices,differentResidueIndices] = find(differenceMatrix);
  for i = 1:length(differentPeakIndices)
    u = find(peakIndices == differentPeakIndices(i));
    v = find(residueIndices == differentResidueIndices(i));
    if (isempty(u))
      %the pruned peak is not unassigned, i.e. it is assigned already.
      assert(isempty(v));
      %the difference is not in one of the unassigned positions. It
      %is at one of the assigned positions which has been pruned.
      impossibleToAssign = 1;
      break;
    end
  end
end
NP                                 = NP           (peakIndices,residueIndices);
ASSIGNTABLE                        = and          (NP,NP);