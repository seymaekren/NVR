function reduceSizeOfTheProblem
fprintf(1, 'restricting the set of peaks\n');
numTestPeaks    = 20;
numTestResidues = 22;
%  peakIndicesOffset = 11;
ROWIN           = ROWIN(1:numTestPeaks);
COLIN           = COLIN(1:numTestResidues);
ASSIGNTABLE     = ASSIGNTABLE(1:numTestPeaks, 1:numTestResidues);
MASTER          = MASTER(1:numTestPeaks, 1:numTestResidues);
HDE             = HDE(1:numTestPeaks, 1:numTestResidues);

[HDE, impossibleToAssign] = myRenormalize(HDE);
assert(~impossibleToAssign);

TP          = TP(1:numTestPeaks, 1:numTestResidues);
[TP, impossibleToAssign] = myRenormalize(TP);
assert(~impossibleToAssign);


CP          = CP(1:numTestPeaks, 1:numTestResidues);
[CP, impossibleToAssign] = myRenormalize(CP);
assert(~impossibleToAssign);


SXCP        = SXCP(1:numTestPeaks, 1:numTestResidues);
[SXCP, impossibleToAssign] = myRenormalize(SXCP);
assert(~impossibleToAssign);


SSCP        = SSCP(1:numTestPeaks, 1:numTestResidues);
[SSCP, impossibleToAssign] = myRenormalize(SSCP);
assert(~impossibleToAssign);

RP1         = RP1(1:numTestPeaks, 1:numTestResidues);
[RP1, impossibleToAssign] = myRenormalize(RP1);
assert(~impossibleToAssign);  

RP2         = RP2(1:numTestPeaks, 1:numTestResidues);
[RP2, impossibleToAssign] = myRenormalize(RP2);
assert(~impossibleToAssign);  

RDC1        = RDC1(1:numTestPeaks);
RDC2        = RDC2(1:numTestPeaks);

peakIDs     = peakIDs(1:numTestPeaks);
RESNUMS     = RESNUMS(1:numTestResidues);
