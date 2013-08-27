%function [HD_SCORE,assignmentAccuracy] = betterHD(peaks,rdcs,HDEXCHANGE, ...
%						  peakIDs, NOES, ...
%						  VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, SHIFTS_Filename,SHIFTX_Filename, useMBM_EM)

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

if (~exist('ROWIN'))
  
  [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, RDC1, RDC2, NTH] = initialize(peaks,rdcs,HDEXCHANGE, ...
						  peakIDs, NOES, ...
						  VECTORS,TYPES, ...
						  RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, SHIFTS_Filename, SHIFTX_Filename);
  
  
  
  keyboard
  numTestPeaks    = 5;
  numTestResidues = 7;
  ROWIN       = ROWIN(1:numTestPeaks);
  COLIN       = COLIN(1:numTestResidues);
  ASSIGNTABLE = ASSIGNTABLE(1:numTestPeaks, 1:numTestResidues);
  MASTER      = MASTER(1:numTestPeaks, 1:numTestResidues);
  HDE         = HDE(1:numTestPeaks, 1:numTestResidues);
  
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
  
  RDC1        = RDC1(1:numTestPeaks);
  RDC2        = RDC2(1:numTestPeaks);

  peakIDs     = peakIDs(1:numTestPeaks);
  RESNUMS     = RESNUMS(1:numTestResidues);
  keyboard
end

[totalNumPeaks,totalNumResidues] = size(MASTER);
numAssignedPeaks                 = 0;
numUnassignedPeaks               = totalNumPeaks;

numVoters           = 5;

%change in variable names:
peakIndices    = ROWIN;
residueIndices = COLIN;
voter          = cell(7,1);
voter          = initializeVoters(voter,SSCP,SXCP,CP,TP,HDE);

%analyzeVoters(voter, numVoters, peakIndices, residueIndices);

%parameters of the program...

%useMBM_EM = 1
debugIncorrectAssignment = 0;


totalNumAssignments = 0;%foundAnAssignment = 0;
%minScore = 0;
overallMASTER = MASTER;
scoreSoFar = 0;

overallMASTER = assign(totalNumAssignments, MASTER, scoreSoFar, voter, ...
		       numVoters, ASSIGNTABLE, peakIndices, residueIndices, NOES, ...
		       ALLDISTS, IALLDISTS, NTH, overallMASTER, ...
		       peakIDs, RESNUMS, RDC1, RDC2, VECTORS);

analyzeOverallAssignments(overallMASTER, peakIDs, RESNUMS);

keyboard

