function score = computeMB_Score(MASTER, MB_ShiftScore,MB_RDC_Score, ...
				 numCS,numRDC, HSQCDATA,ALLDISTS,NTH, ...
				 shiftCoefficient, rdcCoefficient, noeCoefficient)

%shiftCoefficient = 2; rdcCoefficient = 1; noeCoefficient = 5;
score = 0;
numPeaks = size(MASTER,1);

assert (numPeaks == size(MB_ShiftScore,1));
assert (numPeaks == size(MB_RDC_Score ,1));


for peakIndex = 1:numPeaks
  residueIndex = find(MASTER(peakIndex,:));
  assert (length(residueIndex) <= 1);
  if (~isempty(residueIndex))
    score = score + shiftCoefficient*MB_ShiftScore(peakIndex,residueIndex) ...
	    + rdcCoefficient * MB_RDC_Score(peakIndex, residueIndex);
  end
end

%computeMB_NoeScore;
noeScore = 0; numHN_NOES = 0;
if (noeCoefficient > 0)
  [noeScore,numHN_NOES] = computeMB_NoeScore(HSQCDATA,ALLDISTS,NTH, ...
					     MASTER);
end
%noeScore   = 0;
%numHN_NOES = 0;
%score = (score + noeCoefficient* noeScore)/(shiftCoefficient*numCS
%+  rdcCoefficient*numRDC + noeCoefficient*numHN_NOES);
score = (score + noeCoefficient* noeScore)/(shiftCoefficient*numCS +  rdcCoefficient*numRDC + noeCoefficient*numHN_NOES);
%score = score/numRDC ;
%score = noeScore/numHN_NOES;
%keyboard