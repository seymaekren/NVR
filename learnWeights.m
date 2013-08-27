%load ('allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat');

load ('allPeaksAfterInitialAssignments--For1UBI-WithFullyCorrectAlignmentTensorsForMB_RDC_Score.mat');

[I,J] = find(MB_RDC_Score > 1);
if (isempty(I))
  fprintf(1, 'MB_RDC_Score has no score greater than 1.\n');
else
  fprintf(1, 'MB_RDC_Score has some scores greater than 1.\n');
end

[I,J] = find(MB_ShiftScore > 1);
if (isempty(I))
  fprintf(1, 'MB_ShiftScore has no score greater than 1.\n');
else
  fprintf(1, 'MB_ShiftScore has some scores greater than 1.\n');
end

numAdditions = 0; csWeight = 0; rdcWeight = 0; noeWeight = 0;
soleCS_Distinguishes = 0; soleRDC_Distinguishes = 0; soleNOE_Distinguishes ...
    = 0;nooneDistinguishes = 0;

[numPeaks,numResidues] = size(MASTER);

myAlmostCompleteMASTER = MASTER*0;

for peakIndex = 1:numPeaks
  myAlmostCompleteMASTER(peakIndex,peakIndex) =  1;
end
[correctAssignmentNoeScore,numHN_NOEs] = computeMB_NoeScore(HSQCDATA,ALLDISTS,NTH,myAlmostCompleteMASTER);
correctAssignmentNoeScore              = correctAssignmentNoeScore/numHN_NOEs;


csScoreDiff = []; rdcScoreDiff = [];noeScoreDiff = [];
csScoreEntropy = 0; rdcScoreEntropy = 0;

for peakIndex = 1:numPeaks
  
  if (peakIndex > 1)
    myAlmostCompleteMASTER(peakIndex-1,:)           = 0;
    myAlmostCompleteMASTER(peakIndex-1,peakIndex-1) = 1;
  end

  for residueIndex = 1:numResidues

    if (residueIndex == peakIndex)
      continue;
    end

    
    myAlmostCompleteMASTER(peakIndex,:) = 0;
    myAlmostCompleteMASTER(peakIndex,residueIndex) = 1;

    newNoeScore = computeMB_NoeScore(HSQCDATA,ALLDISTS,NTH,myAlmostCompleteMASTER);
    newNoeScore = newNoeScore/numHN_NOEs;

    p_noe       = correctAssignmentNoeScore-newNoeScore;
    
    p_cs        = (MB_ShiftScore(peakIndex,peakIndex)- ...
    		   MB_ShiftScore(peakIndex,residueIndex))/2;
    p_rdc       = (MB_RDC_Score(peakIndex,peakIndex)- ...
    		   MB_RDC_Score(peakIndex,residueIndex))/2;
    
    assert (abs(p_cs) <= 1);
    assert (abs(p_rdc) <= 1);
    assert (abs(p_noe) <= 1);
    
%    fprintf(1, 'p_cs = %f p_rdc = %f p_noe = %f\n', p_cs, p_rdc, p_noe);
    
    
%    p_cs        = (MB_ShiftScore(peakIndex,peakIndex)- ...
%		   MB_ShiftScore(peakIndex,residueIndex))/numCS;
%    p_rdc       = (MB_RDC_Score(peakIndex,peakIndex)- ...
%		   MB_RDC_Score(peakIndex,residueIndex))/numRDC;
    
%    plot(p_cs,p_rdc,'*'); hold on;

    
%    if (p_cs >= p_rdc) & (p_cs >= p_noe) 
     if (p_cs > 0)
       csWeight  = csWeight + 1;
     end

     if (p_cs >0) & (p_rdc <= 0) & (p_noe <= 0)
       soleCS_Distinguishes = soleCS_Distinguishes + 1;
     end
     
     if (p_cs <=0) & (p_rdc > 0) & (p_noe <= 0)
       soleRDC_Distinguishes = soleRDC_Distinguishes + 1;
     end
     
     
     if (p_cs <=0) & (p_rdc <= 0) & (p_noe > 0)
       soleNOE_Distinguishes = soleNOE_Distinguishes + 1;
     end
     
     if (p_cs <=0) & (p_rdc <= 0) & (p_noe <= 0)
       nooneDistinguishes = nooneDistinguishes + 1;
     end

     
     
     if (p_rdc > 0)
       %    elseif (p_rdc>=p_cs) & (p_rdc>=p_noe)
       rdcWeight = rdcWeight + 1;
     end
      %    else
%      assert ((p_noe>=p_rdc) & (p_noe>=p_cs));
    if (p_noe > 0)
      noeWeight = noeWeight + 1;
    end
    
    numAdditions = numAdditions + 1;
    csScoreDiff(numAdditions)  = p_cs;
    rdcScoreDiff(numAdditions) = p_rdc;
    noeScoreDiff(numAdditions) = p_noe;
  end
end

fprintf(1, 'pairwise comparisonlarda sadece csin distinguish ettigi');

fprintf(1, 'case ler %f\n', soleCS_Distinguishes/numAdditions);
fprintf(1, 'rdc: %f noe: %f\n', soleRDC_Distinguishes/numAdditions, ...
	soleNOE_Distinguishes/numAdditions);

fprintf(1, 'no data source distinguishes correctly in %f percent cases.\n',nooneDistinguishes*100/numAdditions)
fprintf(1, 'the cs weight would be %f\n',csWeight/numAdditions);
fprintf(1, 'the rdc weight would be %f\n',rdcWeight/numAdditions);
fprintf(1, 'the noe weight would be %f\n',noeWeight/numAdditions);
figure; 
plot(csScoreDiff,rdcScoreDiff,'*');
xlabel('cs');ylabel('rdc');
figure;
plot(csScoreDiff,noeScoreDiff,'*');
xlabel('cs');ylabel('noe');
figure;
plot(noeScoreDiff,rdcScoreDiff,'*');
xlabel('noe');ylabel('rdc');
