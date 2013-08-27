load ('allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat');

csCoefficient = 0.22; rdcCoefficient = 0.05; noeCoefficient = 0.73;
scores        = (csCoefficient*MB_ShiftScore(1:3,:)+rdcCoefficient* ...
    MB_RDC_Score(1:3,:))/(csCoefficient*numCS+rdcCoefficient*numRDC);

while (1)
  maxScore = max(max(scores));
  if (maxScore == 0)
    break;
  end
  [peakIndex,residueIndex] = find(scores == maxScore);
  if (length(peakIndex) > 1)
    peakIndex = peakIndex(1);
    residueIndex = residueIndex(1);
  end
  if ((peakIndex == 1) & (residueIndex == 1) )
    break;
  end
  scores(peakIndex,residueIndex) = 0;
  fprintf(1, 'max. scoring candidate is peak %d and residue %d, score=%f\n',peakIndex,residueIndex,maxScore);
  fprintf(1, 'MB_ShiftScore for the entry = %f, MB_RDC_Score is %f\n',MB_ShiftScore(peakIndex,residueIndex),MB_RDC_Score(peakIndex,residueIndex));
  %  keyboard
end


