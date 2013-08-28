function [possibleAssignmentsBPG, impossibleToAssign] = NVR_NOE2PROB(possibleAssignmentsBPG,NOES,ALLDISTS,NTH, ROWIN,COLIN);

impossibleToAssign = 0;
possibleAssignmentsBPG=and(possibleAssignmentsBPG,possibleAssignmentsBPG);

for relPeakIndex=1:size(possibleAssignmentsBPG,1)

  absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex = find(NOES(ROWIN(relPeakIndex),:));%set of peaks having an noe with this one
  
  if (length(absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex) == 0 )      
    %if (peak # relPeakIndex has NOE constraints)
    %NOES(ROWIN(relPeakIndex),:) refers
    %to the NOES peak #relPeakIndex has.
    continue;
  end
    
  possibleAbsoluteResidueIndicesForPeak_RelPeakIndex = COLIN(find(possibleAssignmentsBPG(relPeakIndex,:))); %find out who it is could be
          

  for (j=1:length(absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex)) 
    %go through all NOE peaks with peak relPeakIndex
    
    relPeakIndexHavingNOEWithPeak_relPeakIndex=find(ROWIN==absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex(j));
    
    %previously it may have been possible not to find the peak having
    %an NOE with peak relPeakIndex, since that peak could have been
    %assigned. In that case the old HD did not make use of that
    %information. This allowed not to make incorrect pruning in case
    %previous assignment was incorrect, but also did not take advantage
    %of extra data, in case the previous assignment was correct.
    
    if(length(relPeakIndexHavingNOEWithPeak_relPeakIndex) == 0)
      %if peak with relative index
      %relPeakIndexHavingNOEWithPeak_relPeakIndex is assigned, then it
      %may be possible not to find its index in ROWIN.
      continue;
    end
    
    
    relResidueIndices = find(possibleAssignmentsBPG(relPeakIndexHavingNOEWithPeak_relPeakIndex,:));
      
    for k = 1:length(relResidueIndices)
      relResidueIndex = relResidueIndices(k);
      
      assert (possibleAssignmentsBPG(relPeakIndexHavingNOEWithPeak_relPeakIndex,relResidueIndex)~= 0);
      
      absResidueIndex = COLIN(relResidueIndex);
      
      distOfRes_AbsResIdxToPotOfPk_RelPkIdx = ALLDISTS(possibleAbsoluteResidueIndicesForPeak_RelPeakIndex,absResidueIndex);
    end %for k
  end %for j
end%for relPeakIndex






