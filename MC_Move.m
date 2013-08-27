function MASTER = MC_Move(MASTER)


[numPeaks,numResidues] = size(MASTER);

if (numPeaks == numResidues)
  moveType = 1;
else
  moveType=  floor(rand() * 2);
end

assert ((moveType == 0) | (moveType == 1));
if (moveType == 0)
  peakIndex = floor(rand() * numPeaks) + 1;
  assert ((peakIndex >= 1) & (peakIndex <= numPeaks));
  assert ((numResidues >= numPeaks));
  sumOfRows  = sum(MASTER);
  unassignedResidueIndices = find(sumOfRows == 0);
  
  newResidueIndex  = unassignedResidueIndices(floor(rand() * length(unassignedResidueIndices))+1);
  
  
  
  MASTER(peakIndex,:) = 0;
  MASTER(peakIndex,newResidueIndex) = 1;
  
else
  
  peak1Index = floor(rand() * numPeaks) + 1;
  while  (1)
    peak2Index = floor(rand() * numPeaks) + 1;
    if (peak2Index ~= peak1Index)
      break;
    end
  end
  residue1Index = find(MASTER(peak1Index,:));
  residue2Index = find(MASTER(peak2Index,:));
  assert (length(residue1Index) == 1);
  assert (length(residue2Index) == 1);
  
  MASTER (peak1Index,residue1Index) = 0;
  MASTER (peak1Index,residue2Index) = 1;
  MASTER (peak2Index,residue1Index) = 1;
  MASTER (peak2Index,residue2Index) = 0;
  
  
end