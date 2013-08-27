load unambiguousNOEs.txt

numResidues = max(max(unambiguousNOEs));

for resIndex = min(min(unambiguousNOEs)):numResidues
  [rowIndex,colIndex] = find(unambiguousNOEs == resIndex);
  if (length(rowIndex) == 4)
    fprintf(1, 'there are 4 NOEs for residue %d\n', resIndex);
  end
end