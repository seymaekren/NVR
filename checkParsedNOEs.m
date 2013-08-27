%filename   = '1EZA.NOE_List.HnHn';
%filename   = '1EZA.NOE_List.HnHn.C_TerminusRemoved';
filename   = 'unambiguousNOEs.txt';

parsedNOEs = load (filename);

maxIndices = max(parsedNOEs);

numResidues           = max(maxIndices(1),maxIndices(2));

NOEs                  = zeros(numResidues, numResidues);

for i = 1:size(parsedNOEs,1)
  NOEs(parsedNOEs(i,1),parsedNOEs(i,2)) = 1;
  NOEs(parsedNOEs(i,2),parsedNOEs(i,1)) = 1;
end

count = 0;
for i = 1:numResidues
  for j = i+1:numResidues
    if (NOEs(i,j) == 1)
      count  = count + 1;
    end
  end
end
fprintf(1, 'there are %d NOEs.\n',count);