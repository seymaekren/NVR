fullNOE_List      = load('1EZA.NOE_List.HnHn.C_TerminusRemoved');
ambiguousNOE_List = load('ambiguousHSQC_ResidueIndices.txt');

fid = fopen('unambiguousNOEs.txt','w');

for i = 1:size(fullNOE_List)
  residueIndex1 = fullNOE_List(i,1);
  residueIndex2 = fullNOE_List(i,2);
  [relIndex1,relIndex2]  = find(ambiguousNOE_List == residueIndex1);
  if (isempty(relIndex1))
    [relIndex3,relIndex4]  = find(ambiguousNOE_List == residueIndex2);
    if (isempty(relIndex3))
      fprintf(fid, '%d %d\n',residueIndex1,residueIndex2);
    end
  end
end