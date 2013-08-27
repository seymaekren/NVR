function noes = convertResidueIndicesToPeakIndicesInNOE_File(unadjustedNOEs, order);
%changes the format of the NOES file to the new one.
%The original one has many columns in addition to the peak1 and
%peak2 indices.

%order: contains the residue indices as they appear in the modeldata file.
%unadjustedNOEs : contains the residue indices for protons having
%an NOE.

%the returned value contains the peak indices (in contrast to
%residue indices).


[m,n] = size(unadjustedNOEs);
noes = zeros(m,2);
for i = 1:m
  index1    = unadjustedNOEs(i,1);
  index2    = unadjustedNOEs(i,2);
  newIndex1 = find(order == index1);
  newIndex2 = find(order == index2);

  if ((isempty(newIndex1)) | (isempty(newIndex2)))

    fprintf(1, 'could not find the peak indices corresponding');
    fprintf(1, ' to a pair of residues having an NOE.\n');
    fprintf(1, 'the residue indices are %d and %d\n',index1,index2);
    continue;
  end
  
  
  assert (~isempty(newIndex1));
  noes(i,1) = newIndex1;

  assert (~isempty(newIndex2));  
  noes(i,2) = newIndex2;

end
