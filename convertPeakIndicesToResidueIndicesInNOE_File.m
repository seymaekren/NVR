function convertPeakIndicesToResidueIndicesInNOE_File;

%order = load ('order.m.XinGao');
%unadjustedNOEs = load('NOE_List_peakIndices.txt.XinGao'); 
unadjustedNOEs = load('unambiguousNOEs.txt'); 
order = load('order.m.CAM');
outFilename = 'NOEs_withResidueIndices.txt'; 

%order          : contains the residue indices as they appear in the modeldata file.
%unadjustedNOEs : contains the HSQC peak indices of the NOEs.

%writes the absolute residue indices to a file.

%This is written to convert the NOEs obtained directly from raw
%data and containing no info on the residue indices into
%the residue indices so the program readNMR_Data2.m can read them.



[m,n] = size(unadjustedNOEs);
noes  = zeros(m,2);
numUnambiguousNOEs = 0;

for i = 1:m
  index1    = unadjustedNOEs(i,1);
  index2    = unadjustedNOEs(i,2);
  assert ((index1 >= 1) & (index1 <= size(order,1)));
  assert ((index2 >= 1) & (index2 <= size(order,1)));

  newIndex1 = order(index1,2);
  newIndex2 = order(index2,2);

%  if ((isempty(newIndex1)) | (isempty(newIndex2)))
  if ((newIndex1 == -1) | (newIndex2 == -1))
    continue;
  end
  
  numUnambiguousNOEs = numUnambiguousNOEs + 1;  
%  assert (~isempty(newIndex1));
  noes(numUnambiguousNOEs,1) = newIndex1;

%  assert (~isempty(newIndex2));  
  noes(numUnambiguousNOEs,2) = newIndex2;

end


fprintf(1, 'check out %s\n', outFilename);
fid   = fopen(outFilename,'w');
for i = 1:numUnambiguousNOEs
  fprintf(fid, '%d %d\n', noes(i,1),noes(i,2));
end
fclose(fid);