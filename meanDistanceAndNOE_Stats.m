useOrigData       = 1;
[HSQCDATA, NOES]  = readNMR_Data2(useOrigData);
[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
 ignoredHSQCDATA] = loaddata('InputFiles/myinput.m');

 

mu=mean(mean(ALLDISTS));

maxNOE_Distance = 0;
for i = 1:size(NOES,1)
  for j = 1:size(NOES,2)
    if (NOES(i,j) == 1)
      residue1 = i;
      residue2 = j;
      if (ALLDISTS(residue1,residue2) > maxNOE_Distance)
	maxNOE_Distance = ALLDISTS(residue1,residue2);
	%keyboard
      end
    end
  end
end

fprintf(1, 'maxNOE_Distance = %f\n', maxNOE_Distance);
fprintf(1, 'mean mean ALLDISTS = %f\n', mu);