useOrigData = 1;
[HSQCDATA, NOES]  = readNMR_Data2(useOrigData);
[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
 ignoredHSQCDATA] = loaddata('InputFiles/myinput.m');

 
%keyboard
%compute a tollerance for NOE distances
NTH=4.8;
mu=mean(mean(ALLDISTS));
if(mu-12.9>0)
   NTH=NTH+(mu-12.9);
%  NTH=min(NTH,8)
end

keyboard

NTH = NTH + 2.34;

fprintf(1, 'incrementing NTH automatically by 1.5 A\n');
NTH = NTH + 1.5;
NTH = min(NTH,9.33) %this is due to FF2.

fprintf(1, 'NTH = %f\n', NTH);
%fprintf(1, 'testing an NTH 1.45A higher than normal.\n');
%NTH = NTH + 1.45;

%fprintf(1, 'setting NTH deliberately to the following for MBP:\n');
%NTH = 6;

NTH = 5;
%NTH = 7.2;
fprintf(1, 'NTH set to %f deliberately.\n',NTH);
fprintf(1, 'set NTH to a different value here if you want to.\n');
keyboard

maxViolation = 0; maxNOE_Connected_Proton_Proton_Distance = 0;
for i = 1:size(NOES,1)
  for j = 1:size(NOES,2)
    if (NOES(i,j) == 1)
      residue1 = i; %NOTE: THESE ARE NOT ABSOLUTE (RESIDUE) INDICES.
      residue2 = j;
      if (ALLDISTS(residue1,residue2) > NTH)
	fprintf(1, 'warning. peak %d and %d have an NOE but',i,j);
	fprintf(1, ' the corresponding residues have a distance');
	fprintf(1, ' of %f whereas NTH is %f\n', ALLDISTS(residue1,residue2),NTH);
	if (ALLDISTS(residue1,residue2) > maxViolation )
	  maxViolation = ALLDISTS(residue1,residue2);
	end
	%keyboard
      end
      if (ALLDISTS(residue1, residue2) > ...
	  maxNOE_Connected_Proton_Proton_Distance)
	maxNOE_Connected_Proton_Proton_Distance = ALLDISTS(residue1, ...
						  residue2);
      end
    end
  end
end

fprintf(1, 'maxViolationDistance = %f\n', maxViolation);
if (maxViolation > 0)
  fprintf(1, 'need to increase NTH by %f\n', maxViolation-NTH);
end
fprintf(1, 'maxNOE_Connected_Proton_Proton_Distance = %f\n', maxNOE_Connected_Proton_Proton_Distance);