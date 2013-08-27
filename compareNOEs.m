load parsedBackboneNOEs.txt
NOES = load ('~/speech/Workdir/InputFiles/NOES.txt.1UBI');

temp = parsedBackboneNOEs;
parsedBackboneNOEs = NOES;
NOES = temp;

NOE_Matrix1 = zeros(200,200);
for i = 1:length(parsedBackboneNOEs)
  NOE_Matrix1(parsedBackboneNOEs(i,1),parsedBackboneNOEs(i,2)) = 1;
  NOE_Matrix1(parsedBackboneNOEs(i,2),parsedBackboneNOEs(i,1)) = 1;
end


[numNOEs2, dummyVar] = size(NOES);
for i = 1:numNOEs2
  residue1 = NOES(i,1);
  residue2 = NOES(i,2);
 
  if (NOE_Matrix1(residue1,residue2) == 1)
    fprintf(1, 'found the NVR NOE in the MZ noe file.\n');
  else
    fprintf(1, 'did not find the NVR NOE between %d and %d in MZ noe_file.\n',residue1,residue2);
  end
end



