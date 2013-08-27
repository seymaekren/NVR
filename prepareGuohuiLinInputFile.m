function prepareGuohuiLinInputFile(name);

[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread(name,'%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

numDiscard = 0;
while (numDiscard < length(R))
  if (H_CS(length(R)-numDiscard) ~= -999)
    break;
  end
  if (NH_CS(length(R)-numDiscard) ~= -999)
    break;
  end
% $$$   if (RDC1(length(R)-numDiscard) ~= -999)
% $$$     break;
% $$$   end
% $$$   if (RDC2(length(R)-numDiscard) ~= -999) 
% $$$     break;
% $$$   end
  numDiscard = numDiscard + 1;
end

fprintf(1, 'discarding %d peak(s)...\n', numDiscard);
H_CS      = H_CS(1:length(H_CS) - numDiscard);
NH_CS     = NH_CS(1:length(NH_CS) - numDiscard);
RDC1      = RDC1 (1:length(RDC1) - numDiscard);
RDC2      = RDC2 (1:length(RDC2) - numDiscard);
MOLMOL_HB = MOLMOL_HB(1:length(MOLMOL_HB) - numDiscard);


fid = fopen('GuohuiLinInputFile.txt','w');
fprintf(1, 'check out GuohuiLinInputFile.txt\n');

%just write the aa number, aa code (translate it), sse code, then ...
%    coil helix strand, one of them to be 1, the others 0.

for i = 1:length(H_CS)
  
  if(strcmp(T(i),'ALA')==1)
    oneLetterAA_Code='A';
  elseif(strcmp(T(i),'CYS')==1)
    oneLetterAA_Code='C';
  elseif(strcmp(T(i),'ASP')==1)
    oneLetterAA_Code='D';
  elseif(strcmp(T(i),'GLU')==1)
    oneLetterAA_Code='E';
  elseif(strcmp(T(i),'PHE')==1)
    oneLetterAA_Code='F';
  elseif(strcmp(T(i),'GLY')==1)
    oneLetterAA_Code='G';
  elseif(strcmp(T(i),'HIS')==1)
    oneLetterAA_Code='H';
  elseif(strcmp(T(i),'ILE')==1)
    oneLetterAA_Code='I';
  elseif(strcmp(T(i),'LYS')==1)
    oneLetterAA_Code='K';
  elseif(strcmp(T(i),'LEU')==1)
    oneLetterAA_Code='L';
  elseif(strcmp(T(i),'MET')==1)
    oneLetterAA_Code='M';
  elseif(strcmp(T(i),'ASN')==1)
    oneLetterAA_Code='N';
  elseif(strcmp(T(i),'GLN')==1)
    oneLetterAA_Code='Q';
  elseif(strcmp(T(i),'ARG')==1)
    oneLetterAA_Code='R';
  elseif(strcmp(T(i),'SER')==1)
    oneLetterAA_Code='S';
  elseif(strcmp(T(i),'THR')==1)
    oneLetterAA_Code='T';
  elseif(strcmp(T(i),'VAL')==1)
    oneLetterAA_Code='V';
  elseif(strcmp(T(i),'TRP')==1)
    oneLetterAA_Code='W';
  elseif(strcmp(T(i),'TYR')==1)
    oneLetterAA_Code='Y';
  else
    PROBLEM = T(i);   
  end
  
  fprintf(fid, '%d %s %s ', R(i), oneLetterAA_Code, SS{i});
  if (strcmp(SS(i),'C') == 1)
    fprintf(fid, '1 0 0\n');
  elseif (strcmp(SS(i),'H') == 1)
    fprintf(fid, '0 1 0\n');
  elseif (strcmp(SS(i),'B') == 1) 
    fprintf(fid, '0 0 1\n');
  end
end

fclose(fid);

