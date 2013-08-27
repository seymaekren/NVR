function compareSHIFTX_Files(SHIFTX_Filename1, SHIFTX_Filename2)

[rn1 TY SS ha hn1 nf1 cb ca co]= textread(SHIFTX_Filename1,'%f %s %s %f %f %f %f %f %f');

PRED1 = [rn1  hn1 nf1];

[rn2 TY SS ha hn2 nf2 cb ca co]= textread(SHIFTX_Filename2,'%f %s %s %f %f %f %f %f %f');

PRED2 = [rn2  hn2 nf2];

xlabelStr = sprintf('%s',SHIFTX_Filename1);
ylabelStr = sprintf('%s',SHIFTX_Filename2);


figure; 

plot(hn1,hn2,'*');
xlabel(xlabelStr);
ylabel(ylabelStr);
title('comparison of H and N chemical shifts');

figure; 

plot(nf1,nf2,'o');
xlabel(xlabelStr);
ylabel(ylabelStr);
title('comparison of H and N chemical shifts');
