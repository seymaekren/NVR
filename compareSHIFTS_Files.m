function compareSHIFTS_Files(SHIFTS_Filename1, SHIFTS_Filename2)

[rn1 hn1 nf1]= textread(SHIFTS_Filename1,'%f %f %f');

PRED1 = [rn1  hn1 nf1];

[rn2 hn2 nf2]= textread(SHIFTS_Filename2,'%f %f %f');

PRED2 = [rn2  hn2 nf2];

xlabelStr = sprintf('%s',SHIFTS_Filename1);
ylabelStr = sprintf('%s',SHIFTS_Filename2);


figure; 

plot(hn1,hn2,'*');
xlabel(xlabelStr);
ylabel(ylabelStr);
title('comparison of H chemical shifts obtained with SHIFTS');

figure; 

plot(nf1,nf2,'o');
xlabel(xlabelStr);
ylabel(ylabelStr);
title('comparison of N chemical shifts obtained with SHIFTS');
