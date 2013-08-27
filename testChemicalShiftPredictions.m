SHIFTX_Filename = 'SHIFTX/MySHIFTX.7.model6.hSRI';
name            = 'myinput.m.hSRI';
SHIFTS_Filename = 'SHIFTS/MySHIFTS.7.model6.hSRI';

[rn TY SS ha hn nf cb ca co]= textread(SHIFTX_Filename,'%f %s %s %f %f %f %f %f %f');

PRED = [rn  hn nf];


[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread(name,'%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

csToPlot = find(H_CS>-999);

figure;
plot(H_CS(csToPlot), hn(csToPlot), '*');
hold on
plot([min(min(H_CS(csToPlot)),min(hn(csToPlot))),max(max(H_CS(csToPlot)),max(hn(csToPlot)))],[min(min(H_CS(csToPlot)),min(hn(csToPlot))),max(max(H_CS(csToPlot)),max(hn(csToPlot)))]);
figure;
plot(NH_CS(csToPlot), nf(csToPlot), '*');
hold on
plot([min(min(NH_CS(csToPlot)),min(nf(csToPlot))),max(max(NH_CS(csToPlot)),max(nf(csToPlot)))],[min(min(NH_CS(csToPlot)),min(nf(csToPlot))),max(max(NH_CS(csToPlot)),max(nf(csToPlot)))])

[rn hn nf]= textread(SHIFTS_Filename,'%f %f %f ');

PRED = [rn  hn nf];


figure;
plot(H_CS(csToPlot), hn(csToPlot), '*');
hold on
plot([min(min(H_CS(csToPlot)),min(hn(csToPlot))),max(max(H_CS(csToPlot)),max(hn(csToPlot)))],[min(min(H_CS(csToPlot)),min(hn(csToPlot))),max(max(H_CS(csToPlot)),max(hn(csToPlot)))]);
figure;
plot(NH_CS(csToPlot), nf(csToPlot), '*');
hold on
plot([min(min(NH_CS(csToPlot)),min(nf(csToPlot))),max(max(NH_CS(csToPlot)),max(nf(csToPlot)))],[min(min(NH_CS(csToPlot)),min(nf(csToPlot))),max(max(NH_CS(csToPlot)),max(nf(csToPlot)))])
