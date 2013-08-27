function testOneMediumAlignmentTensor(NH_RDCS, CH_RDCS,VECTORS_NH, ...
				      VECTORS_CH);

MASTER = zeros(length(NH_RDCS),size(VECTORS_NH,1));
for i = 1:size(MASTER,1)
  MASTER(i,i) = 1;
end
S1 = updateTen_CH(MASTER,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
%S1 = updateTen(MASTER,NH_RDCS, VECTORS_NH);
backcomputed = zeros(1,size(VECTORS_NH,1));
for(i=1:size(VECTORS_NH,1))
   v = VECTORS_NH(i,:);
   backcomputed(i) = v*S1*v';
end
figure;
rdcsToBePlotted = [];
diff_NH_RDCs = zeros(length(NH_RDCS),1);
fid = fopen('differenceNH_RDCs.txt','w');
fprintf(1, 'check differenceNH_RDCs.txt\n');
for i = 1:length(NH_RDCS)
  if (NH_RDCS(i) ~= -999)
    rdcsToBePlotted = [rdcsToBePlotted i];
%    plot(NH_RDCS(i),backcomputed(i),'*');
%    hold on
    diff_NH_RDCs(i) = backcomputed(i) - NH_RDCS(i);
  else
    fprintf(1, 'skipping %d th RDC.\n',i);
  end
  fprintf(fid, '%f %f\n', i, diff_NH_RDCs(i));
end
fclose(fid);
keyboard


plot(NH_RDCS(rdcsToBePlotted), backcomputed(rdcsToBePlotted),'*');
corrcoef(NH_RDCS(rdcsToBePlotted), backcomputed(rdcsToBePlotted))
hold on, 
%plot([min(NH_RDCS(rdcsToBePlotted),backcomputed(rdcsToBePlotted)),
%max(NH_RDCS(rdcsToBePlotted),backcomputed(rdcsToBePlotted))],[min(backcomputed(rdcsToBePlotted)),max(backcomputed(rdcsToBePlotted))])
plot([min(min(NH_RDCS(rdcsToBePlotted)),min(backcomputed(rdcsToBePlotted))), ...
      max(max(NH_RDCS(rdcsToBePlotted)),max(backcomputed(rdcsToBePlotted)))],[min(min(NH_RDCS(rdcsToBePlotted)),min(backcomputed(rdcsToBePlotted))), ...
      max(max(NH_RDCS(rdcsToBePlotted)),max(backcomputed(rdcsToBePlotted)))]);

title('NH')

%S1 = updateTen(MASTER,CH_RDCS, VECTORS_CH);
backcomputed = zeros(1,size(VECTORS_CH,1));
for(i=1:size(VECTORS_CH,1))
   v = VECTORS_CH(i,:);
   backcomputed(i) = v*S1*v';
end
figure;
rdcsToBePlotted = [];
diff_CH_RDCs = zeros(length(CH_RDCS),1);
fid = fopen('differenceCH_RDCs.txt','w');
fprintf(1, 'check differenceCH_RDCs.txt\n');
for i = 1:length(CH_RDCS)
  if (CH_RDCS(i) ~= -999)
    rdcsToBePlotted = [rdcsToBePlotted i];
    %    plot(CH_RDCS(i),backcomputed(i),'*');
    %    hold on
    diff_CH_RDCs(i) = backcomputed(i) - CH_RDCS(i);
  else
    fprintf(1, 'skipping %d th RDC.\n',i);
  end
  fprintf(fid, '%f %f\n', i, diff_CH_RDCs(i));
end
fclose(fid);
keyboard
plot(CH_RDCS(rdcsToBePlotted), backcomputed(rdcsToBePlotted),'*');
corrcoef(CH_RDCS(rdcsToBePlotted), backcomputed(rdcsToBePlotted))
hold on, 
%plot([min(CH_RDCS(rdcsToBePlotted)), max(CH_RDCS(rdcsToBePlotted))],[min(backcomputed(rdcsToBePlotted)),max(backcomputed(rdcsToBePlotted))])
plot([min(min(CH_RDCS(rdcsToBePlotted)),min(backcomputed(rdcsToBePlotted))), ...
      max(max(CH_RDCS(rdcsToBePlotted)),max(backcomputed(rdcsToBePlotted)))],[min(min(CH_RDCS(rdcsToBePlotted)),min(backcomputed(rdcsToBePlotted))), ...
      max(max(CH_RDCS(rdcsToBePlotted)),max(backcomputed(rdcsToBePlotted)))]);
title('CH')

figure;
plot(CH_RDCS(rdcsToBePlotted)'-backcomputed(rdcsToBePlotted),'*')
keyboard
DO_GLYCINE_RDCs = 0;

if (DO_GLYCINE_RDCs)

  hold on
  %figure
  glycineRDC_Matrix     = load ('NjpData/glycineC-H.tbl.m');
  glycineVectors        = load ('NjpData/glycineC-H_vectors.m');
  
  glyResidueIndices = glycineRDC_Matrix(:,1);
  glyRDCs           = glycineRDC_Matrix(:,2);
  
  backcomputedRDCs  = [];
  
  
  for relResIndex = 1:length(glyResidueIndices)
    glyResidueIndex = glyResidueIndices(relResIndex);
    relVectorIndex = find(glycineVectors(:,1) == glyResidueIndex);
    assert (length(relVectorIndex) == 2);
    backcomputed = 0;
    for i = 1:2
      v = glycineVectors(relVectorIndex(i),2:4);
      backcomputed = backcomputed + v * S1 * v';
    end
    backcomputed = backcomputed/2;
    backcomputedRDCs = [backcomputedRDCs backcomputed];
  end
  
  plot(backcomputedRDCs, glyRDCs, 'r*'); hold on;
  %plot([min(backcomputedRDCs),max(backcomputedRDCs)],[min(glyRDCs),max(glyRDCs)])
  plot([min(min(backcomputedRDCs),min(glyRDCs)), ...
	max(max(backcomputedRDCs),max(glyRDCs))],[min(min(backcomputedRDCs),min(glyRDCs)),max(max(backcomputedRDCs),max(glyRDCs))])
end