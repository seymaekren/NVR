function test2media()

% test1medium()
% Test tensor calculation and RDC backcalculation for two sets of the same type of RDCs in different media.
% Uses updateTen() from original NVR distribution.
% Nick Patrick

dbstop if error

RDC1 = load('data/N-H_medium1.m');
RDC2 = load('data/N-H_medium2.m');
RDC1 = RDC1(:,2);
RDC2 = RDC2(:,2);

VECTORS = load('data/N-H_vectors_1ubi.m');
VECTORS = VECTORS(:,2:4);

ANSWERS = load('answerkey.m'); 
MASTER = zeros(size(ANSWERS,1),size(VECTORS,1));
for(i=1:size(ANSWERS,1))
   MASTER(i,i) = 1;
end

[S1, flag] = updateTen(MASTER,RDC1,VECTORS);
[S2, flag] = updateTen(MASTER,RDC2,VECTORS);

S1
S2

backcomputed1 = zeros(1,size(RDC1,1));
backcomputed2 = zeros(1,size(RDC2,1));

sum1 = 0;
sum2 = 0;
ct1 = 0;
ct2 = 0;

for(i=1:size(ANSWERS))   
   v = VECTORS(i,1:3);
   if(RDC1(i)>-999 && v(1)>-999)
      backcomputed1(i) = v * S1 * v';
   end
   if(RDC2(i)>-999 && v(1)>-999)
      backcomputed2(i) = v * S2 * v';
   end
end

fid1 = fopen('out1.txt', 'w');
for(i=1:size(backcomputed1, 2))
   if(RDC1(i)>-999)
      fprintf(fid1, '%f\t%f\n', RDC1(i), backcomputed1(i));
      sum1 = sum1 + (RDC1(i)-backcomputed1(i)).^2;
      ct1 = ct1 + 1;
   end
end
fclose(fid1);

fid2 = fopen('out2.txt', 'w');
for(i=1:size(backcomputed2, 2))
   if(RDC2(i)>-999)
      fprintf(fid2, '%f\t%f\n', RDC2(i), backcomputed2(i)); 
      sum2 = sum2 + (RDC2(i)-backcomputed2(i)).^2; 
      ct2 = ct2 + 1;
   end 
end
fclose(fid2);

RMSD1 = nthroot(sum1/ct1, 2)
RMSD2 = nthroot(sum2/ct2, 2)

