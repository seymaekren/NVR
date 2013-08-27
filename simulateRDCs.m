VECTORS_NH = load('InputFiles/N-H_vectors.m');
s1 = load ('s1Tensor.txt');
s2 = load ('s2Tensor.txt');

fid1 = fopen('N-H_medium1.m.EIN','w');
fid2 = fopen('N-H_medium2.m.EIN','w');

rdcs1 = [];
rdcs2 = [];

for i = 1:size(VECTORS_NH,1)
  v      = VECTORS_NH(i,2:4);
  resId  = VECTORS_NH(i, 1 );
  noise1 = rand
  noise2 = rand
  rdc1   = v*s1*v' + noise1;
  rdc2   = v*s2*v' + noise2;
  fprintf(fid1, '%d %f\n', resId, rdc1);
  fprintf(fid2, '%d %f\n', resId, rdc2);
  rdcs1  = [rdcs1 rdc1];
  rdcs2  = [rdcs2 rdc2];
end

fclose(fid1);
fclose(fid2);