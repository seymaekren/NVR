%load TOCSY.hSRI
%load TOCSY.poln
%TOCSY = load ('TOCSY.m.gb1.original');
load TOCSY.txt
load order.m
[numEntries, dummyVar] = size(TOCSY);
assert (dummyVar == 4);
fid = fopen('TOCSY.m.EIN','w');
fprintf(1, 'check out TOCSY.m.EIN\n');
for i = 1:numEntries
  residueNumber = TOCSY(i,4);
  peakIndex     = find(order == residueNumber);
  if (~isempty(peakIndex))
    fprintf(fid, '%f %f %d %d\n', TOCSY(i,1),TOCSY(i,2),TOCSY(i,3),peakIndex);
  end
end
fclose(fid);