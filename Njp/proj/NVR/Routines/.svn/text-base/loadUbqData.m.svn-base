function [HSQCDATA, NOES] = loadUbqData(useOrigData)

addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/Matlab/FileProcessing');

if (nargin == 0)
  useOrigData = 1;
end

if (useOrigData)
  filename = sprintf('myinput.m')
  [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, HSQCDATA] = loaddata(filename);

  unadjustedNOEs = load ('NOES.txt');
  order = load ('order.m');
  noes = readOriginalNOE_File(unadjustedNOEs, order);
else
  %use distrib data
  fprintf(1, 'error. please do not use useOrigData == 0\n');
%  assert(0);
  HSQCDATA = load('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/exampledata/hsqcdata.m');
  noes = load('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/exampledata/NOES.m');
end

NOES = zeros(size(HSQCDATA,1));

for(i=1:size(HSQCDATA,1))
   rn=HSQCDATA(i,1); 
   x = find(noes(:,1)==rn);%see if there are any NOEs for this spin system
   for(j=1:length(x))
      rn2=noes(x(j),2);%get the 2nd spin systems
      y=find(HSQCDATA(:,1)==rn2);
      if(length(y)>0)
         NOES(i,y)=1; 
         NOES(y,i)=1;
      end
   end
end
