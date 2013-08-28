function [HSQCDATA, NOES] = readNMR_Data2(useOrigData)

if (nargin == 0)
  useOrigData = 1;
end

if (useOrigData)
  filename = sprintf('InputFiles/myinput.m')
  [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, HSQCDATA] = loaddata(filename);
  NOE_ListWithResidueIndices = load ('InputFiles/NOES.txt');
  order                      = load ('InputFiles/order.m');
  noes                       = convertResidueIndicesToPeakIndicesInNOE_File(NOE_ListWithResidueIndices, order);
  fprintf(1, 'read %d noe constraints\n',size(noes,1));
  %fprintf(1, 'note. these should be peak, not residue indices.\n');
else
  %use distrib data
  fprintf(1, 'error. please do not use useOrigData == 0\n');
  assert(0);
  HSQCDATA = load('/home/home4/apaydin/Mist/NVR/exampledata/hsqcdata.m');
  noes     = load('/home/home4/apaydin/Mist/NVR/exampledata/NOES.m');
end

%here converting the noeList to an NOE matrix.

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
%keyboard
fprintf(1, 'numNOEs = %f\n', 0.5*(sum(sum(NOES))-sum(diag(NOES))));
%keyboard