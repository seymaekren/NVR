function runexample(type)

%runexample: This script runs one of 3 tests on the sample data provided in the distribution. 
%Input:  type: if type == 1, it will run an example of the tensor estimation algorithm
%					on rdc data from ubiquitin using the ubiquitin model 1UBI. It estimates
%					the tensors, shows them, shows the 'real' tensors, and then returns a
%					percentile of accuacy (Higher numbers are better). 
%
%					if type == 2, it will run an example of the NVR assignment algorithm. 
%					if type == 3, it will run an example of the HD homology detection algorithm. 


addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/hd')
addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/nvr')
addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/tenest')
addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/Routines');


dbstop if error
dbstop if warning


useOrigData = 0;

if (useOrigData)
  filename       = sprintf('myinput.m');
  
  [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, HSQCDATA] = loaddata(filename);

  unadjustedNOEs = load ('NOES.txt');
  order          = load ('order.m');
  noes           = adjustNOES(unadjustedNOEs, order);
%  keyboard
else
  %use distrib data

  HSQCDATA = load('hsqcdata.m');
  [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = ...
      loadmodeldata('modeldata.m');
  noes = load('NOES.m');
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



%SHIFTX_Filename   = 'SHIFTX.txt';
%SHIFTS_Filename   = 'PREDICTEDSHIFTS.m';
SHIFTX_Filename    = 'SHIFTX.m';
SHIFTS_Filename    = 'SHIFTS.m';

if (useOrigData)
  fout = fopen('resultsWithOrigData.txt', 'w');
  fprintf(1, 'check out resultsWithOrigData.txt\n'); 
% $$$    versionNumber = 11;
% $$$    filename = sprintf('HD_version1.%d.txt',versionNumber);
% $$$    fout     = fopen  (filename, 'w');
% $$$    
% $$$    fprintf(1, 'check out %s\n',filename);
% $$$   

  %SHIFTS_Filename   = 'predictedshifts.txt';
% $$$     

  %  SHIFTS_Filename   = 'PREDICTEDSHIFTS.m';

% $$$    
% $$$    [assignmentAccuracy,HD_SCORE]=HD_v1_11(HSQCDATA(:,2:3),HSQCDATA(: ...
% $$$ 						  ,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS,SHIFTS_Filename,SHIFTX_Filename)

% $$$   modeIndex = 7; modelIndex = 6;
% $$$   SHIFTX_Filename   = sprintf('SHIFT__FileGeneration/SHIFTX/Parsed/Adjusted/MySHIFTX.%d.model%d.adjusted',modeIndex, modelIndex);
% $$$     
% $$$   SHIFTS_Filename   = sprintf('SHIFT__FileGeneration/SHIFTS/Parsed/Adjusted/MySHIFTS.%d.model%d.adjusted',modeIndex, modelIndex);
% $$$     
% $$$   modelDataFilename = sprintf('ModelDataGeneration/ModelDataFiles/Mode%d.coords%d',  modeIndex, modelIndex);
  

%  SHIFTX_Filename   = '1HEZ_5M_PRO.SHIFTX.parsed.adjusted';
%  SHIFTS_Filename   = '1HEZ_5M_PRO.SHIFTS';
%  modelDataFilename = 'modelDataFile';
  
%  [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata(modelDataFilename);
     
 
 
else
  
  fout = fopen('resultsWithNewDataFormat.txt', 'w');
  fprintf(1, 'check out resultsWithNewDataFormat.txt\n'); 

  
end

 [assignmentAccuracy, HD_SCORE] = HD(HSQCDATA(:,2:3),HSQCDATA(:,4),HSQCDATA(:,5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS,SHIFTS_Filename,SHIFTX_Filename)
 

fprintf(fout, '%f %f\n', assignmentAccuracy, HD_SCORE);
fclose (fout);
keyboard
  




function [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, ...
	 HSQCDATA] = loaddata(name);


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
  if (RDC1(length(R)-numDiscard) ~= -999)
    break;
  end
  if (RDC2(length(R)-numDiscard) ~= -999) 
    break;
  end
  numDiscard = numDiscard + 1;
end

fprintf(1, 'discarding %d peaks...\n', numDiscard);

H_CS      =  H_CS      (1:length(H_CS)      - numDiscard);
NH_CS     =  NH_CS     (1:length(NH_CS)     - numDiscard);
RDC1      =  RDC1      (1:length(RDC1)      - numDiscard);
RDC2      =  RDC2      (1:length(RDC2)      - numDiscard);
MOLMOL_HB =  MOLMOL_HB (1:length(MOLMOL_HB) - numDiscard);

HSQCDATA  = [transpose(1:length(R)-numDiscard), H_CS NH_CS RDC1 RDC2 MOLMOL_HB];

VECTORS   = [X,Y,Z]; 
TYPES     =  T;
RESNUMS   =  R;
SSTRUCT   =  SS;
HBOND     =  HB;

ALLDISTS=zeros(length(HX));

for(i=1:size(ALLDISTS,1))
  
   for(j=1:size(ALLDISTS,1))
     
      ALLDISTS(i,j) = sqrt((HX(i)-HX(j)).^2+(HY(i)-HY(j)).^2+(HZ(i)-HZ(j)).^2);
      IALLDISTS(i,j) = ALLDISTS(i,j); 
      
   end
   
   if(strcmp(SSTRUCT(i), 'C')==1)
      IALLDISTS(i,:) = 5;
   end
end

function noes = adjustNOES(unadjustedNOEs, order);

[m,n] = size(unadjustedNOEs);
noes  = zeros(m,2);

for i = 1:m
  index1    = unadjustedNOEs(i,1);
  index2    = unadjustedNOEs(i,2);
  newIndex1 = find(order == index1);
  newIndex2 = find(order == index2);
  noes(i,1) = newIndex1;
  noes(i,2) = newIndex2;
end

function [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS] = loadmodeldata(name);


[R T X Y Z  SS, HB,HX,HY,HZ] =textread(name,'%f %s  %f %f %f %s %s %f %f %f');

VECTORS = [X,Y,Z]; 
TYPES   = T;
RESNUMS = R;
SSTRUCT = SS;
HBOND   = HB;

ALLDISTS=zeros(length(HX));

for(i=1:size(ALLDISTS,1))
  
   for(j=1:size(ALLDISTS,1))
      ALLDISTS(i,j)  = sqrt((HX(i)-HX(j)).^2+(HY(i)-HY(j)).^2+(HZ(i)-HZ(j)).^2);
      IALLDISTS(i,j) = ALLDISTS(i,j); 
   end   
   
   if(strcmp(SSTRUCT(i), 'C')==1)
      IALLDISTS(i,:) = 5;
   end
end
