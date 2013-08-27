function runexample(type)

%runexample: This script runs one of 3 tests on the sample data provided in the distribution. 
%Input:  type: if type == 1, it will run an example of the tensor estimation algorithm
%					on rdc data from ubiquitin using the ubiquitin model 1UBI. It estimates
%					the tensors, shows them, shows the 'real' tensors, and then returns a
%					percentile of accuacy (Higher numbers are better). 
%
%					if type == 2, it will run an example of the NVR assignment algorithm. 
%					if type == 3, it will run an example of the HD homology detection algorithm. 



HSQCDATA = load('hsqcdata.m');
[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata('modeldata.m');
noes = load('NOES.m');
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


if(type==1)
   %estimate tensors
   [S1]=NVR_TENEST(HSQCDATA(:,4),VECTORS);
   [S2]=NVR_TENEST(HSQCDATA(:,5),VECTORS);
   
   T1 = load('REALTENM1');T1=T1.S;
   T2 = load('REALTENM2');T2=T2.S;
   
   ESTIMATED_TENSOR_1 = S1
   ACTUAL_TENSOR_1 = T1
   
   ESTIMATED_TENSOR_2 = S2
   ACTUAL_TENSOR_2 = T2
   
   fprintf('Percentile of accuracy medium 1: %f\n', NVR_COMP_TEN(S1, T1)*100); 
   fprintf('Percentile of accuracy medium 2: %f\n', NVR_COMP_TEN(S2, T2)*100); 
   
elseif(type==2)
   %run NVR on the ubiquitin data   
   NVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS);
   
elseif(type==3)
   %run HD on two models
   
   
   fprintf('Running Example 1: A Ubiquitin Model (1UBI) on ubiquitin data, bb RMSD = 0.6 Angstroms\n\n');
   
   cd('hd_1ubi');
   %model 1, an actual ubiqutin structure
   HD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS);
   
   cd('../hd_1h8c');
   fprintf('Running Example 2: A Homolog of ubiquitin (1H8C) on ubiquitin data, bb RMSD = 1.8 Angstroms\n\n');
   [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata('modeldata2.m');
   HD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS);
   
   cd('../hd_1esr');
   fprintf('Running Example 3: A Non-Homolog of ubiquitin (1ESR) on ubiquitin data, bb RMSD = 5.9 Angstroms\n\n');
   [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata('modeldata3.m');
   HD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS);
   cd('..');
   
   
   
end




function [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS] = loadmodeldata(name);

%[R T X Y Z  SS, HB,HX,HY,HZ] =textread(name,'%f %s  %f %f %f %s %s %f %f %f\n');
[R T X Y Z  SS, HB,HX,HY,HZ] =textread(name,'%f %s  %f %f %f %s %s %f %f %f');

VECTORS= [X,Y,Z]; 
TYPES=T;
RESNUMS=R;
SSTRUCT=SS;
HBOND = HB;

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
