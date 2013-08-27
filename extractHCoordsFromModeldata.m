function extractHCoordsFromModeldata

%This function was written for debugging purposes.
%The aim is to extract the proton coordinates from the model data
%files in order to compare the proton coordinates of two model data files.

addpath('/home/home4/apaydin/Mist/NVR/hd')
addpath('/home/home4/apaydin/Mist/NVR/nvr')
addpath('/home/home4/apaydin/Mist/NVR/tenest')
addpath('/home/home4/apaydin/Mist/ProteinDynamics/Matlab');

dbstop if error
dbstop if warning


useOrigData = 1;
if (useOrigData)
  filename = sprintf('input.m')
  [R, Hx,Hy,Hz] = extractHCoords(filename);
  fid = fopen('extractedProtons.txt','w');
  for i = 1:length(Hx)
    fprintf(fid, '%d %f %f %f\n', R(i),Hx(i),Hy(i),Hz(i));
  end
  fclose(fid);
end


function [R,HX,HY,HZ] = extractHCoords(name);

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
H_CS = H_CS(1:length(H_CS) - numDiscard);
NH_CS = NH_CS(1:length(NH_CS) - numDiscard);
RDC1 = RDC1 (1:length(RDC1) - numDiscard);
RDC2 = RDC2 (1:length(RDC2) - numDiscard);
MOLMOL_HB = MOLMOL_HB(1:length(MOLMOL_HB) - numDiscard);
HSQCDATA = [transpose(1:length(R)-numDiscard),H_CS NH_CS RDC1 RDC2 MOLMOL_HB];
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
