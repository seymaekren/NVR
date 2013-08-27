function [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS] = loadmodeldata(name);

%copied from Langmead's runexample.m


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


  
