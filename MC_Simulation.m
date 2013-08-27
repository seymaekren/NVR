fprintf(1, 'clearing environment variables...\n');
clear all;

b_running1CMZ                     = 0;


cd InputFiles/
!./cleanSymbolicLinks.sh

if (b_running1CMZ)
  !./symbolicLinksFor1CMZ.sh
else
  !./symbolicLinksFor1UBQ.sh
  cd ..;
  startingMASTER = load ('1UBQ.txt');
  cd InputFiles;
end

!./checkFiles.sh
cd ..
    
    

[numPeaks, numColumns] = size(startingMASTER);
numResidues            = numColumns - 3;
startingMASTER         = startingMASTER(1:numPeaks, 1:numResidues);


dbstop if error
dbstop if warning

  if (~exist('HSQCDATA'))
  
    useOrigData      = 1;
    [HSQCDATA, NOES] = readNMR_Data2(useOrigData);
    
    [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ignoredHSQCDATA] = loaddata('InputFiles/myinput.m');
    
    peaks      = HSQCDATA(:,2:3);
    rdcs       = HSQCDATA(:,4:5);
    HDEXCHANGE = HSQCDATA(:,6) ;
    peakIDs    = HSQCDATA(:,1) ;
  end


NTH=4.8;
mu =mean(mean(ALLDISTS));
if(mu-12.9>0)
   NTH=NTH+(mu-12.9); 
end

NTH = NTH + 2.34; %this is due to 1G6J max. violation distance.

%fprintf(1, 'incrementing NTH automatically by 1.5 A\n');
%NTH = NTH + 1.5;

NTH=min(NTH,9.33);

fprintf(1, 'NTH set to %f\n', NTH);

MC_AroundCorrectAssignment(NOES, ALLDISTS, NTH, startingMASTER);


