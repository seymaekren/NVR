function prepareFilesForNVR
clear all

%---------------
%input filenames
%---------------

resonanceFilename = sprintf('parsedResonances.txt')
%resonanceFilename = sprintf('1dmb_cs_parsed.txt');
%resonanceFilename            = 'parsedResonances.txt.1AAR';
%resonanceFilename = 'Peak-CAM13C_Dieckmann-NHSQC-freq-Refined.dat';
%resonanceFilename = 'HN_chemicalShifts.txt';
%resonanceFilename = 'bmrb5480.txt.parsed';

%  protonCoordinatesFilename = sprintf('2A7O.parsedPDB');
%  protonCoordinatesFilename = sprintf('1UBQH.parsedPDB');
  protonCoordinatesFilename = sprintf('1UBQH.parsedPDB');
%  protonCoordinatesFilename = sprintf('1UD7.parsedPDB');
%  protonCoordinatesFilename = sprintf('FF2Final31_90wSCs306.parsedPDB');
%  protonCoordinatesFilename =  sprintf('withMetMethyl_noTemplate_refined.parsedPDB');
%  protonCoordinatesFilename = sprintf('2E71_firstModel.parsedPDB');
%  protonCoordinatesFilename    = sprintf('1dmb.10.model6.parsedPDB');
%  protonCoordinatesFilename    = '1AARH.parsedPDB';
%  protonCoordinatesFilename    = '2I5O.parsedPDB';
%  protonCoordinatesFilename    = '3GB1.parsedPDB';
%  protonCoordinatesFilename    = '1C05.parsedPDB';
%protonCoordinatesFilename    = '1SY9_model1.parsedPDB';

%protonCoordinatesFilename    = '1ZYM_model1_withHydrogens.parsedPDB';


%rdcInFilename1               = 'nhRdc.m.1AAR.1';
%rdcInFilename2               = 'nhRdc.m.1AAR.2';
%rdcInFilename1               = 'HN.m';
%rdcInFilename1               = 'NH.m';
%rdcInFilename2               = 'CH.m';

%rdcInFilename3              = 'ccaRdc.txt';
%rdcInFilename1               = 'N-H_medium1.m.EIN';
%rdcInFilename2               = 'N-H_medium2.m.EIN';

%----------------
%output filenames
%----------------

%answerkeyFilename    = 'answerkey.m';
%  answerkeyFilename = 'answerkey.m.hSRI';

%rdc1OutFilename      = 'N-H_medium1.m';
%rdc1OutFilename      = 'N-H_medium.m';

%  rdc2OutFilename = 'C-H_medium.m';
%  rdc2OutFilename = 'N-C_medium1.m';
%rdc2OutFilename      = 'N-H_medium2.m';

%rdc3Outfilename      = 'C-Ca_medium1.m';

%orderFilename        = 'order.m';

combinedResonancesAndProtonsFilename          = 'combinedResonancesAndProtonCoordinates.txt';

%-------------------
%computation begins.
%-------------------

[protonName protonAA_Name protonAA_Index protonX protonY protonZ] ...
    = textread(protonCoordinatesFilename, '%s %s %d %f %f %f');
  
%writeAnswerkeyFile (answerkeyFilename, protonAA_Index);

%weightRDCs  = 0;
%writeRdcFile       (rdc1OutFilename, rdcInFilename1, protonName, ...
%		    protonAA_Index, weightRDCs);

%writeRdcFile       (rdc2OutFilename, rdcInFilename2, protonName, ...
%		    protonAA_Index, weightRDCs);
  
%writeRdcFile       (rdc3OutFilename, rdcInFilename3, protonName, ...
%		    protonAA_Index, weightRDCs);

%writeOrderFile     (orderFilename, protonAA_Index);

%fprintf(1, 'wrote order file.\n');

writeCombinedResonancesAndProtonCoords (combinedResonancesAndProtonsFilename, protonAA_Index, protonAA_Name, protonX, protonY, protonZ, resonanceFilename);

function writeCombinedResonancesAndProtonCoords(combinedResonancesAndProtonsFilename, protonAA_Index, protonAA_Name, protonX, protonY, protonZ, resonanceFilename);

[resonanceAA_Index N_CS H_CS] = textread(resonanceFilename,'%d %f %f');
%[resonanceAA_Index H_CS N_CS] = textread(resonanceFilename,'%d %f %f');
%[resonanceAA_Index resonanceAA_Name H_CS N_CS] = textread(resonanceFilename,'%d %s %f %f');
%[someCharacters N_CS H_CS intensity] = textread(resonanceFilename,'%s %f %f %f');

numChemicalShifts = length(N_CS)
numAtoms = length(protonX)
fprintf(1, 'enter return to continue.\n');
keyboard



fprintf(1, 'check out %s\n',combinedResonancesAndProtonsFilename);
  
fid         = fopen(combinedResonancesAndProtonsFilename,'w');

for i = 1:length(protonAA_Name)
  relResonanceIndex = find(resonanceAA_Index == protonAA_Index(i));
  if (isempty(relResonanceIndex))
    fprintf(fid, '%d\t%s\t%d\t%d\t%f\t%f\t%f\n', protonAA_Index(i), protonAA_Name{i},-999,-999,protonX(i),protonY(i),protonZ(i));
  else
  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\n', protonAA_Index(i), protonAA_Name{i},H_CS(relResonanceIndex),N_CS(relResonanceIndex),protonX(i),protonY(i),protonZ(i));
end
%  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\n', protonAA_Index(i), protonAA_Name{i},H_CS(i),N_CS(i),protonX(i),protonY(i),protonZ(i));
end

for i = length(protonAA_Name)+1:numChemicalShifts
  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\n', -999, 'XXX', H_CS(i),N_CS(i),-999,-999,-999);
end

fclose(fid);

function writeAnswerkeyFile (answerkeyFilename, protonAA_Index)
fprintf(1, 'check out %s\n',answerkeyFilename);

fid         = fopen(answerkeyFilename,'w');

if (fid == -1)
  error ('error. cant open output file in writeAnswerkeyFile');
end

for i = 1:length(protonAA_Index)
  fprintf(fid, '%d %d\n',i, protonAA_Index(i));
end
fclose(fid);

%fprintf(1, 'type return to continue.\n');
%keyboard

function writeRdcFile       (rdc1OutFilename, rdcInFilename1, protonName, ...
			     protonAA_Index, weightRDCs);

fprintf(1, 'check out %s\n',rdc1OutFilename);

fid         = fopen(rdc1OutFilename,'w');

nhRdcs      = load (rdcInFilename1);

for i = 1:length(protonName)
  relRDC_Index = find(nhRdcs(:,1) == ...
		      protonAA_Index(i));
  if (isempty(relRDC_Index))
    fprintf(fid, '%d\t%d\n', protonAA_Index(i), -999);
  else
    if (weightRDCs)
      %fprintf(fid, '%d\t%f\n', protonAA_Index(i), chRdcs(relRDC_Index,2)*1.000);      
      %      fprintf(fid, '%d\t%f\n', protonAA_Index(i),
      %      chRdcs(relRDC_Index,2)*.491);
      fprintf(1, 'uncomment code here.\n');
      keyboard
    end
    fprintf(fid, '%d\t%f\n', protonAA_Index(i), nhRdcs(relRDC_Index,2));
  end
end

fclose(fid);



%fprintf(1, 'type return to continue.\n');
%keyboard
  
function writeOrderFile     (orderFilename, protonAA_Index);
fprintf(1, 'check out %s\n',orderFilename);
    
fid         = fopen(orderFilename,'w');
  
for i = 1:length(protonAA_Index)
  fprintf(fid, '%d\n',protonAA_Index(i));
end
fclose(fid);
%fprintf(1, 'type return to continue.\n');
%keyboard
