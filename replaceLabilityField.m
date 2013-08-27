function replaceLabilityField(outputModelDataFilename, ...
			      inputModelDataFilename,...
			      parsedHBondFilename)


%parsedHBondFilename = '1h8c_REF.stride.parsed';
%inputModelDataFilename = '1h8c_REF.modeldata.original';
%outputModelDataFilename = '1h8c_REF.modeldata.original.strideHBonds';

dbstop if error
dbstop if warning
%keyboard


  [R T X Y Z  SS  HX HY HZ] = loadmodeldata(inputModelDataFilename);
  HBondInfo = load (parsedHBondFilename);
  adjustedHBondInfo = adjustHBondInfo(HBondInfo, R);
  writeNewModelDataFile(outputModelDataFilename, ...
			R , T ,  X , Y , Z ,...
		 SS ,  HX , HY ,...
			HZ, adjustedHBondInfo);
  





function  adjustedHBondInfo = adjustHBondInfo(HBondInfo, R);
numEntries = length(R);
adjustedHBondInfo = cell(numEntries,1);
for i = 1:numEntries
  resIndex = R(i);
  pos = find(HBondInfo(:,1) == resIndex);
  if (~isempty(pos))
    adjustedHBondInfo{i} = 'N';
  else
    adjustedHBondInfo{i} = 'Y';
  end
end


function [R,T,X,Y,Z, SS, HX,HY,HZ] = loadmodeldata(name);


[R T X Y Z  SS, HB,HX,HY,HZ] =textread(name,'%f %s  %f %f %f %s %s %f %f %f');




% $$$ function [R, T, RDC1,RDC2,X,Y,Z,H_CS,NH_CS,SS,HX_Data,HX,HY,HZ] = loaddata(name);
% $$$ 
% $$$ [R T RDC1 RDC2 X Y Z H_CS NH_CS SS, MOLMOL_hbond,HX_Data HX,HY,HZ] ...
% $$$     =textread(name,'%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

% $$$ function   writeNewModelDataFile(outputModelDataFilename, ...
% $$$ 				 R, T, RDC1, RDC2, X, Y, Z, H_CS, NH_CS, ...
% $$$ 				 SS, ...
% $$$ 				 HX_Data, HX, HY, HZ, adjustedHBondInfo)

function   writeNewModelDataFile(outputModelDataFilename, ...
				 R, T,  X, Y, Z, 				 SS, ...
				 HX, HY, HZ, adjustedHBondInfo)

fid = fopen(outputModelDataFilename, 'w');

[m] = length(R);
for i = 1:m
  fprintf(fid, '%d %s  %f %f %f %s %s %f %f %f\n', R(i), ...
	  T{i}, X(i), Y(i), Z(i), ...
	  SS{i},adjustedHBondInfo{i}, HX(i),HY(i), HZ(i));
end
fclose(fid);
