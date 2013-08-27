function [NOES, RDC1, RDC2, HDEXCHANGE, peakIDs, HSHIFTS, ...
	 NSHIFTS, HN_SHIFTS, RDCs] = readNMR_Data(useOrigData);

if (nargin == 0)
  useOrigData = 1;
end

if (useOrigData == 0)

  HSQCDATA   = load       ('/home/home4/apaydin/Mist/NVR/exampledata/hsqcdata.m');
  
  noes       = load       ('/home/home4/apaydin/Mist/NVR/exampledata/NOES.m');
  
else
  
  filename       = 'myinput.m';
  [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, HSQCDATA] = loaddata(filename);

  unadjustedNOEs = load ('NOES.txt');
  order          = load ('order.m');
  noes           = processRawNOEs(unadjustedNOEs, order);
  
end


NOES       = adjustNOEs (HSQCDATA,noes);

HN_SHIFTS  = HSQCDATA   (:,2:3);
RDCs       = HSQCDATA   (:,4:5);
HDEXCHANGE = HSQCDATA   (:,6);
peakIDs    = HSQCDATA   (:,1);

HSHIFTS    = HN_SHIFTS  (:,1);
NSHIFTS    = HN_SHIFTS  (:,2);
RDC1       = RDCs       (:,1);
RDC2       = RDCs       (:,2);
