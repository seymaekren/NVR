clear
load CP.ubi.debug.mat
ALLDISTS_ubi_debug = ALLDISTS;
COLIN_ubi_debug = COLIN;
HSHIFTS_ubi_debug = HSHIFTS;
NSHIFTS_ubi_debug = NSHIFTS;
ROWIN_ubi_debug = ROWIN;
TYPES_ubi_debug = TYPES;
ASSIGNTABLE_ubi_debug = ASSIGNTABLE;
CP_ubi_debug = CP;
NOES_ubi_debug = NOES;
NTH_ubi_debug = NTH;
SSTRUCT_ubi_debug = SSTRUCT;
load CP.ubi.mat
ALLDISTS(1,:)-ALLDISTS_ubi_debug(1,:)
COLIN - COLIN_ubi_debug(1,:)
HSHIFTS_ubi_debug(1:70) - HSHIFTS
NSHIFTS_ubi_debug(1:70) - NSHIFTS
ROWIN_ubi_debug(1:70) - ROWIN
TYPES_ubi_debug 
TYPES
ASSIGNTABLE_ubi_debug(1,:) - ASSIGNTABLE(1,:)
CP_ubi_debug(1,:) - CP(1,:)
NOES_ubi_debug(1,:) - NOES(1,:)
NTH_ubi_debug-NTH
SSTRUCT_ubi_debug
SSTRUCT

