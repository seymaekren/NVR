function [voter,numVoters,S1,S2] = updateVoters( voter, ASSIGNTABLE,MASTER, RDC1,RDC2,...
				    VECTORS)

numVoters = 5;

assert(sum(sum(MASTER)) >= 5);

[S1,flag]  = updateTen(MASTER,RDC1,VECTORS); if (flag) return; end
[S2,flag]  = updateTen(MASTER,RDC2,VECTORS); if (flag) return; end

numVoters  = 7;


voter{6} = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
voter{7} = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);