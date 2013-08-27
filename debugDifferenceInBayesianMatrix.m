% $$$ load bayesianMatrix_debug.mat
% $$$ bayesianMatrix_debug = bayesianMatrix;
% $$$ 
% $$$ load bayesianMatrix_debug2.mat
% $$$ bayesianMatrix_debug2 = bayesianMatrix;
% $$$ 


load originalVoters_debug.mat
SSCP_debug  = SSCP;
SXCP_debug  = SXCP;
CP_debug    = CP;
RP1_debug   = RP1;
RP2_debug   = RP2;

load originalVoters_debug2.mat
SSCP_debug2 = SSCP;
SXCP_debug2 = SXCP;
CP_debug2   = CP;
RP1_debug2  = RP1;
RP2_debug2  = RP2;

SSCP_debug(51,20)
SSCP_debug2(49,20)

SSCP_debug(1,:)
SSCP_debug2(1,:)

for i = 1:size(SSCP_debug,1)
  assert (SSCP_debug(i, 1:16) == SSCP_debug2(i,1:16);
end

assert (SSCP_debug(1:21,1:21) == SSCP_debug2(1:21,1:21));
assert (SSCP_debug(23:48,23:48) == SSCP_debug2(22:47,22:47));
% $$$ assert (SSCP_debug(48:
% $$$ 
% $$$ 1:21,1:21
% $$$ 22:47,23:48
% $$$ 48:70,50:72
% $$$ 71,22
% $$$ 72,49
% $$$ debug2,debug
% $$$ 
