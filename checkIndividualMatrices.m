s1 = load ('CP/individualScoringMatrix_CP.txt');
s2 = load ('HDE/individualScoringMatrix_HDE.txt');
s3 = load ('RP1/individualScoringMatrix_RP1.txt');
s4 = load ('RP2/individualScoringMatrix_RP2.txt');
s5 = load ('SSCP/individualScoringMatrix_SSCP.txt');
s6 = load ('SXCP/individualScoringMatrix_SXCP.txt');
s7 = load ('TP/individualScoringMatrix_TP.txt');


s_combined = load ('/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1UBI/WithNVR_TOCSY/WithNTH=8.80/NOEsFromMZ_ubq/TruncatingWithLargerThresholds/WithRDC_SecondRound/WithNH_CH_RDCs/combinedScoringMatrix.txt');

s = s1(2,3) + s2(2,3) + s3(2,3) + s4(2,3) + s5(2,3) + s6(2,3) + ...
    s7(2,3)
s_combined(2,3)