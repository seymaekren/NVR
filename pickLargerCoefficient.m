bmrb_1CMZ = load ('OptimizationFiles/1CMZ/1DK8/maxCoefficients_BMRB.txt')
bmrb_ff2  = load ('OptimizationFiles/ff2/maxCoefficients_BMRB.txt')
bmrb_hSRI = load ('OptimizationFiles/hSRI/maxCoefficients_BMRB.txt')


shifts_1CMZ = load ('OptimizationFiles/1CMZ/1DK8/maxCoefficients_SHIFTS.txt')
shifts_ff2  = load ('OptimizationFiles/ff2/maxCoefficients_SHIFTS.txt')
shifts_hSRI = load ('OptimizationFiles/hSRI/maxCoefficients_SHIFTS.txt')

shiftx_1CMZ = load ('OptimizationFiles/1CMZ/1DK8/maxCoefficients_SHIFTX.txt')
shiftx_ff2  = load ('OptimizationFiles/ff2/maxCoefficients_SHIFTX.txt')
shiftx_hSRI = load ('OptimizationFiles/hSRI/maxCoefficients_SHIFTX.txt')


bmrb = [bmrb_1CMZ; bmrb_ff2; bmrb_hSRI];
shifts = [shifts_1CMZ; shifts_ff2; shifts_hSRI];
shiftx = [shiftx_1CMZ; shiftx_ff2; shiftx_hSRI];

max(bmrb)
max(shifts)
max(shiftx)