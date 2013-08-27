function voter = initialize6Voters(differenceMatrixH_SHIFTX, differenceMatrixN_SHIFTX, differenceMatrixH_SHIFTS, differenceMatrixN_SHIFTS,differenceMatrix_RDC1,differenceMatrix_RDC2)

voter        = cell(6,1);
voter{1}     = differenceMatrixH_SHIFTX;
voter{2}     = differenceMatrixN_SHIFTX;
voter{3}     = differenceMatrixH_SHIFTS;
voter{4}     = differenceMatrixN_SHIFTS;
voter{5}     = differenceMatrix_RDC1;
voter{6}     = differenceMatrix_RDC2;