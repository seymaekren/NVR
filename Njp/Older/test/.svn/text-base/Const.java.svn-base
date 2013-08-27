import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

/**
 * A class defining various global constants for the entire program.
 * 
 * @version 0.1
 * @author Lincong Wang
 * Wrote by Lincong Wang At Bruce's Lab at Dartmouth College
 */
public class Const {
    static final double cst = Math.PI / 180.0; //conversion constant for angle
    static final double eps = 1.0E-10;         //for comparing two float numbers
    //For the definition of these theta angles see our paper
    static final double theta1 = 29.1406 * cst;  //Angles mined from 23 untraHigh resolution X-ray structures.
    static final double theta3 = -20.935 * cst;  //theta3 should be minus c.f. its definition
    static final double theta5 = 26.8237 * cst;
    static final double alpha3 = - theta3;
    static final double theta6 = -0.7518 * cst;
    static final double theta7 = 29.1044 * cst;
    static final double alpha1 = 0.50 * Math.PI - theta7;
    static final double theta9 = -0.0021 * cst;
    static final double theta8 =0.00204369051213495*cst;//an magic angle for computing precisely the coordinates 
    static final double alpha2 =  0.50 * Math.PI - theta5;
    static final double alpha5 = theta1;
    static final double alpha6 = theta6;
    static final double alpha8 = theta9;
    static final double alpha9 = -0.0020436905160222014 * cst;
    static final double[] sevenAngles = new double[]{theta1,theta9,theta3, theta5, theta6, theta7, theta8};
    static final double dN2CA  = 1.458;                    //bond lengths for backbone vectors
    static final double dCA2CO = 1.525;
    static final double dCO2N  = 1.329;
    static final double dN2H   = 1.020;  //1.020 before
    static final double dCO2O  = 1.231;
    static final double dCA2HA = 1.090;
    static final double dCA2CA = 3.8100448409449443;//3.8102585214129507; //the distance between CA(i) and CA(i+1)
    static final double dCA2Nseq = 2.4333047472333475;//; //the distance between CA(i) and N(i+1);
    static final double dCA2NHseq = 2.556836222174780;//; //the distance between CA(i) and NH(i+1);

//     static final double dCO2NH = 2.055; //NH2CO, computed with interAngle = 119.1044;
    static final double dCA2CB = 1.531;
    static final double dCO2NH = Math.sqrt(dCO2N * dCO2N + dN2H * dN2H + 2.0 * dCO2N * dN2H * Math.sin(theta7));
    //the plane angle between CO(i)->NH(i+1) vector and CO(i)->N(i+1) vector
    static final double theta10 = - Math.asin(dN2H * Math.sin(theta7+0.50*Math.PI) / dCO2NH ); // is MINUS
    static final double dN2COHBond   = 2.88;               //the distance across the H-bond
    static final double angleNCACOHA = 118.95951 * cst * 0.975; //angles for computing HA and CO atoms
    static final double angleNCAHA   = -19.4 * cst; // -19.4 * cst; 
    static final double deltaHA   = Math.PI - angleNCACOHA; // the dihedral angle
    static final double thetaHA   = -angleNCAHA;  //19.4 * cst; // the plane angle
    static final double angleCACOO   = 59.329 * cst;  
    static final double angleNCaCoCb = (360.0-122.6767595) * cst; //angles for computing HA and CO atoms
    static final double angleNCaCb   = -20.1866 * cst; // -19.4 * cst; 
    static final double deltaCB   = Math.PI - angleNCaCoCb; // the dihedral angle
    static final double thetaCB    = -angleNCaCb; // the plane angle
    static final double deltaO   = -0.7518 * cst; // the dihedral angle
    static final double thetaO   = 32.50530 * cst; // the plane angle

    //Relative RDC strength
    static final double protonGyroRatio = 26.75;
    static final double carbonGyroRatio = 6.73;
    static final double nitrogenGyroRatio = 2.71;
    static final double noeLowerLimit = 1.90;
    static final double noeUpLimit = 5.01;
    static final double dnaLimit = 3.00;
    static final double danLimit = 3.59122;
    static final double dnnLimit = 4.50;
    static final double noDataIndicator = -999.9; //a special number to indicate that no experiment data available
    static final double noSolutionIndicator = -9999.99; //a special number to indicate that no solutions
    //the relative strength of DD between different nuclei
    //static final double cahaRatio = 1.00;  //use as standard, 
    
    static final double cahaRatio = 1.00;  //use as standard, 
    //static final double nhRatio   = Math.pow(dCA2HA /dN2H, 3) * (nitrogenGyroRatio /carbonGyroRatio);
    static final double nhRatio   = 1.0;//-0.4785;//changed here for eta(based on the scaling factor)
  
    
   // static final double nhRatio   = Math.pow(dCA2HA /dN2H, 3) * (nitrogenGyroRatio /carbonGyroRatio);
    static final double conRatio  = Math.pow(dCA2HA /dCO2N,3) * (nitrogenGyroRatio /protonGyroRatio);
    static final double conhRatio = Math.pow(dCA2HA /dCO2NH,3)* (nitrogenGyroRatio /carbonGyroRatio);
    static final double cacoRatio =1.00; //Math.pow(rN2H / rCA2CO, 3)*(c13GyroRatio * c13GyroRatio / (n15GyroRatio * h1GyroRatio ));

    //Matrices used only for testing the degeneracy and computing orientations. 
    static Matrix mat = new Matrix(3,3);   
    static final double [][] signInv = {{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},{0.0, 0.0, -1.0}};//invert coordinates
    static final double [][] xyInv   = {{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},{0.0, 0.0,  1.0}}; //invert x,y direction
    static final double [][] xzInv   = {{-1.0, 0.0, 0.0},{0.0,  1.0, 0.0},{0.0, 0.0, -1.0}}; //invert x,z direction
    static final double [][] yzInv   = {{ 1.0, 0.0, 0.0},{0.0, -1.0, 0.0},{0.0, 0.0, -1.0}}; //invert y,z direction
    static final double [][] xInv    = {{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},{0.0, 0.0,  1.0}}; //invert x direction
    static final double [][] yInv    = {{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},{0.0, 0.0,  1.0}}; //invert Y direction
    static final double [][] zInv    = {{1.0, 0.0, 0.0},  {0.0, 1.0, 0.0},{0.0, 0.0, -1.0}}; //invert Z direction

    static final Matrix zMat  = new Matrix(zInv);
    static final Matrix xyMat = new Matrix(xyInv);
    static final Matrix xzMat = new Matrix(xzInv);
    static final Matrix yzMat = new Matrix(yzInv);
    static Matrix[] mat4ThreeDirs = new Matrix[]{xyMat, xzMat, yzMat};
    static final Matrix signMat = new Matrix(signInv);
    static final Matrix r1x = mat.rotationMat(theta1, "+x"); 
    static final Matrix r9y = mat.rotationMat(theta9, "+y"); 
    static final Matrix r9y1x = r9y.times(r1x);
    static final Matrix r8z = mat.rotationMat(theta8, "+z"); 
    static final Matrix r3x = mat.rotationMat(theta3, "+x"); 
    static final Matrix r5x = mat.rotationMat(theta5, "+x");
    static final Matrix r6y = mat.rotationMat(theta6, "+y"); 
    static final Matrix r7x = mat.rotationMat(theta7, "+x");
    static final Matrix r10x = mat.rotationMat(theta10, "+x");
    static final Matrix r7x6y5x = r7x.times(r6y.times(r5x));

    static final Matrix r1xAlpha = mat.rotationMat(alpha1, "+x");
    static final Matrix r2xAlpha = mat.rotationMat(alpha2, "+x");  
    static final Matrix r3xAlpha = mat.rotationMat(alpha3, "+x");
    static final Matrix r5xAlpha = mat.rotationMat(alpha5, "+x");
    static final Matrix r6zAlpha = mat.rotationMat(alpha6, "+z"); 
    static final Matrix r8yAlpha = mat.rotationMat(alpha8, "+y"); 
    static final Matrix r9zAlpha = mat.rotationMat(alpha9, "+z"); 
    static final Matrix r8y3xAlpha = r8yAlpha.times(r3xAlpha); 
    static final Matrix r2x6z1x9zAlpha = r2xAlpha.times(r6zAlpha.times(r1xAlpha.times(r9zAlpha))); 

    //Rotation Matrices used to compute HA and O coordinates
    static final Matrix rHA1 = mat.rotationMat(angleNCACOHA,"+y");//theta_{41}=angleNCACOHA=PhiPsiCalPDB(caToCOVec,nToCAVec,caToHAVec);
    static final Matrix rHA2 = mat.rotationMat(angleNCAHA,  "+x");//theta_{42}=angleNCAHA=interAngle(nToCAVec,caToHAVec)-Math.PI/2;
    static final Matrix rHA2HA1 = rHA2.times(rHA1);
    static final Matrix rO   = mat.rotationMat(angleCACOO,  "+x");//angleCACOO = interAngle(caToCOVec,coToOVec);

    static final Matrix matDeltaHA = mat.rotationMat(deltaHA,"+y");//theta_{41}=angleNCACOHA=PhiPsiCalPDB(caToCOVec,nToCAVec,caToHAVec);
    static final Matrix matThetaHA = mat.rotationMat(thetaHA, "+x");//theta_{42}=angleNCAHA=interAngle(nToCAVec,caToHAVec)-Math.PI/2;
    static final Matrix matHA = matThetaHA.times(matDeltaHA);

    //Rotation Matrices used to compute CB
    static final Matrix rCb1 =mat.rotationMat(angleNCaCoCb,"+y");//theta_{41}=angleNCACOHA=PhiPsiCalPDB(caToCOVec,nToCAVec,caToHAVec);
    static final Matrix rCb2 =mat.rotationMat(angleNCaCb,  "+x");//theta_{42}=angleNCAHA=interAngle(nToCAVec,caToHAVec)-Math.PI/2;
    static final Matrix rCb2Cb1 = rCb2.times(rCb1);
    static final Matrix matDeltaCB = mat.rotationMat(deltaCB,"+y");//theta_{41}=angleNCACOHA=PhiPsiCalPDB(caToCOVec,nToCAVec,caToHAVec);
    static final Matrix matThetaCB = mat.rotationMat(thetaCB, "+x");//theta_{42}=angleNCAHA=interAngle(nToCAVec,caToHAVec)-Math.PI/2;
    static final Matrix matDeltaO = mat.rotationMat(deltaO,"+z");
    static final Matrix matThetaO = mat.rotationMat(thetaO, "+x");
    static final Matrix matO = matThetaO.times(matDeltaO);

    //Inverse matrices
    static final Matrix rCb1Inv = rCb1.transpose();
    static final Matrix rCb2Inv = rCb2.transpose();
    static final Matrix rCb2Cb1Inv = rCb1Inv.times(rCb2Inv);
    static final Matrix rHA1Inv = rHA1.transpose();
    static final Matrix rHA2Inv = rHA2.transpose();
    static final Matrix rHA2HA1Inv = rHA1Inv.times(rHA2Inv);
    static final Matrix rOInv   = rO.transpose();    
    static final Matrix matDeltaHAInv = matDeltaHA.transpose();
    static final Matrix matThetaHAInv = matThetaHA.transpose();
    static final Matrix matHAInv = matDeltaHAInv.times(matThetaHAInv);
    static final Matrix matDeltaCBInv = matDeltaCB.transpose();
    static final Matrix matThetaCBInv = matThetaCB.transpose();
    static final Matrix matCBInv = matDeltaCBInv.times(matThetaCBInv);

    static final Matrix matDeltaOInv = matDeltaO.transpose();
    static final Matrix matThetaOInv = matThetaO.transpose();
    static final Matrix matOInv = matDeltaOInv.times(matThetaOInv);

    static final Matrix r1xInv  = r1x.transpose();
    static final Matrix r9yInv  = r9y.transpose();
    static final Matrix r1x9yInv  = r1xInv.times(r9yInv);
    static final Matrix r8zInv  = r8z.transpose();
    static final Matrix r8z1xInv  = r8zInv.times(r1xInv); //Const.r8zInv.times(Const.r1xInv.times(coordCA)));
    static final Matrix r3xInv  = r3x.transpose();
    static final Matrix r1x9y3xInv = r1xInv.times(r9yInv.times(r3xInv));
    static final Matrix r5xInv  = r5x.transpose();
    static final Matrix r6yInv  = r6y.transpose();
    static final Matrix r7xInv  = r7x.transpose();
    static final Matrix r10xInv = r10x.transpose();
    static final Matrix r7x6yInv = r6yInv.times(r7xInv);
    static final Matrix r7x6y5xInv = r5xInv.times(r6yInv.times(r7xInv)); //for computing NH direction
    static final Matrix r6y7x8z1xInv = r7x6yInv.times(r8zInv.times(r1xInv)); 
    static final Matrix r1x8z7x6y5xInv = r7x6y5xInv.times(r8zInv.times(r1xInv)); //for computing NCA direction

    static final Matrix r1xAlphaInv = r1xAlpha.transpose();
    static final Matrix r3xAlphaInv = r3xAlpha.transpose();
    static final Matrix r5xAlphaInv = r5xAlpha.transpose();
    static final Matrix r2xAlphaInv = r2xAlpha.transpose();
    static final Matrix r6zAlphaInv = r6zAlpha.transpose();
    static final Matrix r8yAlphaInv = r8yAlpha.transpose();
    static final Matrix r9zAlphaInv = r9zAlpha.transpose();
    static final Matrix r9z1xAlphaInv = r9zAlphaInv.times(r1xAlphaInv);
    static final Matrix r1x2z3xAlphaInv = r1xAlphaInv.times(r2xAlphaInv.times(r3xAlphaInv));
    static final Matrix matL    = r1xInv.times(r9yInv);

    /** Vectors for convenience and speed **/
    static final double [] ca2coVec = {0.0, 0.0, dCA2CO}; //constants for computing phi,psi from CA direction
    static final double [] dirA  = r3xInv.times(ca2coVec);
    static final double [] co2nVec = {0.0, dCO2N, 0.0};
    static final double [] dirB1 = r5xInv.times(co2nVec);
    static final double [] coordCA = {0.0, dN2CA, 0.0};
    static final double [] dirB2 = r1x8z7x6y5xInv.times(coordCA);
    static final double [] dirB  = {dirB1[0] + dirB2[0],  dirB1[1] + dirB2[1],  dirB1[2] + dirB2[2]};
    static final double  cosPsiC = dirB[0] / Math.sqrt(dirB[0] * dirB[0] + dirB[1] * dirB[1]); //one number later
    static final double  sinPsiC = dirB[1] / Math.sqrt(dirB[0] * dirB[0] + dirB[1] * dirB[1]);
//     double  psiC1 = 0.0;
//     if ( cosPsiC >= 0.0 )
// 	psiC1 = Math.asin(sinPsiC);
//     else if ( cosPsi < 0)
// 	psiC1 = Math.PI-Math.asin(sinPsiC);
    static final double psiC = Math.PI - Math.asin(sinPsiC); //since ( cosPsi < 0)
    /** END ***/
    static final double [] dirCAcnt = r8z1xInv.times(coordCA); //a constant needed for the computing nNhCaCal
    static final double [] coordNH = {0.0, 0.0, -dN2H};       
    static final double [] coordNHCnt = {0.0, 0.0, -1.0};       
    static final double [] matRNH  = r5xInv.times(r6yInv.times(r7xInv.times(coordNHCnt)));
    static final double [] matR = r5xInv.times(r6yInv.times(r7xInv.times(r8zInv.times(r1xInv.times(coordCA)))));
    static final double [] unitZ   = new double[]{0, 0, 1};  //For computing the Matrix $M$ in the MATH file
    static final double [] dirCosCHcnt = rHA1Inv.times(rHA2Inv.times(unitZ));
    static final double [] unitMinusZ  = new double[]{0.0, 0.0, -1.0};
    static final double [] unitN2CA  = {0.0, Math.cos(theta1), Math.sin(theta1)};
    static final double [] dirNHcnt = r7x6yInv.times(coordNH);
    static final double [] dirCosNHcnt = r5xInv.times(r6yInv.times(r7xInv.times(unitMinusZ)));
    static final double [] unitCo2Nh = {0.0, dCO2NH, 0.0};
    static final double [] dirCO2NHcnt = r5xInv.times(r6yInv.times(r10xInv.times(unitCo2Nh)));
    static final double [] coordHA = {0.0, 0.0, dCA2HA};
    static final double [] cntHAVec = rHA2HA1Inv.times(coordHA);
    static final double [] phiCnt = {0.0, Math.sin(alpha5), Math.cos(alpha5)};  //a const vector for computing psi by backwards
    static final double [] psiCnt = r3xAlphaInv.times(matHAInv.times(coordHA));  //a const vector for computing psi by backwards
    static final double [] psiCntUnit = r3xAlphaInv.times(matHAInv.times(unitZ));//a const vector for computing psi by backwards

    // A set of constants for computing \phi and \psi angles from dan, dna
    static final double Cx = dCA2HA * dirCosCHcnt[0]; 
    static final double Cy = dCA2HA * dirCosCHcnt[1];
    static final double Cz = dCA2HA * dirCosCHcnt[2];
    static final double [][] phiMat = r1x9yInv.getArray();
    static final double a1 = Cz * phiMat[0][0] - Cx * phiMat[0][2];
    static final double a2 = Cz * phiMat[1][0] - Cx * phiMat[1][2];
    static final double a3 = Cz * phiMat[2][0] - Cx * phiMat[2][2];
    static final double b1 = Cx * phiMat[0][0] + Cz * phiMat[0][2];
    static final double b2 = Cx * phiMat[1][0] + Cz * phiMat[1][2];
    static final double b3 = Cx * phiMat[2][0] + Cz * phiMat[2][2];
    static final double c1 = Cy * phiMat[0][1];
    static final double c2 = Cy * phiMat[1][1] + dN2CA * Math.cos(theta1);
    static final double c3 = Cy * phiMat[2][1] + dN2CA * Math.sin(theta1) + dN2H;
    static final double e0 = Cx * Cx + Cz * Cz + c1*c1 + c2* c2 + c3 * c3;
    static final double e1 = 2* (b1 * c1 + b2 * c2 + b3 * c3);
    static final double e2 = 2* (a1 * c1 + a2 * c2 + a3 * c3);

    static final double [] caToHAVec = r1x9yInv.times(cntHAVec);
    static final double [] ha ={caToHAVec[0], caToHAVec[1] + dN2CA * Math.cos(theta1), caToHAVec[2] + dN2CA * Math.sin(theta1)};
    static final double [] caToCOVec = r1x9y3xInv.times(ca2coVec);
    static final double [] co ={caToCOVec[0], caToCOVec[1] + dN2CA * Math.cos(theta1), caToCOVec[2] + dN2CA * Math.sin(theta1)};
    static final double [][] psiMat = r1x9y3xInv.getArray();
    static final double Cxx = dirCO2NHcnt[0];
    static final double Cyy = dirCO2NHcnt[1];
    static final double Czz = dirCO2NHcnt[2];
    static final double a11 = Cyy * psiMat[0][0] - Cxx * psiMat[0][1];
    static final double a22 = Cyy * psiMat[1][0] - Cxx * psiMat[1][1];
    static final double a33 = Cyy * psiMat[2][0] - Cxx * psiMat[2][1];
    static final double b11 = -Cxx * psiMat[0][0] - Cyy * psiMat[0][1];
    static final double b22 = -Cxx * psiMat[1][0] - Cyy * psiMat[1][1];
    static final double b33 = -Cxx * psiMat[2][0] - Cyy * psiMat[2][1];
    static final double c11 = Czz * psiMat[0][2] + co[0] - ha[0];
    static final double c22 = Czz * psiMat[1][2] + co[1] - ha[1];
    static final double c33 = Czz * psiMat[2][2] + co[2] - ha[2];
    static final double ee0 = Cxx * Cxx + Cyy * Cyy + c11*c11 + c22* c22 + c33 * c33;
    static final double ee1 = 2* (b11 * c11 + b22 * c22 + b33 * c33);
    static final double ee2 = 2* (a11 * c11 + a22 * c22 + a33 * c33);

    /**
     * A set of constants for computing the side-chain conformation of Proline
     * These constants are just back-computed from Pro19 in 1UBQ,
     */
    static final double dCa2Cb = 1.49;
    static final double dC2C   = 1.54;
    static final double thetaProCB  =  1.3317861674705895;
    static final double thetaProCBd = -2.093056454581993;
    static final Matrix rxCBInv  = mat.rotationMat(thetaProCB,  "-x");
    static final Matrix ryCBdInv = mat.rotationMat(thetaProCBd, "-y");
    static final Matrix rCBInv   = ryCBdInv.times(rxCBInv);

    static final double cbTheta1 = 1.3469915625340805;
    static final double cbTheta2 = 1.1816623139143054;
    static final double cbTheta3 = 1.2351203414182543;
    static final double chi1 = -0.5042748429337971;
    static final double cbTheta2d = 2.112251747644531;
    static final double cbTheta3d = -2.0196971725548947; 
    static final sp3Carbon proCB = new sp3Carbon(cbTheta1,cbTheta2,cbTheta3,chi1,cbTheta2d,cbTheta3d,dCa2Cb,dCA2HA,dCA2HA,false);

    static final double cgTheta1 = 1.3414373700743338; 
    static final double cgTheta2 = 1.2393003417313453;
    static final double cgTheta3 = 1.1826139579031945;
    static final double chi2 =     0.7316244067220934;
    static final double cgTheta2d = 2.023077294407706;
    static final double cgTheta3d = -2.1117884874529365; 
    static final sp3Carbon proCG = new sp3Carbon(cgTheta1,cgTheta2,cgTheta3,chi2,cgTheta2d,cgTheta3d,dC2C,dCA2HA,dCA2HA,false);

    static final double cdTheta1 = 1.3843470838350387;
    static final double cdTheta2 = 1.1883376662609455;
    static final double cdTheta3 = 1.2090038592514238;
    static final double chi3 = -0.6692454232754876;
    static final double cdTheta2d =  2.0937261165122902;
    static final double cdTheta3d = -2.0324030990577624; 
    static final sp3Carbon proCD = new sp3Carbon(cdTheta1,cdTheta2,cdTheta3,chi3,cdTheta2d,cdTheta3d,dC2C,dCA2HA,dCA2HA, true);
    //A set of constants mined from ultra-resolution PDBs from CB to other atoms of the same residue
    static final double cb2HG    = 2.06;
    static final double cb2HDLow = 2.71;
    static final double cb2HDUp  = 3.35;
    static final double cb2HD  = 3.35;
    static final double cb2HE  = 4.00;
    static final double cb2HDPhe = 2.67;
    static final double cb2HEPhe = 4.55;
    static final double cb2HZPhe = 5.00;
    static final double cb2HDTyr = 2.67;
    static final double cb2HETyr = 4.55;

    /**
     * A correction for methyl group
     */
    static final double methylCorrection = 2.5; //0.5;//2.50;changed by zeng
    /** A set of constants for computing the backbone \phi and \psi
     * angles when the orientation of the two peptide planes i, i+2
     * are known.
     */
    static final Matrix rzPiInv = mat.rotationMat(Math.PI, "-z"); 
    static final Matrix rot1 = rzPiInv.times(r7x6y5xInv.times(r8zInv)); 
    static final double[] n2caUnit = {0.0, Math.cos(theta1), Math.sin(theta1)};
    static final double[] Cu = rot1.times(n2caUnit);  //for n2ca vector
    static final double[] Cv = rot1.times(unitMinusZ);//for n2nh vector
    static final Matrix rot1Inv = rot1.transpose();
    static final Matrix matA = rot1.times(r1x9yInv);

    static final double[][] arrLeftHand = new double[][]{{-1, 0, 0}, {0,1,0},{0,0,1}}; //for changing  handedness
    static final Matrix mLeftHand = new Matrix(arrLeftHand);

    static final double phiPro = -65.0 * cst;
    static final Matrix r2yPro =  mat.rotationMat(phiPro, "+y"); 
    static final Matrix r2yInvPro =  r2yPro.transpose(); 

    //The average  Phi/Psi for a typical helix or strand
    static final double phiAveHelix =  -65.3 * cst;  
    static final double psiAveHelix =  -39.4 * cst;
    static final double phiAveBeta  = -120.0 * cst; 
    static final double psiAveBeta  =  138.0 * cst; 

    //The favorable Ramachandran (Ram) regions
    static final double phiHighHelix =  -30.0 * cst;  //helix region
    static final double phiLowHelix  = -100.0 * cst;
    static final double psiHighHelix =  -15.0 * cst;
    static final double psiLowHelix  =  -90.0 * cst;
//     static final double phiHighHelix =  -55.0 * cst;  //helix region
//     static final double phiLowHelix  = -75.0 * cst;
//     static final double psiHighHelix =  -30.0 * cst;
//     static final double psiLowHelix  =  -50.0 * cst;

    static final double phiHighBeta =  -70.0 * cst;  //Beta region
    static final double phiLowBeta  = -170.0 * cst;
    static final double psiHighBeta =  180.0 * cst;
    static final double psiLowBeta  =   80.0 * cst;
    
    
    //added by zeng for ring resonance assignments
    //Phe:
    static final double PheCSHD1Min= 6.02 ;
    static final double PheCSHD1Max= 8.08 ;
    
    static final double PheCSHD2Min=   6.02 ;
    static final double PheCSHD2Max= 8.15 ;   
    
    static final double PheCSHE1Min=  6.02 ;
    static final double PheCSHE1Max= 8.80;
    
    static final double PheCSHE2Min=  6.02  ;
    static final double PheCSHE2Max= 8.80;
    
    static final double PheCSHZMin=  6.02 ;
    static final double PheCSHZMax= 9.50 ;
    
    //His:
    static final double HisCSHD1Min=   6.02   ;
    static final double HisCSHD1Max= 17.20 ;
    
    static final double HisCSHD2Min=   6.02 ;
    static final double HisCSHD2Max= 9.01  ;   
    
    static final double HisCSHE1Min=   6.02   ;
    static final double HisCSHE1Max= 10.26   ;
    
    static final double HisCSHE2Min=  6.62;
    static final double HisCSHE2Max=16.53    ;  
    
    //Trp:
    static final double TrpCSHD1Min=   6.02;
    static final double TrpCSHD1Max= 8.93;  
    
    static final double TrpCSHE1Min=   6.02  ;
    static final double TrpCSHE1Max= 13.29;  
    
    static final double TrpCSHE3Min=   6.02;
    static final double TrpCSHE3Max=  8.98;  
    
    static final double TrpCSHZ2Min=   6.02  ;
    static final double TrpCSHZ2Max=  8.50 ;  
    
    static final double TrpCSHZ3Min=   6.02  ;
    static final double TrpCSHZ3Max=  8.90   ; 
    
    static final double TrpCSHH2Min=    6.02   ;
    static final double TrpCSHH2Max=  10.90  ; 
    
    //Tyr:
    static final double TyrCSHD1Min=   6.02;
    static final double TyrCSHD1Max= 8.53;  
    
    static final double TyrCSHD2Min=   6.02;
    static final double TyrCSHD2Max=10.50;
    
    static final double TyrCSHE1Min=   6.02  ;
    static final double TyrCSHE1Max= 7.86;
    
    static final double TyrCSHE2Min=   6.02;
    static final double TyrCSHE2Max=8.50;
    
    static final double TyrCSHHMin=   6.02;
    static final double TyrCSHHMax=13.75;
    
    ////////////////////////////////////////////////////
    //for simulating more 12C from aliphatic side-chains
    
    //Ile
    static final double IleCSHG12Min=    -2.02;
    static final double IleCSHG12Max= 2.69;
    
    static final double IleCSHG13Min=   -2.04;       
    static final double IleCSHG13Max=  2.99;
    
    //Gln
    static final double GlnCSHG2Min=    0.04;      
    static final double GlnCSHG2Max=   3.66 ;
    
    static final double GlnCSHG3Min=    0.05;      
    static final double GlnCSHG3Max=   3.66 ;
    
    //Lys
    static final double LysCSHG2Min=   -0.82;              
    static final double LysCSHG2Max=   3.01;
    
    static final double LysCSHG3Min=    -0.87;          
    static final double LysCSHG3Max=   2.99 ;
    
    static final double LysCSHE2Min=   1.17;                      
    static final double LysCSHE2Max=   4.32;
    
    static final double LysCSHE3Min=   1.10;                
    static final double LysCSHE3Max=    4.55 ;
    
    
    //////////////////////////////////////////////////////////
    //chemical shifts of ring protons from BMRB, added by zeng
    //Phe4:
    static final double Phe4HD1= 7.11 ;
    static final double Phe4HD2= 7.11 ;
    
    static final double Phe4HE1=  7.28     ;
    static final double Phe4HE2= 7.28     ;   
    
    static final double Phe4HZ= 7.28     ;
    
    //Phe45:
    static final double Phe45HD1= 7.39     ;
    static final double Phe45HD2= 7.39     ;
    
    static final double Phe45HE1=  7.56        ;
    static final double Phe45HE2= 7.56         ;   
    
    static final double Phe45HZ= 7.5         ;
    
    
    //His68:
    static final double His68HD1=  -9999.99; //unavailable from bmrb
    static final double His68HD2=  7.15       ;
    static final double His68HE1= 8.34     ;
    static final double His68HE2=  -9999.99; //unavailable from bmrb
    
    //Tyr59:
    static final double Tyr59HD1=  7.28    ;
    static final double Tyr59HD2=  7.28    ;
    static final double Tyr59HE1=  6.92;
    static final double Tyr59HE2=  6.92;
    /**
     * The configuration of each amino acid: bond information
     * Ling's old naming scheme
     */  /*  
   static final Vector backboneBonds(){
	Vector allBonds = new Vector();
	allBonds.add("H_N");
	allBonds.add("N_H");
	allBonds.add("N_CA");
	allBonds.add("CA_N");
	allBonds.add("CA_HA");
	allBonds.add("HA_CA");
	allBonds.add("CA_C");
	allBonds.add("C_CA");
	allBonds.add("CA_CB");
	allBonds.add("CB_CA");
	allBonds.add("C_O");
	allBonds.add("O_C");
	return allBonds;
    }
   static final Vector glyBonds(){
       Vector allBonds = new Vector();
	allBonds.add("H_N");
	allBonds.add("N_H");
	allBonds.add("N_CA");
	allBonds.add("CA_N");
	allBonds.add("CA_C");
	allBonds.add("C_CA");

	allBonds.add("CA_2HA");
	allBonds.add("2HA_CA");
	allBonds.add("CA_1HA");
	allBonds.add("1HA_CA");
	allBonds.add("C_O");
	allBonds.add("O_C");
	allBonds.add("2HA_1HA");
	allBonds.add("1HA_2HA");
	Collections.sort(allBonds);
	return allBonds;
    }
	
    static final Vector alaBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CB_3HB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("3HB_CB");

	allBonds.add("1HB_2HB");  //We do not compute the internal vdw violation for methyl group
	allBonds.add("1HB_3HB");
	allBonds.add("2HB_3HB");
	allBonds.add("2HB_1HB");
	allBonds.add("3HB_1HB");
	allBonds.add("3HB_2HB");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector argBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CG_CB");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CG_1HG");
	allBonds.add("CG_2HG");
	allBonds.add("1HG_CG");
	allBonds.add("2HG_CG");
	allBonds.add("CD_CG");
	allBonds.add("CG_CD");
	allBonds.add("CD_1HD");
	allBonds.add("CD_2HD");
	allBonds.add("1HD_CD");
	allBonds.add("2HD_CD");
	allBonds.add("CD_NE");
	allBonds.add("NE_CD");
	allBonds.add("HE_NE");
	allBonds.add("NE_HE");
	allBonds.add("CZ_NE");
	allBonds.add("NE_CZ");
	allBonds.add("CZ_NH1");
	allBonds.add("NH1_CZ");
	allBonds.add("CZ_NH2");
	allBonds.add("NH2_CZ");
	allBonds.add("NH1_1HH1");
	allBonds.add("NH1_2HH1");
	allBonds.add("1HH1_NH1");
	allBonds.add("2HH1_NH1");
	allBonds.add("NH2_1HH2");
	allBonds.add("NH2_2HH2");
	allBonds.add("1HH2_NH2");
	allBonds.add("2HH2_NH2");
	//These are inherent vdw violations
	allBonds.add("1HB_2HB");  //do not compute the  vdw violation for methylene group
	allBonds.add("2HB_1HB");  //do not compute the  vdw violation for methylene group
	allBonds.add("1HG_2HG");  //do not compute the  vdw violation for methylene group
	allBonds.add("2HG_1HG");  //do not compute the  vdw violation for methylene group
	allBonds.add("1HD_2HD");  //do not compute the  vdw violation for methylene group
	allBonds.add("2HD_1HD");  //do not compute the  vdw violation for methylene group
	allBonds.add("2HH1_1HH1");
	allBonds.add("1HH1_2HH1");
	allBonds.add("2HH2_1HH2");
	allBonds.add("1HH2_2HH2");
	allBonds.add("CZ_HE");
	allBonds.add("CZ_1HH1");
	allBonds.add("CZ_2HH1");
	allBonds.add("CZ_1HH2");
	allBonds.add("CZ_2HH2");
	allBonds.add("HE_CZ");
	allBonds.add("1HH1_CZ");
	allBonds.add("2HH1_CZ");
	allBonds.add("1HH2_CZ");
	allBonds.add("2HH2_CZ");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector proBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CG_CB");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CG_1HG");
	allBonds.add("CG_2HG");
	allBonds.add("1HG_CG");
	allBonds.add("2HG_CG");
	allBonds.add("CD_CG");
	allBonds.add("CG_CD");
	allBonds.add("CD_1HD");
	allBonds.add("CD_2HD");
	allBonds.add("1HD_CD");
	allBonds.add("2HD_CD");
	allBonds.add("CD_N");
	allBonds.add("N_CD");
	//These are inherent vdw violations
	allBonds.add("1HB_2HB");  
	allBonds.add("2HB_1HB");  
	allBonds.add("1HG_2HG");  
	allBonds.add("2HG_1HG");  
	allBonds.add("1HD_2HD");  
	allBonds.add("2HD_1HD");  
	allBonds.add("CA_CG");  
	allBonds.add("CG_CA");  
	allBonds.add("CB_CD");  
	allBonds.add("CD_CB");  
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector asnBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_OD1");
	allBonds.add("CG_ND2");
	allBonds.add("ND2_1HD2");
	allBonds.add("ND2_2HD2");

	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("OD1_CG");
	allBonds.add("ND2_CG");
	allBonds.add("1HD2_ND2");
	allBonds.add("2HD2_ND2");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector aspBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_OD1");
	allBonds.add("CG_OD2");

	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("OD1_CG");
	allBonds.add("OD2_CG");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector glnBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_CD");
	allBonds.add("CG_1HG");
	allBonds.add("CG_2HG");
	allBonds.add("CD_OE1");
	allBonds.add("CD_NE2");
	allBonds.add("NE2_1HE2");
	allBonds.add("NE2_2HE2");
	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CD_CG");
	allBonds.add("1HG_CG");
	allBonds.add("2HG_CG");
	allBonds.add("OE1_CD");
	allBonds.add("NE2_CD");
	allBonds.add("1HE2_NE2");
	allBonds.add("2HE2_NE2");

	allBonds.add("1HE2_2HE2");
	allBonds.add("2HE2_1HE2");
	allBonds.add("1HG_2HG");
	allBonds.add("2HG_1HG");
	allBonds.add("CD_2HE2");
	allBonds.add("CD_1HE2");
	allBonds.add("2HE2_CD");
	allBonds.add("1HE2_CD");
	allBonds.add("OE1_NE2");
	allBonds.add("NE2_OE1");

	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector gluBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_CD");
	allBonds.add("CG_1HG");
	allBonds.add("CG_2HG");
	allBonds.add("CD_OE1");
	allBonds.add("CD_OE2");
	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CD_CG");
	allBonds.add("1HG_CG");
	allBonds.add("2HG_CG");
	allBonds.add("OE1_CD");
	allBonds.add("OE2_CD");

	allBonds.add("1HG_2HG");
	allBonds.add("2HG_1HG");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector hisBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_ND1");
	allBonds.add("CG_CD2");
	allBonds.add("ND1_HD1");
	allBonds.add("ND1_CE1");
	allBonds.add("CE1_HE1");
	allBonds.add("CE1_NE2");
	allBonds.add("NE2_CD2");
	allBonds.add("CD2_HD2");

	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("ND1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HD1_ND1");
	allBonds.add("CE1_ND1");
	allBonds.add("HE1_CE1");
	allBonds.add("NE2_CE1");
	allBonds.add("CD2_NE2");
	allBonds.add("HD2_CD2");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector ileBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG1");
	allBonds.add("CB_HB");
	allBonds.add("CB_CG2");
	allBonds.add("CG1_CD1");
	allBonds.add("CG1_1HG1");
	allBonds.add("CG1_2HG1");
	allBonds.add("CG2_1HG2");
	allBonds.add("CG2_2HG2");
	allBonds.add("CG2_3HG2");
	allBonds.add("CD1_1HD1");
	allBonds.add("CD1_2HD1");
	allBonds.add("CD1_3HD1");

	allBonds.add("CG1_CB");
	allBonds.add("HB_CB");
	allBonds.add("CG2_CB");
	allBonds.add("CD1_CG1");
	allBonds.add("1HG1_CG1");
	allBonds.add("2HG1_CG1");
	allBonds.add("1HG2_CG2");
	allBonds.add("2HG2_CG2");
	allBonds.add("3HG2_CG2");
	allBonds.add("1HD1_CD1");
	allBonds.add("2HD1_CD1");
	allBonds.add("3HD1_CD1");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector leuBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_CD1");
	allBonds.add("CG_CD2");
	allBonds.add("CG_HG");
	allBonds.add("CD1_1HD1");
	allBonds.add("CD1_2HD1");
	allBonds.add("CD1_3HD1");
	allBonds.add("CD2_1HD2");
	allBonds.add("CD2_2HD2");
	allBonds.add("CD2_3HD2");
	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CD1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HG_CG");
	allBonds.add("1HD1_CD1");
	allBonds.add("2HD1_CD1");
	allBonds.add("3HD1_CD1");
	allBonds.add("1HD2_CD2");
	allBonds.add("2HD2_CD2");
	allBonds.add("3HD2_CD2");

	allBonds.add("1HB_2HB");
	allBonds.add("2HB_1HB");
	allBonds.add("1HD1_2HD1");
	allBonds.add("1HD1_3HD1");
	allBonds.add("2HD1_3HD1");
	allBonds.add("1HD2_2HD2");
	allBonds.add("1HD2_3HD2");
	allBonds.add("2HD2_3HD2");
	allBonds.add("2HD1_1HD1");
	allBonds.add("3HD1_1HD1");
	allBonds.add("3HD1_2HD1");
	allBonds.add("2HD2_1HD2");
	allBonds.add("3HD2_1HD2");
	allBonds.add("3HD2_2HD2");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector metBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_1HG");
	allBonds.add("CG_2HG");
	allBonds.add("CG_SD");
	allBonds.add("SD_CE");
	allBonds.add("CE_1HE");
	allBonds.add("CE_2HE");
	allBonds.add("CE_3HE");

	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("1HG_CG");
	allBonds.add("2HG_CG");
	allBonds.add("SD_CG");
	allBonds.add("CE_SD");
	allBonds.add("1HE_CE");
	allBonds.add("2HE_CE");
	allBonds.add("3HE_CE");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector pheBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_CD1");
	allBonds.add("CG_CD2");
	allBonds.add("CD1_HD1");
	allBonds.add("CD1_CE1");
	allBonds.add("CE1_HE1");
	allBonds.add("CE1_CZ");
	allBonds.add("CZ_HZ");
	allBonds.add("CZ_CE2");
	allBonds.add("CE2_CD2");
	allBonds.add("CE2_HE2");
	allBonds.add("CD2_HD2");

	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CD1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HD1_CD1");
	allBonds.add("CE1_CD1");
	allBonds.add("HE1_CE1");
	allBonds.add("CZ_CE1");
	allBonds.add("HZ_CZ");
	allBonds.add("CE2_CZ");
	allBonds.add("CD2_CE2");
	allBonds.add("HE2_CE2");
	allBonds.add("HD2_CD2");

	allBonds.add("1HB_2HB");
	allBonds.add("2HB_1HB");
	Collections.sort(allBonds);
	return allBonds;
    }
    static final Vector trpBonds(){
    	Vector allBonds = backboneBonds();
    	allBonds.add("CB_CG");
    	allBonds.add("CB_1HB");
    	allBonds.add("CB_2HB");    	
    	allBonds.add("CG_CD1");
    	allBonds.add("CG_CD2");    	
    	allBonds.add("CD2_CE2");
    	allBonds.add("CE2_NE1");
    	allBonds.add("NE1_CD1");    	
    	allBonds.add("CD1_1HD1");
    	allBonds.add("NE1_1HE1");    	
    	allBonds.add("CD2_CE3");
    	allBonds.add("CE3_CZ3");
    	allBonds.add("CZ3_CH2");
    	allBonds.add("CH2_CZ2");
    	allBonds.add("CZ2_CE2");
    	allBonds.add("CE2_CD2");
    	allBonds.add("CE3_3HE");
    	allBonds.add("CZ3_3HZ");
    	allBonds.add("CH2_2HH");
    	allBonds.add("CZ2_2HZ");
    	
    	allBonds.add("CG_CB");
    	allBonds.add("1HB_CB");
    	allBonds.add("2HB_CB");
    	allBonds.add("CD1_CG");
    	allBonds.add("CD2_CG");    	
    	allBonds.add("CE2_CD2");
    	allBonds.add("NE1_CE2");
    	allBonds.add("CD1_NE1");    	
    	allBonds.add("1HD1_CD1");    	
    	allBonds.add("1HE1_NE1");    	
    	allBonds.add("CE3_CD2");
    	allBonds.add("CZ3_CE3");
    	allBonds.add("CH2_CZ3");
    	allBonds.add("CZ2_CH2");
    	allBonds.add("CE2_CZ2");
    	allBonds.add("CD2_CE2");
    	allBonds.add("3HE_CE3");
    	allBonds.add("3HZ_CZ3");
    	allBonds.add("2HH_CH2");
    	allBonds.add("2HZ_CZ2");    	

    	allBonds.add("1HB_2HB");
    	allBonds.add("2HB_1HB");
    	Collections.sort(allBonds);
    	return allBonds;
        }


    static final Vector tyrBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_CD1");
	allBonds.add("CG_CD2");
	allBonds.add("CD1_HD1");
	allBonds.add("CD1_CE1");
	allBonds.add("CE1_HE1");
	allBonds.add("CE1_CZ");
	allBonds.add("CZ_OH");
	allBonds.add("CZ_CE2");
	allBonds.add("CE2_CD2");
	allBonds.add("CE2_HE2");
	allBonds.add("CD2_HD2");
	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CD1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HD1_CD1");
	allBonds.add("CE1_CD1");
	allBonds.add("HE1_CE1");
	allBonds.add("CZ_CE1");
	allBonds.add("OH_CZ");
	allBonds.add("CE2_CZ");
	allBonds.add("CD2_CE2");
	allBonds.add("HE2_CE2");
	allBonds.add("HD2_CD2");

	allBonds.add("1HB_2HB");
	allBonds.add("2HB_1HB");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector lysBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("CG_CD");
	allBonds.add("CG_1HG");
	allBonds.add("CG_2HG");
	allBonds.add("CD_CE");
	allBonds.add("CD_1HD");
	allBonds.add("CD_2HD");
	allBonds.add("CE_NZ");
	allBonds.add("CE_1HE");
	allBonds.add("CE_2HE");
	allBonds.add("NZ_1HZ");
	allBonds.add("NZ_2HZ");
	allBonds.add("NZ_3HZ");

	allBonds.add("CG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("CD_CG");
	allBonds.add("1HG_CG");
	allBonds.add("2HG_CG");
	allBonds.add("CE_CD");
	allBonds.add("1HD_CD");
	allBonds.add("2HD_CD");
	allBonds.add("NZ_CE");
	allBonds.add("1HE_CE");
	allBonds.add("2HE_CE");
	allBonds.add("1HZ_NZ");
	allBonds.add("2HZ_NZ");
	allBonds.add("3HZ_NZ");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector thrBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_OG1");
	allBonds.add("CB_HB");
	allBonds.add("CB_CG2");
	allBonds.add("CG2_1HG2");
	allBonds.add("CG2_2HG2");
	allBonds.add("CG2_3HG2");
	allBonds.add("OG1_CB");
	allBonds.add("HB_CB");
	allBonds.add("CG2_CB");
	allBonds.add("1HG2_CG2");
	allBonds.add("2HG2_CG2");
	allBonds.add("3HG2_CG2");

	allBonds.add("1HG2_2HG2");
	allBonds.add("1HG2_3HG2");
	allBonds.add("2HG2_3HG2");
	allBonds.add("2HG2_1HG2");
	allBonds.add("3HG2_1HG2");
	allBonds.add("3HG2_2HG2");
	Collections.sort(allBonds);
	return allBonds;
    }

    static final Vector valBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG1");
	allBonds.add("CB_HB");
	allBonds.add("CB_CG2");
	allBonds.add("CG1_1HG1");
	allBonds.add("CG1_2HG1");
	allBonds.add("CG1_3HG1");
	allBonds.add("CG2_1HG2");
	allBonds.add("CG2_2HG2");
	allBonds.add("CG2_3HG2");

	allBonds.add("CG1_CB");
	allBonds.add("HB_CB");
	allBonds.add("CG2_CB");
	allBonds.add("1HG1_CG1");
	allBonds.add("2HG1_CG1");
	allBonds.add("3HG1_CG1");
	allBonds.add("1HG2_CG2");
	allBonds.add("2HG2_CG2");
	allBonds.add("3HG2_CG2");
	Collections.sort(allBonds);
	return allBonds;
    }

   static final Vector serBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_OG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("OG_HG");

	allBonds.add("OG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("HG_OG");
	Collections.sort(allBonds);
	return allBonds;
    }

   static final Vector cysBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_SG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("SG_HG");

	allBonds.add("SG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("HG_SG");
	Collections.sort(allBonds);
	return allBonds;
    }*/
//////////////////////////////////////////////////////////
   /**
    * The configuration of each amino acid: bond information
    * BMRB new naming scheme
    */    
  static final Vector backboneBonds(){
	Vector allBonds = new Vector();
	allBonds.add("H_N");
	allBonds.add("N_H");
	allBonds.add("N_CA");
	allBonds.add("CA_N");
	allBonds.add("CA_HA");
	allBonds.add("HA_CA");
	allBonds.add("CA_C");
	allBonds.add("C_CA");
	allBonds.add("CA_CB");
	allBonds.add("CB_CA");
	allBonds.add("C_O");
	allBonds.add("O_C");
	return allBonds;
   }
  static final Vector glyBonds(){
      Vector allBonds = new Vector();
	allBonds.add("H_N");
	allBonds.add("N_H");
	allBonds.add("N_CA");
	allBonds.add("CA_N");
	allBonds.add("CA_C");
	allBonds.add("C_CA");

	allBonds.add("CA_HA3");
	allBonds.add("HA3_CA");
	allBonds.add("CA_HA2");
	allBonds.add("HA2_CA");
	allBonds.add("C_O");
	allBonds.add("O_C");
	allBonds.add("HA2_HA3");
	allBonds.add("HA3_HA2");
	Collections.sort(allBonds);
	return allBonds;
   }
	
   static final Vector alaBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_HB1");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("HB1_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");

	allBonds.add("HB1_HB2");  //We do not compute the internal vdw violation for methyl group
	allBonds.add("HB1_HB3");
	allBonds.add("HB2_HB3");
	allBonds.add("HB2_HB1");
	allBonds.add("HB3_HB1");
	allBonds.add("HB3_HB2");
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector argBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CG_CB");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");	
	allBonds.add("CG_HG2");
	allBonds.add("CG_HG3");
	allBonds.add("HG2_CG");
	allBonds.add("HG3_CG");
	allBonds.add("CD_CG");
	allBonds.add("CG_CD");
	allBonds.add("CD_HD2");
	allBonds.add("CD_HD3");
	allBonds.add("HD2_CD");
	allBonds.add("HD3_CD");
	allBonds.add("CD_NE");
	allBonds.add("NE_CD");
	allBonds.add("HE_NE");
	allBonds.add("NE_HE");
	allBonds.add("CZ_NE");
	allBonds.add("NE_CZ");	
	allBonds.add("CZ_NH1");
	allBonds.add("NH1_CZ");
	allBonds.add("CZ_NH2");
	allBonds.add("NH2_CZ");	
	allBonds.add("NH1_HH11");
	allBonds.add("NH1_HH12");
	allBonds.add("HH11_NH1");
	allBonds.add("HH12_NH1");
	allBonds.add("NH2_HH21");
	allBonds.add("NH2_HH22");
	allBonds.add("HH21_NH2");
	allBonds.add("HH22_NH2");
	//These are inherent vdw violations
	allBonds.add("HB2_HB3");  //do not compute the  vdw violation for methylene group
	allBonds.add("HB3_HB2");  //do not compute the  vdw violation for methylene group
	allBonds.add("HG2_HG3");  //do not compute the  vdw violation for methylene group
	allBonds.add("HG3_HG2");  //do not compute the  vdw violation for methylene group
	allBonds.add("HD2_HD3");  //do not compute the  vdw violation for methylene group
	allBonds.add("HD3_HD2");  //do not compute the  vdw violation for methylene group
	allBonds.add("HH12_HH11");
	allBonds.add("HH11_HH12");
	allBonds.add("HH22_HH21");
	allBonds.add("HH21_HH22");
	allBonds.add("CZ_HE");
	allBonds.add("CZ_HH11");
	allBonds.add("CZ_HH12");
	allBonds.add("CZ_HH21");
	allBonds.add("CZ_HH22");
	allBonds.add("HE_CZ");
	allBonds.add("HH11_CZ");
	allBonds.add("HH12_CZ");
	allBonds.add("HH21_CZ");
	allBonds.add("HH22_CZ");
	
	allBonds.add("NH1_NE");
	allBonds.add("NE_NH1");
	allBonds.add("NH2_NE");
	allBonds.add("NE_NH2");
	
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector proBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CG_CB");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("CG_HG2");
	allBonds.add("CG_HG3");
	allBonds.add("HG2_CG");
	allBonds.add("HG3_CG");
	allBonds.add("CD_CG");
	allBonds.add("CG_CD");
	allBonds.add("CD_HD2");
	allBonds.add("CD_HD3");
	allBonds.add("HD2_CD");
	allBonds.add("HD3_CD");
	allBonds.add("CD_N");
	allBonds.add("N_CD");
	//These are inherent vdw violations
	allBonds.add("HB2_HB3");  
	allBonds.add("HB3_HB2");  
	allBonds.add("HG2_HG3");  
	allBonds.add("HG3_HG2");  
	allBonds.add("HD2_HD3");  
	allBonds.add("HD3_HD2");  
	allBonds.add("CA_CG");  
	allBonds.add("CG_CA");  
	allBonds.add("CB_CD");  
	allBonds.add("CD_CB");  
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector asnBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_OD1");
	allBonds.add("CG_ND2");
	allBonds.add("ND2_HD21");
	allBonds.add("ND2_HD22");

	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("OD1_CG");
	allBonds.add("ND2_CG");
	allBonds.add("HD21_ND2");
	allBonds.add("HD22_ND2");
	
	allBonds.add("HB2_HB3");  //added by zeng
	allBonds.add("HB3_HB2");  
	allBonds.add("HD21_HD22");
	allBonds.add("HD22_HD21");
	
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector aspBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_OD1");
	allBonds.add("CG_OD2");

	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("OD1_CG");
	allBonds.add("OD2_CG");
	
	allBonds.add("OD2_HD2");//added by zeng
	allBonds.add("HD2_OD2");//added by zeng
	allBonds.add("HB2_HB3");  //added by zeng
	allBonds.add("HB3_HB2");  
	
	allBonds.add("OD2_OD1");//added by zeng
	allBonds.add("OD1_OD2");//added by zeng
	
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector glnBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_CD");
	allBonds.add("CG_HG2");
	allBonds.add("CG_HG3");
	allBonds.add("CD_OE1");
	allBonds.add("CD_NE2");
	allBonds.add("NE2_HE21");
	allBonds.add("NE2_HE22");
	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("CD_CG");
	allBonds.add("HG2_CG");
	allBonds.add("HG3_CG");
	allBonds.add("OE1_CD");
	allBonds.add("NE2_CD");
	allBonds.add("HE21_NE2");
	allBonds.add("HE22_NE2");

	allBonds.add("HE21_HE22");
	allBonds.add("HE22_HE21");
	allBonds.add("HG2_HG3");
	allBonds.add("HG3_HG2");
	allBonds.add("CD_2HE2");
	allBonds.add("CD_1HE2");
	allBonds.add("HE21_CD");
	allBonds.add("HE22_CD");
	allBonds.add("OE1_NE2");
	allBonds.add("NE2_OE1");

	allBonds.add("HB2_HB3");  //added by zeng
	allBonds.add("HB3_HB2");  
	allBonds.add("HG2_HG3");
	allBonds.add("HG3_HG2");
	allBonds.add("HE21_HE22");
	allBonds.add("HE22_HE21");
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector gluBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_CD");
	allBonds.add("CG_HG2");
	allBonds.add("CG_HG3");
	allBonds.add("CD_OE1");
	allBonds.add("CD_OE2");
	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("CD_CG");
	allBonds.add("HG2_CG");
	allBonds.add("HG3_CG");
	allBonds.add("OE1_CD");
	allBonds.add("OE2_CD");

	allBonds.add("HB2_HB3");  //added by zeng
	allBonds.add("HB3_HB2");  
	allBonds.add("HG2_HG3");
	allBonds.add("HG3_HG2");
	
	allBonds.add("OE1_OE2");
	allBonds.add("OE2_OE1");
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector hisBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_ND1");
	allBonds.add("CG_CD2");
	allBonds.add("ND1_HD1");
	allBonds.add("ND1_CE1");
	allBonds.add("CE1_HE1");
	allBonds.add("CE1_NE2");
	allBonds.add("NE2_CD2");
	allBonds.add("CD2_HD2");

	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("ND1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HD1_ND1");
	allBonds.add("CE1_ND1");
	allBonds.add("HE1_CE1");
	allBonds.add("NE2_CE1");
	allBonds.add("CD2_NE2");
	allBonds.add("HD2_CD2");
	
	allBonds.add("HB2_HB3");  //added by zeng
	allBonds.add("HB3_HB2");  
	allBonds.add("HE2_NE2");
	allBonds.add("NE2_HE2");
	
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector ileBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG1");
	allBonds.add("CB_HB");
	allBonds.add("CB_CG2");
	allBonds.add("CG1_CD1");
	allBonds.add("CG1_HG12");
	allBonds.add("CG1_HG13");
	allBonds.add("CG2_HG21");
	allBonds.add("CG2_HG22");
	allBonds.add("CG2_HG23");
	allBonds.add("CD1_HD11");
	allBonds.add("CD1_HD12");
	allBonds.add("CD1_HD13");

	allBonds.add("CG1_CB");
	allBonds.add("HB_CB");
	allBonds.add("CG2_CB");
	allBonds.add("CD1_CG1");
	allBonds.add("HG12_CG1");
	allBonds.add("HG13_CG1");
	allBonds.add("HG21_CG2");
	allBonds.add("HG22_CG2");
	allBonds.add("HG23_CG2");
	allBonds.add("HD11_CD1");
	allBonds.add("HD12_CD1");
	allBonds.add("HD13_CD1");
	
	allBonds.add("HG21_HG22");//added by zeng
	allBonds.add("HG21_HG23");
	allBonds.add("HG22_HG21");
	allBonds.add("HG22_HG23");
	allBonds.add("HG23_HG21");
	allBonds.add("HG23_HG22");	
	allBonds.add("HD11_HD12");
	allBonds.add("HD11_HD13");
	allBonds.add("HD12_HD11");
	allBonds.add("HD12_HD13");
	allBonds.add("HD13_HD11");
	allBonds.add("HD13_HD12");
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector leuBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_CD1");
	allBonds.add("CG_CD2");
	allBonds.add("CG_HG");
	allBonds.add("CD1_HD11");
	allBonds.add("CD1_HD12");
	allBonds.add("CD1_HD13");
	allBonds.add("CD2_HD21");
	allBonds.add("CD2_HD22");
	allBonds.add("CD2_HD23");
	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("CD1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HG_CG");
	allBonds.add("HD11_CD1");
	allBonds.add("HD12_CD1");
	allBonds.add("HD13_CD1");
	allBonds.add("HD21_CD2");
	allBonds.add("HD22_CD2");
	allBonds.add("HD23_CD2");

	allBonds.add("HB2_HB3");
	allBonds.add("HB3_HB2");
	allBonds.add("HD11_HD12");
	allBonds.add("HD11_HD13");
	allBonds.add("HD12_HD13");
	allBonds.add("HD21_HD22");
	allBonds.add("HD21_HD23");
	allBonds.add("HD22_HD23");
	allBonds.add("HD12_HD11");
	allBonds.add("HD13_HD11");
	allBonds.add("HD13_HD12");
	allBonds.add("HD22_HD21");
	allBonds.add("HD23_HD21");
	allBonds.add("HD23_HD22");
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector metBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_HG2");
	allBonds.add("CG_HG3");
	allBonds.add("CG_SD");
	allBonds.add("SD_CE");
	allBonds.add("CE_HE1");
	allBonds.add("CE_HE2");
	allBonds.add("CE_HE3");

	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("HG2_CG");
	allBonds.add("HG3_CG");
	allBonds.add("SD_CG");
	allBonds.add("CE_SD");
	allBonds.add("HE1_CE");
	allBonds.add("HE2_CE");
	allBonds.add("HE3_CE");
	
	allBonds.add("HB2_HB3");  //added by zeng
	allBonds.add("HB3_HB2"); 	
	allBonds.add("HG2_HG3");  //added by zeng
	allBonds.add("HG3_HG2");  
	allBonds.add("HE1_HE2");
	allBonds.add("HE1_HE3");
	allBonds.add("HE2_HE1");
	allBonds.add("HE2_HE3");  //added by zeng
	allBonds.add("HE3_HE2");  
	allBonds.add("HE3_HE1"); 
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector pheBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_CD1");
	allBonds.add("CG_CD2");
	allBonds.add("CD1_HD1");
	allBonds.add("CD1_CE1");
	allBonds.add("CE1_HE1");
	allBonds.add("CE1_CZ");
	allBonds.add("CZ_HZ");
	allBonds.add("CZ_CE2");
	allBonds.add("CE2_CD2");
	allBonds.add("CE2_HE2");
	allBonds.add("CD2_HD2");

	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("CD1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HD1_CD1");
	allBonds.add("CE1_CD1");
	allBonds.add("HE1_CE1");
	allBonds.add("CZ_CE1");
	allBonds.add("HZ_CZ");
	allBonds.add("CE2_CZ");
	allBonds.add("CD2_CE2");
	allBonds.add("HE2_CE2");
	allBonds.add("HD2_CD2");

	allBonds.add("HB2_HB3");
	allBonds.add("HB3_HB2");
	Collections.sort(allBonds);
	return allBonds;
   }
   static final Vector trpBonds(){
   	Vector allBonds = backboneBonds();
   	allBonds.add("CB_CG");
   	allBonds.add("CB_HB2");
   	allBonds.add("CB_HB3");    	
   	allBonds.add("CG_CD1");
   	allBonds.add("CG_CD2");    	
   	allBonds.add("CD2_CE2");
   	allBonds.add("CE2_NE1");
   	allBonds.add("NE1_CD1");    	
   	allBonds.add("CD1_HD1");
   	allBonds.add("NE1_HE1");    	
   	allBonds.add("CD2_CE3");
   	allBonds.add("CE3_CZ3");
   	allBonds.add("CZ3_CH2");
   	allBonds.add("CH2_CZ2");
   	allBonds.add("CZ2_CE2");
   	allBonds.add("CE2_CD2");
   	allBonds.add("CE3_HE3");
   	allBonds.add("CZ3_HZ3");
   	allBonds.add("CH2_HH2");
   	allBonds.add("CZ2_HZ2");
   	
   	allBonds.add("CG_CB");
   	allBonds.add("HB2_CB");
   	allBonds.add("HB3_CB");
   	allBonds.add("CD1_CG");
   	allBonds.add("CD2_CG");    	
   	allBonds.add("CE2_CD2");
   	allBonds.add("NE1_CE2");
   	allBonds.add("CD1_NE1");    	
   	allBonds.add("HD1_CD1");    	
   	allBonds.add("HE1_NE1");    	
   	allBonds.add("CE3_CD2");
   	allBonds.add("CZ3_CE3");
   	allBonds.add("CH2_CZ3");
   	allBonds.add("CZ2_CH2");
   	allBonds.add("CE2_CZ2");
   	allBonds.add("CD2_CE2");
   	allBonds.add("HE3_CE3");
   	allBonds.add("HZ3_CZ3");
   	allBonds.add("HH2_CH2");
   	allBonds.add("HZ2_CZ2");    	

   	allBonds.add("HB2_HB3");
   	allBonds.add("HB3_HB2");
   	Collections.sort(allBonds);
   	return allBonds;
       }


   static final Vector tyrBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_CD1");
	allBonds.add("CG_CD2");
	allBonds.add("CD1_HD1");
	allBonds.add("CD1_CE1");
	allBonds.add("CE1_HE1");
	allBonds.add("CE1_CZ");
	allBonds.add("CZ_OH");
	allBonds.add("CZ_CE2");
	allBonds.add("CE2_CD2");
	allBonds.add("CE2_HE2");
	allBonds.add("CD2_HD2");
	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("CD1_CG");
	allBonds.add("CD2_CG");
	allBonds.add("HD1_CD1");
	allBonds.add("CE1_CD1");
	allBonds.add("HE1_CE1");
	allBonds.add("CZ_CE1");
	allBonds.add("OH_CZ");
	allBonds.add("CE2_CZ");
	allBonds.add("CD2_CE2");
	allBonds.add("HE2_CE2");
	allBonds.add("HD2_CD2");
	allBonds.add("OH_HH");
	allBonds.add("HH_OH");
	
	allBonds.add("HB2_HB3");
	allBonds.add("HB3_HB2");
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector lysBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("CG_CD");
	allBonds.add("CG_HG2");
	allBonds.add("CG_HG3");
	allBonds.add("CD_CE");
	allBonds.add("CD_HD2");
	allBonds.add("CD_HD3");
	allBonds.add("CE_NZ");
	allBonds.add("CE_HE2");
	allBonds.add("CE_HE3");
	allBonds.add("NZ_HZ1");
	allBonds.add("NZ_HZ2");
	allBonds.add("NZ_HZ3");

	allBonds.add("CG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("CD_CG");
	allBonds.add("HG2_CG");
	allBonds.add("HG3_CG");
	allBonds.add("CE_CD");
	allBonds.add("HD2_CD");
	allBonds.add("HD3_CD");
	allBonds.add("NZ_CE");
	allBonds.add("HE2_CE");
	allBonds.add("HE3_CE");
	allBonds.add("HZ1_NZ");
	allBonds.add("HZ2_NZ");
	allBonds.add("HZ3_NZ");
	
	allBonds.add("HB2_HB3");//added by zeng
	allBonds.add("HB3_HB2");
	allBonds.add("HG2_HG3");
	allBonds.add("HG3_HG2");
	allBonds.add("HD2_HD3");
	allBonds.add("HD3_HD2");
	allBonds.add("HE2_HE3");
	allBonds.add("HE3_HE2");
	allBonds.add("HZ1_HZ2");
	allBonds.add("HZ1_HZ3");
	allBonds.add("HZ2_HZ1");
	allBonds.add("HZ2_HZ3");
	allBonds.add("HZ3_HZ1");
	allBonds.add("HZ3_HZ2");
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector thrBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_OG1");
	allBonds.add("CB_HB");
	allBonds.add("CB_CG2");
	allBonds.add("CG2_HG21");
	allBonds.add("CG2_HG22");
	allBonds.add("CG2_HG23");
	allBonds.add("OG1_CB");
	allBonds.add("HB_CB");
	allBonds.add("CG2_CB");
	allBonds.add("HG21_CG2");
	allBonds.add("HG22_CG2");
	allBonds.add("HG23_CG2");

	allBonds.add("HG21_HG22");
	allBonds.add("HG21_HG23");
	allBonds.add("HG22_HG23");
	allBonds.add("HG22_HG21");
	allBonds.add("HG23_HG21");
	allBonds.add("HG23_HG22");
	allBonds.add("HG1_OG1");
	allBonds.add("OG1_HG1");
	
	Collections.sort(allBonds);
	return allBonds;
   }

   static final Vector valBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_CG1");
	allBonds.add("CB_HB");
	allBonds.add("CB_CG2");
	allBonds.add("CG1_HG11");
	allBonds.add("CG1_HG12");
	allBonds.add("CG1_HG13");
	allBonds.add("CG2_HG21");
	allBonds.add("CG2_HG22");
	allBonds.add("CG2_HG23");

	allBonds.add("CG1_CB");
	allBonds.add("HB_CB");
	allBonds.add("CG2_CB");
	allBonds.add("HG11_CG1");
	allBonds.add("HG12_CG1");
	allBonds.add("HG13_CG1");
	allBonds.add("HG21_CG2");
	allBonds.add("HG22_CG2");
	allBonds.add("HG23_CG2");
	
	allBonds.add("HG21_HG22");//added by zeng
	allBonds.add("HG21_HG23");
	allBonds.add("HG22_HG23");
	allBonds.add("HG22_HG21");
	allBonds.add("HG23_HG21");
	allBonds.add("HG23_HG22");
	
	allBonds.add("HG11_HG12");//added by zeng
	allBonds.add("HG11_HG13");
	allBonds.add("HG12_HG13");
	allBonds.add("HG12_HG11");
	allBonds.add("HG13_HG11");
	allBonds.add("HG13_HG12");
	Collections.sort(allBonds);
	return allBonds;
   }

  static final Vector serBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_OG");
	allBonds.add("CB_1HB");
	allBonds.add("CB_2HB");
	allBonds.add("OG_HG");

	allBonds.add("OG_CB");
	allBonds.add("1HB_CB");
	allBonds.add("2HB_CB");
	allBonds.add("HG_OG");
	allBonds.add("HB2_HB3");//added by zeng
	allBonds.add("HB3_HB2");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	Collections.sort(allBonds);
	return allBonds;
   }

  static final Vector cysBonds(){
	Vector allBonds = backboneBonds();
	allBonds.add("CB_SG");
	allBonds.add("CB_HB2");
	allBonds.add("CB_HB3");
	allBonds.add("SG_HG");

	allBonds.add("SG_CB");
	allBonds.add("HB2_CB");
	allBonds.add("HB3_CB");
	allBonds.add("HG_SG");
	
	allBonds.add("HB2_HB3");//added by zeng
	allBonds.add("HB3_HB2");
	Collections.sort(allBonds);
	return allBonds;
   }
  
   
   
   
//////////////////////////////////////////////////////////
    static final Map residueBonds(){
    	Map mapBonds = new TreeMap();
    	mapBonds.put("ALA", alaBonds());
    	mapBonds.put("ARG", argBonds());
    	mapBonds.put("ASN", asnBonds());
    	mapBonds.put("ASP", aspBonds());
    	mapBonds.put("GLY", glyBonds());
    	mapBonds.put("CYS", cysBonds());
    	mapBonds.put("GLN", glnBonds());
    	mapBonds.put("GLU", gluBonds());
    	mapBonds.put("HIS", hisBonds());
    	mapBonds.put("ILE", ileBonds());
    	mapBonds.put("LEU", leuBonds());
    	mapBonds.put("LYS", lysBonds());
    	mapBonds.put("MET", metBonds());
    	mapBonds.put("PHE", pheBonds());
    	mapBonds.put("SER", serBonds());
    	mapBonds.put("THR", thrBonds());
    	mapBonds.put("TYR", tyrBonds());
    	mapBonds.put("VAL", valBonds());
    	mapBonds.put("PRO", proBonds());
    	mapBonds.put("TRP", trpBonds());
    	return mapBonds;
        }
    
    static final double vdwH2H = 1.80;  //We need to find the authetic values for these vdw numbers
    static final double vdwC2C = 2.45;  //Some parameters were calibrated with ideal helix
    static final double vdwO2O = 2.50; 
    static final double vdwN2N = 2.50; 
    static final double vdwS2S = 2.50; 
    static final double vdwH2O = 2.00;//1.90;//2.00;
    static final double vdwO2H = 2.00;
    static final double vdwH2C = 2.03;
    static final double vdwC2H = 2.03;
    static final double vdwH2N = 2.00; 
    static final double vdwN2H = 2.00; 
    static final double vdwC2O = 2.40;//2.25;//2.40;
    static final double vdwO2C = 2.40;
    static final double vdwC2N = 2.45;
    static final double vdwN2C = 2.45;
    static final double vdwO2N = 2.25; 
    static final double vdwN2O = 2.25; 

//     static final double vdwH2H = 1.90;  //The values was taken from the Flory book, p254.
//     static final double vdwC2C = 3.00;  //
//     static final double vdwO2O = 2.70; 
//     static final double vdwN2N = 2.60; 
//     static final double vdwS2S = 2.90; //not used and not from Flory
//     static final double vdwH2O = 2.20;
//     static final double vdwO2H = 2.20;
//     static final double vdwH2C = 2.20;
//     static final double vdwC2H = 2.20;
//     static final double vdwH2N = 2.20; 
//     static final double vdwN2H = 2.20; 
//     static final double vdwC2O = 2.70;
//     static final double vdwO2C = 2.70;
//     static final double vdwC2N = 2.80;
//     static final double vdwN2C = 2.80;
//     static final double vdwO2N = 2.60; 
//     static final double vdwN2O = 2.60; 

    static final double vdwS2O = 2.20;
    static final double vdwO2S = 2.20;
    static final double vdwS2H = 2.10;
    static final double vdwH2S = 2.20;
    static final double vdwC2S = 2.30;
    static final double vdwS2C = 2.30; 
    static final double vdwN2S = 2.10; 
    static final double vdwS2N = 2.10; 
    
    static final double hBondH2O = 1.569;
    static final double hBondO2H = 1.569;
    static final Map bondsMap = residueBonds();
}
