This file explains how to prepare the required data for input into NVR
and the associated programs. Then, it explains how to run the various programs. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP # 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1st, prepare a file containing the data one might get from an HSQC, and its variants. 
In particular, prepare a file where each row corresponds to a peak in the HSQC. Each 
row must have the following format

PEAK_ID H_shift	N_SHIFT	RDC_MEDIUM_1  RDC_MEDIUM_2 SLOW_EXCHANGE

The first column is the ID of the peak (using some arbitrary numbering system). The
second column is the amide proton shift, the third is the 
amide nitrogen shift, the fourth column is the rdc recorded for that peak in the first 
aligning medium, and the fifth column is the rdc recorded for that peak in the
second aligning medium. The final column corresponds to the results of an amide
exchange experiment. Put a 1 in the final column if that peak is a slow-exchanger (note: this is correct 28.6.09. slow exchanges-->hydrogen bond in .mr file), 
otherwise, put a 0 in the final column. (IT LOOKS LIKE THIS IS CORRECT, 28/6/09)

Here is an example, taken from an HSQC for human ubiquitin

1 8.942	123.077	-8.7	13.2 0

If you cannot obtain an rdc for a given peak, use the 'dummy' value of 
-999 instead. 

For example, this peak had an rdc only in 1 medium

66 8.593	123.967	3.6	-999  0

and this one had no recordable RDCs

68 8.438	122.163	-999	-999 0 


An example of the HSQC data file is provided in the file hsqcdata.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP # 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Next, prepare an input file corresponding to the model. 

Each row of the corresponds to the backbone amide bond from a PDB file

The format of each row is as follows

RES_NUMBER RES_TYPE X Y Z  SS_TYPE LABILE X_pos Y_pos Z_pos

RES_NUMBER is just the residue number
	 
RES_TYPE is the three letter amino acid code

X, Y and Z are the x,y, and z components of a normalized vector describing the
orientation of the amide bond vector in the PDB frame. You can construct
this vector by subtracting the x,y and z coordinates of the amide nitrogen
from the x,y and z coordinates of the amide proton. Then, you must
normalize the vector. 

SS_TYPE is the secondary structure type of the residue. Put a B if the residue
is in a beta strand, H if it is in an alpha helix and C if it is in random coil. 

LABILE  This is a flag which states whether the amide proton can exchange with 
the water. This can be determined using a program such as MOLMOL.  Put a Y if the 
amide proton is solvent accessible and not involved in a hydrogen bond. Otherwise, 
put a N there. (THIS SHOULD BE REVERSED. 'Y'-->IN AN H-BOND ACCORDING TO MOLMOL. 'N'-->NOT IN AN H-BOND ACCORDING TO MOLMOL.

X_pos, Y_pos and Z_pos are the coordinates of the amide proton taken from the 
PDB file. These will be used to compute distances between protons which are, in turn, 
correlated to Dnns.

Here is an example of an entry in the file taken from the PDB file 1UBI

6 LYS -0.701261	-0.661642 0.265449 B Y -2.82 6.056 -3.992


An example of the model data file is provided in the file modeldata.m.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP # 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Next, prepare a file called NOES.m

Each row of this file should contain two values, 

PEAK_ID_1 PEAK_ID_2

where PEAK_ID_1 and PEAK_ID_2 are the peak id's involved in a Dnn. note, that the ids should 
correspond to the peak ids in the file hsqcdata.m. That is, you need to find a set of dnn's for which 
you can unambiguously determine which HSQC peaks they correspond to. 

So for example if there is a dnn between peaks 1 and 16 your file would have the following entry

1 16

An example of an noe data file is presented in the file NOES.m




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP # 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Next, prepare a file called SHIFTS.m

Each row should have the following format

RES_NUM h_shift n_shift

Where RES_NUM is the residue number, h_shift is the predicted proton shift, n_shift is
the predicted nitrogen shift. 

Here is an example row from the SHIFTS predictions for 1UBI

2 8.02  127.42


Note that the rows in this file should have the same ordering in terms of resiudes as the model data file. 

An example of a SHIFTS data file is presented in the file SHIFTS.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP # 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Next, prepare a file called SHIFTX.m

Each row should have the following format

RES_NUM AA_TYPE SS_TYPE  d1  d2 d3 d4 d5 d6

Where RES_NUM is the residue number, AA_TYPE is the amino acid type (single letter code)
and SS_TYPE is the secondary structure type. The remaining values are chemical shift values. d2 and d3
will be used by NVR to compute assignment probabilities. 

Here is an example row from the SHIFTX predictions for 1UBI

2	Q	B	4.7603	8.6366	122.3549	53.3187	30.6096	173.0614


Note that the rows in this file should have the same ordering in terms of resiudes as the model data file. 

An example of a SHIFTX data file is presented in the file SHIFTS.m



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP # 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Finally, preapre a file called TOCSY.m.  

Each row should have the following format

HN_SHIFT  d2 N_SHIFT PEAK_ID  

Where HN_SHIFT is the amide proton shift, d2 is a side-chain shift, and
N_SHIFT is the amide nitrogen shift

Note that you must cross-reference the side chain peaks to a unique peak in the 
hsqc (from which you get the peak id). 

Here is an example of several side chain resonances from ubiquitin

8.949	1.602	123.129	1
8.96	1.855	123.132	1
8.947	2.222	123.135	1
8.96	2.223	123.132	1
8.961	5.262	123.133	1
8.947	8.95	123.132	1
8.96	8.952	123.133	1


An example of a TOCSY input file is presented in the file TOCSY.m