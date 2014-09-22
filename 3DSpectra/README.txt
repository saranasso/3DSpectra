% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This software is an implementation of the method described in
% "3DSpectra: A 3-Dimensional Quantification Algorithm For LC-MS Labeled Profile Data" 
% Copyright (C) 2012 Sara Nasso
%     
% 3DSpectra is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at
% your option) any later version.
% 
% 3DSpectra is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instructions for the 3DSpectra code released by Sara Nasso under GNU GPL license
% This instructions are meant to enable the user to reproduce submitted results

1.- preliminary step: download mzXML data files from:

https://maspectras.genome.tugraz.at/maspectras/

login as guest (pwd: guest). On the up left corner click on Menu and unfold the MASPECTRAS 1 Data folder and then the AutoQuantificationProfile sub-folder. Download all data in there and save them in a single folder on your file system.

2.- start a Matlab session and type in the command window: open classpath.txt

add at the end of the classpath.txt file your local file paths to the needed jars, as in the example here below:

/Users/saran/work/3DSpectra/jrap_StAX_v5.jar
/Users/saran/work/3DSpectra/maspectras.jar
/Users/saran/work/3DSpectra/massMatlab.jar
/Users/saran/work/3DSpectra/mzRTree.jar
/Users/saran/work/3DSpectra/igbbiojava.jar
/Users/saran/work/3DSpectra/outliers.jar

3.- mex update on new hardware (hist2data) as follows: mex hist2data.c

4.- check additional needed toolboxes (e.g., Bioinformatics, Statistics) are available typing in the command window: ver

5.- cd to 3DSpectra folder

6.- (not needed for evaluation data) create the needed peptide libraries following the instructions given in the readme.txt provided in the ID_files folder 

7.- cd to 3DSpectra folder once more if needed

8.- run the 3DSpectra algorithm by the main.m. 
In this main execution mzRTree for each data file will be created running the script mzRTree_utils.m and data analysis by means of 3DSpectra executed.
If you run into problems with 3DSpectra uncomment "rethrow(ME)" in the catch statement to check out the occurring error.
Finally, result evaluation is executed as well.

9.- analysis can take some hours, depending on the machine  (e.g., ~3 hours on a 2.2 GHz Intel Core i7 with 8 GB 1333 MHz DDR3)

10.- to evaluate results and compare them to MASPECTRAS ones run the main_results.m after cd-ing to the results folder (not needed if main.m is run)

Please notice that in this implementation the following data feature was exploited: samples resolution is constant both along retention times (~1sec) and m/z (0.04Da).
This allowed to exchange retention time (m/z) indexes and  their actual retention time (m/z) values. 
Other data, without constant sampling rates, would need slight modifications to the current implementation. Next release will add this additional feature.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LC/MS data analysis requires extended amounts of memory from the operating system.

If you receive errors related to memory, try the following:
Increase the virtual memory (swap space) for your operating system (with a recommended initial size of 3,069 and a maximum size of 16,368) as described in Memory Usage.
Set the 3 GB switch (32-bit Windows XP only) as described in Memory Usage.
If you receive errors related to Java heap space, increase your Java heap space:
If you have MATLAB version 7.10 (R2010a) or later, see
Java Heap Memory Preferences
If you have MATLAB version 7.9 (R2009b) or earlier, see
http://www.mathworks.com/support/solutions/data/1-18I2C.html

