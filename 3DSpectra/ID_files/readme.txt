Instructions on how to create the peptide libraries.

First of all store your identifications either in a txt or xls file strictly following given schemas. The xls file will only work under windows systems. 

If you stored in a txt file then run main_library.m to create your library (suggested option) otherwise if you stored in a xls file then run main_pepLibrary.m. Both main scripts have some parameters to be set accordingly to your local setup, so please change them properly before running them.

Once the libraries have been created and saved (done automatically by the main script) move the .mat files to the upper folder (i.e., 3DSpectra folder).

PLEASE NOTICE: if you choose to use the txt file and you put in there the MS1 scan number corresponding to the MSn scan where the identification occurred,mzXML data files don't need to be stored locally . Conversely, if you cannot retrieve the MS1 scan number, setting properly the parameters (e.g., idsFileIsStoringMS1scans=0) the code will retrieve it automatically, but all mzXML files need to be stored locally in a single data folder.