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
% clear
load library1li5he
clearvars -except library* x* mz_res* dataDir settings 

%% adjusting to local paths
fileHeaderA=x1zu5A;
fileHeaderA.localFullPath=[dataDir filesep x1zu5A.folderName filesep];
fileHeaderB=x1zu5B;
fileHeaderB.localFullPath=[dataDir filesep x1zu5B.folderName filesep];
fileHeaderC=x1zu5C;
fileHeaderC.localFullPath=[dataDir filesep x1zu5C.folderName filesep];

%% initializing variables & parameters
pep_idx=1;
j=1; 
mzMultiplicationFactorForInt=0;
results.quants=struct('pep_seq',[],'pep_charge',[],'pep_idx',[],'ratios',[],'spettro_pep',[],'spettro_partner',[],'replicate',[],'corr',[]);
results.IDs=0;
resultsA=results;
resultsB=results;
resultsC=results;
%% settings
% settings.npeaks=4;
% settings.diff_mass=6.02; % Da
% settings.tosum=0; % binning
settings.k=round(inv(fileHeaderA.mz_resolution)); % replicates must feature the same mz resolution 
% settings.number_of_gaussians=settings.npeaks+1;
% settings.cut=100;
% settings.ratio='1zu5';
% settings.visualize=0;
% settings.elution_width_left=25;
% settings.elution_width_right=35;
% settings.elution_identification_error_tolerance=10;


pause off
% pause on
%% 3DSpectra execution
settings.replicate='A';
[resultsA]=Spectra3D(resultsA,libraryA,fileHeaderA,settings,pep_idx,j);

settings.replicate='B';
[resultsB]=Spectra3D(resultsB,libraryB,fileHeaderB,settings,pep_idx,j);

settings.replicate='C';
[resultsC]=Spectra3D(resultsC,libraryC,fileHeaderC,settings,pep_idx,j);

results=[resultsA,resultsB,resultsC];
results_quants=[resultsA.quants,resultsB.quants,resultsC.quants];
quants=results.quants;

save(settings.ratio)

	
