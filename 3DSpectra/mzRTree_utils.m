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
clear
close all
% clc

% settings for mzRTree creation, if needed increase your java heap space in
% matlab, instructions on the mathworks website.
number_of_strips=50;
mz_resolution=[];

%% importing the data
mzXMLFilesPath_tmp=[];
dataDir = uigetdir(pwd, 'Select the data files directory');
cd(dataDir)
fileMzXML=struct2cell(dir('*.mzXML'));
mzXMLFilesPath=fileMzXML(1,:);
for i=1:numel(mzXMLFilesPath)
    filePath=[pwd filesep mzXMLFilesPath{i}];
    varName_tmp=strtok(mzXMLFilesPath{i},'.')
    mzRTreePath=[pwd filesep varName_tmp filesep];
    tic
    eval([genvarname(varName_tmp) '=mzRTreeCreation(filePath,mzRTreePath,number_of_strips, mz_resolution);'])
    toc
    dataFile{i}=varName_tmp;
end
save([dataDir filesep 'dataRead.mat'])
	
