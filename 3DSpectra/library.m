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
function [pep_library]=library(filePath,varargin)
% it creates the pep library
% pep_library--> metadata
% filePath--> path;
% varargin --> mzXML file Path. To be used if MS1 scans associated to the identifications are not known and the ids file stores MSn scans only, the mzXML
% file path can be provided and the function will map MSn scans to MS1 ones

if isempty(filePath)
    [file path]=uigetfile('*.*', 'Select the file storing the metadata info');
    filePath=[path file];
end

metadata = importdata(filePath,'\t');
pep_library.peptidesSequence=[metadata.textdata(2:end,1)];
pep_library.charges=[metadata.data(:,1)];
pep_library.massH1=[metadata.data(:,2)];
pep_library.delta=[metadata.data(:,3)];
pep_library.massToSearch=[metadata.data(:,4)];

if isempty(varargin{1})
    pep_library.MS_scan=[metadata.data(:,5)];
else
    %% MS1 scan retrieval
    tic
    [MSall_to_MS1 tot_scan]=ms_ruler(varargin{1});
    pep_library.MS_scan=MSall_to_MS1(ScanNumberOriginal)';
    toc
end
	
