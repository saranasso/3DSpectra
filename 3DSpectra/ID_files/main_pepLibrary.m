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
% PLEASE NOTICE: 
% this script creates libraries starting from xls files formatted as the provided example
% this script calls functions using xlsread that works only under windows environment. You can use main_library instead and format ids information as in the example txt file 
% rename mzXML and xls files properly: for each data file name should be equal to the corresponding ids file and it should end by the replicate suffix.

%% parameters to be set by user
fileNames{1}='1li1he'; % file name without replicate suffix. It helps grouping by replicate afterwards
fileNames{2}='1li2he';
fileNames{3}='1li5he';
fileNames{4}='1li10he';
fileNames{5}='2li1he';
fileNames{6}='5li1he';
fileNames{7}='10li1he';
replicates{1}='A'; % replicate suffix, please list all of them as shown here
replicates{2}='B';
replicates{3}='C';
metadata_folder_path='C:\ICPL data\IDs\';
data_folder_path='C:\ICPL data\profile\mzXML\';

%% building and saving libraries
for iter_ratios=1:length(fileNames)
    ratio=fileNames{iter_ratios}
    for iter_replicates=1:length(replicates)
        replicate=replicates{iter_replicates};
        path_metadata_xls=[metadata_folder_path ratio replicate '.xls'];
        path_data_storage=[data_folder_path ratio replicate '.mzXML']; % PLEASE NOTICE: rename mzXML files properly: name should be equal to the xls IDs file and it should end by the replicate suffix.
        eval(['library' replicate '=pepLibrary(path_metadata_xls,path_data_storage);'])
%         save(['library' ratio replicate])
    end
    save(['library' ratio])
    clearvars -except fileNames iter_ratios replicates
end

	
