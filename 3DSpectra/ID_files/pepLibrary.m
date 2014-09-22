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
function pep_library=pepLibrary(path_metadata_xls,path_data_storage)
% path_xls--> metadata file containing all IDs structured according to a
% known pattern
% pep_library--> peptide's library
% path_metadata_xls= xls file storing ids
% path_data_storage= mzXML data file;

[num text metadata]=xlsread(path_metadata_xls,'IDs');
%% variable modifications for our data
% 1	ICPL (K)	105.021469
% 2	ICPL (N-term N)	105.021454
% 3	ICPL:13C(6) (K)	111.041595
% 4	ICPL:13C(6) (N-term)	111.041595
% 5	Oxidation (M)	15.994919
number_of_peptides=size(metadata,1);
modifications=metadata(2:number_of_peptides,22);
for i=1:size(modifications,1)
    mods=modifications{i};
    if ~isnan(mods)
        num_mod_light(i,1)=numel([findstr(mods,'1'), findstr(mods,'2')]);
        num_mod_light_K(i,1)=numel(findstr(mods,'1'));
        num_mod_light_Nterm(i,1)=numel(findstr(mods,'2'));
        num_mod_heavy(i,1)=numel([findstr(mods,'3'), findstr(mods,'4')]);
        num_mod_heavy_K(i,1)=numel(findstr(mods,'3'));
        num_mod_heavy_Nterm(i,1)=numel(findstr(mods,'4'));
        if num_mod_light(i,1)~=0 || num_mod_heavy(i,1)~=0
            modification_flag(i,1)=1;
            if num_mod_light(i,1)~=0 && num_mod_heavy(i,1)==0
                label(i,:)='light';
            elseif num_mod_heavy(i,1)~=0 && num_mod_light(i,1)==0
                label(i,:)='heavy';
            else
                label(i,:)='noMod';
            end
        else
            modification_flag(i,1)=0;
            label(i,:)='noMod';
        end        
    else
        modification_flag(i,1)=0;
        label(i,:)='noMod';
    end
end

%% keeping only modified peptides
idxs_to_be_kept=find(modification_flag);
pep_library.peptidesSequence=metadata(idxs_to_be_kept+1,19);
pep_library.charges=cell2mat(metadata(idxs_to_be_kept+1,12));
pep_library.massH1=cell2mat(metadata(idxs_to_be_kept+1,13)); % experimental
% pep_library.massH1=metadata(idxs_to_be_kept+1,11); % theoretical
pep_library.delta=cell2mat(metadata(idxs_to_be_kept+1,14));
pep_library.massToSearch=cell2mat(metadata(idxs_to_be_kept+1,10));
ScanNumberOriginal=cell2mat(metadata(idxs_to_be_kept+1,24));
pep_library.proteinID=cell2mat(metadata(idxs_to_be_kept+1,1));
pep_library.score=cell2mat(metadata(idxs_to_be_kept+1,16));
pep_library.expectation=cell2mat(metadata(idxs_to_be_kept+1,17));
% here below no +1 on indexes since they've been already read from the xls file
pep_library.modification_flag=modification_flag(idxs_to_be_kept);
pep_library.num_mod_light=num_mod_light(idxs_to_be_kept);
pep_library.num_mod_light_K=num_mod_light_K(idxs_to_be_kept);
pep_library.num_mod_light_Nterm=num_mod_light_Nterm(idxs_to_be_kept);
pep_library.num_mod_heavy=num_mod_heavy(idxs_to_be_kept);
pep_library.num_mod_heavy_K=num_mod_heavy_K(idxs_to_be_kept);
pep_library.num_mod_heavy_Nterm=num_mod_heavy_Nterm(idxs_to_be_kept);
pep_library.label=label(idxs_to_be_kept,:);

%% MS1 scans retrieval
file_mzXML=path_data_storage;
tic
[MSall_to_MS1 tot_scan]=ms_ruler(file_mzXML);
pep_library.MS_scan=MSall_to_MS1(ScanNumberOriginal)';
toc
	
