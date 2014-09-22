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
function [pep_sequence,flag, index, label]=findPeptide(pep_sequence,pep_idx,pep_library,varargin)
% it look for the current peptide in the pep library
% flag--> true if pep is in library
% index--> library indexes related to the pep
% pep_sequence--> pep sequence searched
% pep_library-->pep library for the dataset
% varargin--> peptide charge

i=1;
index=[];
flag=false;
[pepseq pepzero]=strtok(pep_sequence,'0');
if ~isempty(pepzero)
    pep_heavy=find(pepzero=='2',1);
    pep_light=find(pepzero=='3',1);
    pep_err=find(pepzero=='1',1);
    pep_nz=find(pepzero~='0',1); 
    
    if (isempty(pep_heavy) && isempty(pep_light)) || ~isempty(pep_err) || (~isempty(pep_heavy) && ~isempty(pep_light)) || isempty(pep_nz)
        flag=false;
        index=[];
        disp('NOT LABELED')
        label='NOT LABELED';
        return
    elseif ~isempty(pep_light) && isempty(pep_heavy)
        label='light';
    elseif ~isempty(pep_heavy) && isempty(pep_light)
        label='heavy';
    else
        flag=false;
        index=[];
        disp('?????????????????????????')
        label='NOT LABELED';
        return
    end
    while (pep_library.peptidesSequence{i}(1)~=pep_sequence(1))
        
        i=i+1;
        
    end
    
    while i<=length(pep_library.peptidesSequence) && (pep_library.peptidesSequence{i}(1)==pep_sequence(1))
        
        if isequal(pep_library.peptidesSequence{i},pep_sequence) 
            if ~isempty(varargin)
                pep_charge=pep_library.charges(i);
                if varargin{1}==pep_charge
                    if numel(varargin)==2
                        if isequal(varargin{2},round(pep_library.massH1(i)))
                            flag=true;
                            index=[index i]; 
                        else
                            flag=false;
                        end
                    else
                        flag=true;
                        index=[index i];
                    end
                end
            else
                flag=true;
                index=[index i];
            end
        end
        i=i+1;
    end
    
    
elseif isempty(pepzero)
    
    if ~pep_library.modification_flag(pep_idx)
        flag=false;
        index=[];
        disp('NOT LABELED')
        label='NOT LABELED';
        return
    elseif isequal(pep_library.label(pep_idx,:),'light')
        label='light';
    elseif isequal(pep_library.label(pep_idx,:),'heavy')
        label='heavy';
    else
        flag=false;
        index=[];
        disp('?????????????????????????')
        label='NOT LABELED';
        return
    end
    
    occ_idxs=strmatch(pep_sequence,pep_library.peptidesSequence,'exact');
    if numel(occ_idxs)>0
        flag=true;
        same_charge_occ_idxs=occ_idxs(pep_library.charges(occ_idxs)==pep_library.charges(pep_idx)); 
        [min_expect, min_expect_occ_idx]=min(pep_library.expectation(same_charge_occ_idxs));
        [max_score, max_score_occ_idx]=max(pep_library.score(same_charge_occ_idxs));
        if isequal(min_expect_occ_idx,max_score_occ_idx)
            index=same_charge_occ_idxs(min_expect_occ_idx);
            label=pep_library.label(index,:);
        else
            index=same_charge_occ_idxs(min_expect_occ_idx); 
            label=pep_library.label(index,:);
        end
    end
    
    pepzero=num2str(zeros(numel(pep_sequence),1))';
    pos_Ks=findstr(pep_sequence,'K');
    num_Ks=numel(pos_Ks);
    if isequal(label,'light')
        num_mods=pep_library.num_mod_light(index);
        if num_mods==num_Ks
            pepzero(pos_Ks)='3';
        elseif num_mods>num_Ks
            pepzero(pos_Ks)='3';
            remaining_mod=num_mods-num_Ks;
            free_pos=setxor(pos_Ks,(1:numel(pepzero)));
            for i=1:remaining_mod
                pepzero(free_pos(i))='3';
            end
        elseif num_mods<num_Ks
            pos_Ks_util=pos_Ks(1:num_mods);
            pepzero(pos_Ks_util)='3';
        end
        
    elseif isequal(label,'heavy')
        num_mods=pep_library.num_mod_heavy(index);
        if num_mods==num_Ks
            pepzero(pos_Ks)='2';
        elseif num_mods>num_Ks
            pepzero(pos_Ks)='2';
            remaining_mod=num_mods-num_Ks;
            free_pos=setxor(pos_Ks,(1:numel(pepzero)));
            for i=1:remaining_mod
                pepzero(free_pos(i))='2';
            end
        elseif num_mods<num_Ks
            pos_Ks_util=pos_Ks(1:num_mods);
            pepzero(pos_Ks_util)='2';
        end
    end
    pepzero=['0' pepzero '0'];
    pep_sequence=[pep_sequence pepzero];
end




	
