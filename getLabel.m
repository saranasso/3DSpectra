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
function [pepname pep_label]=getLabel(pep_seq)
% it recognizes the labeling

if ~isempty(pep_seq)
        [pepname pepzero]=strtok(pep_seq,'0');
        pep_heavy=find(pepzero=='2',1);
        pep_light=find(pepzero=='3',1);
        pep_err=find(pepzero=='1',1);
        pep_nz=find(pepzero~='0',1); % redundancy

        if (isempty(pep_heavy) && isempty(pep_light)) || ~isempty(pep_err) || (~isempty(pep_heavy) && ~isempty(pep_light)) || isempty(pep_nz)
            disp('NOT LABELED')
            pep_label='NOT LABELED';
            return
        elseif ~isempty(pep_light) && isempty(pep_heavy)
            pep_label='light';
        elseif ~isempty(pep_heavy) && isempty(pep_light)
            pep_label='heavy';
        else
            disp('????????????????????????????????????')
            pep_label='NOT LABELED';
            return
        end
 else
     disp('Unvalid/empty sequence')
     pepname=[];
     pep_label=[];
%      return
 end	
