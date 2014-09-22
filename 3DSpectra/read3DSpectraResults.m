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
function [names_x3 ratio_peaks_all quants_light quants_heavy charges_x3 replicate label_x3 id_x3]=read3DSpectraResults(results)
% it imports results from 3DSpectra. 

i=1;
for k=1:numel(results) 
    pep_seq=results(k).pep_seq;
    if ~isempty(pep_seq)
        [pepname{i,1} pep_label{i,1}]=getlabel(pep_seq);
        id{i}=results(k).pep_idx;
        names{i,1}=results(k).pep_seq;
        replicate{i,1}=results(k).replicate;
        charges{i}=results(k).pep_charge';
        id_to_check_outliers(i,1)=k;
        if pep_label{i}=='light'
            quants_light(i,1)=sum(results(k).spettro_pep);
            quants_heavy(i,1)=sum(results(k).spettro_partner);
        elseif pep_label{i}=='heavy'
            quants_light(i,1)=sum(results(k).spettro_partner);
            quants_heavy(i,1)=sum(results(k).spettro_pep);
        end
        ratio_peaks_all(i,1)=results(k).ratios(1);
        for pp=1:length(results(k).ratios)
            mean_ratio(i,pp)=results(k).ratios(pp)';
        end
        i=i+1;
    end
end
names_x3=cell(size(names,1),1);
k=1;
for i=1:length(names)
        replicate(k,1)=replicate(i);
        names_x3{k,1}=names{i};
        charges_x3(k,1)=charges(i);
        label_x3{k,1}=pep_label{i};
        id_x3{k,1}=id(i);
        k=i+1;
end	
