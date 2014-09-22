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
function [pep_quant_j quantj_all lightj_all heavyj_all]=readMASPECTRASresults(file_xls_J,Foglio_J,hasReplicates)
% it imports MASPECTRAS results

[numj,textj,rawj]=xlsread(file_xls_J,Foglio_J);

for i=2:size(textj,1)
    pep_quant_juerg{i-1,1}=strtok(textj{i,2}, '');
    pep_quant_juerg{i-1,1}=pep_quant_juerg{i-1,1}(find(pep_quant_juerg{i-1,1}~='.'));
    
    if ~isempty(textj{i,1}) && (strcmp(textj{i,1}(1),'A') || strcmp(textj{i,1}(1),'I'))
        replicate_juerg{i-1,1}=strtok(textj{i,1}(end-3), '');
    end
    
end

for i=1:size(pep_quant_juerg,1)
    pep_quant_j{i,1}=pep_quant_juerg{i,1}(pep_quant_juerg{i,1}~='*');
end

null_idxs=ismember(pep_quant_j,' ');
[pep_quant_j_tmp idxs]=setdiff(pep_quant_j,' ');
pep_quant_j=pep_quant_j(null_idxs==0);

if hasReplicates
    for i=1:length(pep_quant_j)
        pep_quant_j{i,1}=[pep_quant_j{i,1} replicate_juerg{i,1}];
    end
end

heavyj_all=[rawj{2:end,4}]';
idxs_to_keep= intersect(find(~isnan(heavyj_all)),find((heavyj_all~=0)));
heavyj_all=heavyj_all(idxs_to_keep);
quantj_all=[rawj{idxs_to_keep+1,7}]';
lightj_all=[rawj{idxs_to_keep+1,5}]';
pep_quant_j=pep_quant_j(idxs_to_keep);
[pep_quant_j idxs_to_keep idxs_useless]=unique(pep_quant_j,'first','rows');
heavyj_all=heavyj_all(idxs_to_keep);
quantj_all=quantj_all(idxs_to_keep);
lightj_all=lightj_all(idxs_to_keep);
	
