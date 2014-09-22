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
function [results]=getresults(spettro_pep,spettro_partner,pep_seq,pep_charge,pep_idx,label,replicate)
% it computes the final ratio and stores it and additional info in the
% results variable

if strcmp(label,'light')
    ratio=sum(spettro_pep)/sum(spettro_partner);
elseif strcmp(label,'heavy')
    ratio=sum(spettro_partner)/sum(spettro_pep);
end

ratios(1,:)=ratio;
results.pep_seq=pep_seq;
results.pep_charge=pep_charge;
results.pep_idx=pep_idx;
results.ratios=ratios;
results.spettro_pep=spettro_pep;
results.spettro_partner=spettro_partner;
results.replicate=replicate;

	
