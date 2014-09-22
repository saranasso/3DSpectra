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
function [flag partner mzpartner]=findPartner(pep_sequence,pep_library,index,label,diff_mass)
% It returns partner m/z position and sequence and look for it in the pepLibrary returning true if partner is in library


flag=false;

pep_charge=pep_library.charges(index);
peak_mz=(pep_library.massH1(index)+1)/pep_charge;

[pepseq pepzero]=strtok(pep_sequence,'0');



pep_heavy=find(pepzero=='2');
pep_light=find(pepzero=='3');
if ~isempty(pep_light) && isequal(label,'light')
    shift=length(pep_light)*diff_mass/pep_charge;
    pepzero(pep_light)='2';
    partnerzero=pepzero;
elseif ~isempty(pep_heavy) && isequal(label,'heavy')
    shift=-length(pep_heavy)*(diff_mass/pep_charge);
    pepzero(pep_heavy)='3';
    partnerzero=pepzero;
elseif isequal(label,'NOT LABELED')
    flag=false;
    partner=[];
    mzpartner=[];
    shift=[];
    return
end

mzpartner=peak_mz(1)+shift;
partner=[pepseq partnerzero];

putative_idx_Partner=strmatch(partner,pep_library.peptidesSequence,'exact');
index_partner=find(pep_library.charges(putative_idx_Partner), pep_charge);
if ~isempty(index_partner)
    flag=true;
end

	
