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
function pep=jrapDTquery(pep,rti,rtf,mzi,mzf)
% it allows to perform range queries on data read fby means of JRAP

%% range query
mzi_trans=round((mzi-pep.mz_min)./pep.mz_res+1);
mzf_trans=round((mzf-pep.mz_min)./pep.mz_res+1);
rti_trans=round((rti-pep.rtt_min)./pep.rtt_res)+1;
rtf_trans=round((rtf-pep.rtt_min)./pep.rtt_res)+1;
pep.query=pep.counts(mzi_trans:mzf_trans,rti_trans:rtf_trans);


	
