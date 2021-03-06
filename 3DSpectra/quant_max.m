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
function max_Spectra=quant_max(Spectra)
% it averages the top 33% of data to evaluate the max of the signal

for i=1:size(Spectra,1)
    tmp=Spectra(i,:);
    tmp=tmp((tmp~=0));
    N=numel(tmp);
    idx=ceil((33*N)/100); % 33%
    tmp_sorted=sort(tmp,'descend');
    tmp=(tmp_sorted(1:idx));
    max_Spectra(i)=mean(tmp);
end	
