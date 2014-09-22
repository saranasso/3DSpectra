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
function [Spectra]=getSpectra(Chroms,mz_r,isot_distr,peak_distance,k)
% it fits the isotopic distribution model on each spectrum using WLLS

max_peak_idx=mz_r;
mz_data=1:size(Chroms.int,1); %  # rows = # bin m/z
for i=1:size(Chroms.int,2) % for each spectrum (column)
    raw_spectrum=Chroms.int(mz_data,i);
    idx_int=find(raw_spectrum~=0);
    if numel(idx_int)>2
        try
            raw_spectrum=msbackadj(mz_data',raw_spectrum);
        end
    end
    y=mssgolay(mz_data',raw_spectrum,'Span',5,'Degree',4,'Showplot',false); % savitzky golay smoothing
    [peak_distr]=getIsoDistrb(max_peak_idx,isot_distr,peak_distance,(y(mz_data))',k);
    [I1,stdx,mse]=lscov(isot_distr.i_peaks',peak_distr');
    spettro_fitted=isot_distr.i_peaks*I1;
    Spectra(:,i)=spettro_fitted';
end	
