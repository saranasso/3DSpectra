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
function [peak_distr]=getIsoDistrb(max_peak_idx,isot_distr,peak_distance,maxima,k)
% it clusters together data belonging to the same peak

weightedmaxima=maxima;
npeaks=length(isot_distr.i_peaks);
[sorted_peaks idx_sorted_peaks]=sort(isot_distr.prob_peaks,'descend');
coeff=round(peak_distance*k); % distance between peaks
idx_peak1=max_peak_idx;

for j=1:npeaks   
    if j==1 || j==npeaks
        if j==1
            p{j}=(((idx_peak1+(j-1)*coeff)-ceil(coeff/2)):((idx_peak1+(j-1)*coeff)+floor(coeff/2)));
            l(j)=idx_peak1-round((coeff)/2);
            if l(j)<1
               l(j)=1; 
            end
        elseif j==npeaks
            p{j}=(((idx_peak1+(j-1)*coeff)-floor(coeff/2)):((idx_peak1+(j-1)*coeff)+floor(coeff/2)));
            l(j+1)=round(idx_peak1+(j-1)*coeff+(coeff/2));
            if l(j+1)>length(weightedmaxima)
                l(j+1)=length(weightedmaxima);
            end
        end
    else
        p{j}=(((idx_peak1+(j-1)*coeff)-floor(coeff/2)):((idx_peak1+(j-1)*coeff)+floor(coeff/2)));
    end

    if j>1
        l(j)=ceil(p{j}(1));
        if l(j)<1
            l(j)=1;
        elseif l(j)>length(weightedmaxima)
            l(j)=length(weightedmaxima);
        end
        peak_distr(j-1)=sum(weightedmaxima(l(j-1):l(j)-1));
    end
end
peak_distr(npeaks)=sum(weightedmaxima(l(npeaks):((l(npeaks+1)-1))));
	
