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

for i=1:numel(litohe_ratios)
    min_quants_light_ASAPRatio(i)=min(common_quants(i).light_J_all);
    min_quants_light_3DSpectra(i)=min(common_quants(i).light_S_all);
    min_quants_heavy_ASAPRatio(i)=min(common_quants(i).heavy_J_all);
    min_quants_heavy_3DSpectra(i)=min(common_quants(i).heavy_S_all);
    
    max_quants_light_ASAPRatio(i)=max(common_quants(i).light_J_all);
    max_quants_light_3DSpectra(i)=max(common_quants(i).light_S_all);
    max_quants_heavy_ASAPRatio(i)=max(common_quants(i).heavy_J_all);
    max_quants_heavy_3DSpectra(i)=max(common_quants(i).heavy_S_all);
    
    dynamicRange_3DSpectra_light(i)= log10(max_quants_light_3DSpectra(i)/min_quants_light_3DSpectra(i));
    dynamicRange_3DSpectra_heavy(i)= log10(max_quants_heavy_3DSpectra(i)/min_quants_heavy_3DSpectra(i));
    dynamicRange_ASAPRatio_light(i)= log10(max_quants_light_ASAPRatio(i)/min_quants_light_ASAPRatio(i));
    dynamicRange_ASAPRatio_heavy(i)= log10(max_quants_heavy_ASAPRatio(i)/min_quants_heavy_ASAPRatio(i));
    
end

dynamicRangeEstimate=[dynamicRange_3DSpectra_light;dynamicRange_ASAPRatio_light;dynamicRange_3DSpectra_heavy;dynamicRange_ASAPRatio_heavy];
dynamicRange_percentDiff=[(dynamicRange_3DSpectra_light-dynamicRange_ASAPRatio_light)./max(dynamicRange_3DSpectra_light,dynamicRange_ASAPRatio_light)*100;(dynamicRange_3DSpectra_heavy-dynamicRange_ASAPRatio_heavy)./max(dynamicRange_3DSpectra_heavy,dynamicRange_ASAPRatio_heavy)*100];
gain=[(dynamicRange_3DSpectra_light-dynamicRange_ASAPRatio_light)./dynamicRange_ASAPRatio_light*100;(dynamicRange_3DSpectra_heavy-dynamicRange_ASAPRatio_heavy)./dynamicRange_ASAPRatio_heavy*100];
mean_gain=mean(mean(gain));
disp(['The average gain in dynamic range is ' num2str(mean_gain)])
	
