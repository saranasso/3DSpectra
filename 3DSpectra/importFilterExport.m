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
function [names_x3 charges_x3 mean_ratio quants_light quants_heavy ratio_peaks]=importFilterExport(results,file_xls,ratio,outliers_removal,percentiles)
% it imports results from 3DSpectra and applies to them the outlier removal
% as implemented in MASPECTRAS. Thereafter, data are also exported to xls for
% user convenience. The xls will be further updated in next steps.

%% import 3DSpectra results
[names_x3 ratio_peaks_all quants_light quants_heavy charges_x3 replicate label_x3 id_x3]=read3DSpectraResults(results);
method=char(inputname(1));
idxs_ = findstr('_', method);
method(idxs_)=' ';
ratio_peaks=ratio_peaks_all;

%% outlier removal
if outliers_removal
    if isempty(percentiles)
         [ratio_peaks idxs_no_out idxs_out mean_ratio sd_ratio cv_ratio efficiency]=rmOutliers(ratio_peaks);
    else
         [ratio_peaks idxs_no_out idxs_out mean_ratio sd_ratio cv_ratio efficiency]=rmOutliers(ratio_peaks,percentiles(1),percentiles(2));
    end
    names_x3=names_x3(idxs_no_out);
    charges_x3=charges_x3(idxs_no_out);
    label_x3=label_x3(idxs_no_out);
    quants_light=quants_light(idxs_no_out);
    quants_heavy=quants_heavy(idxs_no_out);
    replicate=replicate(idxs_no_out);
else
    mean_ratio=mean(ratio_peaks);
    sd_ratio=std(ratio_peaks);
    cv_ratio=(sd_ratio/mean_ratio)*100;
    efficiency=100; % all kept
end

%% xls writing
xlswrite(file_xls,names_x3,method,'A1');
xlswrite(file_xls,id_x3,method,'B1');
% xlswrite(file_xls,charges_x3,method,'C1');
xlswrite(file_xls,label_x3,method,'D1');
xlswrite(file_xls,quants_light,method,'E1');
xlswrite(file_xls,quants_heavy,method,'F1');
xlswrite(file_xls,ratio_peaks,method,'G1');
xlswrite(file_xls,replicate,method,'H1');
N_ratio=numel(ratio_peaks);
info=[{'Mean'};{'SD'};{'CV %'};{'# pep'};{'efficiency'}];
stats=[mean_ratio;sd_ratio;cv_ratio;N_ratio;efficiency];
xlswrite(file_xls,{'All'},method,'V1');
xlswrite(file_xls,info,method,'V2:V6');
xlswrite(file_xls,{'3DSpectra'},method,'W1');
xlswrite(file_xls,stats,method,'W2:W6');
	
