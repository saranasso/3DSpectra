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
function [pep_scan mz_peak_max_g w_mzi w_mzf w_rti w_rtf ]=peakDetection(Chroms,settings,peak_scan,mz_shift,distr_width,isot_distr,peak_dist)
% it refines the estimate given from the identifications about the position
% of the max peak

k=settings.k;
tolsx=settings.elution_width_left;
toldx=settings.elution_width_right;
lim_rt=settings.elution_identification_error_tolerance;
npeaks=length(isot_distr.prob_peaks);
[sorted_peaks idx_sorted_peaks]=sort(isot_distr.prob_peaks,'descend');
idx_peak_max=idx_sorted_peaks(1);
pep_scan=[];
mz_peak_max_g=[];
numberofscans=size(Chroms.int,2);
flag=true;

% refine window of interest to look for elution time using only 3
% most intense isotopes
k_max=k+(idx_peak_max-1)*peak_dist-1;
lowerBound=k_max-round(peak_dist*(3/2)); % positioning lower bound to one peak before
upperBound=k_max+round(peak_dist*(3/2)); % positioning upper bound to one peak after

time=(peak_scan-lim_rt):(peak_scan+lim_rt);
% for each chromatogram identify putative elution time in the sub-window of
% interest
for i=lowerBound:upperBound
    chrom=Chroms.int(i,:);
    data=msbackadj([1:numberofscans]',chrom','STEPSIZE',50);
    % fit to pick peaks
    try
        y=fit(time',data(time),'gauss1');
    catch
        y=fit(time',data(time),'gauss4');
    end
    [maxValue(i-(lowerBound-1)),idx_maxValue(i-(lowerBound-1))]=max(y(time));
end
peak_scan;
pep_scan_median=median(idx_maxValue);
pep_scan_mode=mode(idx_maxValue);
pep_scan_median=pep_scan_median-1+peak_scan-lim_rt;
pep_scan_mode=pep_scan_mode-1+peak_scan-lim_rt;
if pdist([pep_scan_median; peak_scan])<pdist([pep_scan_mode; peak_scan])
    pep_scan=pep_scan_median;
else
    pep_scan=pep_scan_mode;
end

% for low res data, for high res the ID m/z will be used
w_mzi=k-ceil(peak_dist/2)+round(mz_shift*k); % additional peak_dist/2 as tolerance
if w_mzi<=0
    w_mzi=1;
end
w_mzf=w_mzi+distr_width;
w_rti=pep_scan-tolsx;
w_rtf=pep_scan+toldx;
strip=Chroms.int(w_mzi:w_mzf,pep_scan);
[maxPeak mz_peak_max_l]=max(strip);
mz_peak_max_g=mz_peak_max_l+w_mzi-1;
w_mzi=mz_peak_max_g-ceil((peak_dist)*idx_peak_max+peak_dist/2);  % additional peak_dist as tolerance; first peak maximum at position mzi+1.5*peak_dist, as required by getprior
w_mzf=(w_mzi+(npeaks+2)*peak_dist); % additional peak_dist as tolerance (including compensation for tolerance on wzi)
  
% extreme cases to be better handled in future releases
if w_mzf>size(Chroms.int,1)
    w_mzf=size(Chroms.int,1);
end
if w_mzi<1
    w_mzi=1;
end

%% plot selected window vs whole considered strip
if settings.visualize
     pause on
    idx_peak1=mz_peak_max_g-((idx_peak_max-1)*peak_dist); % absolute position first peak
    idx_peak1_w=mz_peak_max_g-w_mzi-((idx_peak_max-1)*peak_dist);%((peak_dist)*idx_peak_max-peak_dist/2) %mz_peak_max_l-ceil((peak_dist)*idx_peak_max+peak_dist/2);
    int=Chroms.int(w_mzi:w_mzf,w_rti:w_rtf);
    az = 0;
    el = 90;
    hf=figure;
    subplot 211
    surf(int), shading interp
    hold on
    for j=1:npeaks
        stem3(pep_scan-w_rti,idx_peak1_w+(j-1)*peak_dist,maxPeak*isot_distr.i_peaks(j),'g','filled');
    end
    title(['rti=' num2str(w_rti) ' rtf=' num2str(w_rtf) ' mzi=' num2str(w_mzi) ' mzf=' num2str(w_mzf) ' pep scan=' num2str(pep_scan)])
    view(az, el);
    subplot 212
    surf(Chroms.int), shading interp
    hold on
    for j=1:npeaks
        stem3(pep_scan,idx_peak1+(j-1)*peak_dist,maxPeak*isot_distr.i_peaks(j),'g','filled');
    end
    title(['ID peak scan was ' num2str(peak_scan)])
    view(az, el);
    pause(1)
    pause off
    attachPlotToFile(hf,'3DSpectra')
    close all
end

	
