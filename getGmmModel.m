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
function [Chromw gmm_opt tmp_opt s mz_pep_in NlogL noisy_gaussian]=getGmmModel(Chroms,peak1_scan,w_mzi,w_mzf,w_rti,w_rtf,peak_dist,pep_seq,settings,mz_peaks_pep)
% this function fit the GMM on peptide distribution; statistically define a boolean signal mask (in/out) exploiting the GMM and filter the signal to be quantified in next step

cut=settings.cut;
number_of_gaussians=settings.number_of_gaussians;
k=settings.k;

Chromw=Chroms;
Chroms_tmp_int=zeros(size(Chroms.int,1),size(Chroms.int,2));
gmm=cell(1);
tmp=cell(1);
s=cell(1);
AIC=[];
NlogL=[];
BIC=[];

% on window of interest
tmp_chrom=Chroms.int(w_mzi:w_mzf,w_rti:w_rtf);
old_rti=w_rti;
old_rtf=w_rtf;
Chroms.int=tmp_chrom;
%% local coordinates on window
peak1_scan=peak1_scan-w_rti;
w_rtf=w_rtf-w_rti;
w_rti=1;
w_mzf=w_mzf-w_mzi;
w_mzi=1;

%% model fitted with prior info as for mean, sigma and probability for each component.
if ~isempty(number_of_gaussians)
    iterator=number_of_gaussians;
else
    iterator=(1:1:6); % 6 is the max number of isotopes considered if the value is not provided by user
end
j=1;
for i=iterator
    prior=getprior(Chroms,settings,w_mzi,w_mzf,w_rti,w_rtf,peak_dist,i,'MSP_prior',pep_seq,peak1_scan);
    try
        [gmm{j} tmp]=mygmm(Chroms,cut,i,prior);
        AIC(j)=gmm{j}.AIC;
        BIC(j)=gmm{j}.BIC;
        NlogL(j)=gmm{j}.NlogL;
    end
    j=j+1;
end

%% Optimal model choice
% Both the Akaike and Bayes information are negative log-likelihoods for
% the data with penalty terms for the number of estimated parameters.
% They are often used to determine an appropriate number of components
% for a model when the number of components is unspecified.
if ~isempty(AIC)
    [AIC_sort idx]=sort(AIC);
    idxs_sort=find(AIC_sort>0);
    j_ottimo=idx(idxs_sort(1));
    gmm_opt=gmm{j_ottimo};
    tmp_opt=tmp;
else
    gmm_opt=gmm{1};
    tmp_opt=tmp; 
    NlogL=0;
end

[P,nlogl] = posterior(gmm_opt,tmp_opt);
disp('negLogLikelihood')
nlogl
pause(1)
%% fixing also m/z and rt indexes swap with respect to original matrix
idx_pesi(:,1)=tmp_opt(:,2); % indexes m/z
idx_pesi(:,2)=tmp_opt(:,1); % indexes scan (time)
idx_pesi(:,3:(2+number_of_gaussians))=P;
idx_unique=unique(idx_pesi,'rows');
mask=zeros(size(Chroms.int,1),size(Chroms.int,2));

%%  pdf/isopdf estimate
h=(@(x,y)pdf(gmm_opt,[x y])); % pdf function handle
tmp_unique=unique(tmp_opt,'rows');
pdf_est=h(tmp_unique(:,1),tmp_unique(:,2)); % pdf estimate on ions (idxs)
Z=accumarray({tmp_unique(:,2),tmp_unique(:,1)},pdf_est); % matrix storing the pdf estimates indexed by the ion coordinates (idxs)

%% cluster gaussians
idx=cluster(gmm_opt,tmp_unique);
for i=1:number_of_gaussians
    idx_Gcluster{i}=find(idx==i);
    if any(idx_Gcluster{i})
        for it=1:size(idx_Gcluster{i},1)
            c=tmp_unique(idx_Gcluster{i}(it),1);
            r=tmp_unique(idx_Gcluster{i}(it),2);
            gmm_clusters{i}(it,1)=Chroms.int(r,c);
        end   
        gmm_counts(i)=size(idx_Gcluster{i},1);
        gmm_mean_intensities(i)=mean(gmm_clusters{i});
        gmm_median_intensities(i)=median(gmm_clusters{i});
        gmm_maxs(i)=max(gmm_clusters{i});
        gmm_densities(i)=gmm_mean_intensities(i)/gmm_maxs(i);
        gmm_PComponents_densities(i)=gmm_opt.PComponents(i)/gmm_counts(i);
    else
        continue
    end
end
%% gaussian holding the highest number of different ions --> noisy gaussian
[mc most_ions_gauss]=max(gmm_counts);
%% lowest density gaussian --> noisy gaussian
[ld lowest_density_gauss]=min(gmm_densities);
%% highest variance gaussian --> noisy gaussian
[hv highest_variance_gauss]=max(gmm_opt.Sigma,[],3);
highest_variance_gauss_rt=highest_variance_gauss(1,1);
highest_variance_gauss_mz=highest_variance_gauss(2,2);
%% less dense of probability gaussian -->noisy gaussian
[lPCd lowest_PComponents_density_gauss]=min(gmm_PComponents_densities);
%% noisy gaussian
noisy_gaussian=mode([most_ions_gauss,lowest_density_gauss,highest_variance_gauss_rt,highest_variance_gauss_mz,lowest_PComponents_density_gauss]);

%% creating the mask according to 2 conditions: (1-P_noise(i,j)>0.1) && Z(i,j)>0.0001
for kk=1:length(idx_unique)
    i=idx_unique(kk,1);
    j=idx_unique(kk,2);
    P_noise(i,j)=idx_unique(kk,(2+noisy_gaussian));
    
    if  (1-P_noise(i,j)>0.1) && Z(i,j)>settings.pdf_threshold % profile
        mask(i,j)=1;
    else
        mask(i,j)=0;
    end
end

%% first gaussian mean on m/z
mz_peak1=floor(min(gmm_opt.mu(:,2)));
mz_pep_in=k+mz_peak1-1;
masked_Chroms=Chroms.int.*mask;
Chroms_tmp_int(k:(k+w_mzf-w_mzi+1),(1000):1000+(w_rtf-w_rti)+1)=Chroms.int.*mask; %%  1000 arbitrary choice, any greater than the elution width would have fitted
Chromw.int=Chroms_tmp_int;

% % %% VISUALIZATION
if settings.visualize
    az = 0;
    el = 90;
    pause on
    loci_peaks=sort(gmm_opt.mu(setdiff([1:size(gmm_opt.mu,1)],noisy_gaussian),2));
    y_loci_first_peak=loci_peaks(1);
    y_loci=y_loci_first_peak:peak_dist:y_loci_first_peak+peak_dist*(numel(loci_peaks)-1);
    hf=figure;
    subplot 231
    hold on;
    surf(Chroms.int);
    shading interp;
    setYaxis(mz_peaks_pep, old_rti,old_rtf,az,el,y_loci)
    title([{[strtok(pep_seq,'0') '  ']};{'before mask'}])


    subplot 236
    surf(Chroms.int.*mask), shading interp, title('after mask')
    setYaxis(mz_peaks_pep, old_rti,old_rtf,az,el,y_loci)

    subplot 235
    surf(mask), shading interp, title('signal mask')
    setYaxis(mz_peaks_pep, old_rti,old_rtf,az,el,y_loci)

    subplot 232
    surf(Z), shading interp, title('GMM PDF')
    setYaxis(mz_peaks_pep, old_rti,old_rtf,az,el,y_loci)
    subplot 234
    surf(1-P_noise), shading interp, title(['1-P_n_o_i_s_e. Noisy Gaussian: ' num2str(noisy_gaussian)])
    setYaxis(mz_peaks_pep, old_rti,old_rtf,az,el,y_loci)

    subplot 233
    prob_levels=[settings.pdf_threshold,0.0005,0.001];
    [C,h]=contour(Z,prob_levels);
    %     title(['GMM pdf iso levels ' num2str(prob_levels)])
    title('GMM PDF iso-levels')
    setYaxis(mz_peaks_pep, old_rti,old_rtf,az,el,y_loci)

    attachPlotToFile(hf,'3DSpectra')

    pause(1)
    pause off
    close all
end

	
