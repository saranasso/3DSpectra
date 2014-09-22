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
function [prior]=getprior(Chroms,settings,w_mzi,w_mzf,w_rti,w_rtf,peak_dist,number_of_gaussians,method,varargin)
% it computes the prior to be passed to the fit.gmdistribution function

if  strcmp(method, 'data_prior') % deprecated
    prior=zeros(size(Chroms.int,1),size(Chroms.int,2));
    for i=1:number_of_gaussians-1
        prior(w_mzi+(i-1)*peak_dist:(w_mzi+i*peak_dist-1),w_rti:w_rtf)=i;
    end
    prior(1:w_mzi-1,:)=number_of_gaussians;
    prior(w_mzi+i*peak_dist:end,:)=number_of_gaussians;
    prior(:,1:w_rti-1)=number_of_gaussians;
    prior(:,w_rtf+1:end)=number_of_gaussians;

elseif strcmp(method, 'MSP_prior')
    path1=[pwd filesep 'igbbiojava.jar'];
    path2=[pwd filesep 'maspectras.jar'];
    javaaddpath(path1)
    javaaddpath(path2)
    pep_seq=varargin{1};
    isot_distr=distribution(pep_seq,number_of_gaussians);
    prior.PComponents=isot_distr.prob_peaks;
    peak1_mz=w_mzi+floor(peak_dist/2)+floor(peak_dist);
    prior.mu=[];
    prior.Sigma=[];
    peak1_scan=varargin{2};
    for i=1:number_of_gaussians -1
        prior.mu=[prior.mu; peak1_scan peak1_mz+(i-1)*peak_dist]; % gaussian means
        prior.Sigma=cat(3,prior.Sigma,[(settings.elution_width_left+settings.elution_width_right)/6*isot_distr.i_peaks(i) 0;0 peak_dist/2]); % divided by 6 because 99,7% lay in 3*sigma per each tail=6*sigma
    end
    prior.mu=[prior.mu; peak1_scan peak1_mz+(number_of_gaussians-1)*peak_dist];
    prior.Sigma=cat(3,prior.Sigma,[10000 0;0 900]); % appending prior for noisy gaussian
end	
