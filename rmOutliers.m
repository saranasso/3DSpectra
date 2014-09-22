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
function [ratios_no_out idxs_no_out idxs_out r_mean r_std r_CV efficiency]=rmOutliers(ratios,varargin)
%% it performs the outlier removal

if isstruct(ratios)
    ratios=[ratios.ratios];
end
if  size(ratios,2)>1
    ratios=ratios(:);
end
l=length(ratios);

if ~isempty(varargin) && ~isempty(varargin{1}) && ~isempty(varargin{2})
    %% simple
    pl=varargin{1};
    pu=varargin{2};
    y = prctile(ratios,[pl pu]);
    idxs_up=find(ratios<y(2));
    idxs_low=find(ratios>y(1));
    idxs_no_out=intersect(idxs_low,idxs_up);
else
    %% MASPECTRAS version
    ratios_inv = ratios.^(-1);
    idxs=Outliers_Juergen.removeOutliers(ratios,ratios_inv); % idxs is equal to 1 when is an outlier, 0 otherwise
    idxs_no_out=find(idxs==0);
end

%% group inlier values
ratios_no_out=ratios(idxs_no_out);
r_mean=mean(ratios_no_out);
r_std=std(ratios_no_out);
r_CV=r_std/r_mean*100;
efficiency=((numel(ratios_no_out))/l)*100;
%% group outlier values
all_idxs=1:l;
idxs_out=setxor(idxs_no_out,all_idxs);
	
