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
function [int]=plotGmm(gmm,tmp,pep_seq,varargin)
% it helps visualizing gmm

if ~isempty(varargin)
    method=varargin{:};
else
    method='';
end
y=pdf(gmm,tmp);
idx_weights(:,1)=tmp(:,2); % indexes m/z 
idx_weights(:,2)=tmp(:,1); % indexes scan (time)
idx_weights(:,3)=y; % weight from pdf
idx_unique=unique(idx_weights,'rows');
% number_of_gaussians=gmm.NComponents;
int=sum(idx_unique(:,3));
xi=min(idx_weights(:,2));
xf=max(idx_weights(:,2));
yi=min(idx_weights(:,1));
yf=max(idx_weights(:,1));
h=figure;
subplot(2,1,1)
ezsurf(@(x,y)pdf(gmm,[x y]),[xi,xf],[yi,yf])
title([{pep_seq};{['GMM: ' num2str(gmm.NComponents) ' gaussians, EM starting from ' method]}])
subplot(2,1,2)
ezcontour(@(x,y)pdf(gmm,[x y]),[xi,xf],[yi,yf]);
title([{pep_seq};{['GMM: ' num2str(gmm.NComponents) ' gaussians, EM starting from ' method]}])
close all
	
