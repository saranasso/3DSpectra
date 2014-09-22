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
function [obj tmp varargout]=mygmm(Chroms,cut,number_of_gaussians,varargin)

if ~isempty(varargin)
   prior=varargin{1};
else   
    prior=[];
end
dati_filt=Chroms.int./cut;
mz_max=size(dati_filt,1);
data=dati_filt(:);
tot_counts=sum(round(data(data>1)));
if ~isempty(varargin) && isequal(size(Chroms.int),size(prior)) % if data prior
    prior_srot=prior(:);
    s=[];
end

options = statset('Display','final','MaxIter',100);
if ~isempty(varargin) && isequal(size(Chroms.int),size(prior)) % data prior -deprecated 
    disp('NEXT')
    [tmp s]=hist2data(data, mz_max,prior); %from histogram to occurrences
    varargout{1}=s;
    try
    obj = gmdistribution.fit(tmp,number_of_gaussians,'Start',s,'Options',options);
    disp('NLOGL norm')
    obj.NlogL/tot_counts
    pause(1)
    catch
        disp('PRIOR')
        obj = gmdistribution(prior.mu, prior.Sigma, prior.PComponents);
    end
elseif ~isempty(varargin) && ~isequal(size(Chroms.int),size(prior)) %  prior MSP
    [tmp]=hist2data(data, mz_max);
    try
    obj = gmdistribution.fit(tmp,number_of_gaussians,'Start',prior,'Options',options);
    disp('NLOGL norm')
    obj.NlogL/tot_counts
    pause(1)
    catch
        disp('PRIOR')
        obj = gmdistribution(prior.mu, prior.Sigma, prior.PComponents);
    end
else % no available prior
    [tmp]=hist2data(data, mz_max);
    try
    obj = gmdistribution.fit(tmp,number_of_gaussians,'CovType','diagonal','Replicates',10,'SharedCov',false,'Options',options);
    disp('NLOGL norm')
    obj.NlogL/tot_counts
    pause(1)
    catch
        disp('PRIOR')
        obj = gmdistribution(prior.mu, prior.Sigma, prior.PComponents);
    end
end
	
