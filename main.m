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
% Please read the READMEbySara.txt before trying to run this main.
mzRTree_utils.m
clear
close all
clc

[file path]=uigetfile('*.mat', 'Select the matlab file storing the data info'); % created from mzRTree_utils.m : dataRead.mat
load([path file])
dataDir = uigetdir(pwd, 'Select the data files directory'); % directory where data are stored in mzRTree format
clearvars -except x* mz_res* dataDir

settings.npeaks=4;
settings.diff_mass=6.02; % Da
settings.tosum=0; % binning
settings.number_of_gaussians=settings.npeaks+1;
settings.cut=100;
settings.pdf_threshold=0.0001;
settings.ratio='1to1';
settings.visualize=0;
settings.elution_width_left=25;
settings.elution_width_right=35;
settings.elution_identification_error_tolerance=10;

try
    settings.ratio='1to1';
    main1to1
    save 1to1
    pause(0.1)
    clearvars -except x* mz_res* dataDir settings
    pause(0.1)
catch
end

try
    settings.ratio='1to2';
    main1to2
    save 1to2
    pause on
    pause(0.1)
    clearvars -except x* mz_res* dataDir settings
    pause(0.1)
catch
end

try
    settings.ratio='1to5';
    main1to5
    save 1to5
    pause on
    pause(0.1)
    clearvars -except x* mz_res* dataDir settings
    pause(0.1)
catch
end

try
    settings.ratio='1to10';
    main1to10
    save 1to10
    pause on
    pause(0.1)
    clearvars -except x* mz_res* dataDir settings
    pause(0.1)
catch
end


try
    settings.ratio='2to1';
    main2to1
    save 2to1
    pause on
    pause(0.1)
    clearvars -except x* mz_res* dataDir settings
    pause(0.1)
catch
end


try
    settings.ratio='5to1';
    main5to1
    save 5to1
    pause on
    pause(0.1)
    clearvars -except x* mz_res* dataDir settings
    pause(0.1)
catch
end

try
    settings.ratio='10to1';
    main10to1
    save 10to1
    pause on
    pause(0.1)
    clearvars -except x* mz_res* dataDir settings
catch
end

try
    main_results
catch
    disp('Switch to Windows please')
end
