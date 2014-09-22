%  This main is meant to be run after the 3DSpectra algorithm to compare
%  its results to those from MASPECTRAS as stored in the provided results folder.

warning off
%% names of the results .mat files to be loaded
litohe_ratios{1}='1to1';
litohe_ratios{2}='1to2';
litohe_ratios{3}='1to5';
litohe_ratios{4}='1to10';
litohe_ratios{5}='2to1';
litohe_ratios{6}='5to1';
litohe_ratios{7}='10to1';

%% settings for results statistics
outliers_removal=0;  % 0=results including outliers, 1=results excluding outliers
regression_all=1;
percentiles=[];
hasReplicates=1;
expected_ratios=[1,0.5,0.2,0.1,2,5,10];
resultsDir=pwd;
if outliers_removal
    mkdir(pwd,'without outliers');
    cd('without outliers')
else
    mkdir(pwd,'with outliers');
    cd('with outliers')
end

for iter=1:length(litohe_ratios)
    ratio=litohe_ratios{iter}
    try
        load([resultsDir filesep ratio])
    catch
        disp('Did you cd to the results folder?')
        return
    end
    file_xls_J=[resultsDir filesep 'Maspectras' filesep 'outliers_included' filesep ratio];
    results_name=[ratio  '_analyzed'];
    [results_quants]=uniquePeptideRatio(results_quants);
    pause(0.1)
    [names_mst charges_mst mean_ratio_mst quants_light_mst{iter} quants_heavy_mst{iter} ratio_peaks_mst{iter}]=importFilterExport(results_quants,results_name,litohe_ratios{iter},outliers_removal,percentiles);
    [common_quants(iter)]=commons(results_name,results_quants,file_xls_J,'Quantification Comparison',outliers_removal,hasReplicates);
    [stats(iter)]=regression_norm_simm(common_quants(iter),ratio,regression_all);
    %% saving
    save(results_name)
    pause(0.001)
end

% return
%% dynamic range
dynamicRange
save analyzed_results.mat
cd ..