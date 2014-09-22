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
function pep=jrapAccess(filePath,queryBounds,binstep,onlyInfo)
% it allows to access the mzXML file and to perform a range query and
% retrieve information

parser = javaObject('org.systemsbiology.jrap.stax.MSXMLParser', filePath);
% error=sprintf('Error!Right Paramaters are: \n filePath \n filePath,start_scan,end_scan \n filePath,start_scan,end_scan,mzi,mzf');
j=1;

if isempty(queryBounds.startScan) && isempty(queryBounds.endScan)
    
    queryBounds.endScan=parser.getMaxScanNumber();
    queryBounds.MS1=0;
    queryBounds.startScan=1;
    
end

if queryBounds.MS1==0
    
    for scan_iter=queryBounds.startScan:queryBounds.endScan
        current_scan= parser.rap(scan_iter);
        if isempty(current_scan)
            continue
        else
            msLevel=current_scan.getHeader().getMsLevel();
            if msLevel==1
                retentionTime = current_scan.getHeader().getRetentionTime();
                RT{j} = java.lang.Float.parseFloat(retentionTime.substring(2, retentionTime.length() - 1));
                peakList{j,1} = current_scan.getMassIntensityList();
                j=j+1;
            end
        end
    end
    
elseif queryBounds.MS1==1 % scan number is relative to MS1, not absolute to all levels
    
    scan_iter=0;
    scan_counter_MS1=1;
    while scan_counter_MS1<=queryBounds.endScan
        scan_iter=scan_iter+1;
        current_scan= parser.rap(scan_iter);
        if isempty(current_scan)
            continue
        else
            msLevel=current_scan.getHeader().getMsLevel();
            if msLevel==1
                if scan_counter_MS1>=queryBounds.startScan
                    retentionTime = current_scan.getHeader().getRetentionTime();
                    RT{j} = java.lang.Float.parseFloat(retentionTime.substring(2, retentionTime.length() - 1));
                    peakList{j,1} = current_scan.getMassIntensityList();
                    j=j+1;
                end
                scan_counter_MS1=scan_counter_MS1+1;
            end
        end
    end
    
end



%%  rtt indexes calculation
rtts=([RT{:}])';
rtt_diff=diff(rtts);
rtt_min_diff=min(rtt_diff(rtt_diff>0)); % factor for minimum distance
rounded_rtts=round((rtts-rtts(1))./rtt_min_diff+1); %add 1 otherwise 0 as an index, it works as well for rtt_min_diff>1
pep.rtt_min=rtts(1);
pep.rtt_max=rtts(end);
pep.rtts=rtts;
pep.rounded_rtts=rounded_rtts;
pep.rtt_res=rtt_min_diff;
scans=1:size(peakList,1);
pep.scan_min=1;
pep.scan_max=size(peakList,1);
pep.scan_res=1;

%%  m/z indexes calculation
prova=cell2mat(peakList');
mzs=prova(1,:);
mzs_diff=diff(mzs);
mz_min_diff=min(mzs_diff(mzs_diff>0));
mz_min=(min(mzs));
mz_max=(max(mzs));
rounded_mzs=round((mzs-mz_min)./mz_min_diff+1);
pep.mz_min=mz_min;
pep.mz_max=mz_max;
pep.mz_res=mz_min_diff;


if ~onlyInfo
    %%counts to indexes association
    counts=[];
    mzs_rounded=[];
    rtts_round_curr=[];
    scans_curr=[];
    for i=1:length(peakList)
        l=size(peakList{i,1},2);
        rtts_round_curr=[rtts_round_curr,rounded_rtts(i)*ones(1,l)];
        scans_curr=[scans_curr,scans(i)*ones(1,l)];
        clear mzs
    end
    counts=prova(2,:);
    if ~queryBounds.scan_representation
        pep.counts=accumarray({rtts_curr,rounded_mzs},counts,[],[],0,true)';  
    elseif queryBounds.scan_representation
         pep.counts=accumarray({scans_curr,rounded_mzs},counts,[],[],0,true)';
    end
    pep.rounded_rtts=rtts_round_curr;
    pep.scans=scans_curr;
    pep.rounded_mzs=rounded_mzs;
    pep.counts_ints=counts;
    
    %% range query
    if ~isempty(queryBounds.mzi) && ~isempty(queryBounds.mzf)
        mzi_trans=round((queryBounds.mzi-mz_min)./mz_min_diff+1);
        mzf_trans=round((queryBounds.mzf-mz_min)./mz_min_diff+1);
    else
        mzi_trans=1;
        mzf_trans=size(pep.counts,1);
    end
    if ~isempty(queryBounds.rti) && ~isempty(queryBounds.rtf) && ~queryBounds.scan_representation
        rti=queryBounds.rti;
        rtf=queryBounds.rtf;
        rti_trans=round((rti-pep.rtt_min)./pep.rtt_res+1);
        rtf_trans=round((rtf-pep.rtt_min)./pep.rtt_res+1);
        pep.query=pep.counts(mzi_trans:mzf_trans,rti_trans:rtf_trans);
    else
        pep.query=pep.counts(mzi_trans:mzf_trans,:); 
    end
    
    if ~isempty(binstep) % m/z binning
        tosum=binstep;
        tmp=bin(pep.query,tosum);
        pep.query=tmp.int;
    end
end


	
