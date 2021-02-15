function [SEG,CN]=SegmentCN(meanDepth,ploidy,MAbins,minCNLen)
%% Moving average bins
nBins=size(meanDepth,1);
mBins=ceil(MAbins/2);
bin_id={};
for i=1:nBins
	bin_id{i}=max([1,i-mBins]):min([nBins,i+mBins]);
end
%% Round read depth to nearest integer states
COV1=meanDepth*ploidy;
COV0=round(COV1);
%% Segmentation
delta=sum(COV1~=COV0,1);
while max(delta)>0
    for i=1:nBins
        COV1(i,:)=round(mean(COV0(bin_id{i},:),1,'omitnan'));
    end
    delta=sum(COV1~=COV0,1);
    COV0=COV1;
end
CN=COV0;
%% Smoothing
up=find(diff(CN)>0); % left < right copy number changepoint
down=find(diff(CN)<0); % left > right copy number change point
% get rid of intermediate CN states due to moving window average
% there should be no segments (either gain or loss) shorter than mBins
% given the moving average
spacing=min(diff(up));
while spacing<=mBins
    [~,id]=min(diff(up));
    up(id)=round(median([up(id),up(id+1)]));
    up=up([1:id,id+2:end]);
    spacing=min(diff(up));
end
spacing=min(diff(down));
while spacing<=mBins
    [~,id]=min(diff(down));
    down(id)=round(median([down(id),down(id+1)]));
    down=down([1:id,id+2:end]);
    spacing=min(diff(down));
end
% get rid of short segments
bkps=sortrows(dataset([up;down],[ones(size(up));-ones(size(down))],'VarNames',{'changeIdx','change'}));
spacing=min(diff(bkps.changeIdx));
while spacing<minCNLen 
    [~,id]=min(diff(bkps.changeIdx));	
    if id>1 
        leftCN=round(mean(CN(bkps.changeIdx(id-1):bkps.changeIdx(id)),'omitnan'));
    else
        leftCN=round(mean(CN(1:bkps.changeIdx(id)),'omitnan'));
    end
    if id<length(bkps)-1
        rightCN=round(mean(CN(bkps.changeIdx(id+1):bkps.changeIdx(id+2)),'omitnan'));
    else
        rightCN=round(mean(CN(bkps.changeIdx(id+1):end),'omitnan'));
    end
    if leftCN==rightCN
        bkps=bkps([1:id-1,id+2:end],:);	
    else
        bkps=[bkps([1:id-1,id+2:end],:);dataset(round(mean([bkps.changeIdx(id),bkps.changeIdx(id+1)])),(rightCN>leftCN)-(rightCN<leftCN),'VarNames',{'changeIdx','change'})];
        bkps=sortrows(bkps,{'changeIdx'});
    end
    spacing=min(diff(bkps.changeIdx));
end
% segments	
SEG=[];
COV=meanDepth*ploidy;
if length(bkps.changeIdx)>=1
    cov=mean(COV(1:bkps.changeIdx(1)-1),'omitnan');
    SEG=dataset(1,bkps.changeIdx(1),cov,round(cov),'VarNames',{'Start','End','AvgDepth','CN'});
    for j=1:length(bkps)-1	
        cov=mean(COV(bkps.changeIdx(j):bkps.changeIdx(j+1)-1),'omitnan');
        SEG(j+1,:)=dataset(bkps.changeIdx(j),bkps.changeIdx(j+1),cov,round(cov),'VarNames',{'Start','End','AvgDepth','CN'});
    end
    cov=mean(COV(bkps.changeIdx(end):end),'omitnan');
    SEG=[SEG;dataset(bkps.changeIdx(end),length(COV),cov,round(cov),'VarNames',{'Start','End','AvgDepth','CN'})];
else
	SEG=dataset(1,length(COV),mean(COV,'omitnan'),round(mean(COV,'omitnan')),'VarNames',{'Start','End','AvgDepth','CN'});

end	
