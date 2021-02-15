function [bkps,supp]=collectSampleSupportShortSVs(clusterRanges,SuppReads,SampleInfo,insertSize,readLength)
noBkps=4*length(clusterRanges);
CHR1=repmat({'NA'},noBkps,1);
CHR2=CHR1;
POS1=zeros(noBkps,1);
POS2=POS1;
STR1=POS1;
STR2=POS1;
% best mapping quality of all discordant reads
MAQ=zeros(noBkps,1);
% total number of non-primary (split) alignments
NPCount=zeros(noBkps,1);
% total number of support discordant reads
TotalCount=zeros(noBkps,1);
% number of support discordant reads in each sample
Samples=unique(SampleInfo.SampleNames);
SampleSuppCount=zeros(noBkps,length(Samples));
% number of support discordant reads in each sample group
Groups=unique(SampleInfo.SampleGroupID);
%SampleGroupSuppCount=zeros(noBkps,length(Groups));
% number of samples showing discordant support
SuppSamples=zeros(noBkps,1);
% Samples with maximum support
SVSample=repmat({'NA'},noBkps,1);
% Group with maximum support
SVGroup=repmat({'NA'},noBkps,1);
% number of supporting discordant pairs in control samples (SampleGroup=0)
NMCount=zeros(noBkps,1);
% maximum number of supporting reads in a single sample
MaxSampleCount=zeros(noBkps,1);
% maximum sum of supporting reads from a single sample group
GCount=zeros(noBkps,1);

SVidx=0;
supp=[];
for i=1:length(clusterRanges)
	currBkp=clusterRanges(i,:);
	currBkp.pos1_left=currBkp.pos1_left-3*insertSize;currBkp.pos1_right=min(currBkp.pos1_right+3*insertSize,currBkp.pos2_left);
	currBkp.pos2_left=max(currBkp.pos1_right,currBkp.pos2_left-3*insertSize);currBkp.pos2_right=currBkp.pos2_right+3*insertSize;
	suppReads=SuppReads(strcmp(SuppReads.chr1,currBkp.chr1) & strcmp(SuppReads.chr2,currBkp.chr2),:);
	Supp=collectReadsByRange(suppReads,currBkp);	
	str1=currBkp.str1;str2=currBkp.str2;
	Supp=Supp(Supp.str1==str1 & Supp.str2==str2,:);		
	if size(Supp,1)<2
		continue;
	end	
	split=Supp(Supp.split,:);
	nonsplit=Supp(~Supp.split,:);
	if length(nonsplit)>length(split) % More non-split read pairs than split reads
		% eliminate potential bad split alignments
		idx=find((split.pos1>median(nonsplit.pos1) & (str1==1) | split.pos1<median(nonsplit.pos1) & (str1==-1)) & ...
				 (split.pos2>median(nonsplit.pos2) & (str2==1) | split.pos2<median(nonsplit.pos2) & (str2==-1)));
		if length(idx)>0
			pos1=round(median(split.pos1(idx)));
			pos2=round(median(split.pos2(idx)));
		else
			pos1=(str1==1).*max(double(nonsplit.pos1))+(str1==-1).*min(double(nonsplit.pos1));	
			pos2=(str2==1).*max(double(nonsplit.pos2))+(str2==-1).*min(double(nonsplit.pos2));	
		end
	else 
		pos1=round(median(split.pos1));
		pos2=round(median(split.pos2));
	end

	if str1==1
		pos1_left=pos1-6*insertSize;pos1_right=pos1+0.5*readLength;
	else
		pos1_left=pos1-0.5*readLength;pos1_right=pos1+6*insertSize;
	end
	if str2==1
		pos2_left=pos2-6*insertSize;pos2_right=pos2+0.5*readLength;
	else
		pos2_left=pos2-0.5*readLength;pos2_right=pos2+6*insertSize;
	end
	currRange=dataset(pos1_left,pos1_right,pos2_left,pos2_right,'VarNames',{'pos1_left','pos1_right','pos2_left','pos2_right'});
	Reads=collectReadsByRange(Supp,currRange);	
	Reads=Reads(abs(Reads.pos1-pos1)+abs(Reads.pos2-pos2)<=6*insertSize,:);	
	if (length(Reads)<2)
		continue;
	end	
	SVidx=SVidx+1;
	%pos1,pos2,currBkp.chr1,currBkp.chr2,str1,str2
	CHR1{SVidx}=currBkp.chr1;CHR2{SVidx}=currBkp.chr2;STR1(SVidx)=str1;STR2(SVidx)=str2;
	POS1(SVidx)=pos1;POS2(SVidx)=pos2;
	currSupp=sortrows(Reads,{'SampleGroup','readID'});	
	MAQ(SVidx)=max(min([currSupp.maq1,currSupp.maq2],[],2));
	NMCount(SVidx)=length(unique(currSupp.readID(currSupp.SampleGroup==0)));
	NPCount(SVidx)=sum(currSupp.split);
	sampleCounts=zeros(1,length(SampleInfo));
	for si=1:length(SampleInfo)
		singleCount(si)=length(unique(currSupp.readID(currSupp.SampleGroup==si)));
	end
	for si=1:length(Samples)
		sampleCount(si)=sum(singleCount(strcmp(SampleInfo.SampleNames,Samples{si})));
	end
	for gi=1:length(Groups)
		groupCount(gi)=sum(singleCount(strcmp(SampleInfo.SampleGroupID,Groups{gi})));
	end
	TotalCount(SVidx)=sum(sampleCount);	
	SuppSamples(SVidx)=sum(sampleCount>0);
	MaxSampleCount(SVidx)=max(sampleCount);
	SVSample(SVidx)=Samples(find(sampleCount==max(sampleCount),1,'first'));
	SVGroup(SVidx)=Groups(find(groupCount==max(groupCount),1,'first'));
	SampleSuppCount(SVidx,:)=sampleCount;
	supp=[supp;dataset(repmat(SVidx,size(currSupp.chr1)),'VarNames','SVidx'),currSupp];
	if (SVidx==floor(SVidx/1000)*1000)
		fprintf(1,'..%d',SVidx);
	end
end

% Create dataset
bkps=dataset(SVSample,SVGroup,CHR1,POS1,STR1,CHR2,POS2,STR2,...
	TotalCount,NPCount,MAQ,SuppSamples,MaxSampleCount,NMCount,'VarNames',...
	{'Sample','Group','chr1','pos1','str1','chr2','pos2','str2',...
	'TotalCount','SplitCount','maq','suppSamples','TCount','NCount'});
for k=1:length(Samples)
	bkps=[bkps,dataset(SampleSuppCount(:,k),'VarNames',Samples{k})];
end
bkps=bkps(1:SVidx,:);

function Reads1=collectReadsByRange(reads,range)
Reads1=reads(reads.pos1>=range.pos1_left & reads.pos1<=range.pos1_right & reads.pos2>=range.pos2_left & reads.pos2<=range.pos2_right,:);
