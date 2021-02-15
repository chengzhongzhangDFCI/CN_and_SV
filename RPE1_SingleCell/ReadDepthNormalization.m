addpath('/singlecellcenter/RPE-1/Analysis/CN_analysis');
recurrent_bias_coverage_file='SI_30x_MedianCov_10kb.mat';

%% Read depth processing
Depth_folder='ReadDepth';
SampleCov=[];
if ~exist('ReadCount_10kb.mat')
	for i=1:length(SampleInfo.SampleNames)
		Sample=SampleInfo.SampleNames{i};
		fprintf(1,'Processing 10kb read coverage for sample %s..',Sample);
		fid=fopen(fullfile(Depth_folder,[Sample '.10kb_counts.txt']),'r');	
		while fgets(fid,1)=='@'
			fgetl(fid); % header
		end	
		%temp=textscan(fid,'%s\t%d\t%d\t%d\n','headerLines',128770); % 128770 is the number of header lines
		temp=textscan(fid,'%s\t%d\t%d\t%d\n','headerLines',1,'CommentStyle','@');
		cov=dataset(temp{1},temp{2},temp{3},temp{4},'VarNames',{'Chr','Start','End','count'});
		cov=cov(1:find(strcmp(cov.Chr,'chrY'),1,'last'),:); % up till the last bin in Chr. Y
		binCov=nan(length(hg38_10kb_bins),1);
		binCov(hg38_10kb_bins.NonNBase>0)=cov.count;
		SampleNames{i}=Sample;
		SampleCov=[SampleCov,binCov];
		fprintf(1,'done\n');
	end
	save ReadCount_10kb.mat SampleCov SampleInfo
else
	load ReadCount_10kb.mat;
end

%% GC normalization of read depth
binSize=10000;
Bins_10kb=hg38_10kb_bins(1:find(strcmp(hg38_10kb_bins.Chr,'chrX'),1,'last'),:);
load GCbins.mat
GCbinSize=50000;
GCbins=GCbins_50kb;
GC_strata=100*(floor(min(GCbins.avgGC)/binSize*100):0.5:ceil(max(GCbins.avgGC)/binSize*100));
[~,GC_bin_id]=histc(GCbins.avgGC,GC_strata);

normalizedCounts=double(SampleCov(1:length(Bins_10kb),:));
%normalizedCounts(GCbins.avgGC<3150 | GCbins.avgGC>6300)=nan; less than 50 bins
%normalizedCounts(GCbins.NonNBase<0.9*binSize)=nan;
%normalizedCounts=normalizedCounts/(mean(normalizedCounts(~isnan(normalizedCounts))));

lognormalizedCounts=log2(normalizedCounts);
lognormalizedCounts(Bins_10kb.NonNBase<0.9*binSize)=nan;
lognormalizedCounts=lognormalizedCounts-repmat(median(lognormalizedCounts,1,'omitnan'),length(lognormalizedCounts),1);
median_gc_strength=nan(length(GC_strata)-1,size(normalizedCounts,2));	
for GC_id=1:length(GC_strata)-1
	binID=find(GC_bin_id==GC_id);
	counts=lognormalizedCounts(binID,:);
	no_gc_bins=sum(~isnan(counts),1);
	median_gc_strength(GC_id,:)=median(counts,1,'omitnan');
	median_gc_strength(GC_id,no_gc_bins<100)=nan;
end
GCcorr_logCov=lognormalizedCounts;
GCcorr_logCov(GC_bin_id>0,:)=GCcorr_logCov(GC_bin_id>0,:)-median_gc_strength(GC_bin_id(GC_bin_id>0),:);
GCcorr_Cov=2.^GCcorr_logCov;

%% Arm-level coverage
SampleCovCorr=GCcorr_Cov;
for Armi=1:length(ChrArms)
	chrlabel=ChrArms.chr(Armi);	
	bin_id=find(strcmp(Bins_10kb.Chr,chrlabel) & Bins_10kb.Start>=ChrArms.left(Armi) & Bins_10kb.Start<ChrArms.right(Armi));	
	SampleArmMedianCov(Armi,:)=median(SampleCovCorr(bin_id,:),1,'omitnan');
	ArmCN=SampleArmMedianCov(Armi,:);
	if strcmp(ChrArms.label{Armi},'chr10q') % 10q with clonal gain
		meanCN=median(ArmCN(ArmCN>1.25),'omitnan');
	else
		meanCN=median(ArmCN);
	end
	stdCN=std(log(ArmCN(ArmCN>=0.75*meanCN & ArmCN<=1/0.75*meanCN)));
	ArmMedianCov(Armi,:)=dataset(ChrArmsLabels(Armi),meanCN,stdCN,'VarNames',{'Arm','medianCov','std'});	
end
DiploidSampleIdx=find(sum(abs(SampleArmMedianCov-repmat(ArmMedianCov.medianCov,1,length(SampleInfo.SampleNames)))<=0.15,1)>=35); % at most three arms with CNA
save ArmLevelCov.mat ArmMedianCov SampleArmMedianCov DiploidSampleIdx

%% Normalize recurrent bias
load(recurrent_bias_coverage_file);
normCov=SampleCovCorr./repmat(MedianCov_10kb,1,length(SampleInfo.SampleNames));
for Armi=1:length(ChrArms)
	chrlabel=ChrArms.chr{Armi};	
	bin_id=find(strcmp(Bins_10kb.Chr,chrlabel) & Bins_10kb.Start>=ChrArms.left(Armi) & Bins_10kb.Start<ChrArms.right(Armi));	
	SampleArmMedianCov(Armi,:)=median(normCov(bin_id,:),1,'omitnan');
end
normalizedCov=normCov./repmat(median(SampleArmMedianCov,1,'omitnan'),size(normCov,1),1);
save SI_normalizedCov.mat normalizedCov

%% Moving average
SampleCov=normalizedCov;
SampleCov_50kb=nan(size(SampleCov));
SampleCov_250kb=nan(size(SampleCov));
SampleCov_1Mb=nan(size(SampleCov));
for Armi=1:length(ChrArms)
	chrlabel=ChrArms.chr{Armi};	
	bins=find(strcmp(Bins_10kb.Chr,chrlabel));
	for bi=1:length(bins)
		neighborhood_1=max(1,bi-2):1:min(length(bins),bi+2);
		neighborhood_2=max(1,bi-12):1:min(length(bins),bi+12);
		neighborhood_3=max(1,bi-50):1:min(length(bins),bi+50);
		SampleCov_50kb(bins(bi),:)=mean(SampleCov(bins(neighborhood_1),:),1,'omitnan');
		SampleCov_250kb(bins(bi),:)=mean(SampleCov(bins(neighborhood_2),:),1,'omitnan');
		SampleCov_1Mb(bins(bi),:)=mean(SampleCov(bins(neighborhood_3),:),1,'omitnan');
	end
end
save SI_binnedCov.mat SampleCov_*b

