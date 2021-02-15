%% Allelic depth processing
AD_folder='AllelicCov';
refCounts=[];
altCounts=[];
if ~exist('hetSNPAD.mat')
	for i=1:length(SampleInfo.SampleNames)
		Sample=SampleInfo.SampleNames{i};
		AD_file=fullfile(AD_folder,[Sample '.hetSNP.AD.txt']);
		fprintf(1,'Processing allelic depth data from %s..',AD_file);
		fid=fopen(AD_file,'r');	
		temp=textscan(fid,'%*s\t%*d\t%*s\t%*c\t%*c\t%d\t%d\t%*[^\n]\n','HeaderLines',1);
		AD=dataset(temp{1},temp{2},'VarNames',{'ref','alt'});
		refCounts(:,i)=AD.ref;
		altCounts(:,i)=AD.alt;
		fclose(fid);
		fprintf(1,'done\n');
	end
	refCounts=uint16(refCounts);altCounts=uint16(altCounts);
	save hetSNPAD.mat refCounts altCounts
else
	load hetSNPAD.mat
end
hetCov.SampleNames=SampleInfo.SampleNames;
hetCov.RefCounts=refCounts;
hetCov.AltCounts=altCounts;

%% Resorting allelic coverage by Chr name
load HetSNPPos.mat
[~,sortIdx]=sortrows(HetSNP,{'chr','pos'});
hetCov.RefCounts=refCounts(sortIdx,:);
hetCov.AltCounts=altCounts(sortIdx,:);

load Haplotype.mat
%% Filtering SNP sites -- optional
if exist('Haplotype_filtered.mat')
	load Haplotype_filtered.mat;
else
	%gain_seg_idx=logical(strcmp(Haplotype_pos.chr,'chr10') & Haplotype_pos.pos>=60784157);
	%monosomyFilter=logical(~strcmp(Haplotype_pos.AltBase,'*') & Haplotype_phase.disomic_covered>=50 & (abs(Haplotype_phase.disomic_allelicRatio-0.5)<=0.2 | (gain_seg_idx & abs(min(Haplotype_phase.disomic_allelicRatio,1-Haplotype_phase.disomic_allelicRatio)-0.33)<=0.15)));
	%linkedReadsFilter=logical(abs(Haplotype_phase.linkedReads_frac)>=0.25 & Haplotype_phase.linkedReadsCount>=50 & Haplotype_phase.linkedReads_pval<0.01);
%phased_idx=logical((Haplotype_phase.allele_linkedReads~=0 & linkedReadsFilter) | (Haplotype_phase.allele_monosomies~=0 & monosomyFilter));
	Haplotype1=Haplotype_phase.allele_linkedReads;
	Haplotype2=Haplotype_phase.allele_monosomies;
	Haplotype1(~linkageFilter)=0;
	Haplotype2(~singleCellAlleleFilter)=0;
	Haplotype=double(sign(Haplotype1+Haplotype2));
	Haplotype=Haplotype(sortIdx);
	save Haplotype_filtered.mat Haplotype
end

%% Calculate allelic coverage using allelic coverage at filtered SNPs
Ref=double(hetCov.RefCounts);
Alt=double(hetCov.AltCounts);
AlleleA=Ref.*repmat(Haplotype>0,1,length(hetCov.SampleNames))+Alt.*repmat(Haplotype<0,1,length(hetCov.SampleNames));
AlleleB=Ref.*repmat(Haplotype<0,1,length(hetCov.SampleNames))+Alt.*repmat(Haplotype>0,1,length(hetCov.SampleNames));
loci=Haplotype_pos(:,1:2);
AlleleRatio=AlleleA./(AlleleA+AlleleB);

% allelic coverage at 1kb (to avoid over representation at polymorphic regions)
% allelic coverage is best used to estimate allelic ratio (and performing average in the allelic ratio space)

%% 1kb
load SNP_binning_1kb.mat
[bins,IA,IC]=unique(Bin_id);
AlleleRatio_1kb=nan(length(Bins_1kb),length(hetCov.SampleNames));
for bi=1:length(bins)
	id=Bin_SNP_idx{bins(bi)};
	AlleleRatio_1kb(bins(bi),:)=mean(AlleleRatio(id,:),1,'omitnan');
end

% Allele ratio at 10 kb (average over 10 consecutive 1kb bins)
load SNP_binning_10kb.mat
AlleleRatio_10kb=nan(length(Bins_10kb),length(hetCov.SampleNames));
chrs=unique(loci.chr,'stable');
for i=1:length(chrs)
	chrlabel=chrs{i};	
	bins=find(strcmp(Bins_1kb.Chr,chrlabel));
	bins2=find(strcmp(Bins_10kb.Chr,chrlabel));
	ratio=nan(10*length(bins2),length(hetCov.SampleNames));
	ratio(1:length(bins),:)=AlleleRatio_1kb(bins,:);	
	S=zeros(10,length(bins2),length(hetCov.SampleNames));
	for d=1:10
		S(d,:,:)=ratio(d:10:end,:);
	end
	AlleleRatio_10kb(bins2,:)=squeeze(mean(S,1,'omitnan'));	
end

[bins,IA,IC]=unique(Bin_id);
Bin_SNP_count=zeros(size(AlleleRatio_10kb));
for bi=1:length(bins)
	Bin_SNP_count(bins(bi),:)=sum(~isnan(AlleleRatio(Bin_SNP_idx{bins(bi)},:)),1);	
end
save SI_allelicRatio.mat Bins_10kb AlleleRatio_10kb Bin_SNP_count

%% Moving window average
AlleleARatio_50kb=AlleleRatio_10kb; 
AlleleARatio_250kb=AlleleRatio_10kb; 
AlleleARatio_1Mb=AlleleRatio_10kb;
SNPCount=Bin_SNP_count;
SNPCount_50kb=SNPCount;
SNPCount_250kb=SNPCount;
SNPCount_1Mb=SNPCount;

chrs=unique(loci.chr,'stable');
for i=1:length(chrs)
	chrlabel=chrs{i};	
	bins=find(strcmp(Bins_10kb.Chr,chrlabel));
	for bi=1:length(bins)
		neighborhood_1=max(1,bi-2):1:min(length(bins),bi+2);
		neighborhood_2=max(1,bi-12):1:min(length(bins),bi+12);
		neighborhood_3=max(1,bi-50):1:min(length(bins),bi+50);
		AlleleARatio_50kb(bins(bi),:)=mean(AlleleRatio_10kb(bins(neighborhood_1),:),1,'omitnan');
		AlleleARatio_250kb(bins(bi),:)=mean(AlleleRatio_10kb(bins(neighborhood_2),:),1,'omitnan');
		AlleleARatio_1Mb(bins(bi),:)=mean(AlleleRatio_10kb(bins(neighborhood_3),:),1,'omitnan');
		SNPCount_50kb(bins(bi),:)=sum(SNPCount(bins(neighborhood_1),:));	
		SNPCount_250kb(bins(bi),:)=sum(SNPCount(bins(neighborhood_2),:));	
		SNPCount_1Mb(bins(bi),:)=sum(SNPCount(bins(neighborhood_3),:));
	end	
end
SNPCount_50kb=uint16(SNPCount_50kb);
SNPCount_250kb=uint16(SNPCount_250kb);
SNPCount_1Mb=uint16(SNPCount_1Mb);
save SI_BinnedAllelicRatio.mat AlleleARatio_*b* SNPCount_*b

