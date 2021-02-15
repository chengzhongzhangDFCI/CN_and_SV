%% INPUTS
%insertSize=300;
%readLength=150;
%shortRange=2000;
%minCount=2;
%DISCORDANTS;

%% Processing
load SampleInfo.mat
for i=1:length(SampleInfo)
	SampleInfo.SampleNames{i}=regexprep(SampleInfo.SampleNames{i},'-','_');
	SampleInfo.SampleGroupID{i}=regexprep(SampleInfo.SampleGroupID{i},'-','_');
end

SampleGroups=unique(SampleInfo.SampleGroupID);

%% Intra-chromosomal clusters
IntraSupp=DISCORDANTS;

%%
load hg38.mat
CHRs=ChrLabels(~strcmp(ChrLabels,'chrY'));
%CHRs={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'};
INTRA_BKPS=[];
INTRA_SUPP_READS=[];
SUPP_READS_GP=1;
for i=1:length(CHRs)
	ci=CHRs{i};	
	fprintf('\nIntra %s:',ci);
	supp_reads=sortrows(IntraSupp(strcmp(IntraSupp.chr1,ci),:),{'pos1','pos2'});	
	IntraSupp=IntraSupp(~strcmp(IntraSupp.chr1,ci),:);
	temp=[dataset((1:length(supp_reads))','VarNames','idx'),supp_reads];
	temp.pos1=min(supp_reads.pos1,supp_reads.pos2);
	temp.pos2=max(supp_reads.pos1,supp_reads.pos2);
	temp.str1=supp_reads.str1.*int8(supp_reads.pos1<=supp_reads.pos2)+supp_reads.str2.*int8(supp_reads.pos1>supp_reads.pos2);
	temp.str2=supp_reads.str1.*int8(supp_reads.pos1>supp_reads.pos2)+supp_reads.str2.*int8(supp_reads.pos1<=supp_reads.pos2);
	
	% Clustering
    tic,	
	cluster_reads=temp(temp.maq1>=30 & temp.maq2>=30 & temp.SampleGroup>0,:);
	if length(cluster_reads)<=2
		continue;
	end
	clusters=[];
	for j=1:length(SampleGroups)
		curr_sample=SampleGroups{j};
		fprintf(1,'\n  SampleGroup %d: %s\n',j,curr_sample);	
		curr_cluster=cluster_reads(strcmp(SampleInfo.SampleGroupID(cluster_reads.SampleGroup),curr_sample),{'chr1','pos1','str1','chr2','pos2','str2'});
		if length(curr_cluster)==0
			continue;
		end
		curr_cluster.pos1=double(curr_cluster.pos1);curr_cluster.pos2=double(curr_cluster.pos2);
		curr_cluster.Properties.VarNames(1:3)={'chr','pos','str'};
		curr_cluster.chr=zeros(size(curr_cluster.chr));curr_cluster.chr2=zeros(size(curr_cluster.chr));
		curr_cluster=sortrows(curr_cluster,{'pos','pos2'});	
		span=double(max(curr_cluster.pos)-min(curr_cluster.pos));		
		fprintf(1,'%d discordant reads over %d bp, average density %.1f reads per Mb^2',...
		length(curr_cluster),span,length(curr_cluster)/(span/1e6)^2);
		clusters=[clusters;cluster_pairs(curr_cluster,'JumpSize',insertSize,'Threshold',1,'minCount',minCount,'IntraPairs')];
	end
	fprintf(1,'\nTotal: %d clusters\n',length(clusters));
	if length(clusters)==0
		continue;
	else
		clusters=clusters(clusters.count<=100000,:);
		if length(clusters)==0
			continue;
		end
	end	
	clusters.chr=repmat(ci,size(clusters.chr));clusters.chr2=clusters.chr;
	clusters=sortrows(clusters,{'pos_left','pos2_left'});
	clusters.Properties.VarNames={'count'    'chr1'    'pos1_left'    'pos1_right'    'str1'    'chr2'    'pos2_left'    'pos2_right'    'str2'};	
	fprintf(1,'Collecting support from %d reads:',length(supp_reads));	
	[bkps,temp_supp]=collectSampleSupport(clusters,temp,SampleInfo,insertSize,readLength);
	fprintf(1,'Total breakpoint pairs: %d\n',length(bkps));
	if length(bkps)<1
		continue;
	end
	supp=[temp_supp(:,'SVidx'),supp_reads(temp_supp.idx,:)];	
	supp.SVidx=supp.SVidx+length(INTRA_BKPS);
	INTRA_BKPS=[INTRA_BKPS;bkps];
	INTRA_SUPP_READS=[INTRA_SUPP_READS;supp];	
	if length(INTRA_SUPP_READS)<=4000000
		save INTRA_BKPS.mat INTRA_BKPS INTRA_SUPP_READS
	else
		save INTRA_BKPS.mat INTRA_BKPS;
		save(['INTRA_SUPP_READS_grp' int2str(SUPP_READS_GP) '.mat'],'INTRA_SUPP_READS');
		INTRA_SUPP_READS=[];
		SUPP_READS_GP=SUPP_READS_GP+1;
	end
	fprintf(1,'  Done in %.2f sec\n',toc);
end
INTRA_BKPS=[dataset((1:1:length(INTRA_BKPS))','VarNames',{'SVidx'}),INTRA_BKPS];
save INTRA_BKPS.mat INTRA_BKPS INTRA_SUPP_READS
fprintf(1,'..intra-chromosomal complete.\n');
