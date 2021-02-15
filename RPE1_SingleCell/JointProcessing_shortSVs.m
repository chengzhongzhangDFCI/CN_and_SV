%% Processing
load SampleInfo.mat
for i=1:length(SampleInfo)
	SampleInfo.SampleNames{i}=regexprep(SampleInfo.SampleNames{i},'-','_');
	SampleInfo.SampleGroupID{i}=regexprep(SampleInfo.SampleGroupID{i},'-','_');
end

SampleGroups=unique(SampleInfo.SampleGroupID);

insertSize=200;
readLength=150;

%%
TD_BKPS=[];TD_SUPP_READS=[];
DL_BKPS=[];DL_SUPP_READS=[];
FB_BKPS=[];FB_SUPP_READS=[];

%%
load hg38.mat
CHRs=ChrLabels(~strcmp(ChrLabels,'chrY'));
%CHRs={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'};
for i=1:length(CHRs)
	ci=CHRs{i};
    fprintf('\n Short intra %s:\n',ci);
	supp_reads=SHORT_NONINV(strcmp(SHORT_NONINV.chr1,ci),:);
	temp=[dataset((1:length(supp_reads))','VarNames','idx'),supp_reads];	
	temp.pos1=min(supp_reads.pos1,supp_reads.pos2);
	temp.pos2=max(supp_reads.pos1,supp_reads.pos2);
	temp.str1=supp_reads.str1.*int8(supp_reads.pos1<=supp_reads.pos2)+supp_reads.str2.*int8(supp_reads.pos1>supp_reads.pos2);
	temp.str2=supp_reads.str1.*int8(supp_reads.pos1>supp_reads.pos2)+supp_reads.str2.*int8(supp_reads.pos1<=supp_reads.pos2);
	del_reads=temp(temp.str1==1,:); 
	td_reads=temp(temp.str1==-1,:);
	fprintf(1,'  Deletion type:');
	tic,
	del_cluster_reads=del_reads(del_reads.split | (del_reads.pos2-del_reads.pos1)>=2*insertSize,:);
	curr_cluster=del_cluster_reads(:,{'chr1','pos1','str1','chr2','pos2','str2'});
	curr_cluster.pos1=double(curr_cluster.pos1);curr_cluster.pos2=double(curr_cluster.pos2);
	curr_cluster.Properties.VarNames(1:3)={'chr','pos','str'};
	curr_cluster.chr=zeros(size(curr_cluster.chr));curr_cluster.chr2=zeros(size(curr_cluster.chr));
	curr_cluster=sortrows(curr_cluster,{'pos'});
	disc_spacing=floor(double(max(curr_cluster.pos)-min(curr_cluster.pos))/length(curr_cluster));	
	minCount=max(round(11/log10(disc_spacing)),2)-1;
	fprintf(1,'%d discordant reads over %d bp, average spacing %d bp, minimum threshold %d\n',length(curr_cluster),max(curr_cluster.pos)-min(curr_cluster.pos),disc_spacing,minCount);
	clusters=cluster_pairs(curr_cluster,'JumpSize',readLength,'Threshold',1,'minCount',minCount,'IntraPairs');
	fprintf(1,'\nTotal: %d clusters\n',length(clusters));
	if length(clusters)>0
		clusters.chr=repmat(ci,size(clusters.chr));clusters.chr2=clusters.chr;
		clusters=sortrows(clusters,{'pos_left','pos2_left'});
		clusters.Properties.VarNames={'count'    'chr1'    'pos1_left'    'pos1_right'    'str1'    'chr2'    'pos2_left'    'pos2_right'    'str2'};	
		fprintf(1,'Collecting support from %d reads',length(supp_reads));
		[bkps,temp_supp]=collectSampleSupportShortSVs(clusters,temp,SampleInfo,insertSize,readLength);
		supp=[temp_supp(:,'SVidx'),supp_reads(temp_supp.idx,:)];	
		supp.SVidx=supp.SVidx+length(DL_BKPS);
		DL_BKPS=[DL_BKPS;bkps];
		DL_SUPP_READS=[DL_SUPP_READS;supp];	
	end
	fprintf(1,'..done in %.2f sec\n',toc);
	fprintf(1,'  Duplication type:');
	tic,
	td_cluster_reads=td_reads(td_reads.pos2-td_reads.pos1>=readLength,:);
	curr_cluster=td_cluster_reads(:,{'chr1','pos1','str1','chr2','pos2','str2'});
	curr_cluster.pos1=double(curr_cluster.pos1);curr_cluster.pos2=double(curr_cluster.pos2);
	curr_cluster.Properties.VarNames(1:3)={'chr','pos','str'};
	curr_cluster.chr=zeros(size(curr_cluster.chr));curr_cluster.chr2=zeros(size(curr_cluster.chr));
	curr_cluster=sortrows(curr_cluster,{'pos'});
	disc_spacing=floor(double(max(curr_cluster.pos)-min(curr_cluster.pos))/length(curr_cluster));	
	minCount=max(round(11/log10(disc_spacing)),2)-1;
	fprintf(1,'%d discordant reads over %d bp, average spacing %d bp, minimum threshold %d:',length(curr_cluster),max(curr_cluster.pos)-min(curr_cluster.pos),disc_spacing,minCount);
	clusters=cluster_pairs(curr_cluster,'JumpSize',readLength,'Threshold',1,'minCount',minCount,'IntraPairs');
	fprintf(1,'\nTotal: %d clusters\n',length(clusters));
	if length(clusters)>0
		clusters.chr=repmat(ci,size(clusters.chr));clusters.chr2=clusters.chr;
		clusters=sortrows(clusters,{'pos_left','pos2_left'});
		clusters.Properties.VarNames={'count'    'chr1'    'pos1_left'    'pos1_right'    'str1'    'chr2'    'pos2_left'    'pos2_right'    'str2'};	
		fprintf(1,'Collecting support from %d reads',length(supp_reads));
		[bkps,temp_supp]=collectSampleSupportShortSVs(clusters,temp,SampleInfo,insertSize,readLength);
		supp=[temp_supp(:,'SVidx'),supp_reads(temp_supp.idx,:)];	
		supp.SVidx=supp.SVidx+length(TD_BKPS);
		TD_BKPS=[TD_BKPS;bkps];
		TD_SUPP_READS=[TD_SUPP_READS;supp];	
	end
	fprintf(1,'..done in %.2f sec\n',toc);

	fprintf(1,'  Inverted type:');
	supp_reads=SHORT_INV(strcmp(SHORT_INV.chr1,ci),:);	
	temp=[dataset((1:length(supp_reads))','VarNames','idx'),supp_reads];	
	temp.pos1=min(supp_reads.pos1,supp_reads.pos2);
	temp.pos2=max(supp_reads.pos1,supp_reads.pos2);
	temp.str1=supp_reads.str1.*int8(supp_reads.pos1<=supp_reads.pos2)+supp_reads.str2.*int8(supp_reads.pos1>supp_reads.pos2);
	temp.str2=supp_reads.str1.*int8(supp_reads.pos1>supp_reads.pos2)+supp_reads.str2.*int8(supp_reads.pos1<=supp_reads.pos2);
	FF_reads=temp(temp.str1==1,:); 
	RR_reads=temp(temp.str1==-1,:);
	fprintf(1,'  	FF:');
	tic,
	FF_cluster_reads=FF_reads(FF_reads.maq1>=30 & FF_reads.maq2>=30 & ~FF_reads.split,:);
	curr_cluster=FF_cluster_reads(:,{'chr1','pos1','str1','chr2','pos2','str2'});
	curr_cluster.pos1=double(curr_cluster.pos1);curr_cluster.pos2=double(curr_cluster.pos2);
	curr_cluster.Properties.VarNames(1:3)={'chr','pos','str'};
	curr_cluster.chr=zeros(size(curr_cluster.chr));curr_cluster.chr2=zeros(size(curr_cluster.chr));
	curr_cluster=sortrows(curr_cluster,{'pos'});
	disc_spacing=floor(double(max(curr_cluster.pos)-min(curr_cluster.pos))/length(curr_cluster));	
	minCount=max(round(11/log10(disc_spacing)),2)-1;
	fprintf(1,'%d discordant reads over %d bp, average spacing %d bp, minimum threshold %d:',length(curr_cluster),max(curr_cluster.pos)-min(curr_cluster.pos),disc_spacing,minCount);
	clusters=cluster_pairs(curr_cluster,'JumpSize',readLength,'Threshold',1,'minCount',minCount,'IntraPairs');	
	fprintf(1,'    RR:');
	RR_cluster_reads=RR_reads(RR_reads.maq1>=30 & RR_reads.maq2>=30 & ~RR_reads.split,:);
	curr_cluster=RR_cluster_reads(:,{'chr1','pos1','str1','chr2','pos2','str2'});
	curr_cluster.pos1=double(curr_cluster.pos1);curr_cluster.pos2=double(curr_cluster.pos2);
	curr_cluster.Properties.VarNames(1:3)={'chr','pos','str'};
	curr_cluster.chr=zeros(size(curr_cluster.chr));curr_cluster.chr2=zeros(size(curr_cluster.chr));
	curr_cluster=sortrows(curr_cluster,{'pos'});
	disc_spacing=floor(double(max(curr_cluster.pos)-min(curr_cluster.pos))/length(curr_cluster));	
	minCount=max(round(11/log10(disc_spacing)),2)-1;
	fprintf(1,'%d discordant reads over %d bp, average spacing %d bp, minimum threshold %d:',length(curr_cluster),max(curr_cluster.pos)-min(curr_cluster.pos),disc_spacing,minCount);
	clusters=[clusters;cluster_pairs(curr_cluster,'JumpSize',readLength,'Threshold',1,'minCount',minCount,'IntraPairs')];	
	if length(clusters)>0
		clusters.chr=repmat(ci,size(clusters.chr));clusters.chr2=clusters.chr;
		clusters=sortrows(clusters,{'pos_left','pos2_left'});
		clusters.Properties.VarNames={'count'    'chr1'    'pos1_left'    'pos1_right'    'str1'    'chr2'    'pos2_left'    'pos2_right'    'str2'};	
		clusters=clusters(clusters.pos1_right-clusters.pos1_left>=readLength & clusters.pos2_right-clusters.pos2_left>=readLength,:);
		fprintf(1,'\nTotal: %d clusters\n',length(clusters));
		fprintf(1,'Collecting support from %d reads',length(supp_reads));
		[bkps,temp_supp]=collectSampleSupportShortSVs(clusters,temp,SampleInfo,insertSize,readLength);
		supp=[temp_supp(:,'SVidx'),supp_reads(temp_supp.idx,:)];	
		supp.SVidx=supp.SVidx+length(DL_BKPS);
		FB_BKPS=[FB_BKPS;bkps];
		FB_SUPP_READS=[FB_SUPP_READS;supp];	
	end
	fprintf(1,'..done in %.2f sec\n',toc);
	save SHORT_BKPS.mat FB_BKPS TD_BKPS DL_BKPS FB_SUPP_READS TD_SUPP_READS DL_SUPP_READS;
end
FB_BKPS=[dataset((1:1:length(FB_BKPS))','VarNames',{'SVidx'}),FB_BKPS];
TD_BKPS=[dataset((1:1:length(TD_BKPS))','VarNames',{'SVidx'}),TD_BKPS];
DL_BKPS=[dataset((1:1:length(DL_BKPS))','VarNames',{'SVidx'}),DL_BKPS];
save SHORT_BKPS.mat FB_BKPS TD_BKPS DL_BKPS FB_SUPP_READS TD_SUPP_READS DL_SUPP_READS
fprintf(1,'..intra-chromosomal shortSVs complete.\n');
