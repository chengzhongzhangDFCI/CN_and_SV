%% INPUTS
%insertSize=300;
%readLength=150;
%minCount=2;

%% Processing
load SampleInfo.mat
for i=1:length(SampleInfo)
	SampleInfo.SampleNames{i}=regexprep(SampleInfo.SampleNames{i},'-','_');
	SampleInfo.SampleGroupID{i}=regexprep(SampleInfo.SampleGroupID{i},'-','_');
end
SampleGroups=unique(SampleInfo.SampleGroupID);
%% Inter-chromosomal clusters
InterSupp=CHIMERAS;

%%
load hg38.mat
CHRs=ChrLabels(~strcmp(ChrLabels,'chrY'));
%CHRs={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'};
INTER_BKPS=[];
INTER_SUPP_READS=[];
SUPP_READS_GP=1;
for i=1:length(CHRs)
	ci=CHRs{i};	
	currInterSupp=InterSupp(strcmp(InterSupp.chr1,ci) | strcmp(InterSupp.chr2,ci),:);	
	InterSupp=InterSupp(~strcmp(InterSupp.chr1,ci) & ~strcmp(InterSupp.chr2,ci),:);
    chrjs=unique([currInterSupp.chr2(strcmp(currInterSupp.chr1,ci) & strncmp(currInterSupp.chr2,'chr',3));
				  currInterSupp.chr1(strcmp(currInterSupp.chr2,ci) & strncmp(currInterSupp.chr1,'chr',3))]);
	for j=1:length(chrjs)
		cj=chrjs{j};
		supp_reads=currInterSupp(strcmp(currInterSupp.chr2,cj) | strcmp(currInterSupp.chr1,cj),:);
		currInterSupp=currInterSupp(~strcmp(currInterSupp.chr2,cj) & ~strcmp(currInterSupp.chr1,cj),:);
        temp=[dataset((1:length(supp_reads))','VarNames','idx'),supp_reads];
		idx2=strcmp(supp_reads.chr2,ci);
		temp.chr1(idx2)=supp_reads.chr2(idx2);temp.chr2(idx2)=supp_reads.chr1(idx2);
		temp.pos1(idx2)=supp_reads.pos2(idx2);temp.pos2(idx2)=supp_reads.pos1(idx2);
		temp.str1(idx2)=supp_reads.str2(idx2);temp.str2(idx2)=supp_reads.str1(idx2);
		cluster_reads=temp(temp.maq1>=30 & temp.maq2>=30 & temp.SampleGroup>0,:);	
		if length(cluster_reads)<=2
    		continue;
		end	
		% Clustering	
		tic,	
		fprintf(1,'\n%s <=> %s:',ci,cj);
		clusters=[];
		for si=1:length(SampleGroups)
			curr_sample=SampleGroups{si};
			fprintf(1,'\n  SampleGroup %d: %s\n',si,curr_sample);	
			curr_cluster=cluster_reads(strcmp(SampleInfo.SampleGroupID(cluster_reads.SampleGroup),curr_sample),{'chr1','pos1','str1','chr2','pos2','str2'});
			if length(curr_cluster)==0
				continue;
			end
			curr_cluster.pos1=double(curr_cluster.pos1);curr_cluster.pos2=double(curr_cluster.pos2);
			curr_cluster.Properties.VarNames(1:3)={'chr','pos','str'};
			curr_cluster.chr=zeros(size(curr_cluster.chr));curr_cluster.chr2=zeros(size(curr_cluster.chr));
			curr_cluster=sortrows(curr_cluster,{'pos'});	
			span1=max(curr_cluster.pos)-min(curr_cluster.pos);
			span2=max(curr_cluster.pos2)-min(curr_cluster.pos2);
			fprintf(1,'%d discordant reads over %d bp (%s) x %d bp (%s), average density %.1f reads per Mb^2',...
			length(curr_cluster),span1,ci,span2,cj,...
			length(curr_cluster)/double(span1)/double(span2)*1e12);		
			clusters=[clusters;cluster_pairs(curr_cluster,'JumpSize',insertSize,'Threshold',1,'minCount',minCount,'InterPairs')];
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
		clusters.chr=repmat(ci,size(clusters.chr));clusters.chr2=repmat(cj,size(clusters.chr2));
		clusters=sortrows(clusters,{'pos_left','pos2_left'});
		clusters.Properties.VarNames={'count'    'chr1'    'pos1_left'    'pos1_right'    'str1'    'chr2'    'pos2_left'    'pos2_right'    'str2'};	
		fprintf(1,'Collecting support from %d reads:',length(supp_reads));	
		[bkps,temp_supp]=collectSampleSupport(clusters,temp,SampleInfo,insertSize,readLength);
		fprintf(1,'Total breakpoint pairs: %d\n',length(bkps));
		if length(bkps)<1
			continue;
		end
		supp=[temp_supp(:,'SVidx'),supp_reads(temp_supp.idx,:)];
		supp.SVidx=supp.SVidx+length(INTER_BKPS);
		INTER_BKPS=[INTER_BKPS;bkps];
		INTER_SUPP_READS=[INTER_SUPP_READS;supp];	
		if length(INTER_SUPP_READS)<=4000000
			save INTER_BKPS.mat INTER_BKPS INTER_SUPP_READS
		else
			save INTER_BKPS.mat INTER_BKPS;
			save(['INTER_SUPP_READS_grp' int2str(SUPP_READS_GP) '.mat'],'INTER_SUPP_READS');
			INTER_SUPP_READS=[];
			SUPP_READS_GP=SUPP_READS_GP+1;
		end
		fprintf(1,'  Done in %.2f sec\n',toc);	
	end
end
INTER_BKPS=[dataset((1:1:length(INTER_BKPS))','VarNames',{'SVidx'}),INTER_BKPS];
fprintf(1,'..inter-chromosomal events complete.\n');
fprintf(1,'Saving results..');
save INTER_BKPS.mat INTER_BKPS INTER_SUPP_READS
fprintf(1,'..inter-chromosomal complete.\n');
