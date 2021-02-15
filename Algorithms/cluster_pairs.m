function CLUSTERS=cluster_pairs(curr_disc_reads,varargin)

Optional_Args={'JumpSize','ExactSearch','Step','minCount','Threshold','InterPairs'};
Options=parse_options(Optional_Args,varargin);

if ~isempty(Options.JumpSize) && isnumeric(Options.JumpSize) && Options.JumpSize>0
    insertSize=Options.JumpSize;
else
    insertSize=300;
end

if ~isempty(Options.ExactSearch) && strcmp(Options.ExactSearch,'False') ...
        && ~isempty(Options.Step) && isnumeric(Options.Step) && (Options.Step>1)
    exact_search=0;
    step=Options.Step;
else
    exact_search=1;
    step=1;
end

if ~isempty(Options.minCount) && isnumeric(Options.minCount) && Options.minCount>1
    minCount=Options.minCount;
else
    minCount=1;
end

if ~isempty(Options.Threshold) && isnumeric(Options.Threshold) && Options.Threshold>1
    threshold=Options.Threshold;
else
    threshold=1;
end

if ~isempty(Options.InterPairs) && strcmp(Options.InterPairs,'True')
    InterPairs=1;
else
    InterPairs=0;
end

left_chr_list=unique(curr_disc_reads.chr);

clusters.leftChr=[];
clusters.rightChr=[];
clusters.leftRange=[];
clusters.rightRange=[];
clusters.counts=[];

for i=1:1:length(left_chr_list)
    % Process inter-chromosomal pairs for each chromosome
    if InterPairs==1
        right_chr_list=unique(curr_disc_reads.chr2(curr_disc_reads.chr==left_chr_list(i)));
    else
        right_chr_list=unique(curr_disc_reads.chr(curr_disc_reads.chr==left_chr_list(i)));
    end	
    for j=1:1:length(right_chr_list)
        % If process all inter-chromosomal pairs, only process chrA vs. chrB once (chrA<=chrB)
        if (length(left_chr_list)==1) | (InterPairs==0 & right_chr_list(j)==left_chr_list(i)) | (right_chr_list(j)<0) | (InterPairs==1 & right_chr_list(j)>left_chr_list(i))
            %fprintf('Chr.%d and Chr.%d:',left_chr_list(i),right_chr_list(j));
			if InterPairs==0
            	ind=(curr_disc_reads.chr==left_chr_list(i));
			else
				ind=(curr_disc_reads.chr==left_chr_list(i)) & (curr_disc_reads.chr2==right_chr_list(j));
			end 
			left_pos=double(curr_disc_reads.pos(ind));
            right_pos=double(curr_disc_reads.pos2(ind));
            % Search for clusters based on left-read pos
            left_pos_clusters=cluster_segments(left_pos,'JumpSize',insertSize, 'MinCount',minCount);
            counts=[]; 
            if ~isempty(left_pos_clusters.left)                
                subclusters.left=[];
                subclusters.right=[];
                % Search for discordant pairs that indeed form a cluster both on left-read pos and on right-read pos
                for k=1:length(left_pos_clusters.left)
                    % Pick out all reads that fall in the cluster based on left-read pos
                    ind2=(left_pos>=left_pos_clusters.left(k) & left_pos<=left_pos_clusters.right(k));
                    left_cluster_pos=left_pos(ind2);                    
                    right_cluster_pos=right_pos(ind2);
                    % Pick out discordant read pairs chrLeft:pos1--chrRight:pos2 such that there are at least 1+minCount discordant pairs in the neighborhood
                    % chrLeft:[pos1-threshold*insertSize, pos1+threshold*insertSize] and chrRight:[pos2-threshold*insertSize, pos2+threshold*insertSize]
                    for m=1:length(left_cluster_pos)
                        if sum(logical(abs(right_cluster_pos-right_cluster_pos(m))<=insertSize*threshold & abs(left_cluster_pos-left_cluster_pos(m))<=insertSize*threshold))>=1+minCount
                            subclusters.left=[subclusters.left;left_cluster_pos(m)];
                            subclusters.right=[subclusters.right;right_cluster_pos(m)];
                        end
                    end
                end	
                % Merge subclusters to generate non-overlapping clusters of discordant pairs
				if ~isempty(subclusters.left)
					LeftRange=[];
                    RightRange=[];
					[subclusters.left,sort_ind]=sort(subclusters.left);
                    subclusters.right=subclusters.right(sort_ind);
                    % ind4 labels each discordant pair that has not been assigned a cluster
					ind4=logical(subclusters.left>0); 
                    while sum(ind4)>0
                        % Pick the first subcluster/unassigned discordant pair
						ind5=find(ind4>0,1,'first');
                        Left=subclusters.left(ind5);
                        Right=subclusters.right(ind5);
                        % Merge every subcluster that overlaps with the first subcluster
						ind5=find(abs(subclusters.left-Left)<=2*insertSize*threshold & abs(subclusters.right-Right)<=2*insertSize*threshold);
                        ind4(ind5)=0;
                        LeftRange=[LeftRange;min(subclusters.left(ind5)),max(subclusters.left(ind5))];
                        RightRange=[RightRange;min(subclusters.right(ind5)),max(subclusters.right(ind5))];
                    end
                    [~,ind5]=sort(LeftRange(:,2)-LeftRange(:,1),'descend');
                    LeftRange=LeftRange(ind5,:);
                    RightRange=RightRange(ind5,:); 
                    % Merge clusters that are likely to correspond to one cluster
                    l=1;
                    while l<=size(LeftRange,1) 
                        ind6=logical(((LeftRange(:,1)<=LeftRange(l,2)+insertSize & LeftRange(:,1)>=LeftRange(l,1)-insertSize) ...
                                 | (LeftRange(:,2)<=LeftRange(l,2)+insertSize & LeftRange(:,2)>=LeftRange(l,1)-insertSize)) ...
                                 & ((RightRange(:,1)<=RightRange(l,2)+insertSize & RightRange(:,1)>=RightRange(l,1)-insertSize) ...
                                 | (RightRange(:,2)<=RightRange(l,2)+insertSize & RightRange(:,2)>=RightRange(l,1)-insertSize)));
                        if sum(ind6)>1
                            LeftRange(l,:)=[min(LeftRange(ind6,1)),max(LeftRange(ind6,2))];
                            RightRange(l,:)=[min(RightRange(ind6,1)),max(RightRange(ind6,2))];
                            ind6(l)=0;
                            LeftRange=LeftRange(logical(1-ind6),:);
                            RightRange=RightRange(logical(1-ind6),:);
                        else
                            counts(l)=length(find(left_pos>=LeftRange(l,1) & left_pos<=LeftRange(l,2) & right_pos>=RightRange(l,1) & right_pos<=RightRange(l,2)));
                            l=l+1;
                        end
                    end 
					[~,ind]=sort(LeftRange(:,1));
                    LeftRange=LeftRange(ind,:);
                    RightRange=RightRange(ind,:);
                    counts=counts(ind);
					index=logical(counts'>minCount);
                    fprintf('  \t %d clusters.\t',sum(index));
                    clusters.leftRange=[clusters.leftRange; LeftRange(index,:)];
                    clusters.rightRange=[clusters.rightRange; RightRange(index,:)];
                    clusters.counts=[clusters.counts;reshape(counts(index),sum(index),1)];
                    clusters.leftChr=[clusters.leftChr;repmat(left_chr_list(i),sum(index),1)];
                    clusters.rightChr=[clusters.rightChr;repmat(right_chr_list(j),sum(index),1)];
                end
            end
            if isempty(counts)
                fprintf('  \t no clusters.\t');
            end
        end
    end
end

CLUSTERS=[];
if length(clusters.leftChr)==0
	return;
end

str=zeros(size(clusters.leftChr));
str2=str;
counts=str;
for ci=1:length(clusters.counts)    
	if InterPairs
		ind=find(curr_disc_reads.chr==clusters.leftChr(ci) & (curr_disc_reads.chr2==clusters.rightChr(ci)) ...	
			& curr_disc_reads.pos>=clusters.leftRange(ci,1) & curr_disc_reads.pos<=clusters.leftRange(ci,2) ...			            
			& curr_disc_reads.pos2>=clusters.rightRange(ci,1) & curr_disc_reads.pos2<=clusters.rightRange(ci,2));
	else
		ind=find(curr_disc_reads.chr==clusters.leftChr(ci) ...
			& curr_disc_reads.pos>=clusters.leftRange(ci,1) & curr_disc_reads.pos<=clusters.leftRange(ci,2) ...							           
			& curr_disc_reads.pos2>=clusters.rightRange(ci,1) & curr_disc_reads.pos2<=clusters.rightRange(ci,2));								    
	end										   
	STR1=curr_disc_reads.str(ind);
	STR2=curr_disc_reads.str2(ind);
	counts(ci)=length(STR1);
	str(ci)=median(STR1);
	str2(ci)=median(STR2);	
end

CLUSTERS=dataset(counts,clusters.leftChr,clusters.leftRange(:,1),clusters.leftRange(:,2),str,clusters.rightChr,clusters.rightRange(:,1),clusters.rightRange(:,2),str2,...
            'VarNames',{'count','chr','pos_left','pos_right','str1','chr2','pos2_left','pos2_right','str2'});
CLUSTERS=CLUSTERS(CLUSTERS.str1.*CLUSTERS.str2~=0 & counts>minCount,:);

