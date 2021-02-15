INTER_BKPS.SVidx=INTER_BKPS.SVidx+length(INTRA_BKPS);
INTER_SUPP_READS.SVidx=INTER_SUPP_READS.SVidx+length(INTRA_BKPS);

% intra filtering
intra=INTRA_BKPS(INTRA_BKPS.TotalCount-INTRA_BKPS.TCount<=1 & INTRA_BKPS.TCount>=2 & INTRA_BKPS.NCount==0 & INTRA_BKPS.pos2-INTRA_BKPS.pos1>=150000,:);
intra_supp=[];
for i=1:length(intra)
	intra_supp=[intra_supp;INTRA_SUPP_READS(INTRA_SUPP_READS.SVidx==intra.SVidx(i),:)];
end
intra_supp=sortrows(intra_supp,{'SampleGroup','readID','SVidx'});

% inter filtering
inter=INTER_BKPS(INTER_BKPS.TotalCount-INTER_BKPS.TCount<=1 & INTER_BKPS.TCount>=2 & INTER_BKPS.NCount==0,:);
inter_supp=[];
for i=1:length(inter)
	inter_supp=[inter_supp;INTER_SUPP_READS(INTER_SUPP_READS.SVidx==inter.SVidx(i),:)];
end
inter_supp=sortrows(inter_supp,{'SampleGroup','readID','SVidx'});

% rescue additional SVs that may be linked to selected SVs by supporting reads
Supp=[];
SG=unique(intra_supp.SampleGroup,'stable');
supp=[];
for si=1:length(SG)
	A=intra_supp(intra_supp.SampleGroup==SG(si),:);
	B=INTRA_SUPP_READS(INTRA_SUPP_READS.SampleGroup==SG(si),:);
	readids=unique(A.readID,'stable');
	for ri=1:length(readids)
		supp=[supp;B(B.readID==readids(ri),:)];
	end
end
Supp=supp;

SG=unique(inter_supp.SampleGroup,'stable');
supp=[];
for si=1:length(SG)
	A=inter_supp(inter_supp.SampleGroup==SG(si),:);
	B=INTER_SUPP_READS(INTER_SUPP_READS.SampleGroup==SG(si),:);
	readids=unique(A.readID,'stable');
	for ri=1:length(readids)
		supp=[supp;B(B.readID==readids(ri),:)];
	end
end

supp=sortrows([Supp;supp],{'SampleGroup','readID','SVidx'});

% find SV reads supporting single SV events only
supp_single_idx=zeros(1,length(supp));
if supp.readID(2)>supp.readID(1)
	supp_single_idx(1)=1;
end
idx=1+find((diff(double(supp.readID(1:end-1)))~=0 | diff(double(supp.SVidx(1:end-1)))==0) & (diff(double(supp.readID(2:end)))~=0 | diff(double(supp.SVidx(2:end)))==0));
supp_single_idx(idx)=1;
if supp.readID(end)>supp.readID(end-1)
	supp_single_idx(end)=1;
end
supp_single=supp(logical(supp_single_idx),:);
Supp_single=[dataset(zeros(size(supp_single.SVidx)),'VarNames','SVGroup'),supp_single];

% process SV reads supporting multiple SV events
supp_multiple=supp(logical(1-supp_single_idx),:);
supp_multiple=sortrows(supp_multiple,{'SampleGroup','readID','pos1','SVidx'});

% Merge SV events sharing identical SV-supporting reads
[SVlist,iA,iC]=unique(supp_multiple.SVidx);
iden=zeros(length(SVlist),length(SVlist));
for i=1:length(supp_multiple)-1
	if (supp_multiple.readID(i)==supp_multiple.readID(i+1) && strcmp(supp_multiple.chr1(i),supp_multiple.chr1(i+1)) && supp_multiple.pos1(i)==supp_multiple.pos1(i+1))
		iden(iC(i),iC(i+1))=1;
	end
end
iden=iden+iden';

for i=1:length(SVlist)
	redSVs=SVlist(iden(i,:)>0);
	if length(redSVs)==0
		continue;
	end
	redSVs=[SVlist(i);redSVs];
	newid=min(redSVs);
	for j=1:length(redSVs)
		supp_multiple.SVidx(supp_multiple.SVidx==redSVs(j))=newid;
		SVlist(SVlist==redSVs(j))=newid;
	end
end

supp_multiple=sortrows(unique(supp_multiple),{'SampleGroup','readID','SVidx'});
Supp_multiple=[dataset(zeros(size(supp_multiple.SVidx)),'VarNames','SVGroup'),supp_multiple];
AllSupp=[Supp_single;Supp_multiple];

% Group SV events that are linked by supporting reads
[SVlist,iA,iC]=unique(supp_multiple.SVidx);
dist=zeros(length(SVlist),length(SVlist));
for i=1:length(supp_multiple)-1
	if (supp_multiple.readID(i)==supp_multiple.readID(i+1))
		dist(iC(i),iC(i+1))=1;
	end
end
dist=dist+dist';

SVGroup=zeros(size(SVlist));Id=0;
for i=1:length(SVlist)
	connectIdx=[i,find(dist(i,:)>0)];
	if length(connectIdx)==1
		continue;
	end
	groups=SVGroup(connectIdx);
	if max(groups)==0
		Id=Id+1;
		groupId=Id;
	else
		groupId=min(groups(groups>0));
	end	
	for j=1:length(connectIdx)
		if SVGroup(connectIdx(j))==0
			SVGroup(connectIdx(j))=groupId;
		else
			SVGroup(SVGroup==SVGroup(connectIdx(j)))=groupId;
		end
	end
end

for i=1:length(SVlist)
	AllSupp.SVGroup(AllSupp.SVidx==SVlist(i))=SVGroup(i);
end

