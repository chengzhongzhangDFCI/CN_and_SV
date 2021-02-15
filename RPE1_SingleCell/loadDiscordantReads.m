%% Discordant reads processing
insertSize=300;
readLength=150;
shortRange=10000;

fprintf(1,'   Load discordant reads from RPE-1 bulk..');
load /czlab/Data/RPE-1/bulk_standard/RPE-Bulk2_discordantR1.mat;
load /czlab/Data/RPE-1/bulk_standard/RPE-Bulk2_discordantR2.mat;
load /czlab/Data/RPE-1/bulk_standard/RPE-Bulk2_pairInfo.mat;
pair=[PAIR,dataset(logical(R1.fragment==R2.fragment),CHRs(R1.chr),int32(double(R1.Apos)+double(R1.Aspan)),int16(abs(R1.Aspan)),int8((R1.Aspan>0)-(R1.Aspan<0)),R1.maq,R1.Rpos,int16(R1.Rspan),...
				CHRs(R2.chr),int32(double(R2.Apos)),int16(abs(R2.Aspan)),int8((R2.Aspan<0)-(R2.Aspan>0)),R2.maq,R2.Rpos,int16(R2.Rspan),...
		 'VarNames',{'split','chr1','pos1','aspan1','str1','maq1','rpos1','rspan1','chr2','pos2','aspan2','str2','maq2','rpos2','rspan2'})];
pair=pair(~pair.dup,:);
splits=pair(pair.split,:);
splits=splits(diff(splits.readID)~=0 | diff(splits.pos1)>2 | diff(splits.pos2)>2,:); % It's unlikely that the two split parts map to different chromosomes at similar positions;
nonsplits=pair(~pair.split,:);
chimeras=[splits(~strcmp(splits.chr1,splits.chr2),:);nonsplits(~strcmp(nonsplits.chr1,nonsplits.chr2),:)];
discordants=[splits(strcmp(splits.chr1,splits.chr2),:);nonsplits(strcmp(nonsplits.chr1,nonsplits.chr2) & (nonsplits.str1==nonsplits.str2 | abs(nonsplits.pos2-nonsplits.pos1)>readLength),:)];
short_inv=discordants(discordants.str1==discordants.str2 & abs(discordants.pos2-discordants.pos1)<=shortRange,:);
short_noninv=discordants(discordants.str1~=discordants.str2 & abs(discordants.pos2-discordants.pos1)<=shortRange,:);
discordants=discordants(abs(discordants.pos2-discordants.pos1)>shortRange,:);
CHIMERAS=[dataset(repmat(0,size(chimeras.chr1)),'VarNames','SampleGroup'),chimeras];
DISCORDANTS=[dataset(repmat(0,size(discordants.chr1)),'VarNames','SampleGroup'),discordants];
SHORT_INV=[dataset(repmat(0,size(short_inv.chr1)),'VarNames','SampleGroup'),short_inv];
SHORT_NONINV=[dataset(repmat(0,size(short_noninv.chr1)),'VarNames','SampleGroup'),short_noninv];
fprintf(1,'done\n');
discordant_folder='DiscordantReads';
for si=1:length(SampleInfo.SampleNames)
	SampleName=SampleInfo.SampleNames{si};
	output_MAT=true;
	if exist(fullfile(discordant_folder,[SampleName,'_pairInfo.mat'])) 
		fprintf(1,['Loading discordant reads information for sample ' SampleName ' from MAT files..']);
		try	
			load([discordant_folder '/' SampleName,'_discordantR1.mat']);
			load([discordant_folder '/' SampleName,'_discordantR2.mat']);
			load([discordant_folder '/' SampleName,'_pairInfo.mat']);
			fprintf(1,'done\n');
			output_MAT=false;
		catch
			fprintf(1,'MAT files are broken\n');
			output_MAT=true;
		end
	end
	if (output_MAT)
		try
		read_file=fullfile(discordant_folder,[regexprep(SampleName,'SM_','') '.discordant_reads.txt']);
		if ~exist(read_file)
			continue;
		end
		fid=fopen(read_file,'r');
		fprintf(1,['Collecting discordant reads from ' read_file '\n']);
		%read_id chr start end read_start read_end maq duplicate fragment RGindex	
		temp=textscan(fid,'%d\t%s\t%d\t%d\t%d\t%d\t%u8\t%u8\t%u8\t%u8\n','HeaderLines',1);	
		DATA=dataset(convertIntValue(temp{1}),temp{2},convertIntValue(temp{3}),convertIntValue(temp{4}-temp{3}),convertIntValue(temp{5}),convertIntValue(temp{6}-temp{5}),...
			convertIntValue(temp{7}),temp{8},temp{9},temp{10},...
			'VarNames',{'readID','chr','Apos','Aspan','Rpos','Rspan','maq','dup','fragment','rgID'});
		[CHRs,idx_i,idx_j]=unique(DATA.chr);
		DATA.chr=convertIntValue(idx_j);
		idx=find(diff(DATA.readID)==0);
		PAIR=DATA(idx,{'readID','rgID','dup'});
		R1=DATA(idx,{'chr','Apos','Aspan','Rpos','Rspan','maq','fragment'});
		R2=DATA(idx+1,{'chr','Apos','Aspan','Rpos','Rspan','maq','fragment'});
		save([discordant_folder '/' SampleName,'_discordantR1.mat'],'R1','CHRs');
		save([discordant_folder '/' SampleName,'_discordantR2.mat'],'R2','CHRs');
		save([discordant_folder '/' SampleName,'_pairInfo.mat'],'PAIR');
		fprintf(1,'..Successful.\n');
		catch
			fprintf(1,'..Failed!\n');
		end
	end
	fprintf(1,'   Collect discordant reads for clustering..');
	pair=[PAIR,dataset(logical(R1.fragment==R2.fragment),CHRs(R1.chr),int32(double(R1.Apos)+double(R1.Aspan)),int16(abs(R1.Aspan)),int8((R1.Aspan>0)-(R1.Aspan<0)),R1.maq,R1.Rpos,int16(R1.Rspan),...
				CHRs(R2.chr),int32(double(R2.Apos)),int16(abs(R2.Aspan)),int8((R2.Aspan<0)-(R2.Aspan>0)),R2.maq,R2.Rpos,int16(R2.Rspan),...
		 'VarNames',{'split','chr1','pos1','aspan1','str1','maq1','rpos1','rspan1','chr2','pos2','aspan2','str2','maq2','rpos2','rspan2'})];
	pair=pair(~pair.dup,:);
	splits=pair(pair.split,:);
	splits=splits(diff(splits.readID)~=0 | diff(splits.pos1)>2 | diff(splits.pos2)>2,:); % It's unlikely that the two split parts map to different chromosomes at similar positions;	
	nonsplits=pair(~pair.split,:);
	chimeras=[splits(~strcmp(splits.chr1,splits.chr2),:);nonsplits(~strcmp(nonsplits.chr1,nonsplits.chr2),:)];
	discordants=[splits(strcmp(splits.chr1,splits.chr2),:);nonsplits(strcmp(nonsplits.chr1,nonsplits.chr2) & (nonsplits.str1==nonsplits.str2 | abs(nonsplits.pos2-nonsplits.pos1)>readLength),:)];
	short_inv=discordants(discordants.str1==discordants.str2 & abs(discordants.pos2-discordants.pos1)<=shortRange,:);
	short_noninv=discordants(discordants.str1~=discordants.str2 & abs(discordants.pos2-discordants.pos1)<=shortRange,:);
	discordants=discordants(abs(discordants.pos2-discordants.pos1)>shortRange,:);
	CHIMERAS=[CHIMERAS;dataset(repmat(si,size(chimeras.chr1)),'VarNames','SampleGroup'),chimeras];
	DISCORDANTS=[DISCORDANTS;dataset(repmat(si,size(discordants.chr1)),'VarNames','SampleGroup'),discordants];
	SHORT_INV=[SHORT_INV;dataset(repmat(si,size(short_inv.chr1)),'VarNames','SampleGroup'),short_inv];
	SHORT_NONINV=[SHORT_NONINV;dataset(repmat(si,size(short_noninv.chr1)),'VarNames','SampleGroup'),short_noninv];
	fprintf(1,'done\n');
end
fprintf(1,'Complete loading discordant reads\n');
