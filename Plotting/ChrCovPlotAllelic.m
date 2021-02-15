function ChrCovPlotAllelic(chrlabel,Pos,AlleleA,AlleleB)

load hg38.mat

maxLength=max(Pos);
maxScale=log10(double(maxLength)/1e6);
XScale=maxScale*500;
YScale=6.2;

figure,
set(gcf, 'Unit', 'pixels', 'Position', [15,25,XScale+200 50*(YScale+2.5)]),hold on;
AlleleA=AlleleA.*(AlleleA<=2)+(AlleleA>2).*(log(AlleleA+(AlleleA<2))/log(2)+1);
AlleleB=AlleleB.*(AlleleB<=2)+(AlleleB>2).*(log(AlleleB+(AlleleB<2))/log(2)+1);
AlleleA(AlleleA>5)=5;
AlleleB(AlleleB>5)=5;
plot(Pos/1e6,AlleleA,'b.','MarkerEdgeColor',[0,0,0.5],'MarkerSize',12);
plot(Pos/1e6,AlleleB,'r.','MarkerEdgeColor',[0.85,0,0],'MarkerSize',12);	
gBands=GiemsaBands(strcmp(GiemsaBands.Chr,chrlabel),:);
gBands.Start=double(gBands.Start)/1e6;gBands.End=double(gBands.End)/1e6;
gBandColor=gBands.Stain;
gBandColor(strcmp(gBandColor,'acen'))={[0.9,0,0]};
gBandColor(strcmp(gBandColor,'gpos25'))={[0.75,0.75,0.75]};
gBandColor(strcmp(gBandColor,'gpos50'))={[0.5,0.5,0.5]};
gBandColor(strcmp(gBandColor,'gpos75'))={[0.25,0.25,0.25]};
gBandColor(strcmp(gBandColor,'gpos100'))={[0,0,0]};
gBandColor(strcmp(gBandColor,'stalk'))={[0,0.5,0]};
gBandColor(strcmp(gBandColor,'gvar'))={[0.5,0.25,0]};
gBandColor(strcmp(gBandColor,'gneg'))={[1,1,1]};
rectangle('Position',[gBands.Start(1),-0.95,gBands.End(end)-gBands.Start(1),0.6],'LineWidth',1.5);
for gi=1:length(gBands)
	rectangle('Position',[gBands.Start(gi),-0.95,gBands.End(gi)-gBands.Start(gi),0.6],'FaceColor',gBandColor{gi},'EdgeColor','none');
end
box off,
xmax=ceil(maxLength/1e6);
ymax=YScale;
ymin=-0.95;
xmin=-maxScale*2;
set(gca,'Ylim',[-0.95 ymax],'Xlim',[xmin,ceil(xmax/10)*10]);
XTickPos=0:10:floor(xmax);
minorTicks=[];
for i=1:length(XTickPos)-1
	minorTicks=[minorTicks,(XTickPos(i)+1):1:(XTickPos(i+1)-1)];
end
minorTicks=[minorTicks,XTickPos(end)+1:1:xmax];
set(gca,'TickDir','out','FontSize',16,'XTick',XTickPos,'XMinorTick','On');
set(get(gca,'XAxis'),'MinorTickValues',minorTicks);
set(get(gca,'YAxis'),'Visible','off');

xunit=xmax*0.025/maxScale;
plot([-1,-2]*xunit,[0,0],'k-','LineWidth',1.5)
plot([-1,-2]*xunit,[1,1],'k-','LineWidth',1.5)
plot([-1,-2]*xunit,[2,2],'k-','LineWidth',1.5)
plot([-1,-2]*xunit,[3,3],'k-','LineWidth',1.5)
plot([-1,-2]*xunit,[4,4],'k-','LineWidth',1.5)
minor_ticks=[3,5,6,7,9,10];
for i=1:length(minor_ticks)
	plot([-1,-1.5]*xunit,(log(minor_ticks(i))/log(2)+1)*[1,1],'k-','LineWidth',1);
end
plot(-1*xunit*[1,1],[0,log(10)/log(2)+1],'k-','LineWidth',1.5)
plot([-1*xunit maxLength],[0 0],'k:','LineWidth',1);
plot([-1*xunit maxLength],[1 1],'k:','LineWidth',1);
plot([-1*xunit maxLength],[2 2],'k:','LineWidth',1);

set(gca,'Unit','pixels','Position',[100 50 XScale (YScale+1)*50]);
set(gcf,'Unit', 'inches');
Position=get(gcf,'Position');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', Position(3:4), 'PaperPosition', [0,0,Position(3:4)]);
