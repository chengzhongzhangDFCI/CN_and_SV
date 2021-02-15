function plot_RA(RAEvents,Yaxis,height,varargin)
for i=1:length(RAEvents)
    x0=min(RAEvents.pos1(i),RAEvents.pos2(i));
    x1=max(RAEvents.pos1(i),RAEvents.pos2(i));
    if (isnumeric(RAEvents.chr1(i)) && isnumeric(RAEvents.chr2(i)) && RAEvents.chr1(i)==RAEvents.chr2(i)) || strcmp(RAEvents.chr1(i),RAEvents.chr2(i))
        x=linspace(x0,x1,100);
        y=Yaxis+sqrt(1.01-(2*x/(x1-x0)-(x1+x0)/(x1-x0)).^2)*height;
        plot(x,y,varargin{:});
    else
        plot([x0,x0],[Yaxis,Yaxis+height],varargin{:});
    end
end
