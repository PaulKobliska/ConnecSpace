function [ ret ] = CorrelateXYAndRegression( in )
%Useful to Draw Data
%   Detailed explanation goes here
x=in.x;
y=in.y;
x = reshape(x,[],1);
y = reshape(y,[],1);

%[r,p_value] = nancorr(x,y,'Pearson');
ind = 1:length([x]);
kx = find(isnan(x));kx = reshape(kx,[],1);
ky = find(isnan(y));ky = reshape(ky,[],1);
NanData = unique([kx ;ky]);
k = setdiff(ind,NanData);
if numel(k) > 2
p = polyfit(x(ind(k)),y(ind(k)),1);
else
p=[NaN NaN];    %p = polyfit(x,y,1)
end
% ret.r = r;
% ret.p_value = p_value;
% ret.slope = p(1);
% ret.intercept = p(2);
if ~isfield(in,'dotscolor')  ;  in.dotscolor = [0 0 0];end 
if ~isfield(in,'linecolor')  ;  in.linecolor = 'r';end 
if ~isfield(in,'linewidth')  ;  in.linewidth = 3;end 
if ~isfield(in,'linetype')  ;  in.linetype = '-';end 
if ~isfield(in,'xlabel')  ;  in.xlabel = 'x';end 
if ~isfield(in,'ylabel')  ;  in.ylabel = 'y';end 
if ~isfield(in,'title')  ;  in.title = 'Data';end 
if ~isfield(in,'PLOT_ON')  ;  in.PLOT_ON = 1;end 
if ~isfield(in,'xlim')  ;  in.xlim=[min(x) max(x) ] ;end 
if ~isfield(in,'ylim')  ;  in.ylim=[min(y) max(y) ] ;end 
if ~isfield(in,'dispPvalue'); in.dispPvalue = 1;end;
if ~isfield(in,'Test'); in.Test = 'Pearson';end;
if ~isfield(in,'MarkerType'); in.MarkerType = 'o';end;
if ~isfield(in,'legend'); in.legend = 0;end;

[r,p_value] = nancorr2(x,y,in.Test);

ret.r = r;
ret.p_value = p_value;
ret.slope = p(1);
ret.intercept = p(2);

if in.PLOT_ON 

if numel(k) > 2
    
hold on ;
if strcmp(in.MarkerType,'o')
    hdot = scatter(x,y,'filled');
else
    hdot = scatter(x,y);
end
% if in.legend;
% ax_children = get(gca, 'children');
% leg_strings = get(ax_children, 'displayname')
% if strcmp(leg_strings(end),'');
% legend({in.legend_variable});
% else
% legend({leg_strings(end) in.legend_variable});
% end   
% end

set(hdot,'CData',in.dotscolor,'Marker',in.MarkerType) ; 
hline = line([min(x) max(x)],[(min(x)*ret.slope + ret.intercept) (max(x)*ret.slope + ret.intercept)],'color',in.linecolor,'LineWidth',in.linewidth ,'LineStyle',in.linetype) ; 
ret.hdot = hdot;
ret.hline = hline;
end
xlabel(in.xlabel);
ylabel(in.ylabel);

axis square; xlim(in.xlim);
ylim(in.ylim);
set(gca,'TickDir','out','Box','off');
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 12, 'linewidth', 2);
set(gca,'TickLength',[0.02,0.01]);


if in.dispPvalue == 1
if p_value > 0.05    
    title([in.title ', N.S. (' in.Test ')'],'fontsize',12);
elseif p_value >0.01 & p_value < 0.05
    title([in.title ', * (' in.Test ')'],'fontsize',12);
elseif p_value >0.001 & p_value < 0.01
    title([in.title ', ** (' in.Test ')'],'fontsize',12);
else title([in.title ', *** (' in.Test ')'],'fontsize',12);
end

else
if numel(k) > 3;title([in.title]);else  title([in.title ' missing'],'fontsize',12 ); end;
end;
    
    hold off

end
    
end

