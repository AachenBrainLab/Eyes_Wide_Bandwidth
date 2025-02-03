function lines = nvline(x,varargin)

a = gca;
if ishold(a)
    checker = true;
else
    hold(a,'on'); checker = false;
end

for xx = 1 : length(x)
    lines(xx) = plot([x(xx),x(xx)],a.YLim,varargin{:});
end
if ~checker
    hold(a,'off');
end