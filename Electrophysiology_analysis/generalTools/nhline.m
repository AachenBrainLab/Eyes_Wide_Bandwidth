function lines = nhline(x,varargin)

a = gca;
if ishold(a)
    checker = false;
else
    hold(a,'on'); checker = true;
end

for xx = 1 : length(x)
    lines(xx) = plot(a.XLim,[x(xx),x(xx)],varargin{:});
end

if checker
    hold(a,'off');
end
