function plotGamma(as,bs)
% as = [1 1.5 2];
% b = 1;
% bs = b*ones(1,length(as));
% figure;
[styles, colors, symbols] = plotColors;
legendStr = cell(1, length(as));
%styles = {'k:',  'b--', 'r-'};
% textstr={'Iter.=20', 'Iter.=22', 'Iter.=25', 'Iter.=33'};
for i=1:length(as)
    a = as(i); b = bs(i);
    xs = linspace(max(0, a/b - 5* sqrt(a/(b^2))),a/b + 5*sqrt(a/(b^2)), 50);
    model.a = a;
    model.b = b;
    ps = exp(gammaLogprob(model, xs));
    style = [styles{i}, colors(i), symbols(i)];
    plot(xs , ps, style, 'linewidth', 1.5, 'markersize', 1.2);
%     plot(xs , ps, styles{i}, 'color', colors(i), 'linewidth', 1.5);
%     text(a/b,exp(gammaLogprob(model, a/b)), sprintf('a=%2.1f,b=%2.1f', a, b));
    hold on
    legendStr{i} = sprintf('a=%2.1f,\nb=%2.1f', a, b);
    %     axis tight;
end
% legend(legendStr, 'fontsize', 10);
% title('Gamma distributions');
hold off;
end

function [styles, colors, symbols, str] =  plotColors()
% Use plot(x,y,str{i}) to print in i'th style
% We use colors and linestyles, not markers

% This file is from pmtk3.googlecode.com

colors =  ['b' 'r' 'k' 'g' 'c' 'y' 'm' ...
    'r' 'b' 'k' 'g' 'c' 'y' 'm'];
symbols = ['o' 'x' '*' '>' '<' '^' 'v' ...
    '+' 'p' 'h' 's' 'd' 'o' 'x'];
styles = {'-', ':', '-.', '--', '-', ':', '-.', '--', ...
    '-', ':', '-.', '--', '-', ':', '-.', '--'};

for i=1:length(colors)
    %str{i} = sprintf('-%s%s', colors(i), symbols(i));
    str{i} = sprintf('%s%s', colors(i), styles{i});
end

end


function logp = gammaLogprob(arg1, arg2, arg3)
% logp(i) = log p(X(i) | a, b) 
% logp = gammaLogprob(model, X); OR logp = gammaLogprob(a, b, X); 
% a is the shape,
% b is the rate, i.e. 1/scale

% This file is from pmtk3.googlecode.com


if isstruct(arg1)
    model = arg1;
    X = arg2; 
    a = model.a; 
    b = model.b; 
else
    a = arg1; 
    b = arg2; 
    X = arg3; 
end

logZ = gammaln(a) - a.*log(b);
logp = (a-1).*log(X) - b.*X - logZ;
logp = logp(:); 



end
