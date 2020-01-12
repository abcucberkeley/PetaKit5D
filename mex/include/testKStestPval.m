function test_kstest_pval
% Test the relative difference between the p-value calculated in the Matlab
% built-in kstest function and the p-value coming from the ksone function
% in stats.h
%

% Sebastien Besson, July 2011

% Clear the workspace
clear
clc
close all

% Initialize vectors and pvalues vectors
N=10:1:1000;
D=.05:.1:.95;
pvalue_nr=zeros(length(N),length(D));
pvalue_matlab=zeros(length(N),length(D));

% Fill the pvalues (not optimized)
for i=1:numel(D)
    d=D(i);
    for j=1:numel(N)
        n=N(j);
        pvalue_nr(j,i)=qks((sqrt(n) + 0.12 + 0.11/sqrt(n))*d);
        pvalue_matlab(j,i)=getKSval(d, n);
    end
end

% define small and large fonts
tfont = {'FontName', 'Arial', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Arial', 'FontSize', 18};
lfont = {'FontName', 'Arial', 'FontSize', 22};

% Figure 1
% plot the evolution of pvalues as a function of N for various D
f(1)=figure('Position',[200 200 600 600],'PaperPositionMode', 'auto');
h1=semilogx(N,pvalue_nr,'ob','Linewidth',2);
hold on
h2=semilogx(N,pvalue_matlab,'xr','Linewidth',2);
l=legend([h1(1) h2(1)],'stats.h','Matlab built-in','Location', 'NorthEast', tfont{:});
legend(l, 'boxoff');
annotation('textarrow',[.7 .2],[.7 .2],'Linewidth',2,...
    'String','Increasing D',tfont{:});

set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('N',lfont{:});
ylabel('p-value',lfont{:});
box off;

% Select a subset of plot (low D-values)
nD=5;

% Figure 2
% plot the deviation of of pvalues on the interval [0 1]
f(2)=figure('Position',[200 200 600 600],'PaperPositionMode', 'auto');
set(gca,'ColorOrder',hsv(nD),'NextPlot','replaceChildren');
plot(pvalue_matlab(:,1:nD),(pvalue_matlab(:,1:nD)-pvalue_nr(:,1:nD)),'-','Linewidth',2);
l=legend(arrayfun(@(x) ['D = ' num2str(x)],D(1:nD),'Unif',false),'Location', 'SouthWest', tfont{:});
legend(l, 'boxoff'); %

% Set thickness of axes, ticks and assign tick labels
set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('p-value_{Matlab}',lfont{:})
ylabel('p-value_{Matlab}-pvalue_{stats.h}',lfont{:});
box off;

% Figure 3
% plot the relative deviation of of pvalues on the interval [0.05 .1]
f(3)=figure('Position',[200 200 600 600],'PaperPositionMode', 'auto'); % enable resizing
set(gca,'ColorOrder',hsv(nD),'NextPlot','replaceChildren');
plot(pvalue_matlab(:,1:nD),(pvalue_matlab(:,1:nD)-pvalue_nr(:,1:nD))./pvalue_matlab(:,1:nD),'-','Linewidth',2);
l=legend(arrayfun(@(x) ['D = ' num2str(x)],D(1:nD),'Unif',false),'Location', 'SouthWest', tfont{:});
legend(l, 'boxoff'); 

% Set thickness of axes, ticks and assign tick labels
set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('p-value_{Matlab}',lfont{:})
ylabel('(p-value_{Matlab}-p-value_{stats.h})/p-value_{Matlab}',lfont{:});
axis([.005 .1 -0.04 0.02])
box off;

% Save the figures as png
for i=1:3
    figure(f(i));
    print('-dpng','-r300',['Figure' num2str(i) '.png'])
end


function q = qks(z)
% Adapted from stats.h

if z < 0.0
    error('z value for KS dist. must be positive.');
end
if z == 0.0
    q = 1.0;
elseif z < 1.8
    y = exp(-pi^2/(8*z^2));
    q = 1.0 - sqrt(2*pi) / z * (y + y^9 + y^25 + y^49);
else
    x = exp(-2*z^2);
    q = 2.0*(x - x^4 - x^9);
end


function pValue = getKSval(KSstatistic, n)
% Copied from kstest.m

s = n*KSstatistic^2;

% For d values that are in the far tail of the distribution (i.e.
% p-values > .999), the following lines will speed up the computation
% significantly, and provide accuracy up to 7 digits.
if (s > 7.24) ||((s > 3.76) && (n > 99))
    pValue = 2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
    
else
    % Express d as d = (k-h)/n, where k is a +ve integer and 0 < h < 1.
    k = ceil(KSstatistic*n);
    h = k - KSstatistic*n;
    m = 2*k-1;
    
    % Create the H matrix, which describes the CDF, as described in Marsaglia,
    % et al.
    if m > 1
        c = 1./gamma((1:m)' + 1);
        
        r = zeros(1,m);
        r(1) = 1;
        r(2) = 1;
        
        T = toeplitz(c,r);
        
        T(:,1) = T(:,1) - (h.^[1:m]')./gamma((1:m)' + 1);
        
        T(m,:) = fliplr(T(:,1)');
        T(m,1) = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
    else
        T = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
    end
    
    % Scaling before raising the matrix to a power
    if ~isscalar(T)
        lmax = max(eig(T));
        T = (T./lmax)^n;
    else
        lmax = 1;
    end
    
    % Pr(Dn < d) = n!/n * tkk ,  where tkk is the kth element of Tn = T^n.
    % p-value = Pr(Dn > d) = 1-Pr(Dn < d)
    pValue = (1 - exp(gammaln(n+1) + n*log(lmax) - n*log(n)) * T(k,k));
end
