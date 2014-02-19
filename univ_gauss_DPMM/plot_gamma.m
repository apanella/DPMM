function plot_gamma(a, b, max_x)

if nargin<3
    max_x = 10;
end
x = [0:0.01:max_x];
y = gampdf(x, a, 1/b);
plot(x,y);
ylim([0,max(y)]);

end