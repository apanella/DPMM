function plot_IG(a,b)

invgamma_pdf = @(a,b,x) gamma(a)^(-1) * b^a .* x.^(-a-1) .* exp(-b./x);

x = 0:0.01:400;
y = invgamma_pdf(a,b,x);

figure;
plot(x,y);

end
