%test bessel function of first kind

X = 0:0.1:20;
J = zeros(5, 201);

for i = 0:4
    J(i+1,:) = besselj(i,X);
end

plot(X,J,'LineWidth',1.5)
axis([0 20 -0.5 1]);