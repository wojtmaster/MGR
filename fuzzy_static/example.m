clear
u = -2:0.1:2;
y = zeros(size(u));
Y = zeros(size(u));

for i = 1:length(u)
    y(i) = 0.3163*u(i) / (sqrt(0.1+0.9*u(i)^2));
end

R = cell(1,2);
r(1) = -0.3289;
r(2) = 0.3289;
alpha(1) = -5;
alpha(2) = 5;
c(1) = 0;
c(2) = 0;

for i = 1:length(R)
    for j = 1:length(u)
        R{i}(j) = 1 / (1+exp(-alpha(i)*(u(j)-c(i))));
    end
end

for i = 1:length(u)
    w = 0;
    for j = 1:length(R)
        Y(i) = Y(i) + R{j}(i)*r(j);
        w = w + R{j}(i);
    end
    Y(i) = Y(i) / w;
end
figure;
plot(u, y, u, Y);
grid on;
figure;
plot(u, R{1}, u, R{2});
grid on;