% specify order
order = 134;

% classic method
B1 = zeros(1,order+2);
for n=0:order
    B1(n+1) = factorial(2*order-n) / ...
        (factorial(order-n)*factorial(n)*2^(order-n));
end

% ratio method
B2 = zeros(1,order+2);
for n=0:order
    B2(n+1) = factorialratio(2*order-n,order-n) / ...
        (factorial(n)*2^(order-n));
end

% Stirlings approximation
B3 = zeros(1,order+2);
for n=0:order
    B3(n+1) = sqrt(2) * 2^(2*order-n) * ((order-n)/exp(1))^order / ...
        (sqrt(2*pi*order) * (order/exp(1))^order * 2^(order-n));
end
