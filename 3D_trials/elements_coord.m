% coordinates of array elements for 192 elements probe ('L12-3V')

Ne = 192;             
p  = 200e-6;      

x_i = zeros(1, Ne);
y_i = zeros(1, Ne);
z_i = zeros(1, Ne);

for idx = 1:ne
    x_i(idx) = (p/2) * (2*idx - Ne - 1);
    y_i(idx) = 0;
    z_i(idx) = 0;
end
