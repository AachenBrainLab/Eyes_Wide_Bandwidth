function binOut = pC_dec2bin(decOut,n)
% faster version of dec2bin but with binary output instead of string.

[~,e] = log2(max(decOut)); % How many digits do we need to represent the numbers?
binOut = (rem(floor(decOut.*pow2(1-max(n,e):0)),2));