%insert the M-1 0s in the sequence and expand length of d  to len(d)*M
function[out] = sigexpand(d,M)
N = length(d);
out = zeros(M,N);
out(1,:) = d;
out = reshape(out,1,M*N); 