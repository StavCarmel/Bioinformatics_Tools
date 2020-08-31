function [ GC ] = CG( seq )
count=basecount(seq);
N_G = count.G;
N_T = count.T;
N_C = count.C;
N_A = count.A;
GC = ((N_G+N_C)/(N_G+N_C+N_T+N_A))*100;
end


