
tic;
clear

gamma=1;

L=32;
N=L/2;

N_sample = 64;

dt = 0.05;
tmp_tf = 3*L;
ts = 0:dt:tmp_tf;

subA = (1:L/2)';

H_temp = 1*diag(ones(L-1,1),1);
H_temp(L,1) = 1;
H_temp = H_temp+H_temp';
mat_H =  H_temp;

mat_expH=expm(1i*mat_H*dt);
psi_0 = zeros(L,L/2);
for ii=1:L/2
    psi_0(2*ii-1,ii) = 1;
end
psi_0=sparse(psi_0);

sigma_gdt = sqrt(gamma*dt);

S_A = zeros(N_sample,1);

pos_jp1 = circshift(1:L,-1);
pos_jm1 = circshift(1:L,1);
for ss=1:N_sample
 
    rng("shuffle");
    %rng(ss);
    mat_phi = full(psi_0);

    for tt=2:length(ts)

        D=mat_phi*mat_phi';
        ni_t = diag(D);
        Dij_t = diag(D,1);
        Dij_t(L) =  D(L,1);
        mat_D = ni_t - 2*imag(Dij_t) + ni_t(pos_jp1);
        vec_M = (2*mat_D-1)*gamma*dt;

        Mst=(normrnd(0,sigma_gdt,[L,1]));
        vec_M=Mst+vec_M;

        mat_M = 1i*diag(vec_M(1:L-1),1);
        mat_M(L,1) = 1i*vec_M(L);
        mat_M = mat_M - 1i*diag(vec_M(1:L-1),-1);
        mat_M(1,L) = -1i*vec_M(pos_jp1(L));
        mat_M = mat_M + spdiags(vec_M,0,L,L)+ spdiags(vec_M(pos_jm1),0,L,L);

        mat_phi=expm(mat_M)*(mat_expH*mat_phi);

        [mat_phi,~] = qr(mat_phi,'econ');
    end

    D=mat_phi*mat_phi';
    S_A(ss) = entanglement_entroy_D(D,subA);

 
end
file_name = sprintf('QSD_nodal_L%d_gm%.2f_Nsam%d_SA.mat',L,gamma,N_sample);
save(file_name,'L','N','N_sample','gamma','S_A','-v7.3')
%save(file_name,'-v7.3')

toc
