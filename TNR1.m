%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TNR 1
%
% Andrew Goldsborough - 14/03/2016
%
% basic TNR algorithm for the 2D Ising partition function
% for details see arXiv:1509.07484
%
% index convention
%   -1     -1  -3      -1     -1  -2
%    |      |___|       |       \ / 
%   /_\     |___|   -4 -o- -2    o
%   | |     |   |       |       / \
% -2  -3   -2  -4      -3     -4  -3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%input variables
beta = 0.6;
max_level = 10;
chi_u = 6;
chi_w = 6;
chi_v = 6;
chi_y = 6;
UL_max = 5;
loop_max = 100;

%convergence criteria
UL_convergence = 1e-10;
B_convergence = 1e-7;
D_convergence = 1e-7;
A_convergence = 1e-7;

%build Ising transfer matrix. Break the symmetry with small perturbation
epsilon = 0.1;
Q = [(1+epsilon/4)*exp(beta) exp(-beta); exp(-beta) (1-epsilon/4)*exp(beta)];

%define delta
delta3 = zeros(2,2,2);
for i=1:2
    delta3(i,i,i) = 1;
end

%create plaquette operator
A = ncon({delta3,Q,delta3,Q,delta3,Q,delta3,Q},...
    {[-1,8,1],[1,2],[-2,2,3],[3,4],[-3,4,5],[5,6],[-4,6,7],[7,8]});

%normalise
A = A./sqrt(ncon({A,conj(A)},{[1,2,3,4],[1,2,3,4]}));

%find contraction orders using netcon
B_netcon(netcon({[-1,1,2],[2,5,3,6],[-2,3,4],[5,9,7,1],[6,4,8,9],[11,12,7,10],[13,14,8,12],[15,11,16,13],[-4,10,15],[-3,16,14]},0,2,1,1)) = 1:16;
B_order = {[-1,B_netcon(1),B_netcon(2)],...
    [B_netcon(2),B_netcon(5),B_netcon(3),B_netcon(6)],...
    [-2,B_netcon(3),B_netcon(4)],...
    [B_netcon(5),B_netcon(9),B_netcon(7),B_netcon(1)],...
    [B_netcon(6),B_netcon(4),B_netcon(8),B_netcon(9)],...
    [B_netcon(11),B_netcon(12),B_netcon(7),B_netcon(10)],...
    [B_netcon(13),B_netcon(14),B_netcon(8),B_netcon(12)],...
    [B_netcon(15),B_netcon(11),B_netcon(16),B_netcon(13)],...
    [-4,B_netcon(10),B_netcon(15)],...
    [-3,B_netcon(16),B_netcon(14)]};

D_netcon(netcon({[-1,2,1],[1,3,4,2],[-2,3,4]},0,2,1,1)) = 1:4;
D_order = {[-1,D_netcon(2),D_netcon(1)],...
    [D_netcon(1),D_netcon(3),D_netcon(4),D_netcon(2)],...
    [-2,D_netcon(3),D_netcon(4)]};

A_netcon(netcon({[-1,1,2],[4,1,3],[5,3,2],[-4,6],[6,8,4],[7,5,9],[7,-2],[8,11,10],[9,10,12],[-3,11,12]},0,2,1,1)) = 1:12;
A_order = {[-1,A_netcon(1),A_netcon(2)],...
    [A_netcon(4),A_netcon(1),A_netcon(3)],...
    [A_netcon(5),A_netcon(3),A_netcon(2)],...
    [-4,A_netcon(6)],...
    [A_netcon(6),A_netcon(8),A_netcon(4)],...
    [A_netcon(7),A_netcon(5),A_netcon(9)],...
    [A_netcon(7),-2],...
    [A_netcon(8),A_netcon(11),A_netcon(10)],...
    [A_netcon(9),A_netcon(10),A_netcon(12)],...
    [-3,A_netcon(11),A_netcon(12)]};

env_v_l_netcon(netcon({[-2,2,1,3],[17,1,4],[2,5,6,-1],[3,4,7,5],[10,8,6,9],[11,12,7,8],[13,10,14,11],[15,9,13],[16,14,12],[-3,17,16,15]},0,2,1,1)) = 1:17;
env_v_l_order = {[-2,env_v_l_netcon(2),env_v_l_netcon(1),env_v_l_netcon(3)],...
    [env_v_l_netcon(17),env_v_l_netcon(1),env_v_l_netcon(4)],...
    [env_v_l_netcon(2),env_v_l_netcon(5),env_v_l_netcon(6),-1],...
    [env_v_l_netcon(3),env_v_l_netcon(4),env_v_l_netcon(7),env_v_l_netcon(5)],...
    [env_v_l_netcon(10),env_v_l_netcon(8),env_v_l_netcon(6),env_v_l_netcon(9)],...
    [env_v_l_netcon(11),env_v_l_netcon(12),env_v_l_netcon(7),env_v_l_netcon(8)],...
    [env_v_l_netcon(13),env_v_l_netcon(10),env_v_l_netcon(14),env_v_l_netcon(11)],...
    [env_v_l_netcon(15),env_v_l_netcon(9),env_v_l_netcon(13)],...
    [env_v_l_netcon(16),env_v_l_netcon(14),env_v_l_netcon(12)],...
    [-3,env_v_l_netcon(17),env_v_l_netcon(16),env_v_l_netcon(15)]};

env_v_r_netcon(netcon({[1,2,3],[3,4,-1,5],[4,6,7,2],[5,-2,8,6],[11,9,7,10],[12,13,8,9],[14,11,15,12],[16,10,14],[17,15,13],[1,-3,17,16]},0,2,1,1)) = 1:17;
env_v_r_order = {[env_v_r_netcon(1),env_v_r_netcon(2),env_v_r_netcon(3)],...
    [env_v_r_netcon(3),env_v_r_netcon(4),-1,env_v_r_netcon(5)],...
    [env_v_r_netcon(4),env_v_r_netcon(6),env_v_r_netcon(7),env_v_r_netcon(2)],...
    [env_v_r_netcon(5),-2,env_v_r_netcon(8),env_v_r_netcon(6)],...
    [env_v_r_netcon(11),env_v_r_netcon(9),env_v_r_netcon(7),env_v_r_netcon(10)],...
    [env_v_r_netcon(12),env_v_r_netcon(13),env_v_r_netcon(8),env_v_r_netcon(9)],...
    [env_v_r_netcon(14),env_v_r_netcon(11),env_v_r_netcon(15),env_v_r_netcon(12)],...
    [env_v_r_netcon(16),env_v_r_netcon(10),env_v_r_netcon(14)],...
    [env_v_r_netcon(17),env_v_r_netcon(15),env_v_r_netcon(13)],...
    [env_v_r_netcon(1),-3,env_v_r_netcon(17),env_v_r_netcon(16)]};

env_u_netcon(netcon({[1,2,-2],[3,-4,4],[-1,5,6,2],[-3,4,7,5],[10,8,6,9],[11,12,7,8],[13,10,14,11],[15,9,13],[16,14,12],[1,3,16,15]},0,2,1,1)) = 1:16;
env_u_order = {[env_u_netcon(1),env_u_netcon(2),-2],...
    [env_u_netcon(3),-4,env_u_netcon(4)],...
    [-1,env_u_netcon(5),env_u_netcon(6),env_u_netcon(2)],...
    [-3,env_u_netcon(4),env_u_netcon(7),env_u_netcon(5)],...
    [env_u_netcon(10),env_u_netcon(8),env_u_netcon(6),env_u_netcon(9)],...
    [env_u_netcon(11),env_u_netcon(12),env_u_netcon(7),env_u_netcon(8)],...
    [env_u_netcon(13),env_u_netcon(10),env_u_netcon(14),env_u_netcon(11)],...
    [env_u_netcon(15),env_u_netcon(9),env_u_netcon(13)],...
    [env_u_netcon(16),env_u_netcon(14),env_u_netcon(12)],...
    [env_u_netcon(1),env_u_netcon(3),env_u_netcon(16),env_u_netcon(15)]};

env_y_l_netcon(netcon({[-2,1,2,-1],[3,1,2],[-3,3]},0,2,1,1)) = 1:3;
env_y_l_order = {[-2,env_y_l_netcon(1),env_y_l_netcon(2),-1],...
    [env_y_l_netcon(3),env_y_l_netcon(1),env_y_l_netcon(2)],...
    [-3,env_y_l_netcon(3)]};

env_y_r_netcon(netcon({[1,3,2],[2,-1,-2,3],[1,-3]},0,2,1,1)) = 1:3;
env_y_r_order = {[env_y_r_netcon(1),env_y_r_netcon(3),env_y_r_netcon(2)],...
    [env_y_r_netcon(2),-1,-2,env_y_r_netcon(3)],...
    [env_y_r_netcon(1),-3]};

env_w_netcon(netcon({[2,-1,1],[3,1,-2],[4,6,2],[5,3,7],[6,9,8],[7,8,10],[11,4],[13,9,10],[5,12],[-3,12,13,11]},0,2,1,1)) = 1:13;
env_w_order = {[env_w_netcon(2),-1,env_w_netcon(1)],...
    [env_w_netcon(3),env_w_netcon(1),-2],...
    [env_w_netcon(4),env_w_netcon(6),env_w_netcon(2)],...
    [env_w_netcon(5),env_w_netcon(3),env_w_netcon(7)],...
    [env_w_netcon(6),env_w_netcon(9),env_w_netcon(8)],...
    [env_w_netcon(7),env_w_netcon(8),env_w_netcon(10)],...
    [env_w_netcon(11),env_w_netcon(4)],...
    [env_w_netcon(13),env_w_netcon(9),env_w_netcon(10)],...
    [env_w_netcon(5),env_w_netcon(12)],...
    [-3,env_w_netcon(12),env_w_netcon(13),env_w_netcon(11)]};

Z_netcon(netcon({[1,2,3,4],[5,4,6,2],[5,7,6,8],[1,8,3,7]},0,2,1,1)) = 1:8;
Z_order = {[Z_netcon(1),Z_netcon(2),Z_netcon(3),Z_netcon(4)],...
    [Z_netcon(5),Z_netcon(4),Z_netcon(6),Z_netcon(2)],...
    [Z_netcon(5),Z_netcon(7),Z_netcon(6),Z_netcon(8)],...
    [Z_netcon(1),Z_netcon(8),Z_netcon(3),Z_netcon(7)]};

M1_netcon(netcon({[2,1,-1,-5],[3,-7,-3,1],[3,-8,-4,5],[4,5],[2,4,-2,-6]},0,2,1,1)) = 1:5;
M1_order = {[M1_netcon(2),M1_netcon(1),-1,-5],...
    [M1_netcon(3),-7,-3,M1_netcon(1)],...
    [M1_netcon(3),-8,-4,M1_netcon(5)],...
    [M1_netcon(4),M1_netcon(5)],...
    [M1_netcon(2),M1_netcon(4),-2,-6]};

G_u_netcon(netcon({[-1,1,2],[-3,3,4],[2,5,3,6],[5,7,-2,1],[6,4,-4,7]},0,2,1,1)) = 1:7;
G_u_order = {[-1,G_u_netcon(1),G_u_netcon(2)],...
    [-3,G_u_netcon(3),G_u_netcon(4)],...
    [G_u_netcon(2),G_u_netcon(5),G_u_netcon(3),G_u_netcon(6)],...
    [G_u_netcon(5),G_u_netcon(7),-2,G_u_netcon(1)],...
    [G_u_netcon(6),G_u_netcon(4),-4,G_u_netcon(7)]};

G_y_netcon(netcon({[-1,1,2],[4,1,3],[-2,3,2],[5,-4,4],[-3,5]},0,2,1,1)) = 1:5;
G_y_order = {[-1,G_y_netcon(1),G_y_netcon(2)],...
    [G_y_netcon(4),G_y_netcon(1),G_y_netcon(3)],...
    [-2,G_y_netcon(3),G_y_netcon(2)],...
    [G_y_netcon(5),-4,G_y_netcon(4)],...
    [-3,G_y_netcon(5)]};

G_x_netcon(netcon({[-1,1,2],[-2,1,3],[4,3,2],[5,4,-4],[5,-3]},0,2,1,1)) = 1:5;
G_x_order = {[-1,G_x_netcon(1),G_x_netcon(2)],...
    [-2,G_x_netcon(1),G_x_netcon(3)],...
    [G_x_netcon(4),G_x_netcon(3),G_x_netcon(2)],...
    [G_x_netcon(5),G_x_netcon(4),-4],...
    [G_x_netcon(5),-3]};

M_netcon(netcon({[-1,1,-5,3],[1,4,2,5],[-3,2,-7,6],[3,7,11,9],[4,12,5,13,7,9,8,10],[8,6,10,14],[-2,15,-6,11],[15,12,16,13],[-4,16,-8,14]},0,2,1,1)) = 1:16;
M_order = {[-1,M_netcon(1),-5,M_netcon(3)],...
    [M_netcon(1),M_netcon(4),M_netcon(2),M_netcon(5)],...
    [-3,M_netcon(2),-7,M_netcon(6)],...
    [M_netcon(3),M_netcon(7),M_netcon(11),M_netcon(9)],...
    [M_netcon(4),M_netcon(12),M_netcon(5),M_netcon(13),M_netcon(7),M_netcon(9),M_netcon(8),M_netcon(10)],...
    [M_netcon(8),M_netcon(6),M_netcon(10),M_netcon(14)],...
    [-2,M_netcon(15),-6,M_netcon(11)],...
    [M_netcon(15),M_netcon(12),M_netcon(16),M_netcon(13)],...
    [-4,M_netcon(16),-8,M_netcon(14)]};

%observable - magnetism
pauli_z = [1 0;0 -1];
M = ncon({conj(A),conj(A),A,pauli_z,A},M1_order);

%calculate Z in the same way as the magnetisation to keep normalisation in check
MZ = ncon({conj(A),conj(A),A,[1 0;0 1],A},M1_order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:max_level
    
    %start with tensors of random numbers
    u = random('unif',-1,1,[min(chi_u,size(A,1)),size(A,1),min(chi_u,size(A,1)),size(A,1)]);
    w = random('unif',-1,1,[min(chi_w,size(u,3)*size(u,1)),size(u,3),size(u,1)]);
    v_l = random('unif',-1,1,[min(chi_v,size(A,4)*size(u,1)),size(A,4),size(u,1)]);
    v_r = random('unif',-1,1,[min(chi_v,size(A,2)*size(u,1)),size(u,3),size(A,2)]);
    y_l = random('unif',-1,1,[min(chi_y,size(v_l,1)*size(v_l,1)),size(v_l,1),size(v_l,1)]);
    y_r = random('unif',-1,1,[min(chi_y,size(v_r,1)*size(v_r,1)),size(v_r,1),size(v_r,1)]);

    %update v_l, v_r and u first
    for loop = 1:loop_max
        
        %store old B to check convergence
        if loop == 1
            B_old = ncon({v_l,u,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r)},B_order);
        else
            B_old = B./sqrt(ncon({B,conj(B)},{[1,2,3,4],[1,2,3,4]}));
        end
        
        %v_l - by unconstrained linear optimisation (see arXiv:0707.1454v2)
        for UL = 1:UL_max
            old = v_l;
            
            %build environment
            B = ncon({v_l,u,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r)},B_order);
            env_v_l = ncon({u,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r),conj(B)},env_v_l_order);
            
            %optimise by SVD
            [U,~,V] = svd(tfuse(env_v_l,[1,1,-2]),'econ');
            v_l = V*U';
            v_l = tsplit(v_l,2,[size(old,2),size(old,3)]);
            
            %check error - compare with old tensor using difference in the
            %Hilbert-Schmidt norm
            error = sqrt(ncon({v_l,conj(v_l)},{[1,2,3],[1,2,3]})...
                - ncon({v_l,conj(old)},{[1,2,3],[1,2,3]})...
                - ncon({old,conj(v_l)},{[1,2,3],[1,2,3]})...
                + ncon({old,conj(old)},{[1,2,3],[1,2,3]}));
            
            %if error is below the convergence criterion then break UL loop
            if error <= UL_convergence
%                 fprintf('Level: %d, v_l error = %.5e.\n',i,error);
                break
            elseif UL == UL_max
%                 fprintf('Level: %d, v_l error = %.5e. Not converged!\n',i,error);
            end
        end
        
        %v_r
        for UL = 1:UL_max
            old = v_r;
            
            %build environment
            B = ncon({v_l,u,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r)},B_order);
            env_v_r = ncon({v_l,u,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r),conj(B)},env_v_r_order);
            
            %optimise by SVD
            [U,~,V] = svd(tfuse(env_v_r,[1,1,-2]),'econ');
            v_r = V*U';
            v_r = tsplit(v_r,2,[size(old,2),size(old,3)]);
            
            %check error
            error = sqrt(ncon({v_r,conj(v_r)},{[1,2,3],[1,2,3]})...
                - ncon({v_r,conj(old)},{[1,2,3],[1,2,3]})...
                - ncon({old,conj(v_r)},{[1,2,3],[1,2,3]})...
                + ncon({old,conj(old)},{[1,2,3],[1,2,3]}));
            
            %if error is below the convergence criterion then break UL loop
            if error <= UL_convergence
%                 fprintf('Level: %d, v_r error = %.5e.\n',i,error);
                break
            elseif UL == UL_max
%                 fprintf('Level: %d, v_r error = %.5e. Not converged!\n',i,error);
            end
        end
        
        %u
        for UL = 1:UL_max
            old = u;
            
            %build environment
            B = ncon({v_l,u,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r)},B_order);
            env_u = ncon({v_l,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r),conj(B)},env_u_order);
            
            %optimise by SVD
            [U,~,V] = svd(tfuse(tfuse(env_u,[1,-2,1,-3]),[-1,2,2]),'econ');
            u = V*U';
            u = tsplit(u,1,[size(old,1),size(old,3)]);
            u = tsplit(u,3,[size(old,2),size(old,4)]);
            u = permute(u,[1,3,2,4]);
            
            %check error
            error = sqrt(ncon({u,conj(u)},{[1,2,3,4],[1,2,3,4]})...
                - ncon({u,conj(old)},{[1,2,3,4],[1,2,3,4]})...
                - ncon({old,conj(u)},{[1,2,3,4],[1,2,3,4]})...
                + ncon({old,conj(old)},{[1,2,3,4],[1,2,3,4]}));
            
            %if error is below the convergence criterion then break UL loop
            if error <= UL_convergence
%                 fprintf('Level: %d, u error = %.5e.\n',i,error);
                break
            elseif UL == UL_max
%                 fprintf('Level: %d, u error = %.5e. Not converged!\n',i,error);
            end
        end
        
        %update B
        B = ncon({v_l,u,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r)},B_order);
        
        %check error in B
        B_new = B./sqrt(ncon({B,conj(B)},{[1,2,3,4],[1,2,3,4]}));
        
        B_error = sqrt(ncon({B_new,conj(B_new)},{[1,2,3,4],[1,2,3,4]})...
            - ncon({B_new,conj(B_old)},{[1,2,3,4],[1,2,3,4]})...
            - ncon({B_old,conj(B_new)},{[1,2,3,4],[1,2,3,4]})...
            + ncon({B_old,conj(B_old)},{[1,2,3,4],[1,2,3,4]}));
        
        %if error is below the convergence criterion then break the loop
        if B_error <= B_convergence
            fprintf('Level: %d, B error = %.5e.\n',i,B_error);
            break
        elseif loop == loop_max
            fprintf('Level: %d, B error = %.5e. Not converged!\n',i,B_error);
        end
    end
    
    %second, update y_l and y_r, symmetrising before SVD
    for loop = 1:loop_max
        if loop == 1
            D_old = ncon({y_l,B,y_r},D_order);
        else
            D_old = D;
        end
        
        %y_l
        Ddagger = ncon({conj(y_l),conj(B),conj(y_r)},{[-1,1,2],[2,4,3,1],[-2,4,3]});
        env_y_l = ncon({B,y_r,Ddagger},env_y_l_order);
        
        %optimise by SVD, symmertrising for Hermiticity
        [U,~,V] = svd(tfuse(env_y_l,[1,1,-2])+tfuse(permute(conj(env_y_l),[2,1,3]),[1,1,-2]),'econ');
        y_l = V*U';
        y_l = tsplit(y_l,2,[size(B,4),size(B,1)]);
        
        %y_r
        Ddagger = ncon({conj(y_l),conj(B),conj(y_r)},{[-1,1,2],[2,4,3,1],[-2,4,3]});    
        env_y_r = ncon({y_l,B,Ddagger},env_y_r_order);
        
        %optimise by SVD, symmertrising for Hermiticity
        [U,~,V] = svd(tfuse(env_y_r,[1,1,-2])+tfuse(permute(conj(env_y_r),[2,1,3]),[1,1,-2]),'econ');
        y_r = V*U';
        y_r = tsplit(y_r,2,[size(B,2),size(B,3)]);
        
        %update D
        D = ncon({y_l,B,y_r},D_order);

        %check error
        D_error = sqrt(ncon({D,conj(D)},{[1,2],[1,2]})...
            - ncon({D,conj(D_old)},{[1,2],[1,2]})...
            - ncon({D_old,conj(D)},{[1,2],[1,2]})...
            + ncon({D_old,conj(D_old)},{[1,2],[1,2]}));
        
        %if error is below the convergence criterion then break the loop
        if D_error <= D_convergence
            fprintf('Level: %d, D error = %.5e.\n',i,D_error);
            break
        elseif loop == loop_max
            fprintf('Level: %d, D error = %.5e. Not converged!\n',i,D_error);
        end
    end
    
    %finally, update w
    for loop = 1:loop_max
        
        %build environment
        Aprime = ncon({w,conj(v_r),conj(v_l),sqrt(D),y_r,y_l,sqrt(D),v_r,v_l,conj(w)},A_order);
        Aprime = Aprime./sqrt(ncon({Aprime,conj(Aprime)},{[1,2,3,4],[1,2,3,4]}));
        A_old = Aprime;
        env_w = ncon({conj(v_r),conj(v_l),y_r,y_l,v_r,v_l,sqrt(D),conj(w),sqrt(D),conj(Aprime)},env_w_order);
        
        %optimise by SVD
        [U,S,V] = svd(tfuse(env_w,[1,1,-2]),'econ');
        w = V*U';
        w = tsplit(w,2,[size(u,3),size(u,1)]);

        %check error on A
        Aprime = ncon({w,conj(v_r),conj(v_l),sqrt(D),y_r,y_l,sqrt(D),v_r,v_l,conj(w)},A_order);
        Aprime = Aprime./sqrt(ncon({Aprime,conj(Aprime)},{[1,2,3,4],[1,2,3,4]}));
        
        A_error = sqrt(ncon({Aprime,conj(Aprime)},{[1,2,3,4],[1,2,3,4]})...
            - ncon({Aprime,conj(A_old)},{[1,2,3,4],[1,2,3,4]})...
            - ncon({A_old,conj(Aprime)},{[1,2,3,4],[1,2,3,4]})...
            + ncon({A_old,conj(A_old)},{[1,2,3,4],[1,2,3,4]}));

        %if error is below the convergence criterion then break the loop
        if A_error <= A_convergence
            fprintf('Level: %d, A error = %.5e.\n',i,A_error);
            break
        elseif loop == loop_max
            fprintf('Level: %d, A error = %.5e. Not converged!\n',i,A_error);
        end
    end
    
    %bring observable to next level
    G_l = ncon({v_l,conj(v_l)},{[-2,-1,1],[-4,-3,1]});
    G_r = ncon({v_r,conj(v_r)},{[-1,1,-2],[-3,1,-4]});
    G_u = ncon({v_l,v_r,u,A,A},G_u_order);
    G_y = ncon({w,conj(v_r),conj(v_l),y_r,sqrt(D)},G_y_order);
    G_x = ncon({w,conj(v_r),conj(v_l),y_l,sqrt(D)},G_x_order);
    M = ncon({G_y,G_u,G_x,G_r,M,G_l,conj(G_y),conj(G_u),conj(G_x)},M_order); 
    MZ = ncon({G_y,G_u,G_x,G_r,MZ,G_l,conj(G_y),conj(G_u),conj(G_x)},M_order); 
    
    %normalise
    normM = sqrt(ncon({MZ,conj(MZ)},{[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]}));
    MZ = MZ./normM;
    M = M./normM;
    
    %update A for next level
    A = ncon({w,conj(v_r),conj(v_l),sqrt(D),y_r,y_l,sqrt(D),v_r,v_l,conj(w)},A_order);
    A = A./sqrt(ncon({A,conj(A)},{[1,2,3,4],[1,2,3,4]}));
end

%full expectation value
M_exp = ncon({M},{[1,1,2,2,3,4,3,4]})/ncon({MZ},{[1,1,2,2,3,4,3,4]})











toc