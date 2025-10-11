function [x_kkt] = search_of_kkt_point( ...
    Q,d,A,b,Aeq,beq,M,N,x_0,epsilon,met_gqp,tol_gqp,n)

% ======================================================================= %
% Search for x_kkt satisfying x_kkt is a KKT point, 
% and either x_kkt is a vertex of {x: Ax <= b}, or Q|_H_I_x_kkt > 0
%
% INPUT
% Q, d             Parameters in the objective function
% A, b             Parameters in constraint Ax <= b
% Aeq, beq         Parameters in constraint Aeq x = beq
% M, N             Positive semidefinite matrices M and N such that Q = M-N
% x_0              Initial point
% epsilon          Numerical precision control parameters (epsilon > 0)
% met_gqp          Method used in gurobiqp
%                      0: primal simplex
%                      1: dual simplex
%                      2: barrier
% tol_gqp          Accuracy tolerance used in gurobiqp
% n                Size(A,2)
%
% OUTPUT             
% x_kkt            x_kkt
% ======================================================================= %

x_ite = x_0;

 iteration_number = 0; % Iteration number

         
while iteration_number<=1000
     iteration_number = iteration_number +1;
   
     x_hat = x_ite;
     [x_hat_vertex_spd_flag, x_hat_vertex_flag, I_x_hat_c,B_A_I,BH] = check_vertex_spd(Q, A, b, Aeq,x_hat, epsilon);
     if x_hat_vertex_spd_flag
         v_hat=x_hat'*Q*x_hat+2*d'*x_hat;
        if x_hat_vertex_flag == 1  % x_hat is a vertex
            x_tilde = x_hat;
        else
            [y_tilde_I,~] = gurobiqp(BH, B_A_I'*(Q*x_hat+d), A(I_x_hat_c,:)*B_A_I, (b(I_x_hat_c) - A(I_x_hat_c,:)*x_hat),[],[],met_gqp,tol_gqp,size(B_A_I,2));
            x_tilde = x_hat + B_A_I * y_tilde_I;
        end
        v_tilde=x_tilde'*Q*x_tilde+2*d'*x_tilde;
        if abs(v_tilde-v_hat) <= max(1e-8*abs(v_hat),1e-9)
                [x_bar] = gurobiqp(0.5*M,0.5*(d - N*x_hat),A,b,Aeq,beq,met_gqp,tol_gqp,n);
                v_bar=x_bar'*Q*x_bar+2*d'*x_bar;      
                if v_hat-v_bar <= max(1e-8*abs(v_hat),1e-9) 
                    x_kkt=x_hat;
                    break
                end
           
        else
            [x_bar] = gurobiqp(0.5*M, 0.5*(d - N*x_tilde),A,b,Aeq,beq,met_gqp,tol_gqp,n);
            v_bar=x_bar'*Q*x_bar+2*d'*x_bar;
            if v_bar >= v_tilde-max(1e-8*abs(v_tilde),1e-9) 
                x_bar = x_tilde;
            end
           
        end
     else
        [x_bar] = down_hill(Q,d,B_A_I,A,b,x_hat,I_x_hat_c);
    end   
    x_ite = x_bar;
  
end

end