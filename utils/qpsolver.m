function [bestub,best_sol,bestlb,nb_iters]=qpsolver(Q,d,A,b,Aeq,beq,lb,ub,sol,parameters)

%==========================================================================%
% Solve the following quadratic program:
% min  x'*Q*x + 2*d'*x
% s.t. Ax <= b
%      Aeq x = beq
%
% INPUT
% Q, d             Parameters of the quadratic objective function
% A, b             Parameters of the inequality constraint A x <= b
% Aeq, beq         Equality constraints: Aeq x = beq
% lb, ub           Known lower and upper bounds of the objective value
% sol              Initial feasible solution 
% parameters       Structure containing algorithm parameters:
%
% OUTPUT
% bestub           Best upper bound of the objective value found
% best_sol         Corresponding solution vector achieving bestub
% bestlb           Best lower bound of the objective value found
%==========================================================================

met_glp = parameters.metglp;
met_gqp = parameters.metgqp;
max_N= parameters.max_N;

spn = parameters.eps_dcdecomposition;
eps_checkpsd = parameters.eps_checkpsd;
gap_tol =parameters.gap_tol;
tol_mosek = parameters.tol_mosek;
tol_glp = parameters.tol_glp;
tol_gqp = parameters.tol_gqp;
eta=parameters.eta;


bestub= ub;
bestlb= lb;


lb_record=zeros(max_N,1);
ub_record=zeros(max_N,1);
bestub_record=zeros(max_N,1);
bestlb_record=zeros(max_N,1);
cut_lb_record=zeros(max_N,1);
alpha_record=zeros(max_N,1);

time_record_lb=zeros(max_N,1);
time_record_cut_lb=zeros(max_N,1);
time_record_kkt=zeros(max_N,1);
time_record_generate_cut=zeros(max_N,1);


[M,N] = DC_decomposition(Q,spn);

A_bar=A;
b_bar=b;



cut_val=Inf;

n=size(A_bar,2);
[~,tstar0,~]=gurobilp(-ones(n,1),A_bar,b_bar,Aeq,beq,[],[],met_glp,tol_glp);

if abs(tstar0)>1
    tstar0=ceil(abs(tstar0));
    A_bar=A_bar*abs(tstar0);
    Aeq=Aeq*abs(tstar0);
    Q=Q*abs(tstar0^2);
    M=M*abs(tstar0^2);
    N=N*abs(tstar0^2);
    d=d*abs(tstar0);
end
best_sol=sol/tstar0;

for i=1:max_N

    m=size(A_bar,1);
    n=size(A_bar,2);

    tic

    [~,tstar,exitflag]=gurobilp(-ones(n,1),A_bar,b_bar,Aeq,beq,[],[],met_glp,tol_glp);
    if exitflag==-1
        if i==1
            error('The feasible region is empty.\n');
        else
            lb_record(i)=min(cut_lb_record(1:i-1));
            break
        end
    end

    if cut_val<1 || mod(i-1,10)==0
        [lb,sdp_status,S,res]=lower_bound_dnn(Q,d,A_bar,b_bar,Aeq,beq,tol_mosek,m,n);
        if sdp_status==1
            delta=min(eig(S));
             if parameters.verbose==true
                fprintf('compute lower bound, minimal eigenvalue of S=%4.2e\n\n\n',delta);
             end
            lb_record(i)=lb+(1+tstar^2)*min(delta,0);
            u=res.sol.itr.doty;
            U=sMat(u,n+1);
            x_0=U(end,1:end-1)';
        else
            error('DNN lower bound was not solved successfully. Consider increase the mosek tolerance and rety.');
        end
    else
        lb_record(i)=lb_record(i-1);
    end


    bestlb=max(bestlb,lb_record(i));
    bestlb_record(i)=bestlb;

    x0=x_0;
    if max(A_bar*x_0-b_bar) >1e-9 || (~isempty(beq) &&norm(Aeq*x_0-beq)>1e-9)
        x0 = gurobiqp(eye(n),-x_0,A_bar,b_bar,Aeq,beq,met_gqp,tol_gqp,n);
    end

    time_record_lb(i)=toc;

    tic

    save('badexample.mat',"Q","d","A_bar","b_bar","Aeq","beq","M","N","x0","eps_checkpsd","met_gqp","tol_gqp","n");

    [x_kkt] = search_of_kkt_point(Q,d,A_bar,b_bar,Aeq,beq,M,N,x0,eps_checkpsd,met_gqp,tol_gqp,n);

    time_record_kkt(i)=toc;

    barx=x_kkt;


    v=barx'*Q*barx+2*d'*barx;


    if v<bestub
        best_sol=barx;
    end
    bestub=min(v,bestub);
    ub_record(i)=v;
    bestub_record(i)=bestub;


    tic;
    if bestub<=bestlb+abs(bestub)*gap_tol
        break
    end

    nuR=bestub-abs(bestub)*gap_tol;

    nuR2=bestub-eta*abs(bestub)*gap_tol;

    beta=abs(bestub)*gap_tol*0.01;

    if abs(bestub)<gap_tol
       nuR=bestub-gap_tol;
       nuR2=bestub-0.9*gap_tol; 
       beta=min(1e-7,0.01*gap_tol);
    elseif v>bestub
        nuR2=0.99*bestub+0.01*v;
        beta=min(1e-6,0.1*(v-bestub));
    end

 
    tol_mosek_cut=parameters.tol_mosek_cut;
    [c,cut_val,~,S] = generate_cut_dnn(Q,d,A_bar,b_bar,Aeq,beq,m,n,nuR2,barx,x_0,tol_mosek_cut,beta);

    while tol_mosek_cut<=1e-5 && isempty(c)
        tol_mosek_cut=tol_mosek_cut*10;
        fprintf("reducing tol_mosek_cut to %4.2e. \n", tol_mosek_cut);
        [c,cut_val,~,S] = generate_cut_dnn(Q,d,A_bar,b_bar,Aeq,beq,m,n,nuR2,barx,x_0,tol_mosek_cut,beta);    
    end
    if isempty(c)
        error('fail to generate cut with mosek tolerance %4.2e! try to decrease the error tolerance epsilon\n\n\n', tol_mosek_cut)
    else
        delta=min(eig(S));
        if parameters.verbose==true
            fprintf('generate cut with mosek tolerance %4.2e: minimal eigenvalue of S=%4.2e\n\n\n',tol_mosek_cut,delta);
        end
    end


    time_record_generate_cut(i)=toc;

    tic

    A_cut=[A_bar;c'/norm(c)];
    b_cut=[b_bar; (1+c'*barx)/norm(c)];

    [~,tstar_cut,~]=gurobilp(-ones(n,1),A_cut,b_cut,Aeq,beq,[],[],met_glp,tol_glp);

    cut_lb=nuR2+(1+tstar_cut^2)*min(delta,0);

    cut_lb_record(i)=cut_lb;

    if cut_lb>=nuR
        A_bar=[A_bar;-c'/norm(c)];
        b_bar=[b_bar; (-1-c'*barx)/norm(c)];
    else
        alpha=1;
        A_cut=[A_bar;c'/norm(c)];
        b_cut=[b_bar; (alpha+c'*barx)/norm(c)];
        [lb_cut,~,S,~]=lower_bound_dnn(Q,d,A_cut,b_cut,Aeq,beq,tol_mosek,m+1,n);
        delta_cut=min(eig(S));
        if parameters.verbose==true
            fprintf('cut lower bound computing: minimal eigenvalue of S=%4.2e\n\n\n',delta_cut);
        end
        lb_cut2=lb_cut+(1+tstar_cut^2)*min(delta_cut,0);
        if lb_cut2<nuR
            alpha=0.9;
            A_cut=[A_bar;c'/norm(c)];
            b_cut=[b_bar; (alpha+c'*barx)/norm(c)];
            [lb_cut,~,S,~]=lower_bound_dnn(Q,d,A_cut,b_cut,Aeq,beq,tol_mosek,m+1,n);
            delta_cut=min(eig(S));
            if parameters.verbose==true
                fprintf('reduced alpha, cut lower bound recomputing: minimal eigenvalue of S=%4.2e\n\n\n',delta_cut);
            end
            lb_cut2=lb_cut+(1+tstar_cut^2)*min(delta_cut,0);
        end
        cut_lb_record(i)=lb_cut2;
        alpha_record(i)=alpha;
        A_bar=[A_bar;-c'/norm(c)];
        b_bar=[b_bar; (-alpha-c'*barx)/norm(c)];
    end

    
    if parameters.verbose==true
    fprintf('iteration %5d: bestub=%4.8f, nuR2=%4.8f, nuR=%4.8f, delta=%4.8f,current v=%4.8f, cut_lb=%4.10f, best_lb=%4.10f\n\n\n', i, bestub,nuR2,nuR,delta,v,cut_lb_record(i),bestlb);
    end
    time_record_cut_lb(i)=toc;

end

best_sol=best_sol*tstar0;
if i>1
    bestlb=min(min(cut_lb_record(1:i-1)),bestlb);
end
nb_iters=i;
bestub=bestub/parameters.scaling;
bestlb=bestlb/parameters.scaling;
lb_record=lb_record/parameters.scaling;
ub_record=ub_record/parameters.scaling;
bestub_record=bestub_record/parameters.scaling;
bestlb_record=bestlb_record/parameters.scaling;
cut_lb_record=cut_lb_record/parameters.scaling;

lb_record=lb_record(1:i);
ub_record=ub_record(1:i);
bestub_record=bestub_record(1:i);
bestlb_record=bestlb_record(1:i);
cut_lb_record=cut_lb_record(1:i-1);
alpha_record=alpha_record(1:i);
time_record_lb=time_record_lb(1:i);
time_record_cut_lb=time_record_cut_lb(1:i);
time_record_kkt=time_record_kkt(1:i);
time_record_generate_cut=time_record_generate_cut(1:i);


resultsfolder = fullfile(pwd, 'testresults');
if ~exist(resultsfolder, 'dir')
        mkdir(resultsfolder);
end

timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
filename=fullfile(resultsfolder, sprintf('%s-%s.mat', parameters.filename,  timestamp));


save(char(filename),"best_sol","alpha_record","lb_record","time_record_generate_cut","time_record_kkt","time_record_cut_lb","cut_lb_record","time_record_lb","bestlb_record","bestub_record","ub_record")

fprintf("************************************ End of Computation  ****************************\n\n\n\n");




end