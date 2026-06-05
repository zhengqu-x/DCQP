function solve_synthetic_with_dcqp_eliminate(group_name,instance_id)


valid_groups = {'qp_n_0_1', 'qp_n_0_3', 'qp_n_0_9', 'qp_u_0_1', 'qp_u_0_3', 'qp_u_0_9', 'qp_u_25_1'};
if ~ismember(group_name, valid_groups)
    error('Invalid group_name. Must be one of: %s', strjoin(valid_groups, ', '));
end




diary diaryfile-synthetic-dcqp-eliminate.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



datalocation='../data/synthetic/';
myrecord=zeros(20,6);

fprintf("**************************************start to solve ")
fprintf(group_name);
fprintf("****************************\n\n\n\n");

% Check if instance_id is provided
if nargin < 2
    K = 1:20;
else
    K = instance_id;
end

for i=1:length(K)
    k=K(i);


    fprintf("************************************QP Instance %d ****************************\n",k);


    filename=[group_name '_' num2str(k)];

    datafilenamek=[datalocation filename '.mat'];
    filename=[filename '_eliminate'];
    F=load(datafilenamek);

    H=F.H;
    f=F.f;
    A=F.A;
    b=F.b;
    Aeq=F.Aeq;
    beq=F.beq;
    UB=F.UB;
    LB=F.LB;
    if norm(Aeq)==0
        Aeq=[];
        beq=[];
    end

    n = size(A,2);

    Q = H*0.5;
    d = f*0.5;


    A_bar = [A; eye(n); -eye(n)];
    b_bar = [b; UB; -LB];

  
    params = dcqp_default_params();
    params.nb_rounds=10;
    params.filename=filename;
    params.do_scaling=true;
    params.verbose=true;
    %params.eta=0.1;
    %params.mosek_tolerance = 1e-10;

    fprintf('eta=%.2e\n',params.eta);
    
    [Q_tilde, d_tilde, A_tilde, b_tilde, M, x0,obj_constant] = eliminate_equalities(Q, d, A_bar, b_bar, Aeq, beq,n);


    fprintf('obj_constant=%.5f\n',obj_constant);


    [x_opt, fval, info] = dcqp_solve(Q_tilde, d_tilde, A_tilde, b_tilde, [], [], params);
    
    

    %[x_opt, fval, info] = dcqp_solve(Q, d, A_bar, b_bar, Aeq, beq, params);
     
    if info.status=="not solved"
       params.eta=0.8;
       params.gap_tolerance = 1e-5;  
       [x_opt, fval, info] = dcqp_solve(Q_tilde, d_tilde, A_tilde, b_tilde, [], [], params);
       
    end

    x_real=M * x_opt + x0;
    f_real=obj_constant+fval;
    fprintf('optimal value=%.6f\n', f_real);
    myrecord(k,1)=info.gap;
    myrecord(k,2)=max(A_bar*x_real-b_bar);
    if ~isempty(beq)
        myrecord(k,3)=norm(Aeq*x_real-beq);
    end
    myrecord(k,4)=f_real;
    myrecord(k,5)=info.lower_bound;
    myrecord(k,6)=info.time;
end

if length(K)==20 
    save(['summary_results/' group_name '_eliminate.mat'],"myrecord");
end

end
