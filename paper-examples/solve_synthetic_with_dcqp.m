function solve_synthetic_with_dcqp(group_name,instance_id)


valid_groups = {'qp_n_0_1', 'qp_n_0_3', 'qp_n_0_9', 'qp_u_0_1', 'qp_u_0_3', 'qp_u_0_9', 'qp_u_25_1'};
if ~ismember(group_name, valid_groups)
    error('Invalid group_name. Must be one of: %s', strjoin(valid_groups, ', '));
end




diary diaryfile-synthetic-dcqp.txt

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
    params.nb_rounds=100;
    params.filename=filename;
    params.do_scaling=true;
    %params.verbose=true;
    

    [x_opt, fval, info] = dcqp_solve(Q, d, A_bar, b_bar, Aeq, beq, params);
     
    if info.status=="not solved"
       params.eta=0.8;
       [x_opt, fval, info] = dcqp_solve(Q, d, A_bar, b_bar, Aeq, beq, params);
     
    end
    myrecord(k,1)=info.gap;
    myrecord(k,2)=max(A_bar*x_opt-b_bar);
    if ~isempty(beq)
        myrecord(k,3)=norm(Aeq*x_opt-beq);
    end
    myrecord(k,4)=fval;
    myrecord(k,5)=info.lower_bound;
    myrecord(k,6)=info.time;
end

if length(K)==20 
    save(['summary_results/' group_name '.mat'],"myrecord");
end

end