function solve_benchmark_with_dcqp(group_name)


valid_groups = {'qp20_10', 'qp30_15', 'qp40_20', 'qp50_25'};
if ~ismember(group_name, valid_groups)
    error('Invalid group_name. Must be one of: %s', strjoin(valid_groups, ', '));
end

diary diaryfile-benchmark-dcqp.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


datalocation='../data/benchmark/';
myrecord=zeros(16,6);

fprintf("**************************************start to solve ")
fprintf(group_name);
fprintf("****************************\n\n\n\n");


i=0;

for k=1:4

    for k2=1:4
    

    fprintf("************************************QP Instance %d - %d ****************************\n",k,k2);



    filename=[group_name '_' num2str(k) '_' num2str(k2)];

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


    n = size(A,2);

    Q = H*0.5;
    d = f*0.5;


    A_bar = [A; eye(n); -eye(n)];
    b_bar = [b; UB; -LB];


    params = dcqp_default_params();
    params.nb_rounds=1;
    params.filename=filename;
    %params.verbose=true;
    

    [x_opt, fval, info] = dcqp_solve(Q, d, A_bar, b_bar, Aeq, beq, params);
 
    i=i+1;
    myrecord(i,1)=info.gap;
    myrecord(i,2)=max(A_bar*x_opt-b_bar);
    if ~isempty(beq)
        myrecord(i,3)=norm(Aeq*x_opt-beq);
    end
    myrecord(i,4)=fval;
    myrecord(i,5)=info.lower_bound;
    myrecord(i,6)=info.time;


     end
end

save(['summary_results/' group_name '.mat'],"myrecord");

end