

function solve_synthetic_with_gurobi(group_name,instance_id,Tlimit)


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

% Check if instance_id and time limit is provided
if nargin < 2
    K = 1:20;
    Tlimit=3600;
elseif nargin<3
    K = instance_id;
    Tlimit=3600;
else
    K=instance_id;
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


    tol_gqp=1e-4;
    maintimer=tic;
    
    [result] = gurobiqp_with_time_limit(Q,d,A_bar,b_bar,Aeq,beq,tol_gqp,n,Tlimit);

    time_record=toc(maintimer);
    
    bestub=result.objval;
    bestlb=result.objbound;
    bestsol=result.x;

     % Compute final relative gap
    if abs(bestub) > 1e-8
        gap = (bestub - bestlb) / abs(bestub);
    else
        gap = bestub-bestlb;
    end
    
    myrecord(k,1)=gap;
    myrecord(k,2)=max(A_bar*bestsol-b_bar);
    if ~isempty(beq)
        myrecord(k,3)=norm(Aeq*bestsol-beq);
    end
    myrecord(k,4)=bestub;
    myrecord(k,5)=bestlb;
    myrecord(k,6)=time_record;

    resultsfolder = fullfile(pwd, 'testresults');
    if ~exist(resultsfolder, 'dir')
        mkdir(resultsfolder);
    end
    timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
    file_name=fullfile(resultsfolder, sprintf('gurobi_%s-%s.mat', filename,  timestamp));
    
    info=myrecord(k,:);
    save(char(file_name),"bestsol","info");
end

if length(K)==20 
    save(['summary_results/' 'gurobi_' group_name '.mat'],"myrecord");
end

end
