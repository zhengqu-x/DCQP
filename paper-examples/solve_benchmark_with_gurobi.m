function solve_benchmark_with_gurobi(group_name,Tlimit)


valid_groups = {'qp20_10', 'qp30_15', 'qp40_20', 'qp50_25'};
if ~ismember(group_name, valid_groups)
    error('Invalid group_name. Must be one of: %s', strjoin(valid_groups, ', '));
end

diary diaryfile-benchmark-gurobi.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


datalocation='../data/benchmark/';
myrecord=zeros(16,6);


fprintf("**************************************start to solve ")
fprintf(group_name);
fprintf("****************************\n\n\n\n");



% Check if time limit is provided
if nargin < 2
    Tlimit=3600;
end


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


        tol_gqp=1e-4;
        maintimer=tic;

        [result] = gurobiqp_with_time_limit(Q,d,A_bar,b_bar,Aeq,beq,tol_gqp,n,Tlimit);

        time_record=toc(maintimer);

        bestub=result.objval;
        bestlb=result.objbound;
        bestsol=result.x;
        i=i+1;


        % Compute final relative gap
        if abs(bestub) > 1e-8
            gap = (bestub - bestlb) / abs(bestub);
        else
            gap = bestub-bestlb;
        end
        myrecord(i,1)=gap;
        myrecord(i,2)=max(A_bar*bestsol-b_bar);
        if ~isempty(beq)
            myrecord(i,3)=norm(Aeq*bestsol-beq);
        end
        myrecord(i,4)=bestub;
        myrecord(i,5)=bestlb;
        myrecord(i,6)=time_record;



        resultsfolder = fullfile(pwd, 'testresults');
        if ~exist(resultsfolder, 'dir')
            mkdir(resultsfolder);
        end
        timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
        file_name=fullfile(resultsfolder, sprintf('gurobi_%s-%s.mat', filename,  timestamp));

        info=myrecord(k,:);
        save(char(file_name),"bestsol","info");

    end
end

save(['summary_results/' 'gurobi_' group_name '.mat'],"myrecord");

end
