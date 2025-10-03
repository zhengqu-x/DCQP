n=100;

m=50;
meq=0;

s=50;
hs=floor(s/2);

sparsity=0.9;


x0=rand(n,1);
x0=x0/sum(x0);


A=[ones(1,n);sprandn(m,n,sparsity)];
A=full(A);
b=A*x0+0.1*rand(m+1,1);



if meq>0
Aeq=sprandn(meq,n,sparsity);
Aeq=full(Aeq);
beq=Aeq*x0;
else
    Aeq=[];
    beq=[];
end



for i=1:20
    L=sprandn(s,n,sparsity);
    L_1 = L(1:hs,:);
    L_2 = L(hs+1:s,:);
    Q = L_1'*L_1 - L_2'*L_2;
    Q=full(Q);
    d = 0.01*randn(n,1);


    LB=zeros(n,1);
    UB=ones(n,1);
    H=2*Q;
    f=2*d;

    filename=['../data/testrand/normalrandomqp' num2str(n) '_' num2str(m) '_' num2str(meq) '_' num2str(s) '_' num2str(sparsity*10) '_' num2str(i) ];

    save(filename,"LB","UB","H","f","A","b","Aeq","beq");

end


