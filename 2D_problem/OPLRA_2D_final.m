function [z,MAE,REG,Es,Fr1s,Fr2s,Br1r2,Wrf,Pred,Xrf1,Xrf2]=OPLRA_2D_final(fbest1,fbest2,exp,regions1,regions2,Asf,Asf1,Asf2,epsilon,U,U2,l)

features1=size(Asf1,2);
features2=size(Asf2,2);
features=size(Asf,2);
samples=size(Asf1,1);

Fr1s=binvar(regions1,samples);
Fr2s=binvar(regions2,samples);
Br1r2=sdpvar(regions1,regions2,'full');
Pred=sdpvar(regions1,regions2,samples,'full');
Xrf1=sdpvar(regions1,features1);
Xrf2=sdpvar(regions2,features2);
Ers1=sdpvar(regions1,samples);
Ers2=sdpvar(regions2,samples);
Es=sdpvar(samples,1);
Wrf=sdpvar(regions1,regions2,features,'full');
Wrfpos=sdpvar(regions1,regions2,features,'full');
MAE=sdpvar(1);
REG=sdpvar(1);

constraints=[];

constraints=[constraints,MAE==sum(Es)/samples];

for re1=1:regions1
    for re2=1:regions2
        for fe=1:features       
        constraints=[constraints,Wrfpos(re1,re2,fe)>=Wrf(re1,re2,fe)];
        constraints=[constraints,Wrfpos(re1,re2,fe)>=-Wrf(re1,re2,fe)];
        end
    end
end

constraints=[constraints, REG==sum(sum(sum(Wrfpos)))]; 

for re=2:regions1
    constraints=[constraints,Xrf1(re,fbest1)>=Xrf1(re-1,fbest1)+epsilon];
end 

constraints=[constraints,Xrf1(regions1,fbest1)==1];
constraints=[constraints,Xrf1(1,fbest1)>=epsilon];

for re=2:regions2
    constraints=[constraints,Xrf2(re,fbest2)>=Xrf2(re-1,fbest2)+epsilon];
end 

constraints=[constraints,Xrf2(regions2,fbest2)==1];
constraints=[constraints,Xrf2(1,fbest2)>=epsilon];

for i=1:samples
    constraints=[constraints,sum(Fr1s(:,i))==1];
end 

for i=1:samples
    constraints=[constraints,sum(Fr2s(:,i))==1];
end

for re=2:regions1
    for i=1:samples
        constraints=[constraints,Xrf1(re-1,fbest1)+epsilon-U*(1-Fr1s(re,i))<= Asf1(i,fbest1)];
    end
end 

for re=1:regions1-1 
    for i=1:samples
        constraints=[constraints,Asf1(i,fbest1)<=Xrf1(re,fbest1)-epsilon+U*(1-Fr1s(re,i))];
    end 
end 

for re=2:regions2
    for i=1:samples
        constraints=[constraints,Xrf2(re-1,fbest2)+epsilon-U*(1-Fr2s(re,i))<= Asf2(i,fbest2)];
    end
end 

for re=1:regions2-1 
    for i=1:samples
        constraints=[constraints,Asf2(i,fbest2)<=Xrf2(re,fbest2)-epsilon+U*(1-Fr2s(re,i))];
    end 
end 

for re1=1:regions1
    for re2=1:regions2
        for i=1:samples
         constraints=[constraints,Pred(re1,re2,i)==sum(Asf(i,:).*Wrf(re1,re2,:))+Br1r2(re1,re2)];
        end 
    end 
end  


for re1=1:regions1
   for re2=1:regions2
     for i=1:samples
            constraints=[constraints,Es(i,1)>=0];
            %constraints=[constraints,Ers1(re1,i)<=U2*(Fr1s(re1,i))];
            %constraints=[constraints,Ers2(re2,i)<=U2*(Fr2s(re2,i))];
            %constraints=[constraints,Ers1(re1,i)>=exp(i,1)-Pred(re1,re2,i)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            %constraints=[constraints,Ers1(re1,i)>=Pred(re1,re2,i)-exp(i,1)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            %constraints=[constraints,Ers2(re2,i)>=exp(i,1)-Pred(re1,re2,i)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            %constraints=[constraints,Ers2(re2,i)>=Pred(re1,re2,i)-exp(i,1)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            %constraints=[constraints,Es(i,1)>=Ers1(re1,i)];
            %constraints=[constraints,Es(i,1)>=Ers2(re2,i)];  
            constraints=[constraints,Es(i,1)>=exp(i,1)-Pred(re1,re2,i)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            constraints=[constraints,Es(i,1)>=Pred(re1,re2,i)-exp(i,1)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
       
      end
   end
end

z=MAE+l*REG;

options=sdpsettings('solver','mosek'); %,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',50,'mosek.MSK_DPAR_MIO_MAX_TIME',50);
%options = sdpsettings('solver', 'gurobi');

readAcross=optimize(constraints,z,options);

z=value(z);
MAE=value(MAE);
REG=value(REG);
Es=value(Es); 
Pred=value(Pred);
Xrf1=value(Xrf1);
Xrf2=value(Xrf2);
Wrf=value(Wrf);
Br1r2=value(Br1r2);
Fr1s=value(Fr1s);
Fr2s=value(Fr2s);



