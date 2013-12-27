function [f_app_Train_Data, f_app_Test_Data] =LP_iter_chm(Train_Data, Test_Data, f, my_eps)


Sz_Train = size(Train_Data,1);
Sz_Test = size(Test_Data,1);
f_SZ = size(f,2);

dist1 = squareform(pdist([Train_Data;Test_Data]));

D = (dist1.^2)./my_eps; 
ker = exp(-D);


ker = ker(:,1:Sz_Train);


fz = find(ker < 10.^-10);
ker(fz)=0;



for i=1:size(ker,1)
    if(sum(ker(i,:)) > 0)
      ker(i,:) = ker(i,:) / sum(ker(i,:));
    end;
end;




f_app_Data = zeros(Sz_Train+Sz_Test,f_SZ);

for i=1:Sz_Train+Sz_Test
    for j=1:Sz_Train
         f_app_Data(i,:) = f_app_Data(i,:) + ker(i,j)*f(j,:);
    end;
end;



f_app_Train_Data = f_app_Data(1:Sz_Train,:);
f_app_Test_Data = f_app_Data(Sz_Train+1:end,:);