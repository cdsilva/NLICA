



function [f_app_Train_Data, f_app_Test_Data, deltas, my_err, my_err_max, correct_app] = LP_chm(Train_Data, Test_Data, f, my_err_thresh,  my_eps)


SzTR = size(Train_Data,1);
SzTS = size(Test_Data,1);

my_k = size(f,2);

F_TRAIN = f(1:SzTR,:);
F_TO_APP = F_TRAIN;
F_KER = zeros(SzTR+SzTS,1);



[APP_RES_TR, APP_RES_TS] = LP_iter_chm(Train_Data, Test_Data, F_TO_APP, my_eps);


%%
MyErrZeroInd = zeros(SzTR,1);
correct_app = zeros(SzTR,1);
for ER=1:SzTR
    MyErr(ER,1) = sum(abs(F_TO_APP(ER,:) - APP_RES_TR(ER,:)).^2)./my_k;
    if(MyErr(ER,1) < my_err_thresh)
        correct_app(ER) = 1;
    end;
end;

my_err(:,1) = mean(MyErr(:,1));
my_err_max(:,1) = max(MyErr(:,1));


f_app_Train_Data{1} = APP_RES_TR;
f_app_Test_Data{1} = APP_RES_TS;


F_TO_APP = F_TO_APP-APP_RES_TR;
F_KER = [f_app_Train_Data{1};f_app_Test_Data{1}];

deltas{1} = F_TO_APP;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


i=1;


while(my_err_max(:,i) > my_err_thresh) 

    
 i = i+1;    
 

my_eps = my_eps./2;

[APP_RES_TR, APP_RES_TS] = LP_iter_chm(Train_Data, Test_Data, F_TO_APP, my_eps);



for ER=1:SzTR
    MyErr(ER,i) = sum(abs(F_TO_APP(ER,:) - APP_RES_TR(ER,:)).^2)./my_k;
    if(MyErr(ER,i) < my_err_thresh)
        if(correct_app(ER)==0)
            correct_app(ER) = i;
        end;
    end;
end;

my_err(:,i) = mean(MyErr(:,i));
my_err_max(:,i) = max(MyErr(:,i));



f_app_Train_Data{i}= f_app_Train_Data{i-1} + APP_RES_TR;
f_app_Test_Data{i} = f_app_Test_Data{i-1} + APP_RES_TS;


F_TO_APP = F_TO_APP-APP_RES_TR;


deltas{i} = F_TO_APP;

F_KER = [APP_RES_TR;APP_RES_TS];



if(i>50)
    if( mean(my_err(i))./mean(my_err(i-1)) > 0.5)
        n=9;
        break;
    end;
end;

end;

