%from all seizure files,extract an interval 5 sec before and after seizure
%prepare the seizure data from the training data

clear all;
clc;
new_var=[];
%Seizure Data files
sfiles=['chb01_03.edf';'chb01_04.edf';'chb01_15.edf';'chb01_16.edf';'chb01_18.edf';'chb01_21.edf';'chb01_26.edf'];
%Seizure Time Intervals
stime=[2996 3036;1467 1494;1732 1772;1015 1066;1720 1810;327 420;1862 1963];
SeizureData=[];
NonSeizureData=[];

for i=1:1   %For seizure data files
    [sinfo,sdata]=edfread(sfiles(i,:)); % read File No. i in sdata
    i
    [r,c]=size(sdata);
    for j=1:r
        for k=(stime(i,1)-300):(stime(i,2)+300) % 5 sec before and after seizure
            data=sdata(j,(k*256:((k+1)*256-1))); %1 sec interval
            % calculate the d3 a4 d4 wavelet coefficients
            [a1,~] = dwt(data,'db1'); 
            [a2,~] = dwt(a1,'db1');
            [a3,d3]= dwt(a2,'db1');
            [a4,d4]= dwt(a3,'db1');
            %calculate mean,var,energy,entropy,max,min for each
            vect1=[];
            vect2=[];
            vect3=[];
            vect1=[mean(d3);var(d3);d3*d3';(log(eps+d3.^2))*d3';max(d3);min(d3)];
            vect2=[mean(d4);var(d4);d4*d4';(log(eps+d4.^2))*d4';max(d4);min(d4)];
            vect3=[mean(a4);var(a4);a4*a4';(log(eps+a4.^2))*a4';max(a4);min(a4)];
            temp=[vect1;vect2;vect3;iqr(data);mad(data)];
            if(k<=stime(i,2) && k>=stime(i,1)) % if data belongs to Seizure set
                SeizureData=[SeizureData,temp];
            else
                NonSeizureData=[NonSeizureData,temp];
            end
            %each Seizure set contains 20 feature of interval
        end
    end
end
SeizureData=SeizureData';
NonSeizureData=NonSeizureData';

for i=2:length(stime)   %For seizure data files
    SeizureData=SeizureData';
    NonSeizureData=NonSeizureData';
    [sinfo,sdata]=edfread(sfiles(i,:)); % read File No. i in sdata
    i
    [r,c]=size(sdata);
    for j=1:r
        for k=(stime(i,1)-300):(stime(i,2)+300) % 5 sec before and after seizure
            data=sdata(j,(k*256:((k+1)*256-1))); %1 sec interval
            % calculate the d3 a4 d4 wavelet coefficients
            [a1,~] = dwt(data,'db1'); 
            [a2,~] = dwt(a1,'db1');
            [a3,d3]= dwt(a2,'db1');
            [a4,d4]= dwt(a3,'db1');
            %calculate mean,var,energy,entropy,max,min for each
            vect1=[];
            vect2=[];
            vect3=[];
            vect1=[mean(d3);var(d3);d3*d3';(log(eps+d3.^2))*d3';max(d3);min(d3)];
            vect2=[mean(d4);var(d4);d4*d4';(log(eps+d4.^2))*d4';max(d4);min(d4)];
            vect3=[mean(a4);var(a4);a4*a4';(log(eps+a4.^2))*a4';max(a4);min(a4)];
            temp=[vect1;vect2;vect3;iqr(data);mad(data)];
            if(k<=stime(i,2) && k>=stime(i,1)) % if data belongs to Seizure set
                SeizureData=[SeizureData,temp];
            else
                NonSeizureData=[NonSeizureData,temp];
            end
            %each Seizure set contains 20 feature of interval
        end
    end
SeizureData=SeizureData';
NonSeizureData=NonSeizureData';
    
    %Steps to b followed
%train the svm with whole seizure data
%and same length of non-seizure data
%find the classifier using radial basis function
%100% accuracy for whole Non-Seizure Files 
incorr_s=[];%incorr in seiizure files
mat1=[];mat2=[];incorr_ns=[];%incorr in non seiizure files
cnt1=0;cnt2=0;
for ch=1:23              %for channel attached to a patient
    s=[];ns=[];cl=[];
    mat1=[];mat2=[];
    %indic=randi(length(SeizureData),1,60);
    indic=1:length(SeizureData);
    s=SeizureData(indic,:);
    ns=NonSeizureData(indic,:);
    var1=[ones(length(s),1);zeros(length(ns),1)];
    clear cl;
    cl= svmtrain( var1 , [s;ns] ,'-t 0'); %svmtrain(label,data,'svmoptions')
    %-t 2(radial basis function) shows best results compared to other 2 kernels
    i1=svmpredict(ones(length(s),1),s,cl);
    if(length(find(i1==1))/length(i1)==1)
        cnt1=cnt1+1;
    end
    %i2=svmpredict(randi([0:1],length(ns),1),ns,cl);
    i2=svmpredict(zeros(length(ns),1),ns,cl);
    if(length(find(i2==0))/length(i2)==1)
        cnt2=cnt2+1;
    end
    for(k=1:length(i1)) %data is seiz only if 10 sec interval is classified seiz activity
       if(k>=10)
        if(length(find(i1((k-9):k)==1))>=4)
            mat1=[mat1;1];
        else
            mat1=[mat1;0];
        end
       else
            mat1=[mat1;svmpredict(1,s(k,:),cl)];
       end
    end
        incorr_s=[incorr_s;length(find(mat1==1))/length(s)];
        for(k=1:length(i2))
            if(k>=10)
                if(length(find(i2((k-9):k)==1))>4)
                    mat2=[mat2;1];
                else
                    mat2=[mat2;0];
                end
            else
                mat2=[mat2;svmpredict(0,ns(k,:),cl)];
            end
        end
        incorr_ns=[incorr_ns;length(find(mat2==0))/length(ns)];
    end
    new_var=[new var;length(find(i1==1))/length(i1)*100,length(find(i2==0))/length(i2)*100,length(find(mat1==1))/length(mat1)*100,length(find(mat2==0))/length(mat2)*100];
end
save('Outputfile','new_var')