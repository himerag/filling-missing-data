function status=merge_chlo_inegi_paper(file1_,file2_,file3_,file4_,file5_,p)
status=1;


%% venia el for
vacios=[];
file1=file1_;
file2=file2_;
file3=file3_;
file4=file4_;
control=zeros(4,1);
if exist(file1,'file')==2
    data1=ncread(file1,'chlor_a');    
    data1(data1==0)=nan;
    data1(data1<0)=nan;
    dw1=isnan(data1);
    control(1)=sum(~dw1(:));
    % ventana movil
    dw=dw1;
    data=data1;    
    h_=ones(3);
    h_(2,2)=0;
    Y = filter2(h_,~isnan(data));
    Y(~dw)=0;
    puntos=find(Y>2);
    base=zeros(size(puntos));
    siz=size(data);
    data=double(data);
    vacios(end+1)=numel(puntos);
    aaa=tic;
    parfor i=1:length(puntos)
        [I,J]=ind2sub(siz,puntos(i));
        a=[I,J;I-1,J;I+1,J;...
            I,J-1;I-1,J-1;I+1,J-1;...
            I,J+1;I-1,J+1;I+1,J+1;...
            ];
        ind1=a(:,1)<1|a(:,1)>siz(1);
        ind2=a(:,2)<1|a(:,2)>siz(2);
        a(ind1|ind2,:)=[];
        tmp=data(sub2ind(siz,a(:,1),a(:,2)));
        if sum(~isnan(tmp))>2
            base(i)=nanmean(tmp);
        else
            base(i)=nan;
        end
    end
    
    fprintf('el archivo %s procesado en:\n%f[s]\n',file1,toc(aaa))
    fprintf('Porcentaje llenado: %0.2f%%\n',(sum(~isnan(base))/control(1))*100) 
    data(puntos)=base;
    data1=data;
    dw1=isnan(data);
    
else
    data1=[];
    dw1=[];
end
if exist(file2,'file')==2
    data2=ncread(file2,'chlor_a');
    data2(data2==0)=nan;
    data2(data2<0)=nan;
    dw2=isnan(data2);
    control(2)=sum(~dw(:));    
    dw=dw2;
    data=data2;
    h_=ones(3);
    h_(2,2)=0;
    Y = filter2(h_,~isnan(data));
    Y(~dw)=0;
    puntos=find(Y>2);
    base=zeros(size(puntos));
    siz=size(data);
    data=double(data);
    aaa=tic;
    parfor i=1:length(puntos)
        [I,J]=ind2sub(siz,puntos(i));
        a=[I,J;I-1,J;I+1,J;...
            I,J-1;I-1,J-1;I+1,J-1;...
            I,J+1;I-1,J+1;I+1,J+1;...
            ];
        ind1=a(:,1)<1|a(:,1)>siz(1);
        ind2=a(:,2)<1|a(:,2)>siz(2);
        a(ind1|ind2,:)=[];
        tmp=data(sub2ind(siz,a(:,1),a(:,2)));
        if sum(~isnan(tmp))>2
            base(i)=nanmean(tmp);
        else
            base(i)=nan;
        end
    end
    fprintf('%s\nProcesado en: %f[s]\n',file2,toc(aaa))
    fprintf('Porcentaje llenado: %0.2f%%\n',(sum(~isnan(base))/control(2))*100)
    data(puntos)=base;
    data2=data;
    dw2=isnan(data2);
    
else
    data2=[];
    dw2=[];
end
if exist(file3,'file')==2
    data3=ncread(file3,'chlor_a');
    data3(data3==0)=nan;
    data3(data3<0)=nan;
    dw3=isnan(data3);
    control(3)=sum(~(dw3(:)));    
    dw=dw3;
    data=data3;
    h_=ones(3);
    h_(2,2)=0;
    Y = filter2(h_,~isnan(data));
    Y(~dw)=0;
    puntos=find(Y>2);
    base=zeros(size(puntos));
    siz=size(data);
    data=double(data);
    vacios(end+1)=numel(puntos);
    aaa=tic;
    parfor i=1:length(puntos)
        [I,J]=ind2sub(siz,puntos(i));
        a=[I,J;I-1,J;I+1,J;...
            I,J-1;I-1,J-1;I+1,J-1;...
            I,J+1;I-1,J+1;I+1,J+1;...
            ];
        ind1=a(:,1)<1|a(:,1)>siz(1);
        ind2=a(:,2)<1|a(:,2)>siz(2);
        a(ind1|ind2,:)=[];
        tmp=data(sub2ind(siz,a(:,1),a(:,2)));
        if sum(~isnan(tmp))>2
            base(i)=nanmean(tmp);
        else
            base(i)=nan;
        end
    end
    fprintf('%s\nProcesado en: %f[s]\n',file3,toc(aaa))
    fprintf('Porcentaje llenado: %0.2f%%\n',(sum(~isnan(base))/control(3))*100)
    data(puntos)=base;
    data3=data;
    dw3=isnan(data);
else
    data3=[];
    dw3=[];
end
if exist(file4,'file')==2
    data4=ncread(file4,'chlor_a');
    data4(data4==0)=nan;
    data4(data4<0)=nan;
    dw4=isnan(data4);
    control(4)=sum(~dw4(:));
    dw=dw4;
    data=data4;
    h_=ones(3);
    h_(2,2)=0;
    Y = filter2(h_,~isnan(data));
    Y(~dw)=0;
    puntos=find(Y>2);
    base=zeros(size(puntos));
    siz=size(data);
    data=double(data);
    vacios(end+1)=numel(puntos);
    aaa=tic;
    parfor i=1:length(puntos)
        [I,J]=ind2sub(siz,puntos(i));
        a=[I,J;I-1,J;I+1,J;...
            I,J-1;I-1,J-1;I+1,J-1;...
            I,J+1;I-1,J+1;I+1,J+1;...
            ];
        ind1=a(:,1)<1|a(:,1)>siz(1);
        ind2=a(:,2)<1|a(:,2)>siz(2);
        a(ind1|ind2,:)=[];
        tmp=data(sub2ind(siz,a(:,1),a(:,2)));
        if sum(~isnan(tmp))>2
            base(i)=nanmean(tmp);
        else
            base(i)=nan;
        end
    end
    fprintf('%s\nProcesado en:%f[s]\n',file3,toc(aaa))
    fprintf('Porcentaje llenado: %0.2f%%\n',(sum(~isnan(base))/control(3))*100)
    data(puntos)=base;
    data4=data;
    dw4=isnan(data);
else
    data4=[];
    dw4=[];
end
%%
if ~any(control)
    retrun
end

[~,dummy]=sort(control,'descend');
eval(['file_s=file' num2str(dummy(1)) ';'])

if length(dummy)<2    
    return
end
eval(['data=data' num2str(dummy(1)) ';'])
eval(['dw=dw' num2str(dummy(1)) ';'])
%Calcula el RMSE de las 4 images
rmse_=inf(size(dummy));
for j=1:length(dummy)
    eval(['dummy1=~(dw|dw' num2str(dummy(j)) ');']);
    tmpx=data(dummy1);
    eval(['tmpy=data' num2str(dummy(j)) '(dummy1);']);
    rmse_(j)=(mean((tmpx-tmpy).^2))^0.5;
end

[~,I]=sort(rmse_);
dummy=dummy(I);
dw=isnan(data);
% base0=zeros(size(data));
% base0(~dw)=1;
% dataB=data;
% dataB(isnan(dataB))=0;
%%
dummy(1)=[];


if numel(dummy)>=1
    i=1;    
    eval(['datat=data' num2str(dummy(i)) ';' ]) 
    error_=(data-datat)./ datat;
    eval(['dw_=dw' num2str(dummy(i)) ';'])
    dataN=datat;
    eval(['dwt_=dw' num2str(dummy(i)) ';'])
    dummy1=~(dw|dwt_);
    dw_(~dw)=true;
    puntos=find(~dw_);
    base=zeros(size(puntos));
    siz=size(data);
    aaa=tic;
    parfor ii=1:length(puntos)
        dummy3=nan(4,4);
        [I,J]=ind2sub(siz,puntos(ii));
        %I=2195
        %J=1689
        
        a=find(dummy1(I,1:J));
        if ~isempty(a)
            J1=a(end);
            I1=I;
            dummy3(1,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(1,1)=1;
            dummy3(1,2)=1;
            dummy3(1,3)=error_(I1,J1);
        end
        a=find(dummy1(I,J+1:siz(2)));
        if ~isempty(a)
            J1=a(1)+J;
            I1=I;
            dummy3(2,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(2,1)=1;
            dummy3(2,2)=1;
            dummy3(2,3)=error_(I1,J1);
        end
        
        a=find(dummy1(1:I,J));
        if ~isempty(a)
            J1=J;
            I1=a(end);
            dummy3(3,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(3,1)=1;
            dummy3(3,2)=1;
            dummy3(3,3)=error_(I1,J1);
        end
        
        a=find(dummy1(I+1:siz(1),J));
        if ~isempty(a)
            J1=J;
            I1=a(1)+I;
            dummy3(4,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(4,1)=1;
            dummy3(4,2)=1;
            dummy3(4,3)=error_(I1,J1);
        end
        dummy3(isnan(dummy3(:,1)),:)=[];
        ee=sum(dummy3(:,3)./(dummy3(:,4).^2))/sum(1./(dummy3(:,4).^2));
        base(ii)=dataN(I,J)+(ee*dataN(I,J)/100);
    end
    fprintf('IDW data%d: %f[s]\n',dummy(i),toc(aaa))
    data(puntos)=base;
    dw=isnan(data);
end
if numel(dummy)>=2
    i=2;
    eval(['datat=data' num2str(dummy(i)) ';' ]) 
    error_=(data-datat)./ datat;
    eval(['dw_=dw' num2str(dummy(i)) ';'])
    dataN=datat;
    eval(['dwt_=dw' num2str(dummy(i)) ';'])
    dummy1=~(dw|dwt_);
    dw_(~dw)=true;

    puntos=find(~dw_);
    base=zeros(size(puntos));
    siz=size(data);
    aaa=tic;
    parfor ii=1:length(puntos)
        dummy3=nan(4,4);
        [I,J]=ind2sub(siz,puntos(ii));
        a=find(dummy1(I,1:J));
        if ~isempty(a)
            J1=a(end);
            I1=I;
            dummy3(1,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(1,1)=1;
            dummy3(1,2)=1;
            dummy3(1,3)=error_(I1,J1);
        end
        a=find(dummy1(I,J+1:siz(2)));
        if ~isempty(a)
            J1=a(1)+J;
            I1=I;
            dummy3(2,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(2,1)=1;
            dummy3(2,2)=1;
            dummy3(2,3)=error_(I1,J1);
        end
        
        a=find(dummy1(1:I,J));
        if ~isempty(a)
            J1=J;
            I1=a(end);
            dummy3(3,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(3,1)=1;
            dummy3(3,2)=1;
            dummy3(3,3)=error_(I1,J1);
        end
        
        a=find(dummy1(I+1:siz(1),J));
        if ~isempty(a)
            J1=J;
            I1=a(1)+I;
            dummy3(4,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(4,1)=1;
            dummy3(4,2)=1;
            dummy3(4,3)=error_(I1,J1);
        end
        dummy3(isnan(dummy3(:,1)),:)=[];
        ee=sum(dummy3(:,3)./(dummy3(:,4).^2))/sum(1./(dummy3(:,4).^2));
        base(ii)=dataN(I,J)+(ee*dataN(I,J)/100);
    end
    fprintf('IDW data%d: %f[s]\n',dummy(i),toc(aaa))
    data(puntos)=base;
    dw=isnan(data);
end
if numel(dummy)>=3
    i=3;
    eval(['datat=data' num2str(dummy(i)) ';' ]) 
    error_=(data-datat)./ datat;
    eval(['dw_=dw' num2str(dummy(i)) ';'])
    dataN=datat;
    eval(['dwt_=dw' num2str(dummy(i)) ';'])
    dummy1=~(dw|dwt_);
    dw_(~dw)=true;

    puntos=find(~dw_);
    base=zeros(size(puntos));
    siz=size(data);
    aaa=tic;
    parfor ii=1:length(puntos)
        dummy3=nan(4,4);
        [I,J]=ind2sub(siz,puntos(ii));
        a=find(dummy1(I,1:J));
        if ~isempty(a)
            J1=a(end);
            I1=I;
            dummy3(1,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(1,1)=1;
            dummy3(1,2)=1;
            dummy3(1,3)=error_(I1,J1);
        end
        a=find(dummy1(I,J+1:siz(2)));
        if ~isempty(a)
            J1=a(1)+J;
            I1=I;
            dummy3(2,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(2,1)=1;
            dummy3(2,2)=1;
            dummy3(2,3)=error_(I1,J1);
        end
        
        a=find(dummy1(1:I,J));
        if ~isempty(a)
            J1=J;
            I1=a(end);
            dummy3(3,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(3,1)=1;
            dummy3(3,2)=1;
            dummy3(3,3)=error_(I1,J1);
        end
        
        a=find(dummy1(I+1:siz(1),J));
        if ~isempty(a)
            J1=J;
            I1=a(1)+I;
            dummy3(4,4)=((J1-J)^2+(I1-I)^2)^0.5;
            dummy3(4,1)=1;
            dummy3(4,2)=1;
            dummy3(4,3)=error_(I1,J1);
        end
        dummy3(isnan(dummy3(:,1)),:)=[];
        ee=sum(dummy3(:,3)./(dummy3(:,4).^2))/sum(1./(dummy3(:,4).^2));
        base(ii)=dataN(I,J)+(ee*dataN(I,J)/100);
    end
    fprintf('IDW data%d: %f[s]\n',dummy(i),toc(aaa))
    data(puntos)=base;
    dw=isnan(data);
end
info_=ncinfo(file1);
ncwriteschema(file5_,info_)
ncwrite(file5_,'chlor_a',data);
data(~isnan(data))=1;
ncwrite(file5_,'chlor_a_count',data);
status=0;

