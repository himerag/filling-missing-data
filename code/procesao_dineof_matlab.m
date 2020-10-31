% Parametros
base_='/home/terra/Descargas/paper/datos/base_dineof/';
basep_='/home/terra/Descargas/paper/datos/base_dineof_p/';
baser1_='/home/terra/Descargas/paper/datos/base_dineof_r1/';
path_dineof='/home/terra/Descargas/paper/Codigo/dineof-3.0/bin/Linux/dineof-3.0-x64-linux';
segmentos=[4000,6000];

lista=dir([base_ '*.nc']);
numimages=length(lista);
info_=ncinfo([base_ lista(1).name],'chlor_a');
dim1=round(linspace(1,info_.Size(1),ceil(info_.Size(1)/segmentos(1))+1));
cuantos1=diff(dim1)+1;
cuantos1(end)=info_.Size(1)-sum(cuantos1(1:end-1));
ini1=cumsum(cuantos1)+1;
ini1=[1 ini1(1:end-1)];
dim2=round(linspace(1,info_.Size(2),ceil(info_.Size(2)/segmentos(2))+1));
cuantos2=diff(dim2)+1;
cuantos2(end)=info_.Size(2)-sum(cuantos2(1:end-1));
ini2=cumsum(cuantos2)+1;
ini2=[1 ini2(1:end-1)];
x1=repmat(ini1,numel(ini2),1);
cuantos1=repmat(cuantos1,numel(ini2),1);
x2=repmat(ini2',numel(ini1),1);
cuantos2=repmat(cuantos2',numel(ini1),1);
control=[x1(:) x2(:) cuantos1(:) cuantos2(:)];
i=1;
j=1;
packini=struct;
packini.data=[];
packini.mask=[];
packini.time=[];
packini.alpha=0.01;
packini.numit=3;
packini.nev=5;
packini.neini=1;
packini.ncv=10;
packini.tol=1.0e-8;
packini.nitemax=300;
packini.toliter=1.0e-3;
packini.rec=0;
packini.eof = 1;
packini.norm = 0;
packini.Output=baser1_;
packini.results=[];
packini.seed=243435;
packini.cloud_size=100;
for i=1:size(control,1)
    %hacemos el cubo
    base1=zeros(control(i,3),control(i,4),numimages);
    for j=1:numimages
        tmp=ncread([base_ lista(j).name],'chlor_a',[control(i,1),control(i,2)],[control(i,3),control(i,4)]);
        base1(:,:,j)=tmp;
    end
    losnan=all(isnan(base1),3);
    file1=[basep_ sprintf('data_dineof%03d.gher',i)];
    file2=[basep_ sprintf('mask_dineof%03d.gher',i)];
    file3=[basep_ sprintf('time_dineof%03d.gher',i)];
    file4=[basep_ sprintf('settings_dineof%03d.init',i)];
    file5=[baser1_ sprintf('outfile_dineof%03d.gher',i)];
    file6=[baser1_ sprintf('cmdout_dineof%03d.init',i)];
    flag1 = gwrite(file1,base1);
    flag2 = gwrite(file2,double(~losnan));
    flag3 = gwrite(file3,1:size(base1,3));
    fid=fopen(file4,'w');
    packini.data=file1;
    packini.mask=file2;
    packini.time=file3;
    packini.results=file5;
    for item=fieldnames(packini)'
        switch item{1}
            case {'data' 'mask' 'results'}
                fprintf(fid,[item{1} ' = [''%s'']\n'],packini.(item{1}));
            case {'time' 'Output'}
                fprintf(fid,[item{1} ' = ''%s''\n'],packini.(item{1}));
            case {'numit' 'nev' 'neini' 'ncv' 'nitemax' 'rec' 'eof' 'norm' 'seed' 'cloud_size'}
                fprintf(fid,[item{1} ' = %d\n'],packini.(item{1}));
            otherwise
                fprintf(fid,[item{1} ' = %0.1e\n'],packini.(item{1}));
        end
    end
    fclose(fid);
    [status,cmdout]=system([path_dineof ' ' file4],'-echo');
    fid=fopen(file6,'w');
    fprintf(fid,'%s\n',cmdout);
    fclose(fid);
end

for j=1:numimages
    info=ncinfo([base_ lista(j).name]);
    ncwriteschema([baser1_ ['Dineof_' lista(j).name]],info)
end
for i=1:size(control,1)
    file5=[baser1_ sprintf('outfile_dineof%03d.gher',i)];
    base1=gread(file5);
    for j=1:numimages
        ncwrite([baser1_ ['Dineof_' lista(j).name]],'chlor_a',base1(:,:,j),[control(i,1),control(i,2)]);    
        ncwrite([baser1_ ['Dineof_' lista(j).name]],'chlor_a_count',int32(~isnan(base1(:,:,j))),[control(i,1),control(i,2)]);    
    end
end