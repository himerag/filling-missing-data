year_=[2018 2019];
api='https://oceandata.sci.gsfc.nasa.gov/ob/getfile/';
sensores={'A%d%03d%d%03d.L3m_8D_CHL_chlor_a_4km.nc',...
    'T%d%03d%d%03d.L3m_8D_CHL_chlor_a_4km.nc',...
    'V%d%03d%d%03d.L3m_8D_SNPP_CHL_chlor_a_4km.nc',...
    'V%d%03d%d%03d.L3m_8D_JPSS1_CHL_chlor_a_4km.nc'};
numsensor=numel(sensores);
snap_graf='/home/gedeon/Documents/procesar_paper_himer/mosaic_esau1.xml';
pathp1='/home/gedeon/Documents/procesar_paper_himer/';
gpt_path='/opt/snap/bin/gpt';
p = gcp;
for i=1:numel(year_)
    p2=(datenum(year_(i),1,1):8:datenum(year_(i),12,31))-datenum(year_(i),1,1);
    p1=p2+1;
    p2(1)=[];
    p2(end+1)=datenum(year_(i),12,31)-datenum(year_(i),1,1)+1;
    p3=repmat(year_(i),1,numel(p2));   
    for j=1:size(p1,1)
        aa=tic;
        control=false(numsensor,1);
        archivos=cell(numsensor,1);
        filename_out=sprintf('M%d%03d%d%03d.L3m_8D_CHL_chlor_a_1km.nc',p3(j),p1(j),p3(j),p2(j));
        if exist([pathp1 filename_out],'file')
            continue
        end
        for s=1:numsensor
            filename=sprintf(sensores{s},p3(j),p1(j),p3(j),p2(j));
            url=[api filename];
            if exist([pathp1 filename],'file')
                control(s)=true;
                continue
            else
                if exist([pathp1 'P1_' filename],'file')
                    control(s)=true;
                    archivos{s}=[pathp1 'P1_' filename];
                    continue
                else
                    a=tic;
                    [status,cmdout]=system(['wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --content-disposition https://oceandata.sci.gsfc.nasa.gov/ob/getfile/' filename]);
                    b=toc(a);
                    if status
                        fprintf('%s: Error descarga NASA\n',filename)
                        continue
                    end
                    fprintf('%s: Descarga NASA en %0.3f[s]\n',filename,b)
                    a=tic;
                    [status,cmdout]=system([gpt_path ' ' snap_graf ' -e -t ' [pathp1 'P1_' filename] ' -f NetCDF4-BEAM ' [pathp1  filename] ]);
                    b=toc(a);
                    if status
                        fprintf('%s: Error Mosaixco GPT\n',['P1_' filename])
                        continue
                    end
                    control(s)=true;
                    archivos{s}=[pathp1 'P1_' filename];
                    fprintf('%s: Mosaico GPT en %0.3f[s]\n',filename,b)
                    delete([pathp1 filename])
                end
            end 
        end
        if ~all(control)
            continue
        end
        status=merge_chlo_inegi_paper(archivos{1},archivos{2},archivos{3},archivos{4},[pathp1 filename_out],p);
        if status
            continue
        end
        for ii=1:numsensor
            delete(archivos{ii})
        end
        fprintf('%s: Datos procesados en: %0.3f[s]\n',filename_out,toc(aa))
    end
end
delete(p)