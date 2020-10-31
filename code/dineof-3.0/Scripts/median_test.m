function dataout = median_test(dataini,mask,box)

M = size(dataini,1);
N = size(dataini,2);
Z = size(dataini,3);

dataout=NaN*zeros(M,N,Z);
disp(['median test: calculating... Please be patient'])
for i=box:M-box
  for j=box:N-box
     if mask(i,j)==1;
  %      for k=1:Z          
          temp = dataini(i-box+1:i+box,j-box+1:j+box,:);
          mid_size=size(temp,1)*size(temp,2);
          temp = reshape(temp,[mid_size Z]);
          med=nanmedian(temp,1);
          med_mat = repmat(med,[mid_size,1]);
          mad = 1.4826 * nanmedian(abs(temp-med_mat),1);
          dataout(i,j,:)=(abs(squeeze(dataini(i,j,:))'-med))./mad;
  %      end
     end
  end
%  disp(['median test: calculating row ' num2str(i) ' out of ' num2str(M-box) '. Please be patient'])
end


          