function temp=read_matrix(file,nz,nx)

fid=fopen(file,'r');
temp=fread(fid,[nz,nx],'float');
fclose(fid);