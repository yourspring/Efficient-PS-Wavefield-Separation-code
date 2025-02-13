function [vx,vz,vxp,vzp,vxs,vzs] = ani_decomposition_wavefront_phase(vx1,vz1,vp,vs,epsilon,delta,dz,dx)%��������ʵ�ָ������Բ����ֽ��һ�ײ�ǰ��λ��
%function [vx,vz,vxp,vzp,vxs,vzs] = ani_decomposition_wavefront_phase(vx1,vz1,vp,vs,epsilon,delta,dz,dx,beta,num)
%��������ֽ��Ĳ�������
%һ�ײ�ǰ��λ��ʵ�ָ������Բ����ֽ�
%%betaΪsor�������ɳ����ӣ�num��ʾ���ɷ�������������

[NZ, NX] = size(vx1);%��ȡ���벨���ĳߴ�
%Ƶ����λ��ͨ������shift_x��shift_z��ʵ�ֵģ���Щ�����ǻ��ڲ���ʸ���Ϳռ䲽������ģ�����ż���ߴ�����飬�ر������ĵ����λ���ӣ�ʹ���Ϊʵ��
%-shift wavefield---
[kz,kx]=generate_wavenumber(NZ,NX,dz,dx);%���ɲ���ʸ��
shift_x=exp(1i*(dx/2).*kx);
shift_z=exp(1i*(-dz/2).*kz);%����Ƶ����λ����
%����Ƶ����λ�����е������������NZ��NXΪż��ʱ��
if mod(NX,2)==0
    shift_x(:,NX/2+1) = real(shift_x(:,NX/2+1));
end
if mod(NZ,2)==0
    shift_z(NZ/2+1,:) = real(shift_z(NZ/2+1,:));
end
%����ע�͵Ĳ��֣����ڼ��㸴������ʸ��
% jkx=1i.*kx;
% jkz=1i.*kz;
% if mod(NX,2)==0
%     jkx(:,NX/2+1) = real(jkx(:,NX/2+1));
% end
% if mod(NZ,2)==0
%     jkz(NZ/2+1,:) = real(jkz(NZ/2+1,:));
% end

vx=vx1;
vz=ifft(shift_z.*fft(vz1,[],1),[],1);%��vz����һϵ��Ƶ��λ�ƺ��渵��Ҷ�任����  ��vz1��z������п����渵��Ҷ�任
vz=ifft(shift_x.*fft(vz,[],2),[],2);  %��vz�Ƶ������������棬֮�������Ĳ���󵼷��벨��
% nt = 1501; % ʱ���������
% dt = 1e-3; % ʱ�䲽��
ux=vx1;
uz=vz;


% tic;
%-----------------
r1=(1+2*epsilon).*vp.*vp - vs.*vs;
r2=sqrt(((1+2*delta).*vp.*vp - vs.*vs).*(vp.*vp - vs.*vs));
r3=vp.*vp - vs.*vs;
r4=2*(delta-epsilon).*vp.*vp.*(vp.*vp - vs.*vs);

% compute direction for vx
gz_vx=derivate1_fd8(vx,1,dz); % d vx/dz ������vx�ֱ���x��z������ݶȣ��������ϵı仯�ʣ�
gx_vx=derivate1_fd8(vx,2,dx); % d vx/dx 

nxs_vx = zeros(NZ,NX);%��ʼ����ر����ľ������ڴ���
nzs_vx = zeros(NZ,NX);
rx_vx = zeros(NZ,NX);
rz = zeros(NZ,NX);

rr_vx = gx_vx.*gx_vx + gz_vx.*gz_vx;%�ݶ�������ģ��ƽ��
%����rr_vx��ֵ�����������������nxs_vx nzs_vx
for ix=1:NX
    for iz=1:NZ
        if rr_vx(iz,ix)==0.0
            nxs_vx(iz,ix) = 0.0;
            nzs_vx(iz,ix) = 0.0;
        else
            nxs_vx(iz,ix) = gx_vx(iz,ix)*gx_vx(iz,ix)/(gx_vx(iz,ix)*gx_vx(iz,ix) + gz_vx(iz,ix)*gz_vx(iz,ix)); % nx*nx
            nzs_vx(iz,ix) = gz_vx(iz,ix)*gz_vx(iz,ix)/(gx_vx(iz,ix)*gx_vx(iz,ix) + gz_vx(iz,ix)*gz_vx(iz,ix)); % nz*nz
        end%���rr_vx=0�������������Ϊ�㣬��Ȼʹ���ݶ������ķ���������õ�ķ������
        
    end
end
nxs_vx = imgaussfilt(nxs_vx,5);%֮���ڶ������������nxs_vx nzs_vxӦ�ø�˹�˲����Ӷ�ƽ����Щ������ȥ������
nzs_vx = imgaussfilt(nzs_vx,5);

for ix=1:NX%����ٴ�ͨ��ѭ������ÿ�������
    for iz=1:NZ
        if rr_vx(iz,ix)==0.0%����rr_vx������ֵΪ0.0��λ�ã�rx_vx��r1��ƽ��
            rx_vx(iz,ix) = r1(iz,ix)*r1(iz,ix);
        else
            rx_vx(iz,ix) = (r1(iz,ix)+r4(iz,ix)*nzs_vx(iz,ix)/(r1(iz,ix)*nxs_vx(iz,ix) + r3(iz,ix)*nzs_vx(iz,ix)))^2;%��Ȼ�ø�ʽ����
        end
        rz(iz,ix) = r2(iz,ix)*r2(iz,ix);%ͬʱrz��ֵ������Ϊr2��ƽ��

    end
end
% 
% a=ones(NZ,NX); b=rz./rx_vx;%��ʼ��������ab��ֵ
% 
% [wx]=possion2dsolver6_sor(a,b,vx,dz,dx,beta,num);%���ò��ɺ�����������ĳ����wx
% 
% 
% wx(isnan(wx))=0;%����wx�е�NaNֵ���������滻Ϊ0
% % % % 

%% compute direction for vz
gz_vz=derivate1_fd8(vz,1,dz); % d vz/dz ת�����vz�ķ���������ݶ�
gx_vz=derivate1_fd8(vz,2,dx); % d vz/dx

nxs_vz = zeros(NZ,NX);%��ʼ�����ڴ洢��������ľ���
nzs_vz = zeros(NZ,NX);
rx_vz = zeros(NZ,NX);%��ʼ�����ں����������
rr_vz = gx_vz.*gx_vz + gz_vz.*gz_vz;%�ݶ�ģ��ƽ��

for ix=1:NX%����ѭ��ͬ��
    for iz=1:NZ
        if rr_vz(iz,ix)==0.0
            nxs_vz(iz,ix) = 0.0;
            nzs_vz(iz,ix) = 0.0;
        else
            nxs_vz(iz,ix) = gx_vz(iz,ix)*gx_vz(iz,ix)/(gx_vz(iz,ix)*gx_vz(iz,ix) + gz_vz(iz,ix)*gz_vz(iz,ix)); % nx*nx
            nzs_vz(iz,ix) = gz_vz(iz,ix)*gz_vz(iz,ix)/(gx_vz(iz,ix)*gx_vz(iz,ix) + gz_vz(iz,ix)*gz_vz(iz,ix)); % nz*nz
        end
    end
end
nxs_vz = imgaussfilt(nxs_vz,5);
nzs_vz = imgaussfilt(nzs_vz,5);


for ix=1:NX
    for iz=1:NZ
        if rr_vz(iz,ix)==0.0
            rx_vz(iz,ix) = r1(iz,ix)*r1(iz,ix);
        else
            rx_vz(iz,ix) = (r1(iz,ix)+r4(iz,ix)*nzs_vz(iz,ix)/(r1(iz,ix)*nxs_vz(iz,ix) + r3(iz,ix)*nzs_vz(iz,ix)))^2;
        end
        rz(iz,ix) = r2(iz,ix)*r2(iz,ix);

    end
   
end

% 
alpha=1e-10;
  vxdx = derivate1_fd8(vx,2,dx);
vzdz = derivate1_fd8(vx,2,dz);
  r_vx=r2./(r1+r4.*nzs_vx./(r1.*nxs_vx + r3.*nzs_vx + alpha));
  r_vz=r2./(r1+r4.*nzs_vz./(r1.*nxs_vz + r3.*nzs_vz + alpha));
  co=(vxdx.*vxdx)./((vxdx.*vxdx)+(r_vz.*r_vz.*(vzdz.*vzdz)))+alpha;
  si=1-co;
  u1px=ux.*(vp.*vp).*(1+2*epsilon.*si);
  u1pz=uz.*(vp.*vp).*(1+2*epsilon.*si);
  u1pxdxx = derivate2_fd8(u1px,2,dx);
u1pzdzz = derivate2_fd8(u1pz,1,dz);
u1pxdx = derivate1_fd8(u1px,2,dx);
u1pxdxz = derivate1_fd8(u1pxdx,1,dz);
u1pzdx = derivate1_fd8(u1pz,2,dx);
u1pzdxz = derivate1_fd8(u1pzdx,1,dz);
u1pxp= u1pxdxx + r_vz.*u1pzdxz;
u1pzp = r_vx.*u1pxdxz + r_vz.^2.*u1pzdzz;
  u2px=ux.*(vp.*vp).*(1+2*epsilon.*si).*co;
  u2pz=uz.*(vp.*vp).*(1+2*epsilon.*si).*co;

  u2pxdxx = derivate2_fd8(u2px,2,dx);
u2pzdzz = derivate2_fd8(u2pz,1,dz);
u2pxdx = derivate1_fd8(u2px,2,dx);
u2pxdxz = derivate1_fd8(u2pxdx,1,dz);
u2pzdx = derivate1_fd8(u2pz,2,dx);
u2pzdxz = derivate1_fd8(u2pzdx,1,dz);
u2pxp = u2pxdxx + r_vz.*u2pzdxz;
u2pzp = r_vx.*u2pxdxz + r_vz.^2.*u2pzdzz;
   u3px=-u2px;
  u3pz=-u2pz;
   u3pxdxx = derivate2_fd8(u3px,2,dx);
u3pzdzz = derivate2_fd8(u3pz,1,dz);
u3pxdx = derivate1_fd8(u3px,2,dx);
u3pxdxz = derivate1_fd8(u3pxdx,1,dz);
u3pzdx = derivate1_fd8(u3pz,2,dx);
u3pzdxz = derivate1_fd8(u3pzdx,1,dz);
u3pxp = u3pxdxx + r_vz.*u3pzdxz;
u3pzp = r_vx.*u3pxdxz + r_vz.^2.*u3pzdzz;
 u1sx=ux.*(vs.*vs);
 u1sz=uz.*(vs.*vs);

   u1sxdzz = derivate2_fd8(u1sx,1,dz);
u1szdxx= derivate2_fd8(u1sz,2,dx);
u1sxdx = derivate1_fd8(u1sx,2,dx);
u1sxdxz = derivate1_fd8(u1sxdx,1,dz);
u1szdx = derivate1_fd8(u1sz,2,dx);
u1szdxz = derivate1_fd8(u1szdx,1,dz);
u1sxp = r_vx.^2.*u1sxdzz-r_vz.*u1szdxz;
u1szp = -r_vx.*u1sxdxz + u1szdxx;
 u2sx=ux.*(vs.*vs).*co;
   u2sz=uz.*(vs.*vs).*co;

%     ����u2sx��u2sz�ĵ���
    u2sxdx = derivate1_fd8(u2sx, 2, dx);
    u2szdxx = derivate2_fd8(u2sz, 2, dx);
    u2szdx = derivate1_fd8(u2sz, 2, dx);
    u2sxdzz = derivate2_fd8(u2sx, 1, dz);
    
%     ���㽻�浼��
    u2sxdxz = derivate1_fd8(u2sxdx, 1, dz);
    u2szdxz = derivate1_fd8(u2szdx, 1, dz);
    
%     ����u2sx��u2sz��������
    u2sxp = r_vx.^2.*u2sxdzz-r_vz.*u2szdxz;
u2szp = -r_vx.*u2sxdxz + u2szdxx;
    
%     ����u3sx��u3sz��������u2sx��u2sz���෴����
    u3sx = -u2sx;
    u3sz = -u2sz;
    

%     ����u3sx��u3sz�ĵ���
    u3sxdx = derivate1_fd8(u3sx, 2, dx);
    u3sxdzz = derivate2_fd8(u3sx, 1, dz);
    u3szdxx = derivate2_fd8(u3sz, 2, dx);
    u3sxdxz = derivate1_fd8(u3sxdx, 1, dz);
    u3szdx = derivate1_fd8(u3sz, 2, dx);
    u3szdxz = derivate1_fd8(u3szdx, 1, dz);
    
%     ����u3sx��u3sz����������ע�⣺�������u2sxp_temp��u3sxp_temp����֮ǰ����õĻ�����Ϊ�����ṩ�ģ����ڱ����ظ����㣩
%     ������ǣ�����Ҫ���¼������ǣ�������Ϊ�˼�ʾ�������Ǽ��������Ѿ����ڡ�
%     ʵ���ϣ���������Ҫ�Ƴ���Щ��������ǰ��ļ�����ֱ�ӻ�ȡu2sxp��u3sxp�ȡ�
   u3sxp = r_vx.^2.*u3sxdzz-r_vz.*u3szdxz;
u3szp = -r_vx.*u3sxdxz + u3szdxx;
%     u2sxp_temp= u2sxp;
%     u2szp_temp=u2szp;
%      u3sxp_temp=u3sxp;
%      u3szp_temp=u3szp ;
%     �������յ��ٶȷ���
    common_term = ((4 .* epsilon) .* vp.^2) ./ (vp.^2 - vs.^2)+alpha;
    delta_term = ((2. * delta) .* vp.^2)./ (vp.^2 - vs.^2)+alpha;
    
    vxp = -(u1pxp + common_term .* u2pxp - delta_term.* u3pxp+ alpha);
    vzp = -(u1pzp + common_term .* u2pzp - delta_term .* u3pzp+ alpha);
    vxs = -(u1sxp + common_term .* u2sxp - delta_term .* u3sxp+ alpha); % ʹ����ʱ����
    vzs = -(u1szp + common_term .* u2szp - delta_term .* u3szp+ alpha); % ʹ����ʱ����
    vxp(isnan(vxp)) = 0;
    vzp(isnan(vzp)) = 0;
    vxs(isnan(vxs)) = 0;
    vzs(isnan(vzs)) = 0;
end

