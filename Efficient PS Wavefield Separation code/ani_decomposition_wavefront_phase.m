function [vx,vz,vxp,vzp,vxs,vzs] = ani_decomposition_wavefront_phase(vx1,vz1,vp,vs,epsilon,delta,dz,dx)%函数用于实现各向异性波场分解的一阶波前相位法
%function [vx,vz,vxp,vzp,vxs,vzs] = ani_decomposition_wavefront_phase(vx1,vz1,vp,vs,epsilon,delta,dz,dx,beta,num)
%函数输出分解后的波场分量
%一阶波前相位法实现各向异性波场分解
%%beta为sor迭代超松弛因子，num表示泊松方程求解迭代次数

[NZ, NX] = size(vx1);%获取输入波场的尺寸
%频域移位是通过乘以shift_x和shift_z来实现的，这些因子是基于波数矢量和空间步长计算的，对于偶数尺寸的数组，特别处理中心点的移位因子，使其变为实数
%-shift wavefield---
[kz,kx]=generate_wavenumber(NZ,NX,dz,dx);%生成波数矢量
shift_x=exp(1i*(dx/2).*kx);
shift_z=exp(1i*(-dz/2).*kz);%计算频域移位因子
%处理频域移位因子中的特殊情况（当NZ或NX为偶数时）
if mod(NX,2)==0
    shift_x(:,NX/2+1) = real(shift_x(:,NX/2+1));
end
if mod(NZ,2)==0
    shift_z(NZ/2+1,:) = real(shift_z(NZ/2+1,:));
end
%下面注释的部分，用于计算复数波数矢量
% jkx=1i.*kx;
% jkz=1i.*kz;
% if mod(NX,2)==0
%     jkx(:,NX/2+1) = real(jkx(:,NX/2+1));
% end
% if mod(NZ,2)==0
%     jkz(NZ/2+1,:) = real(jkz(NZ/2+1,:));
% end

vx=vx1;
vz=ifft(shift_z.*fft(vz1,[],1),[],1);%对vz做了一系列频域位移和逆傅里叶变换处理  对vz1在z方向进行快速逆傅里叶变换
vz=ifft(shift_x.*fft(vz,[],2),[],2);  %将vz移到整数网格上面，之后用中心差分求导分离波场
% nt = 1501; % 时间采样点数
% dt = 1e-3; % 时间步长
ux=vx1;
uz=vz;


% tic;
%-----------------
r1=(1+2*epsilon).*vp.*vp - vs.*vs;
r2=sqrt(((1+2*delta).*vp.*vp - vs.*vs).*(vp.*vp - vs.*vs));
r3=vp.*vp - vs.*vs;
r4=2*(delta-epsilon).*vp.*vp.*(vp.*vp - vs.*vs);

% compute direction for vx
gz_vx=derivate1_fd8(vx,1,dz); % d vx/dz 计算了vx分别在x，z方向的梯度（两方向上的变化率）
gx_vx=derivate1_fd8(vx,2,dx); % d vx/dx 

nxs_vx = zeros(NZ,NX);%初始化相关变量的矩阵用于储存
nzs_vx = zeros(NZ,NX);
rx_vx = zeros(NZ,NX);
rz = zeros(NZ,NX);

rr_vx = gx_vx.*gx_vx + gz_vx.*gz_vx;%梯度向量的模的平方
%根据rr_vx的值来计算两个方向分量nxs_vx nzs_vx
for ix=1:NX
    for iz=1:NZ
        if rr_vx(iz,ix)==0.0
            nxs_vx(iz,ix) = 0.0;
            nzs_vx(iz,ix) = 0.0;
        else
            nxs_vx(iz,ix) = gx_vx(iz,ix)*gx_vx(iz,ix)/(gx_vx(iz,ix)*gx_vx(iz,ix) + gz_vx(iz,ix)*gz_vx(iz,ix)); % nx*nx
            nzs_vx(iz,ix) = gz_vx(iz,ix)*gz_vx(iz,ix)/(gx_vx(iz,ix)*gx_vx(iz,ix) + gz_vx(iz,ix)*gz_vx(iz,ix)); % nz*nz
        end%如果rr_vx=0则两个方向分量为零，不然使用梯度向量的分量来计算该点的方向分量
        
    end
end
nxs_vx = imgaussfilt(nxs_vx,5);%之后在对两个方向分量nxs_vx nzs_vx应用高斯滤波，从而平滑这些分量并去除噪声
nzs_vx = imgaussfilt(nzs_vx,5);

for ix=1:NX%最后再次通过循环遍历每个网格点
    for iz=1:NZ
        if rr_vx(iz,ix)==0.0%对于rr_vx矩阵中值为0.0的位置，rx_vx设r1的平方
            rx_vx(iz,ix) = r1(iz,ix)*r1(iz,ix);
        else
            rx_vx(iz,ix) = (r1(iz,ix)+r4(iz,ix)*nzs_vx(iz,ix)/(r1(iz,ix)*nxs_vx(iz,ix) + r3(iz,ix)*nzs_vx(iz,ix)))^2;%不然用该式计算
        end
        rz(iz,ix) = r2(iz,ix)*r2(iz,ix);%同时rz的值被设置为r2的平方

    end
end
% 
% a=ones(NZ,NX); b=rz./rx_vx;%初始化输入量ab的值
% 
% [wx]=possion2dsolver6_sor(a,b,vx,dz,dx,beta,num);%调用泊松函数，来计算某个量wx
% 
% 
% wx(isnan(wx))=0;%处理wx中的NaN值。将他们替换为0
% % % % 

%% compute direction for vz
gz_vz=derivate1_fd8(vz,1,dz); % d vz/dz 转向计算vz的方向分量的梯度
gx_vz=derivate1_fd8(vz,2,dx); % d vz/dx

nxs_vz = zeros(NZ,NX);%初始化用于存储方向分量的矩阵
nzs_vz = zeros(NZ,NX);
rx_vz = zeros(NZ,NX);%初始化用于后续计算矩阵
rr_vz = gx_vz.*gx_vz + gz_vz.*gz_vz;%梯度模的平方

for ix=1:NX%后续循环同上
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

%     计算u2sx和u2sz的导数
    u2sxdx = derivate1_fd8(u2sx, 2, dx);
    u2szdxx = derivate2_fd8(u2sz, 2, dx);
    u2szdx = derivate1_fd8(u2sz, 2, dx);
    u2sxdzz = derivate2_fd8(u2sx, 1, dz);
    
%     计算交叉导数
    u2sxdxz = derivate1_fd8(u2sxdx, 1, dz);
    u2szdxz = derivate1_fd8(u2szdx, 1, dz);
    
%     计算u2sx和u2sz的物理量
    u2sxp = r_vx.^2.*u2sxdzz-r_vz.*u2szdxz;
u2szp = -r_vx.*u2sxdxz + u2szdxx;
    
%     计算u3sx和u3sz（它们是u2sx和u2sz的相反数）
    u3sx = -u2sx;
    u3sz = -u2sz;
    

%     计算u3sx和u3sz的导数
    u3sxdx = derivate1_fd8(u3sx, 2, dx);
    u3sxdzz = derivate2_fd8(u3sx, 1, dz);
    u3szdxx = derivate2_fd8(u3sz, 2, dx);
    u3sxdxz = derivate1_fd8(u3sxdx, 1, dz);
    u3szdx = derivate1_fd8(u3sz, 2, dx);
    u3szdxz = derivate1_fd8(u3szdx, 1, dz);
    
%     计算u3sx和u3sz的物理量（注意：这里假设u2sxp_temp和u3sxp_temp等是之前计算好的或者作为输入提供的，用于避免重复计算）
%     如果不是，则需要重新计算它们，但这里为了简化示例，我们假设它们已经存在。
%     实际上，您可能需要移除这些参数并从前面的计算中直接获取u2sxp和u3sxp等。
   u3sxp = r_vx.^2.*u3sxdzz-r_vz.*u3szdxz;
u3szp = -r_vx.*u3sxdxz + u3szdxx;
%     u2sxp_temp= u2sxp;
%     u2szp_temp=u2szp;
%      u3sxp_temp=u3sxp;
%      u3szp_temp=u3szp ;
%     计算最终的速度分量
    common_term = ((4 .* epsilon) .* vp.^2) ./ (vp.^2 - vs.^2)+alpha;
    delta_term = ((2. * delta) .* vp.^2)./ (vp.^2 - vs.^2)+alpha;
    
    vxp = -(u1pxp + common_term .* u2pxp - delta_term.* u3pxp+ alpha);
    vzp = -(u1pzp + common_term .* u2pzp - delta_term .* u3pzp+ alpha);
    vxs = -(u1sxp + common_term .* u2sxp - delta_term .* u3sxp+ alpha); % 使用临时变量
    vzs = -(u1szp + common_term .* u2szp - delta_term .* u3szp+ alpha); % 使用临时变量
    vxp(isnan(vxp)) = 0;
    vzp(isnan(vzp)) = 0;
    vxs(isnan(vxs)) = 0;
    vzs(isnan(vzs)) = 0;
end

