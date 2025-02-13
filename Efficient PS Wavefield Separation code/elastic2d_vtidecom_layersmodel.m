
%波前相位方向的各向异性波场分解
clear;clc,

%% 模拟参数
tic;
nt = 1501; % 时间采样点数
dt = 1e-3; % 时间步长
fpeak = 13; % 子波主频
stype = 1; % 加载震源类型
t = (0:nt - 1) * dt;

%建立速度模型
nx=1000;
nz=300;
% nx=600;
% nz=600;
nlayer=30; dx=10;dz=10;
src_z = 0; % 震源激发位置索引
src_x = nx/2; % 震源激发位置索引
temp=read_matrix('.\layer model\vp_300x1000_10m.bin',nz,nx);
vp=imgaussfilt(temp,5);%对纵波速度矩阵进行高斯滤波
temp=read_matrix('.\layer model\vs_300x1000_10m.bin',nz,nx);
vs=imgaussfilt(temp,5);%对横波速度矩阵进行高斯滤波
temp=read_matrix('.\layer model\rho_300x1000_10m.bin',nz,nx);
rho=imgaussfilt(temp,5);%对密度矩阵进行高斯滤波
temp=read_matrix('.\layer model\epsilon_300x1000_10m.bin',nz,nx);
epsilon=imgaussfilt(temp,5);%对epsilon矩阵进行高斯滤波
temp=read_matrix('.\layer model\delta_300x1000_10m.bin',nz,nx);
delta=imgaussfilt(temp,5);%对delta矩阵进行高斯滤波
% temp=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\delta2_0.1.bin',nz,nx);
% delta=imgaussfilt(temp,5);%对delta矩阵进行高斯滤波
% temp=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\vp_0.1.bin',nz,nx);
% vp=imgaussfilt(temp,5);%对纵波速度矩阵进行高斯滤波
% temp=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\vs2_0.1.bin',nz,nx);
% vs=imgaussfilt(temp,5);%对横波速度矩阵进行高斯滤波
% temp=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\epsilon_0.1.bin',nz,nx);
% epsilon=imgaussfilt(temp,5);%对epsilon矩阵进行高斯滤波
% temp=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\rop_0.1.bin',nz,nx);
% rho=imgaussfilt(temp,5);%对密度矩阵进行高斯滤波




%
% load('Hess_2DA_VTI.mat')
% nx=800;nz=300; nlayer=30; dx=10;dz=10;

% nt = 4001; % 时间采样点数
% dx = 10; % X方向空间步长
% dz = 10; % Z方向空间步长
% dt = 1e-3; % 时间步长
% nx=2206;nz=915;
% src_z = ceil(nz/2); % 震源激发位置索引
% src_x = ceil(nx/2); % 震源激发位置索引
% fpeak = 10; % 子波主频15
% stype = 1; % 加载震源类型
% %nlayer = 30; % 吸收边界厚度
% t = (0:nt - 1) * dt;
% % [nz,nx]=size(vp);
% % %

% nlayer=30;
% % dx=10;dz=10;
%  src_z =0; %ceil(nz/2); % 震源激发位置索引
% src_x = ceil(nx/2); % 震源激发位置索引

% temp=read_matrix("D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Hess model\hess_vp_10m.bin",915,2206);
% vp=imgaussfilt(temp,5);
% temp=read_matrix("D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Hess model\hess_vs_10m.bin",915,2206);
% vs=imgaussfilt(temp,5);
% temp=read_matrix("D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Hess model\hess_rho_10m.bin",915,2206);
% rho=imgaussfilt(temp,5);
% temp=read_matrix("D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Hess model\hess_epsilon_10m.bin",915,2206);
% epsilon=imgaussfilt(temp,5);
% temp=read_matrix("D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Hess model\hess_delta_10m.bin",915,2206);
% delta=imgaussfilt(temp,5);


%%
% %Marmousi-I模型
%   nlayer=30; dx=10;dz=10;
%  nx=800;nz=300;
% src_z = 1; % 震源激发位置索引
% src_x = ceil(nx/2); % 震源激发位置索引
% vp=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Marmousi\vp_marmousi-i_306_1611_10m.bin',306,1611);
% vp=vp(51:300,501:1300);
% vp=imgaussfilt(vp,10);
% 
% vs=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\\Marmousi\vs_marmousi-i_306_1611_10m.bin',306,1611);
% vs=vs(51:300,501:1300);
% vs=imgaussfilt(vs,10);
% 
% rho=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\\Marmousi\rho_marmousi-i_306_1611_10m.bin',306,1611);
% rho=rho(51:300,501:1300);
% rho=imgaussfilt(rho,10);
% 
% 
% epsilon=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\\Marmousi\epsilon_marmousi-i_306_1611_10m.bin',306,1611);
% epsilon=epsilon(51:300,501:1300);
% epsilon=imgaussfilt(epsilon,10);
% 
% delta=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\\Marmousi\delta_marmousi-i_306_1611_10m.bin',306,1611);
% delta=delta(51:300,501:1300);
% delta=imgaussfilt(delta,10);


%% 模型边界扩展
NZ = nz + 2 * nlayer;
NX = nx + 2 * nlayer;

[vp] = modpad2d(vp, nlayer, NZ, NX);
[vs] = modpad2d(vs, nlayer, NZ, NX);
[rho] = modpad2d(rho, nlayer, NZ, NX);
[epsilon] = modpad2d(epsilon, nlayer, NZ, NX);
[delta] = modpad2d(delta, nlayer, NZ, NX);


%% 震源处理
% 雷克子波
t_delay = 1.5 / fpeak; % 子波延迟时间
src = (1 - 2 * (pi * fpeak * (t - t_delay)).^2) .* exp(-(pi * fpeak * (t - t_delay)).^2); % 震源子波
dt=0.001;
fmax=1/(2*dt);
df=fmax/round(nt/2);
f=df:df:fmax;
w=2*pi*f;
src_fft=fft(src);
src_fft2=src_fft(1:round(nt/2))./(-w.*w);
src2=real(ifft(src_fft2,nt));


src_z = src_z + nlayer;
src_x = src_x + nlayer;

%% 差分系数
N = 5;
invdx = 1 / dx;
invdz = 1 / dz;
% 经典差分系数
% c1 = 1.211242675807458;
% c2 = -0.089721679689456;
% c3 = 0.013842773437802;
% c4 = -0.001765659877271;
% c5 = 1.186794704887013e-04;
% 优化差分系数
c1 = 1.236425;
c2 = -0.1081130;
c3 = 0.02339911;
c4 = -0.5061550e-2;
c5 = 0.7054313e-3;




%% 弹性参数
invrhox = zeros(NZ, NX);
invrhoz = zeros(NZ, NX);
c11 = rho.*(1+2.*epsilon).*vp.*vp;
c13 = rho.*sqrt((vp.*vp-vs.*vs).*((1+2.*delta).*vp.*vp-vs.*vs))-rho.*vs.*vs;
c33 = rho.*vp.*vp;
c44 = rho.*vs.*vs;


for ix = 1:NX
    for iz = 1:NZ
        if ix < NX
            invrhox(iz, ix) = 2.0 / (rho(iz, ix) + rho(iz, ix+1));
        end
        if iz < NZ
            invrhoz(iz, ix) = 2.0 / (rho(iz, ix) + rho(iz+1, ix));
        end
        if iz < NZ && ix < NX
            c44(iz, ix) = 4 / (1.0 / c44(iz, ix) + ...
                1.0 / c44(iz, ix+1) + ...
                1.0 / c44(iz+1, ix) + ...
                1.0 / c44(iz+1, ix+1));
        end
    end
end
%% 定义需要存储的波场快照和地震记录
% Vx=zeros(nz,nx,nt);Vz=zeros(nz,nx,nt);Vxp=zeros(nz,nx,nt);
% Vzp=zeros(nz,nx,nt);Vxs=zeros(nz,nx,nt);Vzs=zeros(nz,nx,nt);

Srvx=zeros(nt,nx);Srvz=zeros(nt,nx);Srvxp=zeros(nt,nx);
Srvzp=zeros(nt,nx);Srvxs=zeros(nt,nx);Srvzs=zeros(nt,nx);

%% 吸收边界
% EAL边界条件
% [D] = eal2d(vp,nlayer,dx,dz,NX,NZ);
% 以下三行为MEAL边界条件，如果使用MEAL，则注释上一行，激活下8行。
[D, Beta] = meal2d(vp, nlayer, dx, dz, fpeak, NX, NZ);
D = D ./ Beta;
invrhox = invrhox ./ Beta;
invrhoz = invrhoz ./ Beta;
c11 = c11 ./ Beta;
c13 = c13 ./ Beta;
c33 = c33 ./ Beta;
c44 = c44 ./ Beta;
coef1 = (1 - (D * dt) / 2) ./ (1 + D * dt / 2);
coef2 = dt ./ (1 + D * dt / 2);

%% 初始条件
vx = zeros(NZ, NX);
vz = zeros(NZ, NX);
tau_xx = zeros(NZ, NX);
tau_zz = zeros(NZ, NX);
tau_xz = zeros(NZ, NX);



for it = 1:nt
    disp(it);
    % 更新质点振动速度场
    for ix = N + 1:NX - N
        for iz = N + 1:NZ - N
            dtau_xxdx = c1 * (tau_xx(iz, ix+1) - tau_xx(iz, ix)) + ...
                c2 * (tau_xx(iz, ix+2) - tau_xx(iz, ix-1)) + ...
                c3 * (tau_xx(iz, ix+3) - tau_xx(iz, ix-2)) + ...
                c4 * (tau_xx(iz, ix+4) - tau_xx(iz, ix-3)) + ...
                c5 * (tau_xx(iz, ix+5) - tau_xx(iz, ix-4));
            dtau_zzdz = c1 * (tau_zz(iz+1, ix) - tau_zz(iz, ix)) + ...
                c2 * (tau_zz(iz+2, ix) - tau_zz(iz-1, ix)) + ...
                c3 * (tau_zz(iz+3, ix) - tau_zz(iz-2, ix)) + ...
                c4 * (tau_zz(iz+4, ix) - tau_zz(iz-3, ix)) + ...
                c5 * (tau_zz(iz+5, ix) - tau_zz(iz-4, ix));
            dtau_xzdx = c1 * (tau_xz(iz, ix) - tau_xz(iz, ix-1)) + ...
                c2 * (tau_xz(iz, ix+1) - tau_xz(iz, ix-2)) + ...
                c3 * (tau_xz(iz, ix+2) - tau_xz(iz, ix-3)) + ...
                c4 * (tau_xz(iz, ix+3) - tau_xz(iz, ix-4)) + ...
                c5 * (tau_xz(iz, ix+4) - tau_xz(iz, ix-5));
            dtau_xzdz = c1 * (tau_xz(iz, ix) - tau_xz(iz-1, ix)) + ...
                c2 * (tau_xz(iz+1, ix) - tau_xz(iz-2, ix)) + ...
                c3 * (tau_xz(iz+2, ix) - tau_xz(iz-3, ix)) + ...
                c4 * (tau_xz(iz+3, ix) - tau_xz(iz-4, ix)) + ...
                c5 * (tau_xz(iz+4, ix) - tau_xz(iz-5, ix));
            dtau_xxdx = dtau_xxdx * invdx;
            dtau_zzdz = dtau_zzdz * invdz;
            dtau_xzdx = dtau_xzdx * invdx;
            dtau_xzdz = dtau_xzdz * invdz;
            vx(iz, ix) = coef1(iz, ix) * vx(iz, ix) + coef2(iz, ix) * ...
                invrhox(iz, ix) * (dtau_xxdx + dtau_xzdz);
            vz(iz, ix) = coef1(iz, ix) * vz(iz, ix) + coef2(iz, ix) * ...
                invrhoz(iz, ix) * (dtau_xzdx + dtau_zzdz);
        end
    end
    % 震源加载
    vz(src_z, src_x) = vz(src_z, src_x) + src(it) * dt * invdx * invdz;
    %更新应力场
    for ix = N + 1:NX - N
        for iz = N + 1:NZ - N
            dvxdx = c1 * (vx(iz, ix) - vx(iz, ix-1)) + ...
                c2 * (vx(iz, ix+1) - vx(iz, ix-2)) + ...
                c3 * (vx(iz, ix+2) - vx(iz, ix-3)) + ...
                c4 * (vx(iz, ix+3) - vx(iz, ix-4)) + ...
                c5 * (vx(iz, ix+4) - vx(iz, ix-5));
            dvxdz = c1 * (vx(iz+1, ix) - vx(iz, ix)) + ...
                c2 * (vx(iz+2, ix) - vx(iz-1, ix)) + ...
                c3 * (vx(iz+3, ix) - vx(iz-2, ix)) + ...
                c4 * (vx(iz+4, ix) - vx(iz-3, ix)) + ...
                c5 * (vx(iz+5, ix) - vx(iz-4, ix));
            dvzdx = c1 * (vz(iz, ix+1) - vz(iz, ix)) + ...
                c2 * (vz(iz, ix+2) - vz(iz, ix-1)) + ...
                c3 * (vz(iz, ix+3) - vz(iz, ix-2)) + ...
                c4 * (vz(iz, ix+4) - vz(iz, ix-3)) + ...
                c5 * (vz(iz, ix+5) - vz(iz, ix-4));
            dvzdz = c1 * (vz(iz, ix) - vz(iz-1, ix)) + ...
                c2 * (vz(iz+1, ix) - vz(iz-2, ix)) + ...
                c3 * (vz(iz+2, ix) - vz(iz-3, ix)) + ...
                c4 * (vz(iz+3, ix) - vz(iz-4, ix)) + ...
                c5 * (vz(iz+4, ix) - vz(iz-5, ix));
            dvxdx = dvxdx * invdx;
            dvxdz = dvxdz * invdz;
            dvzdx = dvzdx * invdx;
            dvzdz = dvzdz * invdz;
            tau_xx(iz, ix) = coef1(iz, ix) * tau_xx(iz, ix) + coef2(iz, ix) * ...
                (c11(iz, ix) * dvxdx + c13(iz, ix) * dvzdz);
            tau_zz(iz, ix) = coef1(iz, ix) * tau_zz(iz, ix) + coef2(iz, ix) * ...
                (c13(iz, ix) * dvxdx + c33(iz, ix) * dvzdz);
            tau_xz(iz, ix) = coef1(iz, ix) * tau_xz(iz, ix) + coef2(iz, ix) * ...
                c44(iz, ix) * (dvzdx + dvxdz);
        end
    end


    beta=1.9;num=60;

    if(mod(it-1,20)==0)

    [vx1,vz1,vxp,vzp,vxs,vzs] = ani_decomposition_wavefront_phase(vx,vz,vp,vs,epsilon,delta,dz,dx);

     save(sprintf('./ZhiqingyouHessx/Hess_time%d.mat', it),'vx1','vz1','vxp','vzp','vxs','vzs');       

    end

%     Vx(:,:,it)=vx1(nlayer+1:NZ-nlayer,nlayer+1:NX-nlayer);
%     Vz(:,:,it)=vz1(nlayer+1:NZ-nlayer,nlayer+1:NX-nlayer);
%     Vxp(:,:,it)=vxp(nlayer+1:NZ-nlayer,nlayer+1:NX-nlayer);
%     
%     Vzp(:,:,it)=vzp(nlayer+1:NZ-nlayer,nlayer+1:NX-nlayer);
%     Vxs(:,:,it)=vxs(nlayer+1:NZ-nlayer,nlayer+1:NX-nlayer);
%     Vzs(:,:,it)=vzs(nlayer+1:NZ-nlayer,nlayer+1:NX-nlayer);

   Srvx(it,:)=vx1(nlayer+1,nlayer+1:NX-nlayer); 
Srvz(it,:)=vz1(nlayer+1,nlayer+1:NX-nlayer); 
Srvxp(it,:)=vxp(nlayer+1,nlayer+1:NX-nlayer); 
Srvzp(it,:)=vzp(nlayer+1,nlayer+1:NX-nlayer); 
Srvxs(it,:)=vxs(nlayer+1,nlayer+1:NX-nlayer); 
Srvzs(it,:)=vzs(nlayer+1,nlayer+1:NX-nlayer);
% 检查GPU设备
% gpuDev = gpuDevice; disp(gpuDev);
% % 创建两个大的随机矩阵并转移到GPU
% A = gpuArray.rand(1000, 1000, 'single');
% B = gpuArray.rand(1000, 1000, 'single');
% % 使用tic和toc测量时间，但注意GPU计算的异步性 
% tic;
% C = A * B; % 在GPU上执行矩阵乘法 
% wait(gpuDevice);
% % 确保所有GPU操作完成（对于某些MATLAB版本可能不需要）
% elapsedTime = toc;
% % 显示结果 
% disp(['GPU计算时间: ', num2str(elapsedTime), ' 秒']); 
% % 如果需要将结果传回
% % CPU C_cpu = gather(C);

end
toc;

