
%波前相位方向的各向异性波场分解
clear;clc,

%% 模拟参数
nx = 1000; % X方向采样点数
nz = 300; % Z方向采样点数
nt = 601; % 时间采样点数
dx = 10; % X方向空间步长
dz = 10; % Z方向空间步长
dt = 1e-3; % 时间步长
src_z = 300; % 震源激发位置索引
src_x = 300; % 震源激发位置索引
fpeak = 13; % 子波主频15
stype = 1; % 加载震源类型
nlayer = 30; % 吸收边界厚度
t = (0:nt - 1) * dt;



%% 建立速度模型
vp = zeros(nz, nx);
vs = zeros(nz, nx);
rho = zeros(nz, nx);
epsilon = zeros(nz, nx);
delta = zeros(nz, nx);

for ix = 1:nx
    for iz = 1:nz
        vp(iz, ix) = 3000;
        vs(iz, ix) = 1500;
        rho(iz, ix) = 1000;
        epsilon(iz, ix) =0.2;
        delta(iz, ix) =0.1;
    end
end

% nx=2206;nz=915; nlayer=30; dx=10;dz=10;
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

% 
% 
%Marmousi-I模型
%  nx=800;nz=300; nlayer=30; dx=10;dz=10;
% 
% 
% vp=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Marmousi\vp_marmousi-i_306_1611_10m.bin',306,1611);
% vp=vp(1:300,501:1300);
% vp=imgaussfilt(vp,10);
% vs=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Marmousi\vs_marmousi-i_306_1611_10m.bin',306,1611);
% vs=vs(1:300,501:1300);
% vs=imgaussfilt(vs,10);
% 
% rho=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Marmousi\rho_marmousi-i_306_1611_10m.bin',306,1611);
% rho=rho(1:300,501:1300);
% rho=imgaussfilt(rho,10);
% 
% epsilon=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Marmousi\epsilon_marmousi-i_306_1611_10m.bin',306,1611);
% epsilon=epsilon(1:300,501:1300);
% epsilon=imgaussfilt(epsilon,10);
% % 
% delta=read_matrix('D:\数学地物场波模拟\ps波解耦技术相关\波前相位方向波场分解_智1\波前相位方向波场分解_智\Marmousi\delta_marmousi-i_306_1611_10m.bin',306,1611);
% delta=delta(1:300,501:1300);
% delta=imgaussfilt(delta,10);

% 
% %% 模型边界扩展
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
    % 更新应力场
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

    if(mod(it-1,100)==0)

        [vx1,vz1,vxp,vzp,vxs,vzs] = ani_decomposition_wavefront_phase(vx,vz,vp,vs,epsilon,delta,dz,dx,beta,num);

    end
% 

%     disp(it);
%     % 波场快照
%     if(mod(it,100)==0)
% 
%         [colormax0,colormin0] = perclip([vx1;vz1;vxs;vxp;vzs;vzp],0.004);
% 
%         x = (0:nx - 1) * dx;
%         z = (0:nz - 1) * dz;
% 
%         figure(1);imagesc(x,z,vx);shading interp;sucolormap(gca,0);
%         clim([colormin0, colormax0]);
%         set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
%         set(gca,'FontName','Times New Roman','FontSize',10);
%         xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
%         ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
%         title('Horizontal','FontName','Times New Roman','FontSize',11);
% 
%         figure(2);imagesc(x,z,vxp);shading interp;sucolormap(gca,0);
%         clim([colormin0, colormax0]);
%         set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
%         set(gca,'FontName','Times New Roman','FontSize',10);
%         xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
%         ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
%         title('P Horizontal','FontName','Times New Roman','FontSize',11);
% 
% 
%         figure(3);imagesc(x,z,vxs);shading interp;sucolormap(gca,0);
%         clim([colormin0, colormax0]);
%         set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
%         set(gca,'FontName','Times New Roman','FontSize',10);
%         xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
%         ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
%         title('S Horizontal','FontName','Times New Roman','FontSize',11);
% 
%         %%
%         figure(4);imagesc(x,z,vz);shading interp;sucolormap(gca,0);
%         clim([colormin0, colormax0]);
%         set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
%         set(gca,'FontName','Times New Roman','FontSize',10);
%         xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
%         ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
%         title('Vertical','FontName','Times New Roman','FontSize',11);
% 
%         figure(5);imagesc(x,z,vzp);shading interp;sucolormap(gca,0);
%         clim([colormin0, colormax0]);
%         set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
%         set(gca,'FontName','Times New Roman','FontSize',10);
%         xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
%         ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
%         title('P Vertical','FontName','Times New Roman','FontSize',11);
% 
% 
%         figure(6);imagesc(x,z,vzs);shading interp;sucolormap(gca,0);
%         clim([colormin0, colormax0]);
%         set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
%         set(gca,'FontName','Times New Roman','FontSize',10);
%         xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
%         ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
%         title('S Vertical','FontName','Times New Roman','FontSize',11);
% 
%     end
end