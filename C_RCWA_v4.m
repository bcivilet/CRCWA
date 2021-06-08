%% C-RCWA Method version 0.4
% Author: Benjamin Civiletti
% Last edited 2/23/21
% Changelog:


%clear all; close all

%% Physical constants
% Only edit the wavelength lambda0 and incidence-angle theta

epsilon0 = 8.854*10^-12;
mu0 = 4*pi*10^-7;
eta0 = sqrt(mu0/epsilon0);
lambda0 = 600; 
k0 = 2*pi/lambda0;  
omega = k0/eta0;
theta = 0;        

%% Parameters
% M       - Fourier mode parameter
% N       - total number of layers
% L       - period of grating
% h       - discretization parameter in x2
% x1s     - discretization parameter in x1
% lthick  - (N x 1) array of layer thicknesses
% ep      - (N x 1) cell array of relative epsilon
% g       - function handle for grating 
% gprime  - function handle for the derivative
% sampleP - number of subintervals for simpson's rule (even)
% p       - polarization state (either 'p' or 's')
% method  - chooses the solver (either 'C','R' or 'CR')

M = 10; 
m = 2*M+1;  
N = 2;
L = 500;
h = 1; 
x2s = h*ones(1,N); 
x1s = 1;
lthick = [700,700];
ep{1} = 1+1i*10^-9;
ep{2} = -15+1i*4;
%ep{2} = 1+1i*10^-9;
g = @(x1) 100*cos(2*pi*x1/L); 
gprime = @(x1) -200*pi/L*sin(2*pi*x1/L);
sampleP = 500;
p = 'p';
method = 'C';

%triangular grating
%g = @(x1) 4*((x1>0 & x1<250).*(50-x1/5) + (x1<=0 & x1>=-250).*(x1/5+50));
%gprime = @(x1) 4*((x1>0 & x1<250).*(-1/5) + (x1<=0 & x1>=-250).*(1/5));


%% Define the coordinate transform

H = 700;
S = @(x2) (x2-H>0).*(1 - 3/(H^2).*(x2-H).^2+(2/H^3).*(x2-H).^3)+...
    (x2-H<=0).*(1 - 3/(H^2).*(x2-H).^2-(2/H^3).*(x2-H).^3);

Sprime = @(x2) -(x2-H>0).*(-6/(H^2).*(x2-H)+(6/H^3).*(x2-H).^2)-...
    (x2-H<=0).*(-6/(H^2).*(x2-H)-(6/H^3).*(x2-H).^2);

detDG = @(x1,x2) Sprime(x2).*g(x1)+1;

a11 = @(x1,x2) abs(detDG(x1,x2));
a21 = @(x1,x2) -sign(detDG(x1,x2)).*S(x2).*gprime(x1);
a12 = @(x1,x2) a21(x1,x2);
a22 = @(x1,x2) ((S(x2)*gprime(x1)).^2+1)./abs(detDG(x1,x2));



%% Calculate grid points

[x1grid,x2grid,x2bound] = gridC(N,L,lthick,x1s,x2s);
x2gridn = size(x2grid,2);
x1gridn = size(x1grid,2);

%% Plot the grating profile

figure(5)
plot(x1grid,g(x1grid),'LineWidth',2)
axis equal

gR = @(x1) -g(x1) + H;

%% Define the boundary conditions

k = k0*sin(theta) + 2*pi*(-M:M)/L;
K = diag(k);
beta = sqrt(k0^2-k.^2);

O = zeros(m,m);
I = eye(m,m);

Ye_p = [O -diag(beta/k0) ; I  O];
Ye_m = [O  diag(beta/k0) ; I  O];
Yh_p = [-diag(beta/k0) O ; O -I];
Yh_m = [diag(beta/k0)  O ; O -I];


%% Set the polarization state

A = zeros(4*M+2,1);
if p == 's'     %s polarization
    A(M+1) = 1; 
elseif p == 'p' %p polarization
    A(3*M+2) = 1; 
end

%% Preallocate some cells for the RCWA algorithm

Zs = cell(1,x2gridn+1);
Ts = cell(1,x2gridn+1);
Wus = cell(1,x2gridn);
Dus = cell(1,x2gridn);
f = cell(1,x2gridn+1);

Zs{x2gridn+1} = [Ye_p;Yh_p]; 

%% Start the RCWA algorithm 

for i = x2gridn:-1:1
    %finds the current layer
    check = (x2grid(i)-x2bound > 0); layer = find(check,1,'last');
    if isa(ep{layer},'function_handle')
        epsilon = ep{layer}; 
        T = ToeplitzM(m,epsilon,L,sampleP);
    else
        epsilon = ones(1,sampleP+1)*ep{layer};
        T = eye(m,m)*ep{layer};
    end
    
    if method == 'C' %C-RCWA method
        
        %define epsilon and mu
        ep11 = @(x1) ep{layer}*a11(x1,x2grid(i));
        ep21 = @(x1) ep{layer}*a21(x1,x2grid(i));
        ep12 = @(x1) ep{layer}*a12(x1,x2grid(i));
        ep22 = @(x1) ep{layer}*a22(x1,x2grid(i));
        ep33 = @(x1) ep{layer}*a11(x1,x2grid(i));
        
        mu11 = @(x1) a11(x1,x2grid(i));
        mu21 = @(x1) a21(x1,x2grid(i));
        mu12 = @(x1) a12(x1,x2grid(i));
        mu22 = @(x1) a22(x1,x2grid(i));
        mu33 = @(x1) a11(x1,x2grid(i));
        
        Tep11 = ToeplitzM(m,ep11,L,sampleP);
        Tep21 = ToeplitzM(m,ep21,L,sampleP);
        Tep12 = Tep21;
        Tep22 = ToeplitzM(m,ep22,L,sampleP);
        Tep33 = ToeplitzM(m,ep33,L,sampleP);
        
        Tmu11 = (ep{layer}^-1)*Tep11;
        Tmu21 = (ep{layer}^-1)*Tep21;
        Tmu12 = Tmu21;
        Tmu22 = (ep{layer}^-1)*Tep22;
        Tmu33 = (ep{layer}^-1)*Tep33;
        
        P11 = K*(Tep22\Tep21);
        P14 = k0*Tmu33-(1/k0)*K*(Tep22\K);
        P22 = Tmu12*(Tmu22\K);
        P23 = -k0*(Tmu11-Tmu12*(Tmu22\Tmu21));
        P32 = -k0*Tep33+(1/k0)*K*(Tmu22\K);
        P33 = K*(Tmu22\Tmu21);
        P41 = -k0*(Tep12*(Tep22\Tep21)-Tep11);
        P44 = Tep12*(Tep22\K);
        
        
        P = [P11 O O P14;...
            O P22 P23 O;...
            O P32 P33 O;...
            P41 O O P44];
        
    elseif method == 'R' %RCWA method
        
        ind = find(gR(x1grid) < x2grid(i));
        
        if size(ind,2) == size(x1grid,2) %below grating
            ep11 = @(x1) ep{layer};
            ep22 = @(x1) ep11(x1);
            ep33 = @(x1) ep11(x1);
        elseif (size(ind,2) ~= size(x1grid,2)) && isempty(ind) == 0 %g
            x1l = x1grid(ind(1));
            x1r = x1grid(ind(end));
            ep11 = @(x1) (x1 < x1l).*ep{1} + (x1 > x1r).*ep{1} +...
                (x1 >= x1l & x1 <= x1r).*ep{2};
            ep22 = @(x1) ep11(x1);
            ep33 = @(x1) ep11(x1);
        elseif isempty(ind) == 1 %above grating
            ep11 = @(x1) ep{layer};
            ep22 = @(x1) ep11(x1);
            ep33 = @(x1) ep11(x1);
        end
        
        Tep11 = ToeplitzM(m,ep11,L,sampleP);
        Tep21 = O;
        Tep12 = O;
        Tep22 = Tep11;
        Tep33 = Tep11;
        
        Tmu11 = I;
        Tmu21 = O;
        Tmu12 = O;
        Tmu22 = I;
        Tmu33 = I;
        
        P14 = k0*Tmu33-(1/k0)*K*(Tep22\K);
        P23 = -k0*(Tmu11-Tmu12*(Tmu22\Tmu21));
        P32 = -k0*Tep33+(1/k0)*K*(Tmu22\K);
        P41 = -k0*(Tep12*(Tep22\Tep21)-Tep11);
         
        P = [O O O P14;...
             O O P23 O;...
             O P32 O O;...
             P41 O O O];
        
    end
    
    if i == 1 || i==x2gridn
        Delta = h;
    else
        Delta = x2grid(i)-x2grid(i-1);
    end
    
    [G,D] = eig(P);
    Ddiag = diag(D);
    [~,ind] = sort(imag(Ddiag),'descend');
    Ddiag = Ddiag(ind);
    D = diag(Ddiag);
    G = G(:,ind);
    W = G\Zs{i+1};
    Wu = W(1:4*M+2,:);
    Wus{i} = Wu;
    Wl = W(4*M+3:end,:);
    Du = D(1:4*M+2,1:4*M+2);
    Dus{i} = Du;
    Dl = D(4*M+3:end,4*M+3:end);
    expmat1 = expm(-1i*Delta*Dl);
    expmat2 = expm(1i*Delta*Du);
    Zs{i} = G*[eye(4*M+2,4*M+2); expmat1*Wl*(Wu\expmat2)];
end

%solve for T0
Z0u = Zs{1}(1:4*M+2,:);
Z0l = Zs{1}(4*M+3:end,:);
T0R = ([Z0u -Ye_m;Z0l -Yh_m]\[Ye_p;Yh_p])*A;
Ts{1} = T0R(1:4*M+2,:);

for i=1:x2gridn
    if i==1 || i==x2gridn
        Delta = h;
    else
        Delta = x2grid(i)-x2grid(i-1);
    end
    Ts{i+1}=(Wus{i}\expm(1i*Delta*Dus{i}))*Ts{i};
end

%solve for f

for i=1:x2gridn+1
    f{i} = Zs{i}*Ts{i};
end

%reconstruct the solution and plot
basis = @(x1,n) exp(1i*n*x1);


Ex = zeros(x1gridn,x2gridn);
Ey = zeros(x1gridn,x2gridn);
Hx = zeros(x1gridn,x2gridn);
Hy = zeros(x1gridn,x2gridn);

if method == 'C'
    Exhat = zeros(x1gridn,x2gridn);
    Eyhat = zeros(x1gridn,x2gridn);
    Hxhat = zeros(x1gridn,x2gridn);
    Hyhat = zeros(x1gridn,x2gridn);
else
end

if method == 'C'
    for i=1:x2gridn
        for j=1:x1gridn
            Exhat(j,i) = basis(x1grid(j),k)*f{i+1}(1:2*M+1,1);
            Eyhat(j,i) = basis(x1grid(j),k)*f{i+1}(2*M+2:4*M+2,1);
            Hxhat(j,i) = basis(x1grid(j),k)*f{i+1}(4*M+3:6*M+3,1);
            Hyhat(j,i) = basis(x1grid(j),k)*f{i+1}(6*M+4:end,1);
        end
    end
    

for i = 1:x2gridn
    for j = 1:x1gridn
        [x1hat,x2hat,index]=invmap(x1grid(j),x2grid(i),g,S,x2grid,H);
        Ex(j,i) = Exhat(j,index);
        Ey(j,i) = Eyhat(j,index);
        Hx(j,i) = Hxhat(j,index);
        Hy(j,i) = Hyhat(j,index);    
    end
end
elseif method == 'R'
    for i=1:x2gridn
        for j=1:x1gridn
            Ex(j,i) = basis(x1grid(j),k)*f{i+1}(1:2*M+1,1);
            Ey(j,i) = basis(x1grid(j),k)*f{i+1}(2*M+2:4*M+2,1);
            Hx(j,i) = basis(x1grid(j),k)*f{i+1}(4*M+3:6*M+3,1);
            Hy(j,i) = basis(x1grid(j),k)*f{i+1}(6*M+4:end,1);
        end
    end
end
x2plotgrid = linspace(-H,H,2*H/h);
if p == 'p' && method == 'C'
    
    figure(1)
    tiledlayout(1,2);
    nexttile
    set(0,'defaultAxesFontSize',20)
    pcolor(x1grid/L,x2plotgrid/H,fliplr(abs(Hyhat))');
    cb = colorbar('east','YAxisLocation','Right')
    %cb.Position =cb.Position + [2 0 0 0];
    %set(hc,'YAxisLocation','Right')
    %set(gca,'YDir','reverse');
    axis equal;xlim([-1/2,1/2]);ylim([-1,1])
    shading interp
    colormap jet
    hold on
    plot(x1grid/L,0*ones(x1gridn),'Color','white','LineWidth',2)
    ylabel('$\mathbf{\hat{x}_{2}}/H$','interpreter','latex')
    xlabel('$\mathbf{\hat{x}_{1}}/ \Lambda$','interpreter','latex')
    title('$\mathbf{|\hat{H}_{3}|}$','interpreter','latex',...
        'FontWeight','bold')
    xticks(-1/2:1/2:1/2)
    yticks(-1:1/2:1)
    
    nexttile
    pcolor(x1grid/L,x2plotgrid/H,fliplr(abs(Hy))');
    colorbar('east','YAxisLocation','Right')
    %colorbar
    %set(gca,'YDir','reverse');
    axis equal;xlim([-1/2,1/2]);ylim([-1,1])
    shading interp
    colormap jet
    hold on
    plot(x1grid/L,g(x1grid)/H,'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{2}}/H$','interpreter','latex')
    xlabel('$\mathbf{x_{1}}/ \Lambda$','interpreter','latex')
    title('$\mathbf{\left|H_{3} \right|}$','interpreter','latex')
    xticks(-1/2:1/2:1/2)
    yticks(-1:1/2:1)
    
    x0=600;
    y0=600;
    width=550;
    height=450;
    set(gcf,'position',[x0,y0,width,height])
    
    figure(2)
    tiledlayout(1,2);
    nexttile
    pcolor(x1grid,x2plotgrid,fliplr(abs(Exhat))');
    %colorbar('east','YAxisLocation','Right')
    colorbar
    %set(gca,'YDir','reverse');
    axis equal;xlim([-1/2,1/2]);ylim([-1,1])
    shading interp
    colormap jet
    hold on
    plot(x1grid,0*ones(x1gridn),'Color','white','LineWidth',2)
    ylabel('$\mathbf{\hat{x}_{2}}/H$','interpreter','latex')
    xlabel('$\mathbf{\hat{x}_{1}}/ \Lambda$','interpreter','latex')
    title('$\mathbf{|\hat{E}_{1} |}$','interpreter','latex')
    xticks(-1:1/2:1)
    yticks(-1:1/2:1)
    
    nexttile
    pcolor(x1grid,x2plotgrid,fliplr(abs(Ex))');
    %colorbar('east','YAxisLocation','Right')
    colorbar
    %set(gca,'YDir','reverse');
    axis equal;xlim([-1/2,1/2]);ylim([-1,1])
    shading interp
    colormap jet
    hold on
    plot(x1grid,g(x1grid),'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{2}}/H$','interpreter','latex')
    xlabel('$\mathbf{x_{1}}/ \Lambda$','interpreter','latex')
    title('$\mathbf{\left|E_{1} \right|}$','interpreter','latex')
    xticks(-1:1/2:1)
    yticks(-1:1/2:1)
    
    x0=600;
    y0=600;
    width=550;
    height=450;
    set(gcf,'position',[x0,y0,width,height])
    
    %cb.Layout.Tile = 'east';
elseif p == 's' && method == 'C'
    
    figure(1)
    tiledlayout(2,1);
    nexttile
    set(0,'defaultAxesFontSize',20)
    pcolor(x2grid,x1grid,abs(Eyhat));
    colorbar('east','YAxisLocation','Right')
    %set(gca,'YDir','reverse');
    axis equal;ylim([-L/2,L/2]);xlim([0,2*H])
    %set(gca,'YDir','reverse');
    shading interp
    colormap jet
    hold on
    plot(H*ones(x1gridn),x1grid,'Color','white','LineWidth',2)
    xlabel('$\mathbf{\hat{x}_{2}}$','interpreter','latex')
    ylabel('$\mathbf{\hat{x}_{1}}$','interpreter','latex')
    title('$\left|\hat{E}_{y} \right|$','interpreter','latex')
    
    nexttile
    pcolor(x2grid,x1grid,abs(Ey));
    colorbar('east','YAxisLocation','Right')
    %set(gca,'YDir','reverse');
    axis equal;ylim([-L/2,L/2]);xlim([0,2*H])
    %set(gca,'YDir','reverse');
    shading interp
    colormap jet
    hold on
    plot(-g(x1grid)+H,x1grid,'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{1}}$','interpreter','latex')
    xlabel('$\mathbf{x_{2}}$','interpreter','latex')
    title('$\left|E_{y} \right|$','interpreter','latex')
    
    x0=600;
    y0=600;
    width=530;
    height=450;
    set(gcf,'position',[x0,y0,width,height])
    
    figure(2)
    tiledlayout(2,1);
    nexttile
    pcolor(x2grid,x1grid,abs(Hxhat));
    colorbar('east','YAxisLocation','Right')
    %set(gca,'YDir','reverse');
    axis equal;ylim([-L/2,L/2]);xlim([0,2*H])
    %set(gca,'YDir','reverse');
    shading interp
    colormap jet
    hold on
    plot(H*ones(x1gridn),x1grid,'Color','white','LineWidth',2)
    ylabel('$\mathbf{\hat{x}_{1}}$','interpreter','latex')
    xlabel('$\mathbf{\hat{x}_{2}}$','interpreter','latex')
    title('$\left|\hat{H}_{x} \right|$','interpreter','latex')
    
    nexttile
    pcolor(x2grid,x1grid,abs(Hx));
    colorbar('east','YAxisLocation','Right')
    %set(gca,'YDir','reverse');
    axis equal;xlim([-L/2,L/2]);ylim([0,2*H])
    %set(gca,'YDir','reverse');
    shading interp
    colormap jet
    hold on
    plot(-g(x1grid)+H,x1grid,'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{1}}$','interpreter','latex')
    xlabel('$\mathbf{x_{2}}$','interpreter','latex')
    title('$\left|H_{x} \right|$','interpreter','latex')
    
    x0=600;
    y0=600;
    width=530;
    height=450;
    set(gcf,'position',[x0,y0,width,height])
elseif p == 'p' && method == 'R'
    
    figure(3)
    t = tiledlayout(2,1);
    nexttile
    set(0,'defaultAxesFontSize',20)
    pcolor(x2grid,x1grid,abs(Hy));
    cb1 = colorbar('east','YAxisLocation','Right')
    %set(cb1,'Position',[0.8748 0.6053 0.0151 0.3019]);
    %set(gca,'YDir','reverse');
    axis equal;ylim([-L/2,L/2]);xlim([0,2*H])
    shading interp
    colormap jet
    hold on
    plot(-g(x1grid)+H,x1grid,'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{1}}$','interpreter','latex')
    xlabel('$\mathbf{x_{2}}$','interpreter','latex')
    title('$\left|{H}_{x_{3}} \right|$','interpreter','latex')
      
    nexttile
    pcolor(x2grid,x1grid,abs(Ex));
    cb2 = colorbar('east','YAxisLocation','Right');
    %set(cb2,'Position',[0.8748 0.6053 0.0151 0.3019]);
    %set(gca,'YDir','reverse');
    axis equal;ylim([-L/2,L/2]);xlim([0,2*H])
    shading interp
    hold on
    plot(-g(x1grid)+H,x1grid,'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{1}}$','interpreter','latex')
    xlabel('$\mathbf{x_{2}}$','interpreter','latex')
    title('$\left|E_{x_{1}}^{\phantom{5}} \right|$','interpreter','latex')
    
    x0=600;
    y0=600;
    width=530;
    height=450;
    set(gcf,'position',[x0,y0,width,height])
elseif p == 's' && method == 'R'
    
    figure(3)
    t = tiledlayout(2,1);
    nexttile
    set(0,'defaultAxesFontSize',20)
    colorbar('east','YAxisLocation','Right')
    pcolor(x2grid,x1grid,abs(Ey));
    set(gca,'YDir','reverse');
    axis equal;ylim([-L/2,L/2]);xlim([0,2*H])
    set(gca,'YDir','reverse');
    shading interp
    colormap jet
    hold on
    plot(-g(x1grid)+H,x1grid,'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{1}}$','interpreter','latex')
    xlabel('$\mathbf{x_{2}}$','interpreter','latex')
    title('$\left|E_{y} \right|$','interpreter','latex')
    
    nexttile
    pcolor(x2grid,x1grid,abs(Hx));
    colorbar('east','YAxisLocation','Right')
    %set(gca,'YDir','reverse');
    axis equal;ylim([-L/2,L/2]);xlim([0,2*H])
    %set(gca,'YDir','reverse');
    shading interp
    colormap jet
    hold on
    plot(-g(x1grid)+H,x1grid,'Color','white','LineWidth',2)
    ylabel('$\mathbf{x_{1}}$','interpreter','latex')
    xlabel('$\mathbf{x_{2}}$','interpreter','latex')
    title('$\left|H_{x} \right|$','interpreter','latex')
    
    x0=600;
    y0=600;
    width=530;
    height=450;
    set(gcf,'position',[x0,y0,width,height])
end

%clearvars -except f Zs Ts Wus Dus Delta T0R L M m basis x1gridn x2gridn...
    %Hy Ey x1grid x2grid k Hyhat Hy x1grid g Hyhatrcwa Eyhat Ex Exhat

function [T] = ToeplitzM(m,f,L,sampleP)
%Toeplitz Calculates the mxm Toeplitz matrix from sampled data f or 
%                    function handle f
%            Inputs:  m - odd integer for index 
%                     f - function handle or (1 x sampleP+1) array 
%                     L - width of device
%                    
%            Outputs: T - (m x m) array for the Toeplitz matrix of f
%
fhandle = isa(f,'function_handle');
n = 2*m-1;
fftvec = zeros(1,n);
switch fhandle
    case 1
        for i=1:n
            ind = (n-1)/2;
            j = -ind:ind;
            fftvec(i) = (1/L)*simpson(@(x) f(x)...
                .*exp(-2*pi*1i*j(i)*x/L),-L/2,L/2,sampleP);
        end
    case 0
       for i=1:n
            ind = (n-1)/2;
            j = -ind:ind;
            x = linspace(-L/2,L/2,sampleP+1);
            basis = exp(-2*pi*1i*j(i)*x/L);
            fsample = f.*basis;
            fftvec(i) = (1/L)*simpson(fsample,-L/2,L/2,sampleP);
        end 
end
fftvec0pos = fftvec(m:end);
fftvec0neg = flip(fftvec(1:m));
T = toeplitz(fftvec0pos,fftvec0neg);
end

function [x1grid,x2grid,x2bound] = gridC(N,L,lthick,x1s,x2s)
%gridC Calculate the grid points
%   Inputs:         N        - total number of layers
%                   x1s      - discretization size
%                   x2s      - discretization size in each layer
%                   lthick   - thicknesses of the layers
%                   H        - half the thickness in x2-direction
%                   h        - discretization parameter in x2-direction

x1n = L/x1s;
x1grid = linspace(-L/2 + x1s/2, L/2 - x1s/2,x1n);
x2slice = lthick./x2s; %number of slices per layer
x2grid = [];
x2bound = [];
lthick = [0 lthick];

for i = 1:N
    x2cslice = x2slice(i); %number of slices in ith layer
    sbounds = [sum(lthick(1:i)) sum(lthick(1:i+1))]; %endpoints of layer i
    x2grid =  [x2grid linspace(sbounds(1) + x2s(i)/2,sbounds(2) - ...
        x2s(i)/2,x2cslice)];
    x2bound = [x2bound sbounds(1)];
end
end

function [approxInt] = simpson(f,a,b,n)
%Simpson Composite simpson's method for approximating the integral of f(x)
%   Inputs:  f - function handle for f(x) or sampled function
%            a - left end-point of interval
%            b - right end-point of interval
%            n - number of subintervals (has to be even)
fhandle = isa(f,'function_handle');
h = (b-a)/n;
switch fhandle
    case 1
        x = linspace(a,b,n+1);
        y = f(x);
    case 0
        y = f;
end
approxInt=(h/3)*(y(1)+y(n+1)+4*sum(y(2:2:n))+2*sum(y(3:2:n-1)));
end

function [x1hat, x2hat,index] = invmap(x1,x2,g,S,x2grid,H)
    x1hat = x1;
    C = g(x1);
    F = @(x) C*S(x)-x+x2;
    x2hat = bisection(F,0,2*H,10^-6);
    [~,index] = min(abs(x2hat-x2grid));
end

function [r,n] = bisection(f,a,b,TOL)
%   bisection: computes approximate solution to f(x)=0
%   Inputs:  f   - function handle for f(x)
%            a   - left end-point of interval
%            b   - right end-point of interval
%            TOL - tolerance
%   Outputs: r - the approximate root
%            n - number of iterations to get within tolerance
n = 0;
if f(a)*f(b) >= 0
    error('There is no root of f(x) in [a,b]'); %stop
end

while (b-a)/2 > TOL
    c = (a + b)/2;
    if f(c) == 0
        break;
    end
    
    if f(a)*f(c) < 0 %root is in [a,c]
        b = c;
    else             %root is in [c,b]
        a = c;
    end
    n = n + 1;
end
r = (a + b)/2;
if f(r) >10*TOL
    disp('warning')
end
end

