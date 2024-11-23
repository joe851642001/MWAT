load('b_d5.mat') % b_in = b(:,1) in steady state

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Define Parameters %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdas = 1.064; % seed wavelength [um]
lambdap = 0.975; % pump wavelength[um]
lslp = lambdas/lambdap; % quantum defect  
M = size(b_in,1); % M, total # of modes 
nco = 1.5; % refractive index of core
w = 20; % core radius [um]
W = 200; % cladding radius [um] 
Tb = 0; % temperature at conductive boundaries
L = 1e6; % fiber length = L [um]
Pseed = sum(abs(b_in).^2)*nco*(2*w)^2 ; % seed power [W]
Psat = 100 ; % Gain saturation power [W]
G0 = log(sqrt(exp(3.2*log(9))))/1e6 ; % Field gain [1/um]
tUnit = 6 ; % total simulation time for normal time step [ms]
tFactor = 1 ; % enlarge the time step: dt = tFactor * dt0 
% ----------------------------------------- 
k0 = 2*pi/lambdas ; % total wavenumber
bamp = abs(b_in) ; % field amplitude at input
bph = angle(b_in) ; % field phase at input
%%% Perturbed input
  factor = 1 ; pert = 4 ; % peak-peak power variation: factor*10^(-pert)
  fpert = 1e3 ; % frequency of perturbation [Hz]
   rng('shuffle') ; % to avoid repeated rand() values 
   Pphase = 2*pi*rand(M,1) ; % modulation offsets among modes % Tech Noise Correction
%   Pphase = 2*pi*[0.2;0.9;0.4;0.7;0.3]
    PpertRatio = (bamp/norm(bamp)).^2 ;  % Pert Correction 2021-11-03
%  PpertRatio = [1;0] ;               % Pert Correction 2021-11-03
%%%
%------------------------------------------
b0 = sqrt(Pseed/nco)*(bamp/norm(bamp)).*exp(1i*bph)/(2*w); % field at input, divided by fiber diameter % I/n = E^2 index effect corrected
beta = sqrt((2*pi*nco/lambdas)^2 - (linspace(1,M,M)*pi/(2*w)).^2); % propagation constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Define Transverse (x) Steps & psi %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rcore = 1*M*20+1; % # of x steps in core
x = linspace(-w,w,Rcore).'; % x ticks in core [um]
dx = 2*w/(Rcore-1); % step size in x
wdx = w/dx ;
%%% ------------- Mode profiles in x ---------------- 
psi = zeros(Rcore,M); 
psiN = zeros(Rcore,M); 
%%%%% ----- Normalized so that sum(abs(psiN).^2) = 1
for m = 1:M
    if mod(m,2)==0
        psiN(:,m)=sin(m*pi*x/(2*w));
         psiN(1:floor(Rcore/2),m)=-(fliplr(psiN(ceil(Rcore/2+1):Rcore,m).')).'; % psi(x) = -psi(-x) to reduce numerical error
    else
        psiN(:,m)=cos(m*pi*x/(2*w));
        psiN(1:floor(Rcore/2),m)=(fliplr(psiN(ceil(Rcore/2)+1:Rcore,m).')).'; % psi(x) = psi(-x) to reduce numerical error
    end
    psiN(:,m) = psiN(:,m)/sqrt(sum(abs(psiN(:,m)).^2)); % normalization
end
%%%%% ----- Normalized so that sum(abs(psi).^2*dx) = 2*w (core size in y)
for m = 1:M
    if mod(m,2)==0
        psi(:,m)=sqrt(2)*sin(m*pi*x/(2*w));
         psi(1:floor(Rcore/2),m) = -(fliplr(psi(ceil(Rcore/2+1):Rcore,m).')).'; % psi(x) = -psi(-x) to reduce numerical error
    else
        psi(:,m)=sqrt(2)*cos(m*pi*x/(2*w));
        psi(1:floor(Rcore/2),m)=(fliplr(psi(ceil(Rcore/2)+1:Rcore,m).')).'; % psi(x) = psi(-x) to reduce numerical error
    end
%     psi(:,m) = psi(:,m)/sqrt(sum(abs(psi(:,m)).^2)); % normalization
end
% norm_trapz_psi2 = trapz(x,abs(psi(:,1)).^2); % other normalization approach (not used)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Define Longitudinal (z) Steps %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dz = 10*(2*pi/beta(M)); % step size in z for light [um] 
% dz = 2*pi/abs(beta(1)-beta(M))/40; % step size in z for light [um] 
J = round(L/dz); % J, total # of z steps
dzT = 2*pi/abs(beta(M-1)-beta(M))/20; % step size in z for temperature [um]
N = round(dzT/dz) ; % coarsening: 1 step in T = N steps in light
dzT = dz*N ; % Temperature step size
J = J - mod(J,N); % Number of optical steps
Jred = J/N-1; % Number of temperature steps
z = dz*linspace(0,J,J+1); % z ticks
ik0dz = 1i*k0*dz;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Mode-Overlap Matrix (psi_m*psi_n) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:M
    for n = 1:M
      Dbeta(m,n) = exp(1i*(beta(n)-beta(m))*dz);
    end
end
DbetaCol = reshape(Dbeta,[M^2,1]) ;
for m = 1:M
    for n = m:M
      PsimPsin(:,m+n) = psiN(:,m).*psiN(:,n);
    end
end
for m = 1:M
    for n = 1:M
      PsimPsinFull(m,n,:) = Dbeta(m,n)*(psiN(:,m).*psiN(:,n)).';
    end
end
PsimPsinCol = reshape(PsimPsinFull,[M^2,length(x)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Linear Propagation Matrix (HDHG) %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HD = diag(exp(1i*beta*dz).*ones(1,M)); 
% ----------Gain matrix parameters---------
HG = zeros(M);
gamma = 0;
b = cat(2,b0,zeros(M,Jred+1)); 
G = G0;
HG = diag(ones(1,M))*exp(G0*dz);
HDHG = HD*HG;
for m = 1:M
    for n = 1:M
      PsimPsinFull(m,n,:) = Dbeta(m,n)*(psiN(:,m).*psiN(:,n)).';
    end
end
PsimPsinCol = reshape(PsimPsinFull,[M^2,length(x)]);
for m = 1:length(x)
    HDHGPsi2(:,:,m) = HD*HG*PsimPsinFull(:,:,m);
end
HDHGPsi2Col = reshape(HDHGPsi2,[M^2,length(x)]);

%%%%% When Gain Saturation is ON % HG for Gain Saturation
for m = 1:length(x)
    HDPsi2(:,:,m) = HD*PsimPsinFull(:,:,m);
end
HDPsi2Col = reshape(HDPsi2,[M^2,length(x)]);
for m = 1:M
    for n = 1:M
        if m == n
         PsimPsinFull_Diag(m,n,:) = Dbeta(m,n)*(psiN(:,m).*psiN(:,n)).';
        else
         PsimPsinFull_Diag(m,n,:) = 0 ; 
        end
    end
end
PsimPsinCol_Diag = reshape(PsimPsinFull_Diag,[M^2,length(x)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Define Time (t) Steps %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 1.4/1e6 ; % thermal conductivity [W/(K um)]
rhoC = 1.67e6/1e18 ; % volumetric heat capacity [J/(K(um)^3)], rho: density, C: specific heat capacity
Nclad = round((W - w)/dx); % # of x steps on one side of cladding
W = ((Rcore)/2+Nclad)*dx; % redefine cladding radius
X = linspace(-W,W,Rcore+2*Nclad).'; % transverse ticks from air/clad to clad/air 
nx = Rcore+2*Nclad; % total # of transverse steps (X) in fiber 
dt = tFactor*rhoC*dx^2/(2*K)  % time step size dt [s] that satisfies CFL condition 
t_tot = tFactor*tUnit*1e-3; % total simulation time [s] 
P = round(t_tot/dt); % total # of time steps 
tplot = linspace(0,dt*P,P); % time ticks for plots
ax = zeros(Rcore,J+1); % empty intensity vector
%%% Perturbed input
P0 = factor*10^(-pert)*PpertRatio ; %(bamp/norm(bamp)).^2 ; % amplitude of perturbation in power % Tech Noise Correction % Pert Correction 2021-11-03 % Pert phase 2023-05-14
PI = 0.5*P0.*sin(2*pi*repmat(tplot,M,1)*fpert + pi*repmat(Pphase,1,P)) ; % Tech Noise Correction
bamp_pert = sqrt(repmat((bamp/norm(bamp)).^2,1,P) + PI) ;
b0 = sqrt(Pseed/nco)*bamp_pert.*exp(1i*bph)/(2*w); % field at input, divided by fiber diameter % I/n = E^2 index effect corrected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Define Temperature (T) Matrix Parameters %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% N defined previously in z-step section 
Nhalf = floor(N/2); % half step size for coarsening use
Znum_red = floor(J/N); %  total # of z steps in T
Qcore_red = zeros(length(x),Znum_red+1); % empty vector of reduced heat source
for m = 1:Znum_red-1 
       z_red(m+1) = z(N*m+1); % generating reduced z ticks for T
end
z_red(1) = z(1); z_red(end+1) = z(end); % generating reduced z ticks for T
[GridXT,GridZT] = ndgrid(1:length(x),z_red); % ndgrid of (x,z) in T (for griddedInterpolant)
[GridXL,GridZL] = ndgrid(1:length(x),z(1:J+1)); % ndgrid of (x,z) in light (for griddedInterpolant)
[GridXred,GridZred] = ndgrid(1:length(x),z_red); % ndgrid of (x,z) in T (for griddedInterpolant)
[GridXQ,GridZQ] = ndgrid(1:length(x),z(1:J+1)); % ndgrid of (x,z) in light (for griddedInterpolant)
eta = 3.5e-5/(2*nco); % Thermo-optic coeff [1/K] % Index effect corrected
hq = 1000/1e12 ; % Convection coeff [W/(K(um)^2)] (not used)
%%%%% notes: K, rho, C are defined in t-step section, K/rhoC [um^2/s],Q/rhoC [K(um)^2/J] 
%%%%% Constant coeff of tridiagonal system 
c = dt*K/(rhoC*dx^2) ; %Superdiagonal: coeff of T(r+1) % b = c ; %Subdiagonal: coeff of T(r-1)
a = 2*(1+c) ; %Main diagonal: coeff of T(r)
cminus = 2*(1-c);
Qcoef = 2*dt/rhoC; % coefficient for the Q term It should be 2*dt/rhoC
coreUp = Nclad+1 ; coreLow = Nclad+Rcore ; cladLow = length(X)-1;
Q = zeros(length(X),Znum_red+1) ;
d = NaN(1,nx-2); % for CN use
cVec = -c*ones(nx-2,1);
aVec = a*ones(nx-2,1);
TriD = spdiags([cVec aVec cVec],-1:1,nx-2,nx-2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Recording-Vector Generation %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Perturb = zeros(length(x),J+1);
HDGTcol = zeros(M^2,J+1);
HDGT = zeros(M,M,J+1);
H0 = diag(ones(1,M)); 
Nplusone = N+1; 
ax = zeros(Rcore,Jred+2); 
Is = Psat/(2*w)^2 ; % Gain Saturation % Pseed = sumsqr(abs(b0*2*w))*nco
Znumber = 100 ; % MultiZ Multiple Z
Zspacing = floor(Jred/Znumber) ; % MultiZ Multiple Z
% Zspacing = 20 ; % MultiZ Multiple Z
b_rec = NaN(M,Znumber+1,floor(P/100)); % MultiZ Multiple Z
T_rec = NaN(1,Jred+2,floor(P/100)); % empty vector of fiber-end T 
xT_rec = floor(Rcore/4); 

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% S I M U L A T I O N %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 load('DnT.mat')

tic
for idxt = 2:P+1 % time loop

    %%% Perturbed input
    b(:,1) = b0(:,idxt-1) ;
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%% Thermo-Optic matrix %%%
        Perturb = exp(ik0dz*Dn); % perturb calc by overlap integral
% %        % HTnew = PsimPsinCol*Perturb ;  % HT gen (HDHG*HT)      
% %        % HT = reshape(HTnew,[M,M,J+1]); % HT gen (HDHG*HT)
%         HDGTcol = HDHGPsi2Col*Perturb ;  % HDGT gen        
%         HDGT = reshape(HDGTcol,[M,M,J+1]); % HDGT gen

    %%% Gain Saturation %%%
%         HDTcol = HDPsi2Col*Perturb ; HDT = reshape(HDTcol,[M,M,J+1]); % HDGT gen
        HDT = reshape(HDPsi2Col*Perturb,[M,M,J+1]); % HDGT gen
         Gsat = G0./(1+Ipre/Is) ; % HG for Gain Saturation
         GT = Gsat(:,1:(J/(Jred+1)):end) ; % HG for Gain Saturation, gain for heat source
         Gnew = exp(Gsat*dz) ; % HG for Gain Saturation
%          HGnew = PsimPsinCol*Gnew ; HG = reshape(HGnew,[M,M,J+1]); % HG for Gain Saturation
         HG = reshape(PsimPsinCol_Diag*Gnew,[M,M,J+1]); % HG for Gain Saturation
    %%%%%
      
%%%%%% for idxz = 2:J+1 % z loop old
for idxz = 0:Jred 

        H = H0;
    for nstep = 2:Nplusone
%         H = HDGT(:,:,idxz*N+nstep)*H;
        H = HG(:,:,idxz*N+nstep)*HDT(:,:,idxz*N+nstep)*H; % Gain Saturation
    end
       b(:,idxz+2) = H*b(:,idxz+1);

end
    %%% Update intensity distribution %%%
    ax = (abs(psi*b).^2)*nco ; % intensity [W/um^2] Corrected % I/n = E^2 index effect corrected

    %%% Update temperature and index %%% 
%    Qcore = 2*G*(lslp - 1)*ax ; % Corrected % Power gain corrected
    Qcore = 2*(lslp - 1)*GT.*ax ; % Corrected % Power gain corrected % HG for gain saturation
    Q(coreUp:coreLow,:) = Qcore; % coarsening by sampling 
    dVector = cminus*T(2:nx-1,:)+c*(T(3:nx,:)+T(1:nx-2,:))+Qcoef*Q(2:nx-1,:);
    T(2:cladLow,:) = TriD\dVector ;

%%%%%%%
    F = griddedInterpolant(GridXT,GridZT,T(Nclad+1:Nclad+Rcore,:)); % 2D interpolation of T distribution by griddedInterpolant with defalt setting
    T_interp = F(GridXL,GridZL);
    Dn(:,:) = eta*(T_interp-Tb); % refractive index distribution
%%%%%%%

%%%%%%% 
    F = griddedInterpolant(GridXT,GridZT,ax); % HG Gain Saturation
    Ipre = F(GridXL,GridZL); % HG Gain Saturation
%%%%%%% 

    %%% Monitoring the progress %%%
    if ( mod(idxt,100) == 0)
%        disp([num2str(100*idxt/P) ' % at t = ' num2str(toc) ' s']);
        disp([num2str(toc/(3600*idxt/P) - toc/3600),' hours away.']);
        %%% Recording temperature and fields %%%
        T_rec(:,:,idxt/100) = T(Nclad+floor(Rcore/4),:); % temperature at the end of fiber  
        b_rec(:,:,idxt/100) = b(:,(end-Znumber*Zspacing):Zspacing:end); % MultiZ Multiple Z
    end

end
toc  
% save('DnT.mat','Dn','T','Ipre') %%% ctrl+T to continue next sim from this point. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Save Simulation Results %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gnum = round(exp(G0*1e6));
fnum = round(factor) ;
pertMode = find(PpertRatio) ; % Pert Correction 2021-11-03
if size(pertMode,1) > 1       % Pert Correction 2021-11-03
  pertMode = 0 ;              % Pert Correction 2021-11-03
end                           % Pert Correction 2021-11-03
FrontName = ['PertGSMZ_FocusD5_ref_',num2str(round(Pseed)),'_PM',num2str(M),'g',num2str(gnum),'_',num2str(fnum),'E-', num2str(pert),'_m',num2str(pertMode),'_f',num2str(fpert)]; % Pert Correction 2021-11-03
SearchName = [FrontName,'*.mat'];
FilesHere = dir(SearchName); 
CurrStep = size(FilesHere,1); 
FileName = [FrontName,'_', num2str(CurrStep+1),'.mat'];
%%% Perturbed input
save(FileName,'Ipre','factor','pert','fpert','Pphase','b0','ax','b','b_rec','bamp','bph','T','Dn','dt','dx','dz','dzT','t_tot','tplot','w','W','x','X','z','z_red','-v7.3')
% save(FileName,'T_rec','xT_rec','Ipre','factor','pert','fpert','Pphase','b0','ax','b','b_rec','bamp','bph','T','Dn','dt','dx','dz','dzT','t_tot','tplot','w','W','x','X','z','z_red','-v7.3')
