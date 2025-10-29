clear
clc
format long

%% visualization setting for boundary evolution
is_viz = true;

% Number of sampling points
N = [400];

%% initial condition
for rr = 1:length(N)
%first generate 25 points
pts=[...
2.0, 7.8
3.0, 7.5
3.5, 6.8
3.3, 5.8
5.0, 5.3
5.2, 4.3
4.5, 4.0
3.5, 4.0
3.5, 3.0
5.0, 1.0
3.5, 0.5
3.0, 1.0
2.0, 1.5
1.0, 1.0
0.5, 0.5
-1.0, 1.0
0.5, 3.0
0.5, 4.0
-0.5, 4.0
-1.2, 4.3
-1.0, 5.3
0.7, 5.8
0.4, 6.8
1.0, 7.5
2.0, 7.8];
%make the curve more smooth
x_data = (pts(:,1)-2.0);
y_data = 1.5*(pts(:,2));
x_data(x_data<0) = x_data(x_data<0) - 0.2;
x_data(x_data>0) = x_data(x_data>0) + 0.2;
% Use cubic spline interpolation to obtain the initial curve which consists
% of N points and make all of them equally spaced
t_data = 1:length(x_data);
t_query = linspace(t_data(1),t_data(end),N+1);
xx = spline(t_data,x_data,t_query); 
yy = spline(t_data,y_data,t_query); 
scurve = [xx',yy'];
new_curve = zeros(size(scurve));
dx = diff(scurve(:,1));
dy = diff(scurve(:,2));
ds = sqrt(dx.^2 + dy.^2);
deltas = [0; cumsum(ds)];
target_arclength = linspace(0, deltas(end), N+1);
new_curve(:,1) = interp1(deltas, scurve(:,1),target_arclength,'spline');
new_curve(:,2) = interp1(deltas, scurve(:,2), target_arclength, 'spline');
xxx = new_curve(1:N,1);
yyy = new_curve(1:N,2);
%shrink the curve
xxx = xxx./4;
yyy = yyy./4;

curve0 = [xxx,yyy];
curve = curve0;
theta0 = [0:2*pi/N:2*pi-2*pi/N]';
% Here tvec0 is used to align the vectors obtained from local
% svd so that the tangent vectors are counterclockwise
tvec0 = [-sin(theta0) cos(theta0)];
nvec0 = [-cos(theta0) -sin(theta0)];
[tvec,~,~] = localsvd(curve,21,1,tvec0,nvec0);
for pp = 2:N
     if (tvec(:,pp-1)'*tvec(:,pp)<0)
        tvec(:,pp) = -tvec(:,pp);
     end
end
tvec = tvec';

nvec(:,1) = tvec(:,2);
nvec(:,2) = -tvec(:,1);

told = tvec;
nold = nvec;

%time step size
dt = 1.0e-5;
%start time
time = 0.;
%count the number of resampling procedure
resample = 0;
%end time
time_end = 2.1;
%total time steps
T = time_end/dt;
%count the time step
tdt = 0;

while time <= time_end+dt
    [A,b,t_gmls,n_gmls,rightpt,alpha,center_x,center_y] = velocity(curve,N,told,nold);
    told = t_gmls';
    nold = -n_gmls';

    % velocity in normal direction
    d_curve = -A\b;
    vn = n_gmls'.*d_curve;
    %max distance and min distance from boundary points to centroid
    max_radius(tdt+1,1) = max(sqrt((curve(:,1)-center_x).^2 + (curve(:,2)-center_y).^2));
    min_radius(tdt+1,1) = min(sqrt((curve(:,1)-center_x).^2 + (curve(:,2)-center_y).^2));
    curve_old = curve;
    %keep record of boundary at sepcific time to plot them later
    if tdt == 1000
    curve_t1 = curve;
    elseif tdt == 6000
        curve_t2 = curve;
    elseif tdt == 2.1e05
        curve_t3 = curve;
    end

    %% update boundary after dt
    %forward euler
    curve = curve + dt.*(n_gmls'.*d_curve);
   
    %% Resampling by cubic spline interpolation if necessary to prevent points from collapsing
    scurve = curve;
    new_curve = zeros(size(scurve));
    dx = diff(scurve(:,1));
    dy = diff(scurve(:,2));
    ds = sqrt(dx.^2 + dy.^2);
    deltas = [0; cumsum(ds)];
    target_arclength = linspace(0, deltas(end), N);
    new_curve(:,1) = interp1(deltas, scurve(:,1),target_arclength,'spline');
    new_curve(:,2) = interp1(deltas, scurve(:,2), target_arclength, 'spline');
    curve = new_curve;
    resample = resample + 1;
  
   if is_viz
        figure (1)%plot boundary evolution and its initial boundary
        plot([curve0(:,1);curve0(1,1)], [curve0(:,2);curve0(1,2)],'c','LineWidth',2,'MarkerSize',5)
        hold on
        plot(curve_old(:,1), curve_old(:,2),'k--','LineWidth',2,'MarkerSize',5);
        hold off
        xlim([-1 1])
        ylim([0 3.0])
        legend('Initial Boundary','Current Boundary','Location','northeast','interpreter','latex','FontSize',10)
        title(['Boundary Moton at t=',num2str(tdt)])
        pause(0.1)
   end
   time = time + dt;
    tdt = tdt + 1;
    end
end
figure (2)%plot the boundary at specific times
plot([curve0(:,1);curve0(1,1)], [curve0(:,2);curve0(1,2)],'k-','LineWidth',2,'MarkerSize',5);
hold on
plot(curve_t1(:,1), curve_t1(:,2),'b--','LineWidth',2,'MarkerSize',5);
hold on
plot(curve_t2(:,1), curve_t2(:,2),'m--','LineWidth',2,'MarkerSize',5);
hold on
plot(curve_t3(:,1), curve_t3(:,2),'r-.','LineWidth',2,'MarkerSize',5);
hold off
legend('$t=0$','$t=0.01$','$t=0.06$','$t=2.1$','Location','northeast','interpreter','latex','FontSize',14)
figure (3)%plot the max distance and min distance every 100 time steps
plot([0:100:2.1e05]'*1.0e-05,max_radius(1:100:length(max_radius)),'b','LineWidth',2,'MarkerSize',5);
hold on
plot([0:100:2.1e05]'*1.0e-05,min_radius(1:100:length(min_radius)),'k','LineWidth',2,'MarkerSize',5);
hold off
xlim([0 2.1])
xlabel('Time')
ylabel('Distance to the center','interpreter','latex')
legend('maximum distance','minimum distance','FontSize',14)

function [A,b,t_gmls,n_gmls,rightpt,alpha,center_x,center_y] = velocity(curve,N,told,nold)
rightpt = zeros(N,1);
leftpt = zeros(N,1);
A = zeros(N,N);
b_int = zeros(N,1);
b_sint = zeros(N,1);
x = curve(:,1);
y = curve(:,2);
% number of k-nearest points
if mod(ceil(sqrt(N)),2)==0
    k = ceil(sqrt(N))+1;
else
    k = ceil(sqrt(N));
end
epsilon = 1.0e-12;
trial=0;
%approximate tangent vector and normal vector by local svd
[t_rough,n_rough,inds] = localsvd(curve,k,1,told,nold);
n_rough(1,:) = -t_rough(2,:);
n_rough(2,:) = t_rough(1,:);
%approximate tangent vector and normal vector by gmls
[t_gmls6, alpha6] = localgmls(curve,inds,1,t_rough,n_rough,6);
t_gmls = t_gmls6;
n_gmls = zeros(size(t_gmls));
n_gmls(1,:) = -t_gmls(2,:);
n_gmls(2,:) = t_gmls(1,:);
for pp = 1:N
     if (n_rough(:,pp)'*n_gmls(:,pp)<0)
        n_gmls(:,pp) = -n_gmls(:,pp);
    end
end
% iterate gmls process to make \alpha_1 close to 0
while max(abs(alpha6(1,:)'))>epsilon
    [t_gmls6, alpha6] = localgmls(curve,inds,1,t_gmls,n_gmls,6);
    t_gmls = t_gmls6;
n_gmls = zeros(size(t_gmls));
n_gmls(1,:) = -t_gmls(2,:);
n_gmls(2,:) = t_gmls(1,:);
for pp = 1:N
     if (n_rough(:,pp)'*n_gmls(:,pp)<0)
        n_gmls(:,pp) = -n_gmls(:,pp);
    end
end
trial = trial+1;
end
[~, alpha6] = localgmls(curve,inds,1,t_gmls,n_gmls,6);
alpha = alpha6;
alpha(1,:) = 0; 

%approximate curvature at each point
cur6_gmls =-2*alpha6(2,:)'./((1+alpha6(1,:)'.^2).^(3/2));
cur_gmls = cur6_gmls;

%integration bound \delta s_i and \delta s_{-i}
for qq = 1:N-1
  rightpt(qq) = t_gmls(:,qq)'*[x(qq+1)-x(qq);y(qq+1)-y(qq)];
end
rightpt(N) = t_gmls(:,N)'*[x(1)-x(N);y(1)-y(N)];

for qq = 2:N
    leftpt(qq) = t_gmls(:,qq)'*[x(qq-1)-x(qq);y(qq-1)-y(qq)];
end
leftpt(1) = t_gmls(:,1)'*[x(N)-x(1);y(N)-y(1)];

%% polynomial
poly_matrix = [rightpt rightpt.^2 rightpt.^3 rightpt.^4 rightpt.^5 rightpt.^6];
poly = diag(poly_matrix*alpha,0);
%p^2(s)/s^2
fovers_lmatrix = [ones(N,1) leftpt leftpt.^2 leftpt.^3 leftpt.^4 leftpt.^5];
fovers_leftpt = diag(fovers_lmatrix*alpha,0);
fovers_rmatrix = [ones(N,1) rightpt rightpt.^2 rightpt.^3 rightpt.^4 rightpt.^5];
fovers_rightpt = diag(fovers_rmatrix*alpha,0);

%% polynomial_prime
%sqrt(1+(p'(0))^2)=1 as alpha1=0
fprime0 = sqrt(1 + (alpha(1,:))'.^2);
%sqrt(1+(p'(delta s))^2)
poly_rprime_matrix = [ones(N,1) 2*rightpt 3*rightpt.^2 4*rightpt.^3 5*rightpt.^4 6*rightpt.^5];
poly_rprime = diag(poly_rprime_matrix*alpha,0);
fprime_rightpt = sqrt(1 + poly_rprime.^2);
%sqrt(1+(p'(-delta s))^2)
poly_lprime_matrix = [ones(N,1) 2*leftpt 3*leftpt.^2 4*leftpt.^3 5*leftpt.^4 6*leftpt.^5];
poly_lprime = diag(poly_lprime_matrix*alpha,0);
fprime_leftpt = sqrt(1 + poly_lprime.^2);

%centroid
int_x = (curve(:,1) + [curve(2:end,1);curve(1,1)].*fprime_rightpt).*rightpt./2;
int_s = (ones(N,1) + fprime_rightpt).*rightpt./2;
int_y = (curve(:,2) + [curve(2:end,2);curve(1,2)].*fprime_rightpt).*rightpt./2;
center_x = sum(int_x)/sum(int_s);
center_y = sum(int_y)/sum(int_s);

%% matrix A
for j=1:N
    %% Distance matrix and its norm
    D_matrix = [x(j),y(j)] - curve;
    norm_matrix = vecnorm(D_matrix,2,2);
   
%% Matrix A
norm_left = -1/(2*pi).*log(norm_matrix);
norm_right = [norm_left(2:end);norm_left(1)];
G_left = norm_left.*rightpt;
G_mid_right = norm_right.*fprime_rightpt.*rightpt;
G_right = [G_mid_right(end);G_mid_right(1:end-1)];
A(j,:) = (G_left + G_right)'./2;

if j==1
    A(j,N) = G_right(N)/2;
else
    A(j,j-1) = G_right(j-1)/2;
end
if j==N
    A(j,1) = G_left(1)/2;
else
    A(j,j+1) = G_left(j+1)/2;
end
fovers_lterm = - 1/(4*pi)*log(1+(fovers_leftpt(j))^2)*fprime_leftpt(j)*abs(leftpt(j))/2;
fovers_rterm = - 1/(4*pi)*log(1+(fovers_rightpt(j))^2)*fprime_rightpt(j)*rightpt(j)/2;
%singular terms
int_beta0 = -1/(2*pi)*(-leftpt(j)*log(-leftpt(j))+leftpt(j)+rightpt(j)*log(rightpt(j))-rightpt(j));
int_beta1 = -1/(2*pi)*(-1/2*leftpt(j)^2*log(-leftpt(j))+1/4*(leftpt(j))^2+1/2*rightpt(j)^2*log(rightpt(j))-1/4*rightpt(j)^2);
int_beta2 = -1/(2*pi)*(1/3*(-leftpt(j))^3*log(-leftpt(j))-1/9*(-leftpt(j))^3+1/3*rightpt(j)^3*log(rightpt(j))-1/9*rightpt(j)^3);

%psi_i(s) = beta0 + beta1*s + beta2*s^2
A(j,j) = int_beta0 - int_beta1*(rightpt(j)+leftpt(j))/(rightpt(j)*leftpt(j)) + int_beta2/(rightpt(j)*leftpt(j));
if j==1
    A(j,N) = A(j,N) + fovers_lterm - int_beta1*rightpt(j)^2*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    - int_beta2*rightpt(j)*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j)));
else
    A(j,j-1) = A(j,j-1) + fovers_lterm - int_beta1*rightpt(j)^2*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    - int_beta2*rightpt(j)*fprime_leftpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j)));
end
if j==N
    A(j,1) = A(j,1) + fovers_rterm + int_beta1*leftpt(j)^2*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    + int_beta2*leftpt(j)*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j))) ;
else
    A(j,j+1) = A(j,j+1) + fovers_rterm + int_beta1*leftpt(j)^2*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(leftpt(j)-rightpt(j))) ...
    + int_beta2*leftpt(j)*fprime_rightpt(j)/(rightpt(j)*leftpt(j)*(rightpt(j)-leftpt(j))) ;
end

%% b
G_prime = 1/(2*pi)*D_matrix./(norm_matrix.^2);
G_prime(j,:) = 0;

%trapezoidal rule
norm_left_prime = G_prime;
norm_right_prime = [norm_left_prime(2:end,:);norm_left_prime(1,:)];
G_left_prime = norm_left_prime.*rightpt.*cur_gmls;
if j==1
    G_left_prime(N,:) = norm_left_prime(N,:)*fprime_leftpt(1)*abs(leftpt(1))*cur_gmls(N);
else
    G_left_prime(j-1,:) = norm_left_prime(j-1,:)*fprime_leftpt(j)*abs(leftpt(j))*cur_gmls(j-1);
end
G_right_prime = norm_right_prime.*fprime_rightpt.*rightpt.*[cur_gmls(2:end);cur_gmls(1)];
b_int(j) = 0.5*sum(dot(n_gmls,(G_left_prime'+[G_right_prime(end,:);G_right_prime(1:end-1,:)]'),1))+1/(2*pi)*cur_gmls(j)*alpha(2,j)*rightpt(j)/2 + 1/(2*pi)*cur_gmls(j)*alpha(2,j)*abs(leftpt(j))/2;

%simpson rule
norm_lprime = [norm_left_prime(end,:);norm_left_prime(1:end-1,:)];
norm_mprime = norm_left_prime;
norm_rprime = norm_right_prime;
G_lprime = norm_lprime.*abs(rightpt-leftpt).*fprime_leftpt.*[cur_gmls(end);cur_gmls(1:end-1)];
G_mprime = norm_mprime.*abs(rightpt-leftpt).*cur_gmls;
G_rprime = norm_rprime.*abs(rightpt-leftpt).*fprime_rightpt.*[cur_gmls(2:end);cur_gmls(1)];
if j==1
    b_sint(j) = 1/6*sum(dot(n_gmls,([G_lprime(2:end,:);G_lprime(1,:)]' + 4.*G_mprime' + [G_rprime(end,:);G_rprime(1:end-1,:)]'),1))+1/(2*pi)*cur_gmls(j)*alpha(2,j)*(4/6*abs(rightpt(j)-leftpt(j))+1/6*fprime_rightpt(N)*abs(rightpt(N)-leftpt(N))+1/6*fprime_leftpt(j+1)*abs(rightpt(j+1)-leftpt(j+1)));
elseif j==N
    b_sint(j) = 1/6*sum(dot(n_gmls,([G_lprime(2:end,:);G_lprime(1,:)]' + 4.*G_mprime' + [G_rprime(end,:);G_rprime(1:end-1,:)]'),1))+1/(2*pi)*cur_gmls(j)*alpha(2,j)*(4/6*abs(rightpt(j)-leftpt(j))+1/6*fprime_rightpt(j-1)*abs(rightpt(j-1)-leftpt(j-1))+1/6*fprime_leftpt(1)*abs(rightpt(1)-leftpt(1)));
else
    b_sint(j) = 1/6*sum(dot(n_gmls,([G_lprime(2:end,:);G_lprime(1,:)]' + 4.*G_mprime' + [G_rprime(end,:);G_rprime(1:end-1,:)]'),1))+1/(2*pi)*cur_gmls(j)*alpha(2,j)*(4/6*abs(rightpt(j)-leftpt(j))+1/6*fprime_rightpt(j-1)*abs(rightpt(j-1)-leftpt(j-1))+1/6*fprime_leftpt(j+1)*abs(rightpt(j+1)-leftpt(j+1)));
end
end
%trapezoidal rule
% b = b_int + cur_gmls./2;
%simpson rule
b = b_sint./2 + cur_gmls./2;
end

% build multi-indices for multivariate polynomials if necessary
function index =  generatemultiindex(N,dim)
% input
% N : max degree of polynomials
% dim : dimension

P = nchoosek(N+dim,N);
index = zeros(dim,P);

Ntotal = (N+1)^dim;
allindex = zeros(dim,Ntotal);

for i=1:dim
    nskip = (N+1)^(dim-i);
    for k=1:Ntotal/nskip/(N+1)
        for j=1:N+1
            allindex(i,(k-1)*nskip*(N+1)+(j-1)*nskip+1:(k-1)*nskip*(N+1)+j*nskip) = j-1;
        end
    end  
end

index1 = find(sum(allindex)<=N);
index = allindex(:,index1);
end

% find k-nearest points (can use the built-in function (knn) directly)
function [ds,inds] = knnCPU(R,Q,k,maxMem)
% Find the k nearest neighbors for each element of the query set Q among
% the points in the reference set R
% Q is Nxd where N is the number of query points and d is the dimension
% R is Mxd where M is the number of reference points and d is the dimension
% k is the number of nearest neighbors desired
% maxMem is the number of GB of ram knnCPU can use
% ds is Nxk and contains the distances to the k nearest neighbors.
% inds is Nxk and contains the indices of the k nearest neighbors.
    if (nargin<4)
        maxMem = 2;
    end

    M = size(R,1);
    N = size(Q,1);    
    
    maxArray = (maxMem*2500)^2;
    blockSize = floor(maxArray/M);
    blocks = floor(N/blockSize);

    ds = zeros(N,k);
    inds = zeros(N,k);

    Nr = sum(R.^2,2);
    Nq = sum(Q.^2,2);

    for b = 1:blocks
        dtemp = -2*R*Q((b-1)*blockSize + 1:b*blockSize,:)';
        dtemp = bsxfun(@plus,dtemp,Nr);
        dtemp = bsxfun(@plus,dtemp,Nq((b-1)*blockSize + 1:b*blockSize)');
        [dst,indst] = sort(dtemp,1);
        ds((b-1)*blockSize + 1:b*blockSize,:) = dst(1:k,:)';
        inds((b-1)*blockSize + 1:b*blockSize,:) = indst(1:k,:)';
    end

    if (blocks*blockSize < N)
        dtemp = -2*R*Q(blocks*blockSize + 1:N,:)';
        dtemp = bsxfun(@plus,dtemp,Nr);
        dtemp = bsxfun(@plus,dtemp,Nq(blocks*blockSize + 1:N)');
        [dst,indst] = sort(dtemp,1);
        ds(blocks*blockSize + 1:N,:) = dst(1:k,:)';
        inds(blocks*blockSize + 1:N,:) = indst(1:k,:)';
    end

    ds = real(sqrt(ds));
end

% GMLS approximation of tangent vectors
function [t_tilde,alpha] = localgmls(x,inds,d,t_rough,n_rough,ell)
N = size(x,1);
n = size(x,2);
k0 = size(inds,2);

if d > 1
    index =  generatemultiindex(ell,d); % index is d*term    
    index = index(:,2:end);
else
    index = 1:ell;
end
term = size(index,2);

ind_alin = zeros(d,1);
for jj = 1:d
    ind_alin(jj) = find((sum(index,1)==1)&(index(jj,:)==1));
end

for pp=1:N
    dx = x(inds(pp,:),:)';
    dx = dx - repmat(x(pp,:)',1,k0); 

     lcx = zeros(k0,n);
     lcx = dx'*[t_rough(:,pp) n_rough(:,pp)];
                
     %%% Phi = [1 yy yy.^2 yy.^3 ...]
     Phi = ones(k0,term); % k*term
        
     %%% constant has 1, deg 1 has d, deg 2 has (d+1)*d/2
     for ss = 1:term
        for rr = 1:d
            Phi(:,ss) = Phi(:,ss).*lcx(:,rr).^index(rr,ss);
        end
     end

     % coefficients in polynomial by least squares
     a1 =(Phi'*Phi)\(Phi'*lcx(:,n));
     alpha(:,pp) = a1;
        
     t_tilde(:,pp) = t_rough(:,pp) + a1(ind_alin)*n_rough(:,pp);
     t_tilde(:,pp) = t_tilde(:,pp)/vecnorm(t_tilde(:,pp));
    end
end

% Get approximate tangent vectors and normal vectors by local SVD method
function [t_rough,n_rough,inds] = localsvd(x,k0,d,tvec,nvec)
[~,inds] = knnCPU(x,x,k0);
N = size(x,1);
n = size(x,2);
   
for pp = 1:N
    dx = x(inds(pp,:),:)';
    dx = dx - repmat(x(pp,:)',1,k0); 
    [U,~,~] = svd(dx); 
    
    t_rough(:,pp) = U(:,d);        
    n_rough(:,pp) = U(:,d+1:n);

    if (tvec(pp,:)*t_rough(:,pp)<0)
        t_rough(:,pp) = -t_rough(:,pp);
    end

    if (nvec(pp,:)*n_rough(:,pp)<0)
        n_rough(:,pp) = -n_rough(:,pp);
    end
end
end
