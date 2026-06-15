%This code solves the nonlinear Poisson Boltzmann equation for a Kirkwood
%sphere with a charge at its center, and calculates the electrostatic free
%energy. 

%Version of May 2022 by  
%Aidan Winiewicz, Sylvia Amihere, Shan Zhao 

%Contact us: 
%Dr. Shan Zhao, 
%Department of Mathematics
%The University of Alabama, USA
%Email: szhao@ua.edu

%Please cite the following when using this code in a publication:
%S. Amihere, W. Geng, and S. Zhao, Benchmarking electrostatic free energy 
%of the nonlinear Poisson-Boltzmann model for the Kirkwood sphere, 
%Communications in Information & Systems, in press, (2022).


function  nonlinear_npb_energy %solves the nonlinear Poisson Boltzmann equation
prompt1 = 'What is the radius of the molecule you wish to consider? \n';
prompt2 = 'What is the dielectric constant of the molecule? \n';
prompt3 = 'What is the dielectric constant of the solvent? \n';
prompt4 = 'What is the ionic strengh in your solvent? \n';
prompt5 = 'What is the total charge of your molecule \n';

%user may select xstart (which is the same as the molecule radius), epsa, epsw, bulk_strength (which will result in a change
%in xinf), and charge
while 1
    xstart = input(prompt1);
    if xstart > 0
        break;
    else
        disp('Pick a positive molecule radius');
    end
end
while 1
    epsa = input(prompt2);
    if epsa > 0
        break;
    else
        disp('Pick a positive dielectric constant for the molecule region');
    end
end
while 1
    epsw = input(prompt3);
    if epsw > 0
        break;
    else
        disp('Pick a positive dielectric constant for the solvent');
    end
end
while 1
    bulk_strength = input(prompt4);
    if (.1<= bulk_strength) && (bulk_strength <= 10)
        break;
    else
        disp('Pick a bulk strength between [.1,10]');
    end
end

%epsa = 1;
%epsw = 80;
%xstart = 2;
%bulk_strength = 1;
charge = input(prompt5);

n = 12000;      %Domain and mesh
if bulk_strength <= 7
    xend = 38.58410733 + 7.789115646*xstart - 5.420714356*bulk_strength;
else
    xend = 10;
end
dx = (xend-xstart)/(n-1);


kappa = 8.430325455*bulk_strength/epsw; %PDE coefficients
lambda = sqrt(kappa);
ec2_kbt = (332.0716/0.5961574);

nxk = 4; %number of fictitious points in the MIB scheme
L = 8; %order of method
norder = 2; %2nd derivative
norder1 = 1; %1st derivative

phi = -(ec2_kbt*charge)/(epsw*(xstart^2)); %boundary condition on the left
psi = 0; %boundary condition on the right

xi = zeros(n,1);
for i = 1:n
    xi(i) = xstart + (i-1)*dx;
end

%FD coefficients for 2nd and 1st derivatives
vde = lag4(xstart,xend,n,nxk,norder);
vde2 = lag4(xstart,xend,n,nxk,norder1);

%weights in leftidm and rightidm
m0 = 1;
emat = leftidm(dx,m0,L);
remat = rightidm(dx,m0,L);
if nxk<=L
    m=nxk;
    emat_hdm = lefthdm(dx,m0,L,m,emat);
    remat_hdm = righthdm(dx,m0,L,m,remat);
end

%righthand side vector
rhs_vector = RHS_Vector(emat_hdm,remat_hdm,vde,vde2,xi,nxk,phi,psi,n,L);

guessU = zeros(n,1);
finalU = nonlinear_solver(guessU,rhs_vector,vde,vde2,emat_hdm,remat_hdm,n,L,nxk,xi,lambda);

% compute free energy
%potential_a = finalU(1);
%fprintf('u(a) is equal to %17.15f \n', potential_a);
UNC1 = finalU(1) - ((ec2_kbt*charge)/(epsa*xstart));%C1 = U^{+}(a) - (Aq/(e^{-}a))
numerical_energy = (0.5)*(charge)*(UNC1)*(0.5961574);%!(1/2)*q*C1*KBT
%fprintf('Numerical energy is %18.15f \n', numerical_energy);
next_step = next_step_generator(n,xstart,epsa,epsw,bulk_strength,charge);
main_error = abs(next_step - numerical_energy);
avg_error = (next_step + numerical_energy)*1/2;
fprintf('The NPB energy of the Kirkwook sphere is %18.15f Kcal/mol \n', avg_error);
fprintf('The numerical error is estimated to be %d \n',main_error);
end

%LaGrange kernel up to 4th order

function vde = lag4(xstart,xend,n,nxk,nod)
%inputs:xend is endpoint of computational domain
% xstart is beginning of domain
% n is number of gridpoints (includes endpts)
% nxk is points in each direction from center
% nod is derivative order
%vde is output matrix

vde = zeros(nod,2*nxk+1); %initialize a matrix to hold the approximation
%weights of each derivate approximation
NS = 2*nxk+1;   %length of stencil/number of gridpoints
x = zeros(NS,1); %initialize vector which will assume specific gridpoints
DX = (xend - xstart)/(n-1);
ZFD = 0; %center of the approximation
for i = 1:NS
    x(i) = -nxk*DX+(i-1)*DX;
end
c = weights(x,ZFD,nod);

for i = 2:(nod+1)
    for j = 1:(2*nxk+1)
        vde(i-1,j) = c(j,i); %vde is the transpose of c, and vde does not
        %include the first row of c.
    end
end
end
function emat = leftidm(DX,m,L)

emat = zeros(m,L+2);
NS = L+1+m; %length of each FD stencil
ND = 1; %Neumann BC is first derivative
w = zeros(NS,ND+1);
ME = L+2; %unknown weights for single FP
BE = zeros(ME,1);
x = zeros(NS,1);

ZFD = 0; %finite difference weights
for i = 1:NS
    x(i) = -m*DX+(i-1)*DX;
end
c = weights(x,ZFD,ND);

for i = 1:(ND+1)
    for j = 1:(NS)
        
        w(j,i) = c(j,i);
    end
end

%bdy condition
% u^{(nd)} = phi

BE(ME) = 1; %bdy weight
IW = ND+1;

for i = 1:(L+1)
    BE(i) = BE(i)-w(m+i,IW);
end

for i = 1:ME %solve
    BE(i) = BE(i)/w(1,IW);
end

for i = 1:m
    for j = 1:ME
        
        emat(i,j) = BE((i-1)*ME + j);
        
    end
end
end
%implicit derivative matching v(b)=0
function remat = rightidm(DX, m, L)

remat = zeros(m,L+1);
ND = 0; %Dirichlet BC
ME = L+1; %unknown weights for single FP
NS = L+1; %length of each FD stencil
XEE = zeros(NS,1);
w = zeros(NS,ND+1);

%initialization
ZFD = 0;
for i=1:(NS+1)
    XEE(i) = -L*DX+(i-1)*DX;
end

XE = [XEE(1:NS-m);XEE(NS-m+2:end)];

c = weights(XE,ZFD,ND);

for i=1:(ND+1)
    for j=1:NS
        w(j,i) = c(j,i);
    end
end
% Boundary Condition
%Dirichlet BC ==> v(b) = 0
BE = zeros(ME,1);
BE(ME) = 1; %bdy weight
IW = ND+1;
for i =1:L
    BE(i)=BE(i)-w(L+1-i,IW);
end
for i=1:ME %solve
    BE(i) = BE(i)/w(NS,IW);
end

for i = 1:m
    for j = 1:ME
        remat(i,j) = BE((i-1)*ME+j);
    end
end
end

%hierachical derivative matching for more than 1 fictitious pt
function ematL = lefthdm(DX,m0,L,m,emat)
ME=L+2;		%unknown weights for single FP
ND=1;		%consider maximum third order BC
ematold = zeros(m0,ME);

for i=1:m0
    for j=1:ME
        ematold(i,j)= emat(i,j);
    end
end

for mnew = (m0+1):m	% From master FP to all FP
    %increase one more FP at every step
    
    NS=L+1+mnew;	% length of each FD stencil
    XE = zeros(NS,1);
    ematnew = zeros(mnew, ME);
    ematL = zeros(m,ME);
    BE = zeros(ME);
    w = zeros(NS,ND+1);
    
    % inherit the center part.
    for i=1:mnew-1
        for j=1:ME
            ematnew(i,j)= ematold(i,j);
        end
    end
    
    ZFD=0;		%finite difference weights
    for i=1:NS
        XE(i)= -mnew*DX+(i-1)*DX;
    end
    c = weights(XE,ZFD,ND);
    for i=1:ND+1
        for j=1:NS
            w(j,i)= c(j,i);
        end
    end
    
    
    %###################### Boundary condition ###############
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %       v'(a) = -A*q/(e^{+}*a^2)
    %C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BE(ME)=1; %boundary weight
    
    IW=ND+1;		%first order derivative (for Neumann BCs)
    for i=1:L+1
        BE(i)=BE(i)-w(mnew+i,IW);
    end
    
    for i=1:mnew-1
        for j=1:ME
            BE(j)= BE(j)-w(mnew+1-i,IW)*ematnew(i,j);
        end
    end
    
    for i=1:ME		%solve
        BE(i)=BE(i)/w(1,IW);
    end
    
    for i=1:ME
        ematnew(mnew,i)= BE(i);
    end
    
    %	Be ready for next enlargement or computation.
    for i=1:mnew
        for j=1:ME
            ematold(i,j)=ematnew(i,j);
        end
    end
    
end		% MNEW=M0+1:M

for i=1:m
    for j=1:ME
        ematL(i,j)= ematold(i,j);
    end
end
end

%hierarchical derivative matching
function rematL = righthdm(DX, m0,L,m,remat)
rematold = zeros(m0,L+1);

ME=L+1;		%unknown weights for single FP
ND=0;	%consider maximum third order BC
for i=1:m0
    for j=1:ME
        rematold(i,j)= remat(i,j);
    end
end


for mnew=m0+1:m	%From master FP to all FP
    %increase one more FP at every step
    
    NS=L+mnew;	%length of each FD stencil
    w = zeros(NS,ND+1);
    XEE = zeros(NS+1,1);
    rematnew = zeros(mnew,ME);
    BE= zeros(ME,1);
    rematL = zeros(m,ME);
    
    %inherit the center part.
    for i=1:mnew-1
        for j=1:ME
            rematnew(i,j)= rematold(i,j);
        end
    end
    
    ZFD=0;	%finite difference weights
    for i=1:NS+1
        XEE(i)=-L*DX+(i-1)*DX;
    end
    XE = [XEE(1:(NS-mnew));XEE((NS-mnew)+2:end)];
    c = weights(XE,ZFD,ND);
    
    for i=1:ND+1
        for j=1:NS
            w(j,i)=c(j,i);
        end
    end
    
    
    %###################### Boundary condition ###############
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 	    Dirichlet BC
    %        v(b) = 0
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    BE(ME)=1;
    
    IW=ND+1;		%zero order derivative
    for i=1:L
        BE(i)=BE(i)-w(L+1-i,IW);
    end
    
    for i=1:mnew-1
        for j=1:ME
            BE(j)= BE(j)-w(L+i,IW)*rematnew(i,j);
        end
    end
    
    for i=1:ME	%solve
        BE(i)=BE(i)/w(NS,IW);
    end
    
    for i=1:ME
        rematnew(mnew,i)= BE(i);
    end
    
    %Be ready for next enlargement or computation.
    for i=1:mnew
        for j=1:ME
            rematold(i,j)=rematnew(i,j);
        end
    end
end	%MNEW=M0+1,M

for  i=1:m
    for j=1:ME
        rematL(i,j)= rematold(i,j);
    end
end
end

function rhs_vector = RHS_Vector(emat,remat,vde,vde2,xi,nxk,phi,psi,n,L)

%right hand side vector
rhs_vector = zeros(n,1);

%left boundary
for i = 1:(n-1)
    for j = 1:(nxk-(i-1))
        IP = nxk-(i+j-2);
        rhs_vector(i) = rhs_vector(i) - (vde(2,j) + (2/xi(i))*vde2(1,j))*emat(IP,L+2)*phi;
    end
end

%right boundary
for i = 1:n-1
    for j = n-i+1:nxk
        IP = i+j-n;
        rhs_vector(i) = rhs_vector(i) - (vde(2,j+nxk+1) + (2/xi(i))*vde2(1,j+nxk+1))*remat(IP,L+1)*psi;
    end
end

%right boundary explicitly enforced
for i= n:n
    rhs_vector(i) = psi;
end
end

%this function uses the Inexact Newton Method to obtain the solution vector
function finalU = nonlinear_solver(guessU,rhs_vector,vde,vde2,emat,remat,n,L,nxk,xi,lambda)
TOL = 1*(10^-10); %tolerance for inexact newton. Tolerance smaller than
%10^-12 may not converge !!
errINNER = 1; %initial error
%iDisplay = 1 ; %output error at each Newton iteration (1: Yes; 0: No)
nITR = 0 ; %number of iterations for Newton

while errINNER >= TOL
    nITR = nITR + 1;
    nonlinearMAT = nonlinear_matrix_builder(vde,vde2, emat, remat,n,L,nxk,xi,lambda,guessU);
    vstep = seek_residue(n,lambda,guessU,nonlinearMAT,rhs_vector);
    alpha = 1; %alpha = 1 has fastest convergence (typically about 5 iterations)
    errINNER = max(abs(vstep));
    finalU = guessU + alpha*vstep; %increment the guess to approach the solution
    guessU = finalU; %set your approximation as your guess for the next step
    %     fprintf('finalU(n): %f \n',finalU(n)); % You can verify here that your
    % last step produces a solution vector that satisfies the RHS homogeneous Dirichlet
    % boundary condition
end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This subroutine computes the RHS in the 1st step of
% the Inexact Newton's Method, that is,
% -F(u_{n}) = f - L*u_{n} - N(u_{n})
% L is the coefficient matrix (that is the subroutine PAMATRIX)
% N(u_{n}) = -lambda^2*sinh(u_{n})
% AU = L*u_{n}
% U = INITIAL GUESS = u_{n}
% F = MODIFIED RHS VECTOR
% RESIDUE = -F(u_{n}) = f - L*u_{n} - N(u_{n})
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function vstep = seek_residue(n,lambda,guessU,nonlinearMAT,rhs_vector)
rhsau = zeros(n,1);
residue = zeros(n,1);
for i = 1:n-1
    rhsau(i) = -(lambda^2)*cosh(guessU(i))*guessU(i) + (lambda^2)*sinh(guessU(i));
end

[i,j,s] = find(nonlinearMAT);
nonlinearMAT = sparse(i,j,s);
AU = nonlinearMAT*guessU;


for i= 1:n
    residue(i) = rhs_vector(i)-AU(i)+rhsau(i);
end
vstep = nonlinearMAT \ residue;
end
function nonlinearMAT = nonlinear_matrix_builder(vde,vde2, emat, remat,n,L,nxk,xi,lambda,guessU)
%vde is vector of 2nd derivative weights
%vde2 is vector of first derivative weights

nonlinearMAT = zeros(n,n);

for i = 1:(n-1)
    for k = max(-i+nxk+2,1):min(n+nxk+1-i,2*nxk+1)
        nonlinearMAT(i,i+(k-(nxk+1))) = nonlinearMAT(i,i+(k-(nxk+1))) +vde(2,k)+(2/xi(i))*vde2(1,k);
    end
end

% %left boundary
for i = 1:(n-1)
    for j = 1:(nxk-(i-1))
        IP = nxk-(i+j-2);
        for k =1:L+1
            nonlinearMAT(i,k) = nonlinearMAT(i,k)+ (vde(2,j) + (2/xi(i))*vde2(1,j))*emat(IP,k);
        end
    end
end

%main diagonal involving lambda^2*u_i
for i = 1:(n-1)
    nonlinearMAT(i,i) = nonlinearMAT(i,i)-lambda^2*cosh(guessU(i));
end

%right bdy folding
for i = 1:n-1
    for j = n-i+1:nxk
        IP = i+j-n;
        for k = 1:L
            nonlinearMAT(i,n-k) =  nonlinearMAT(i,n-k)+(vde(2,j+nxk+1) + (2/xi(i))*vde2(1,j+nxk+1))*remat(IP,k);
        end
    end
end

%explicitly enforce right boundary condition
for i = n:n
    nonlinearMAT(i,i) = nonlinearMAT(i,i) + 1;
end

end



function c = weights(x,z,m)
%Input Parameters
% z location where approximations are to be accurate,
% x(1:nd) grid point locations, found in x(1:n)
% n = length of x
% x(1:nd) and c(1:nd,1:m+1)
% m highest derivative for which weights are sought,
% Output Parameter
% c(1:nd,1:m+1) weights at grid locations x(1:n) for derivatives
% of order 1:m+1, found in c(1:n,1:m+1)
%%SA: the weights are arranged column wise.
%%SA: for example, the weigths for FD approximation of order 2 will be
%%found in the last column of c


n = length(x);
nd  = n;
c = zeros(nd,m+1);
c1 = 1;
c4 = x(1)-z;

c(1,1) = 1;
for i = 1:n-1% in the fortran code, this goes from 1 to 2 if m=2. This will also go from 1 to 2 since n=3
    ix = i+1;% Since, we want to stick to the fortran code, we have to add 1 to i so we can get the first index of x as 1
    % and not zero for c4 = x(ix)-z.
    mn = min(i,m);%
    c2 = 1;
    c5 = c4;
    c4 = x(ix)-z;%
    for j = 0:i-1%
        jx = j+1;% See line 22 for similar explantion for  c3 = x(ix)-x(jx)
        c3 = x(ix)-x(jx);%
        c2 = c2*c3;
        if j == i-1%
            for k = mn:-1:1%
                kx = k+1;% See line 22 for similar explantion to line 35
                c(ix,kx) = c1*(k*c(ix-1,kx-1)-c5*c(ix-1,kx))/c2;%
            end
            c(ix,1) = -c1*c5*c(ix-1,1)/c2;%
        end
        for k = mn:-1:1%
            kx = k+1;%
            c(jx,kx) = (c4*c(jx,kx)-k*c(jx,kx-1))/c3;%
        end
        c(jx,1) = c4*c(jx,1)/c3;%
    end
    c1 = c2;
end
%c = c';%
%c was equated to c' in order to arrange the weights row-wise and not column-wise
end

function next_step = next_step_generator(n,xstart,epsa,epsw,bulk_strength,charge)
if bulk_strength <= 2
    xinf = 65;
else
    xinf = 20;
end
xend = xinf*xstart;
dx = (xend-xstart)/(n);
kappa = 8.430325455*bulk_strength/epsw; %changes depending on the ionic strength
lambda = sqrt(kappa);
ec2_kbt = (332.0716/0.5961574);
nxk = 4; %number of fictitious points
L = 8;
norder = 2; %2nd derivative
norder1 = 1; %1st derivative
phi = -(ec2_kbt*charge)/(epsw*(xstart^2)); %boundary condition on the left
psi = 0; %boundary condition on the right
xi = zeros(n+1,1);
for i = 1:n+1
    xi(i) = xstart + (i-1)*dx;
end

%FD coefficients for 2nd and 1st derivatives
vde = lag4(xstart,xend,n+1,nxk,norder);
vde2 = lag4(xstart,xend,n+1,nxk,norder1);

%weights in leftidm and rightidm
m0 = 1;
emat = leftidm(dx,m0,L);
remat = rightidm(dx,m0,L);
if nxk<=L
    m=nxk;
    emat_hdm = lefthdm(dx,m0,L,m,emat);
    remat_hdm = righthdm(dx,m0,L,m,remat);
end

%righthand side vector
rhs_vector = RHS_Vector(emat_hdm,remat_hdm,vde,vde2,xi,nxk,phi,psi,n+1,L);

guessU = zeros(n+1,1);
finalU = nonlinear_solver(guessU,rhs_vector,vde,vde2,emat_hdm,remat_hdm,n+1,L,nxk,xi,lambda);
% compute free energy
%potential_a = finalU(1);
%fprintf('u(a) is equal to %17.15f \n', potential_a);
UNC1 = finalU(1) - ((ec2_kbt*charge)/(epsa*xstart));%C1 = U^{+}(a) - (Aq/(e^{-}a))
numerical_energy1 = (0.5)*(charge)*(UNC1)*(0.5961574);%!(1/2)*q*C1*KBT
%fprintf('Numerical energy is %18.15f \n', numerical_energy1);
next_step = numerical_energy1;
end


