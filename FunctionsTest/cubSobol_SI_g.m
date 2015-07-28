function [q,int,out_param] = cubSobol_SI_g(varargin)
%cubSobol_SI_g Quasi-Monte Carlo method using Sobol' cubature over the
%d-dimensional region to compute the Sobol Index within a specified generalized error
%tolerance with guarantees under Walsh-Fourier coefficients cone decay
%assumptions
%
%   [q,out_param] = CUBSOBOL_SI_G(f,hyperbox,u) estimates the Sobol Index of f
%   over dimension u and region described by hyperbox, and with an error
%   guaranteed not to be greater than a specific generalized error tolerance,
%   tolfun:=max(abstol,reltol*| SI(f) |). Input f is a function handle. f should
%   accept an n x d matrix input, where d is the dimension and n is the 
%   number of points being evaluated simultaneously. The input hyperbox is
%   a 2 x d matrix, where the first row corresponds to the lower limits 
%   and the second row corresponds to the upper limits of the integral.
%   Given the construction of Sobol' sequences, d must be a positive 
%   integer with 1<=d<=1111.
%
%   q = CUBSOBOL_SI_G(f,hyperbox,u,measure,abstol,reltol)
%   estimates the Sobol Index of f over dimension u and hyperbox. The answer
%   is given within the generalized error tolerance tolfun. All parameters
%   should be input in the order specified above. If an input is not specified,
%   the default value is used. Note that if an input is not specified,
%   the remaining tail cannot be specified either. Inputs f and hyperbox 
%   are required. The other optional inputs are in the correct order:
%   measure,abstol,reltol,mmin,mmax,fudge,toltype and
%   theta.
%
%   q = CUBSOBOL_SI_G(f,hyperbox,u,'measure',measure,'abstol',abstol,'reltol',reltol)
%   estimates the Sobol Index of f over dimension u and hyperbox. The answer
%   is given within the generalized error tolerance tolfun. All the field-value
%   pairs are optional and can be supplied in any order. If an input is not
%   specified, the default value is used.
%
%   q = CUBSOBOL_SI_G(f,hyperbox,u,in_param) estimates the Sobol Index of f 
%   over dimension u and hyperbox. The answer is given within the
%   generalized error tolerance tolfun.
% 
%   Input Arguments
%
%     f --- the integrand whose input should be a matrix n x d where n is
%     the number of data points and d the dimension, which cannot be
%     greater than 1111. By default f is f=@ x.^2.
%
%     hyperbox --- the integration region defined by its bounds. It must be
%     a 2 x d matrix, where the first row corresponds to the lower limits 
%     and the second row corresponds to the upper limits of the integral.
%     The default value is [0;1].
%
%     u --- Integer that defines the dimensions to be checked. It needs to
%     be different than 0 and the absolute value of u needs to be smaller
%     or equal to the dimension of the function. Negative value of u means
%     all the dimensions except u. By default is 1.
%
%     in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%     measure of a uniformly distributed random variable in the hyperbox
%     or normally distributed with covariance matrix I_d. The only possible
%     values are 'uniform' or 'normal'. For 'uniform', the hyperbox must be
%     a finite volume while for 'normal', the hyperbox can only be defined as 
%     (-Inf,Inf)^d. By default it is 'uniform'.
%
%     in_param.abstol --- the absolute error tolerance, abstol>=0. By 
%     default it is 1e-4.
%
%     in_param.reltol --- the relative error tolerance, which should be
%     in [0,1]. Default value is 1e-2.
% 
%   Optional Input Arguments
% 
%     in_param.mmin --- the minimum number of points to start is 2^mmin.
%     The cone condition on the Fourier coefficients decay requires a
%     minimum number of points to start. The advice is to consider at least
%     mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
%     default it is 10.
%
%     in_param.mmax --- the maximum budget is 2^mmax. By construction of
%     the Sobol' generator, mmax is a positive integer such that
%     mmin<=mmax<=53. The default value is 24.
%
%     in_param.fudge --- the positive function multiplying the finite 
%     sum of Fast Walsh Fourier coefficients specified in the cone of functions.
%     This input is a function handle. The fudge should accept an array of
%     nonnegative integers being evaluated simultaneously. For more
%     technical information about this parameter, refer to the references.
%     By default it is @(m) 5*2.^-m.
%
%     in_param.toltype --- this is the generalized tolerance function.
%     There are two choices, 'max' which takes
%     max(abstol,reltol*| integral(f) | ) and 'comb' which is the linear combination
%     theta*abstol+(1-theta)*reltol*| integral(f) | . Theta is another 
%     parameter to be specified with 'comb'(see below). For pure absolute
%     error, either choose 'max' and set reltol = 0 or choose 'comb' and set
%     theta = 1. For pure relative error, either choose 'max' and set 
%     abstol = 0 or choose 'comb' and set theta = 0. Note that with 'max',
%     the user can not input abstol = reltol = 0 and with 'comb', if theta = 1
%     abstol con not be 0 while if theta = 0, reltol can not be 0.
%     By default toltype is 'max'.
% 
%     in_param.theta --- this input is parametrizing the toltype 
%     'comb'. Thus, it is only active when the toltype
%     chosen is 'comb'. It establishes the linear combination weight
%     between the absolute and relative tolerances
%     theta*abstol+(1-theta)*reltol*| integral(f) |. Note that for theta = 1, 
%     we have pure absolute tolerance while for theta = 0, we have pure 
%     relative tolerance. By default, theta=1.
%
%   Output Arguments
%
%     q --- the estimated value of the Sobol Index.
%
%     out_param.d --- dimension over which the algorithm integrated.
%
%     out_param.n --- number of Sobol' points used for computing the
%     integral of f.
%
%     out_param.bound_err --- predicted bound on the error based on the cone
%     condition. If the function lies in the cone, the real error will be
%     smaller than generalized tolerance.
%
%     out_param.time --- time elapsed in seconds when calling cubSobol_SI_g.
%
%     out_param.exitflag --- this is a binary vector stating whether
%     warning flags arise. These flags tell about which conditions make the
%     final result certainly not guaranteed. One flag is considered arisen
%     when its value is 1. The following list explains the flags in the
%     respective vector order:
%
%                       1 : If reaching overbudget. It states whether
%                       the max budget is attained without reaching the
%                       guaranteed error tolerance.
%      
%                       2 : If the function lies outside the cone. In
%                       this case, results are not guaranteed. For more
%                       information about the cone definition, check the
%                       article mentioned below.
% 
%  Guarantee
% This algorithm computes the Sobol Index of real valued functions in [0,1)^d
% with a prescribed generalized error tolerance. The Walsh-Fourier
% coefficients of the integrand are assumed to be absolutely convergent. If
% the algorithm terminates without warning messages, the output is given
% with guarantees under the assumption that the integrand lies inside a
% cone of functions. The guarantee is based on the decay rate of the
% Walsh-Fourier coefficients. For more details on how the cone is defined,
% please refer to the references below.
% 
%  Examples
% 
% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2:
% 
% >> f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)]; 
% >> q = cubSobol_SI_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
% >> check = abs(exactsol-q) < 1e-5
% check = 1
% 
% 
% Example 2:
% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2
% in the interval R^3 where x1, x2 and x3 are normally distributed:
% 
% >> f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
% >> q = cubSobol_SI_g(f,hyperbox,'normal',1e-3,1e-3); exactsol = 1;
% >> check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max')
% check = 1
% 
% 
% Example 3: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [-1,2)^2:
% 
% >> f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
% >> q = cubSobol_SI_g(f,hyperbox,'uniform',1e-3,1e-2); exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;
% >> check = abs(exactsol-q) < gail.tolfun(1e-3,1e-2,1,exactsol,'max')
% check = 1
%
%
% Example 4: 
% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.
% 
% >> f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); hyperbox = [-inf(1,1);inf(1,1)];
% >> q = cubSobol_SI_g(f,hyperbox,'normal',1e-4,1e-2); price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);
% >> check = abs(price-q) < gail.tolfun(1e-4,1e-2,1,price,'max')
% check = 1
%
%
% Example 5:
% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the interval
% [0,1)^5 with pure absolute error 1e-5.
% 
% >> f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
% >> q = cubSobol_SI_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
% >> check = abs(exactsol-q) < 1e-5
% check = 1
%
%
%   See also CUBLATTICE_G, CUBMC_G, MEANMC_G, MEANMCBER_G, INTEGRAL_G
% 
%  References
%
%   [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama, "Reliable adaptive
%   cubature using digital sequences," 2014. Submitted for publication:
%   arXiv:1410.8615.
%
%   [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   GAIL: Guaranteed Automatic Integration Library (Version 2.1)
%   [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://code.google.com/p/gail/ 
%
%   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
%   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
%   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
%   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
%   Workshop On Sustainable Software for Science: Practice And Experiences
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%

tic
%% Initial important cone factors and Check-initialize parameters
r_lag = 4; %distance between coefficients summed and those computed
[f,hyperbox,out_param] = cubSobol_SI_g_param(r_lag,varargin{:});
u = out_param.u;
l_star = out_param.mmin - r_lag; % Minimum gathering of points for the sums of DFWT
if u > 0
    S = @(b,c) 1 - c(:,1)./(b(:,2)-c(:,3).^2);
else
    S = @(b,c) b(:,1)./(c(:,2)-b(:,3).^2);
end

if strcmp(out_param.measure,'normal')
   f=@(x) f(gail.stdnorminv(x));
elseif strcmp(out_param.measure,'uniform')
   Cnorm = prod(hyperbox(2,:)-hyperbox(1,:));
   f=@(x) Cnorm*f(bsxfun(@plus,hyperbox(1,:),bsxfun(@times,(hyperbox(2,:)-hyperbox(1,:)),x))); % a + (b-a)x = u
end

%% Main algorithm
sobstr = sobolset(2*out_param.d); %generate a Sobol' sequence
sobstr = scramble(sobstr,'MatousekAffineOwen'); %scramble it
Stilde = zeros(out_param.mmax-out_param.mmin+1,3); %initialize sum of DFWT terms
CStilde_low = -inf(out_param.mmax-l_star+1,3); %initialize various sums of DFWT terms for necessary conditions
CStilde_up = inf(out_param.mmax-l_star+1,3); %initialize various sums of DFWT terms for necessary conditions
errest = zeros(out_param.mmax-out_param.mmin+1,3); %initialize error estimates
appxinteg = zeros(out_param.mmax-out_param.mmin+1,1); %initialize approximations to integral
exit_len = 2;
out_param.exit = false(1,exit_len); %we start the algorithm with all warning flags down

%% Initial points and FWT
out_param.n = 2^out_param.mmin; %total number of points to start with
n0 = out_param.n; %initial number of points
xpts = sobstr(1:n0,1:2*out_param.d); %grab Sobol' points
y(:,2) = f(xpts(:,1:out_param.d)).^2; %evaluate integrands
y(:,3) = f(xpts(:,1:out_param.d)); %evaluate integrands
if u > 0
    y(:,1) = y(:,2) - y(:,3).*f([xpts(:,out_param.d+1:out_param.d+u-1) xpts(:,u) xpts(:,out_param.d+u+1:2*out_param.d)]);
else
    y(:,1) = y(:,2) - y(:,3).*f([xpts(:,1:+abs(u)-1) xpts(:,out_param.d+abs(u)) xpts(:,abs(u)+1:out_param.d)]);
end
yval = y;

%% Compute initial FWT
nllstart = int64(2^(out_param.mmin-r_lag-1));
for k = 1:size(y,2)
    for l=0:out_param.mmin-1
       nl=2^l;
       nmminlm1=2^(out_param.mmin-l-1);
       ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
       evenval=y(ptind,k);
       oddval=y(~ptind,k);
       y(ptind,k)=(evenval+oddval)/2;
       y(~ptind,k)=(evenval-oddval)/2;
    end
    %y now contains the FWT coefficients 

    %% Create kappanumap implicitly from the data
    kappanumap(:,k) = (1:out_param.n)'; %initialize map
    for l=out_param.mmin-1:-1:1
       nl=2^l;
       oldone=abs(y(kappanumap(2:nl,k),k)); %earlier values of kappa, don't touch first one
       newone=abs(y(kappanumap(nl+2:2*nl,k),k)); %later values of kappa, 
       flip=find(newone>oldone); %which in the pair are the larger ones
       if ~isempty(flip)
           flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
           flipall=flipall(:);
           temp=kappanumap(nl+1+flipall,k); %then flip 
           kappanumap(nl+1+flipall,k)=kappanumap(1+flipall,k); %them
           kappanumap(1+flipall,k)=temp; %around
       end
    end

    %% Compute Stilde
    Stilde(1,k) = sum(abs(y(kappanumap(nllstart+1:2*nllstart,k),k)));
    err_bound_k(1,k) = out_param.fudge(out_param.mmin)*Stilde(1,k);
    est_int_k(1,k) = mean(yval(:,k));
end

q = 1/2*(S(est_int_k+err_bound_k , est_int_k-err_bound_k) + S(est_int_k-err_bound_k , est_int_k+err_bound_k));
out_param.bound_err = 1/2*(S(est_int_k+err_bound_k , est_int_k-err_bound_k) - S(est_int_k-err_bound_k , est_int_k+err_bound_k));
errest(1) = out_param.bound_err;
int = est_int_k(1,3);

% Necessary conditions for all three integrals
for k = 1:size(y,2)
    for l = l_star:out_param.mmin % Storing the information for the necessary conditions
        C_low = (1+out_param.fudge(out_param.mmin-l))/(1+2*out_param.fudge(out_param.mmin-l));
        C_up = (1+out_param.fudge(out_param.mmin-l));
        CStilde_low(l-l_star+1,k) = C_low*sum(abs(y(kappanumap(2^(l-1)+1:2^l,k),k)));
        CStilde_up(l-l_star+1,k) = C_up*sum(abs(y(kappanumap(2^(l-1)+1:2^l,k),k)));
    end
    if any(CStilde_low(:,k) > CStilde_up(:,k))
       out_param.exit(2) = true;
    end
end

% Check the end of the algorithm
deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-errest(1)),...
    out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+errest(1)),out_param.toltype));
deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-errest(1)),...
    out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+errest(1)),out_param.toltype));

is_done = false;
q=q+deltaminus;
appxinteg(1) = q;
out_param.time=toc;
if out_param.bound_err <= deltaplus
   is_done = true;
elseif out_param.mmin == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
   out_param.exit(1) = true;
end


%% Loop over m
for m=out_param.mmin+1:out_param.mmax
    if is_done,
       break;
    end
    out_param.n=2^m;
    mnext=m-1;
    nnext=2^mnext;
    xnext=sobstr(n0+(1:nnext),1:2*out_param.d); 
    n0=n0+nnext;
    ynext(1:nnext,2) = f(xnext(:,1:out_param.d)).^2; %evaluate integrands
    ynext(1:nnext,3) = f(xnext(:,1:out_param.d)); %evaluate integrands
    if u > 0
        ynext(1:nnext,1) = ynext(:,2) - ynext(:,3).*f([xnext(:,out_param.d+1:out_param.d+u-1) xnext(:,u) xnext(:,out_param.d+u+1:2*out_param.d)]);
    else
        ynext(1:nnext,1) = ynext(:,2) - ynext(:,3).*f([xnext(:,1:+abs(u)-1) xnext(:,out_param.d+abs(u)) xnext(:,abs(u)+1:out_param.d)]);
    end
    yval=[yval; ynext];

   %% Compute initial FWT on next points
   nllstart=int64(2^(m-r_lag-1));
   meff=m-out_param.mmin+1;
   ny = size(y,1);
   for k = 1:size(y,2)
       for l=0:mnext-1
          nl=2^l;
          nmminlm1=2^(mnext-l-1);
          ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
          evenval=ynext(ptind,k);
          oddval=ynext(~ptind,k);
          ynext(ptind,k)=(evenval+oddval)/2;
          ynext(~ptind,k)=(evenval-oddval)/2;
       end

       %% Compute FWT on all points
       y(1:2*ny,k)=[y(1:ny,k);ynext(:,k)];
       nl=2^mnext;
       ptind=[true(nl,1); false(nl,1)];
       evenval=y(ptind,k);
       oddval=y(~ptind,k);
       y(ptind,k)=(evenval+oddval)/2;
       y(~ptind,k)=(evenval-oddval)/2;

       %% Update kappanumap
       kappanumap(1:2*ny,k) = [kappanumap(1:ny,k) ; 2^(m-1)+kappanumap(1:ny,k)]; %initialize map
       for l=m-1:-1:m-r_lag
          nl=2^l;
          oldone=abs(y(kappanumap(2:nl,k),k)); %earlier values of kappa, don't touch first one
          newone=abs(y(kappanumap(nl+2:2*nl,k),k)); %later values of kappa, 
          flip=find(newone>oldone);
          if ~isempty(flip)
              flipall=bsxfun(@plus,flip,0:2^(l+1):2^m-1);
              flipall=flipall(:);
              temp=kappanumap(nl+1+flipall,k);
              kappanumap(nl+1+flipall,k)=kappanumap(1+flipall,k);
              kappanumap(1+flipall,k)=temp;
          end
       end

       %% Compute Stilde
       Stilde(meff,k)=sum(abs(y(kappanumap(nllstart+1:2*nllstart,k),k)));
       err_bound_k(meff,k) = out_param.fudge(out_param.mmin)*Stilde(meff,k);
       est_int_k(meff,k) = mean(yval(:,k));
   end
   q = 1/2*(S(est_int_k(meff,:)+err_bound_k(meff,:) , est_int_k(meff,:)-err_bound_k(meff,:)) + S(est_int_k(meff,:)-err_bound_k(meff,:) , est_int_k(meff,:)+err_bound_k(meff,:)));
   out_param.bound_err = 1/2*(S(est_int_k(meff,:)+err_bound_k(meff,:) , est_int_k(meff,:)-err_bound_k(meff,:)) - S(est_int_k(meff,:)-err_bound_k(meff,:) , est_int_k(meff,:)+err_bound_k(meff,:)));
   errest(meff) = out_param.bound_err;
   int = est_int_k(meff,3);

   % Necessary conditions
   for k = 1:size(y,2)
       for l = l_star:m % Storing the information for the necessary conditions
            C_low = (1+out_param.fudge(m-l))/(1+2*out_param.fudge(m-l));
            C_up = (1+out_param.fudge(m-l));
            CStilde_low(l-l_star+1,k) = max(CStilde_low(l-l_star+1,k),C_low*sum(abs(y(kappanumap(2^(l-1)+1:2^l,k),k))));
            CStilde_up(l-l_star+1,k) = min(CStilde_up(l-l_star+1,k),C_up*sum(abs(y(kappanumap(2^(l-1)+1:2^l,k),k))));
       end
       if any(CStilde_low(:,k) > CStilde_up(:,k))
           out_param.exit(2) = true;
       end
   end
   
   % Check the end of the algorithm
    deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-errest(meff)),...
        out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+errest(meff)),out_param.toltype));
    deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-errest(meff)),...
        out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+errest(meff)),out_param.toltype));
   
   q=q+deltaminus;
   appxinteg(meff)=q;
   out_param.time=toc;
   if out_param.bound_err <= deltaplus
      is_done = true;
   elseif m == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
      out_param.exit(1) = true;
   end
end

% Decode the exit structure
exit_str=2.^(0:exit_len-1).*out_param.exit;
exit_str(out_param.exit==0)=[];
if numel(exit_str)==0;
    out_param.exitflag=0;
else
    out_param.exitflag=exit_str;
end
out_param = rmfield(out_param,'exit');

out_param.time=toc;
end


%% Parsing for the input of cubSobol_SI_g
function [f,hyperbox, out_param] = cubSobol_SI_g_param(r_lag,varargin)

% Default parameter values
default.u = 1;
default.hyperbox = [zeros(1,1);ones(1,1)];% default hyperbox
default.measure  = 'uniform';
default.abstol  = 1e-4;
default.reltol  = 1e-2;
default.mmin  = 10;
default.mmax  = 24;
default.fudge = @(m) 5*2.^-m;
default.toltype  = 'max';
default.theta  = 1;

if numel(varargin)<3
    help cubSobol_SI_g
    warning('GAIL:cubSobol_SI_g:fdnotgiven',...
        'At least, dimension u, function f and hyperbox need to be specified. Example for f(x) = x^2 and u = 1:')
    f = @(x) x.^2;
    out_param.f=f;
    out_param.u = default.u;
    hyperbox = default.hyperbox;
else
    f = varargin{1};
    if ~gail.isfcn(f)
        warning('GAIL:cubSobol_SI_g:fnotfcn',...
            'The given input f was not a function. Example for f(x) = x^2 and u = 1:')
        f = @(x) x.^2;
        out_param.f=f;
        out_param.u = default.u;
        hyperbox = default.hyperbox;
    else
        out_param.f=f;
        hyperbox = varargin{2};
        if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(size(hyperbox,2)<1111)
            warning('GAIL:cubSobol_SI_g:hyperbox_error1',...
                'The hyperbox must be a real matrix of size 2xd where d can not be greater than 1111. Example for f(x) = x^2 and u = 1:')
            f = @(x) x.^2;
            out_param.f=f;
            out_param.u = default.u;
            hyperbox = default.hyperbox;
        else
            u = varargin{3};
            if ~(ceil(u)==u) || abs(u) > size(hyperbox,2) || u == 0
                warning('GAIL:cubSobol_SI_g:u_error1',...
                    'Dimension u must be a an integer such that -d <= u <= d and u != 0. Example for f(x) = x^2 and u = 1:')
                f = @(x) x.^2;
                out_param.f=f;
                out_param.u = default.u;
                hyperbox = default.hyperbox;
            end
        end
    end
end

validvarargin=numel(varargin)>3;
if validvarargin
    in4=varargin(4:end);
    for j=1:numel(varargin)-3
    validvarargin=validvarargin && (isnumeric(in4{j}) ...
        || ischar(in4{j}) || isstruct(in4{j}) || gail.isfcn(in4{j}));
    end
    if ~validvarargin
        warning('GAIL:cubSobol_SI_g:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in4=varargin{4};
end

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
  f_addParamVal = @addParameter;
else
  f_addParamVal = @addParamValue;
end

if ~validvarargin
    out_param.measure = default.measure;
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;  
    out_param.fudge = default.fudge;
    out_param.toltype = default.toltype;
    out_param.theta = default.theta;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'hyperbox',@isnumeric);
    addRequired(p,'u',@isnumeric);
    if isnumeric(in4) || ischar(in4)
        addOptional(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@gail.isfcn);
        addOptional(p,'toltype',default.toltype,...
            @(x) any(validatestring(x, {'max','comb'})));
        addOptional(p,'theta',default.theta,@isnumeric);
    else
        if isstruct(in4) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'reltol',default.reltol,@isnumeric);
        f_addParamVal(p,'mmin',default.mmin,@isnumeric);
        f_addParamVal(p,'mmax',default.mmax,@isnumeric);
        f_addParamVal(p,'fudge',default.fudge,@gail.isfcn);
        f_addParamVal(p,'toltype',default.toltype,...
            @(x) any(validatestring(x, {'max','comb'})));
        f_addParamVal(p,'theta',default.theta,@isnumeric);
    end
    parse(p,f,hyperbox,u,varargin{4:end})
    out_param = p.Results;
end

out_param.d = size(hyperbox,2);

fdgyes = 0; % We store how many functions are in varargin. There can only
            % two functions as input, the function f and the fudge factor.
for j = 1:size(varargin,2)
    fdgyes = gail.isfcn(varargin{j})+fdgyes;
end
if fdgyes < 2 % No fudge factor given as input
    default.fudge = @(m) 5*2.^-(m/d);
end

%hyperbox should be 2 x dimension
if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(out_param.d<1111)
    warning('GAIL:cubSobol_SI_g:hyperbox_error2',...
        'The hyperbox must be a real matrix of size 2 x d where d can not be greater than 1111. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    hyperbox = default.hyperbox;
end

% Force measure to be uniform or normal only
if ~(strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') )
    warning('GAIL:cubSobol_SI_g:notmeasure',['The measure can only be uniform or normal.' ...
            ' Using default measure ' num2str(default.measure)])
    out_param.measure = default.measure;
end

% Force absolute tolerance greater than 0
if (out_param.abstol < 0 )
    warning('GAIL:cubSobol_SI_g:abstolnonpos',['Absolute tolerance cannot be negative.' ...
            ' Using default absolute tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force relative tolerance greater than 0 and smaller than 1
if (out_param.reltol < 0) || (out_param.reltol > 1)
    warning('GAIL:cubSobol_SI_g:reltolnonunit',['Relative tolerance should be chosen in [0,1].' ...
            ' Using default relative tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('GAIL:cubSobol_SI_g:lowmmin',['The minimum starting exponent ' ...
            'should be an integer greater than 0 and smaller or equal than the maxium.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force mmin to be integer greater than r_lag (so that l_star=mmin-r_lag>=0)
if out_param.mmin < r_lag
    warning('GAIL:cubSobol_SI_g:lowmminrlag',['The minimum starting exponent ' ...
            'should be at least ' num2str(r_lag) '.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 54
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=53)
    warning('GAIL:cubSobol_SI_g:wrongmmax',['The maximum exponent for the budget should be an integer biger than mmin and smaller than 54.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('GAIL:cubSobol_SI_g:fudgenonpos',['The fudge factor should be a positive function.' ...
            ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end

% Force toltype to be max or comb
if ~(strcmp(out_param.toltype,'max') || strcmp(out_param.toltype,'comb') )
    warning('GAIL:cubSobol_SI_g:nottoltype',['The error type can only be max or comb.' ...
            ' Using default toltype ' num2str(default.toltype)])
    out_param.toltype = default.toltype;
end

% Force theta to be in [0,1]
if (out_param.theta < 0) || (out_param.theta > 1)
    warning('GAIL:cubSobol_SI_g:thetanonunit',['Theta should be chosen in [0,1].' ...
            ' Using default theta ' num2str(default.theta)])
    out_param.theta = default.theta;
end

% Checking on pure absolute/relative error
if (out_param.abstol==0) && (out_param.reltol==0)
    warning('GAIL:cubSobol_SI_g:tolzeros',['Absolute and relative error tolerances can not be simultaniusly 0.' ...
            ' Using default absolute tolerance ' num2str(default.abstol) ' and relative tolerance ' num2str(default.reltol)])
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
end
if (strcmp(out_param.toltype,'comb')) && (out_param.theta==1) && (out_param.abstol==0)
    warning('GAIL:cubSobol_SI_g:abstolzero',['When choosing toltype comb, if theta=1 then abstol>0.' ...
            ' Using default absolute tolerance ' num2str(default.abstol) ])
    out_param.abstol = default.abstol;
end
if (strcmp(out_param.toltype,'comb')) && (out_param.theta==0) && (out_param.reltol==0)
    warning('GAIL:cubSobol_SI_g:reltolzero',['When choosing toltype comb, if theta=0 then reltol>0.' ...
            ' Using default relative tolerance ' num2str(default.reltol) ])
    out_param.reltol = default.reltol;
end

% Checking on the hyperbox given the measure
if (strcmp(out_param.measure,'uniform')) && ~all(all(isfinite(hyperbox)))
    warning('GAIL:cubSobol_SI_g:hyperboxnotfinite',['If uniform measure, hyperbox must be of finite volume.' ...
            ' Using default hyperbox:'])
    disp([zeros(1,out_param.d);ones(1,out_param.d)])
    hyperbox = [zeros(1,out_param.d);ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (sum(sum(isfinite(hyperbox)))>0)
    warning('GAIL:cubSobol_SI_g:hyperboxfinite',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (any(hyperbox(1,:)==hyperbox(2,:)) || any(hyperbox(1,:)>hyperbox(2,:)))
    warning('GAIL:cubSobol_SI_g:hyperboxnormalwrong',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end

end
