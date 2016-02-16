clearvars

abstol = 1e-4;
reltol = 0; % Pure absolute tolerance
mmax = 22; % I adjust that not to run out of memory. It can go up to 54. Type help cubSobol_SI_g for more information.

%% Ishigami
disp('Running Ishigami example ...')
d = 3;
f = @(x) sin(x(:,1)).*(1+1/10*(x(:,3).^4))+7*sin(x(:,2)).^2;
hyperbox = pi*[-ones(1,d) ; ones(1,d)];

exitflag = [];
exitflagt = [];
for j = 1:d
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflag = [exitflag out_param.exitflag];
    SI(1,j) = q;
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflagt = [exitflagt out_param.exitflag];
    SI(2,j) = q;
end
disp(SI)
if any(exitflag > 0) || any(exitflagt > 0)
    warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
end
    



%% Bratley et al.
disp('Running Bratley et al. example ...')
d = 6;
f = @(x) sum(bsxfun(@times, cumprod(x,2), (-1).^(1:d)),2);
hyperbox = [zeros(1,d) ; ones(1,d)];

exitflag = [];
exitflagt = [];
tic
for j = 1:d
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflag = [exitflag out_param.exitflag];
    SI(1,j) = q;
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflagt = [exitflagt out_param.exitflag];
    SI(2,j) = q;
end
toc
disp(SI)
if any(exitflag > 0) || any(exitflagt > 0)
    warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
end
    



%% Sobol' g-function
disp('Running Sobol g-function example ...')
d = 6;
a = [0 1/2 3 9 99 99];
f = @(x) prod(bsxfun(@(x,a) (abs(4*x-2)+a)./(1+a), x , a),2);
hyperbox = [zeros(1,d) ; ones(1,d)];

exitflag = [];
exitflagt = [];
for j = 1:d
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflag = [exitflag out_param.exitflag];
    SI(1,j) = q;
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflagt = [exitflagt out_param.exitflag];
    SI(2,j) = q;
end
disp(SI)
if any(exitflag > 0) || any(exitflagt > 0)
    warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
end
    



%% Morokoff and Caflish
disp('Running Morokoff and Caflish example ...')
d = 6;
f = @(x) (1+1/d)^d*prod(x,2).^(1/d);
hyperbox = [zeros(1,d) ; ones(1,d)];

exitflag = [];
exitflagt = [];
for j = 1:d
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflag = [exitflag out_param.exitflag];
    SI(1,j) = q;
    [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    exitflagt = [exitflagt out_param.exitflag];
    SI(2,j) = q;
end
disp(SI)
if any(exitflag > 0) || any(exitflagt > 0)
    warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
end