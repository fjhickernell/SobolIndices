clearvars

abstol = 5*1e-3;
pre = 2; % Precision to display according to 1*1e-3 (3 digits)
reltol = 0e-2; % Pure absolute tolerance
mmin = 9;
mmax = 24; % I adjust that not to run out of memory. It can go up to 54. Type help cubSobol_SI_g for more information.
wronga = 0; % Wrong estimates with all indices method
wrong = 0; % Wrong estimates using the one by one algo
s = .0; %time in seconds used delay each function evaluation
samples = 100;

fudge = @(m,d) 10*2.^-(1.*m);

%% Bratley et al.
disp('Running Bratley et al. example ...')
d = 6;
f = @(x) sum(bsxfun(@times, cumprod(x,2), (-1).^(1:d)),2) +  pause_t(s);
hyperbox = [-zeros(1,d) ; ones(1,d)];

R = [.6528636616, .1791303924, 0.3701041165e-1, 0.1332374820e-1, 0.1480416466e-2, 0.1480416466e-2];
SI = zeros(1, d);
SI_n = SI;
error = SI;
reliability = SI;
error_all = [];
for k = 1:samples
    [q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
    SI = SI + q;
    SI_n = SI_n + out_param.n;
    error = error + abs(R-q);
    error_all = [error_all; abs(R-q)];
    reliability = reliability + 1*(abs(R-q) > max(abstol,R*reltol));
end
SI = SI/samples; SI_n = SI_n/samples; reliability = 1 - reliability/samples; error = error/samples;
round(SI, pre, 'significant')
round(SI_n)
round(error, pre, 'significant')
round(reliability*100)

% csvwrite('bratley_replicated.csv', error_all)



%% Sobol' g-function
disp('Running Sobol g-function example ...')
d = 6;
a = [0 1/2 3 9 99 99];
f = @(x) prod(bsxfun(@(x,a) (abs(4*x-2)+a)./(1+a), x , a),2) +  pause_t(s);
hyperbox = [zeros(1,d) ; ones(1,d)];

R = [.5867811893, .2607916397, 0.3667382378e-1, 0.5867811312e-2, 0.5867753221e-4, 0.5867753221e-4];
SI = zeros(1, d);
SI_n = SI;
error = SI;
reliability = SI;
error_all = [];
for k = 1:samples
    [q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
    SI = SI + q;
    SI_n = SI_n + out_param.n;
    error = error + abs(R-q);
    error_all = [error_all; abs(R-q)];
    reliability = reliability + 1*(abs(R-q) > max(abstol,R*reltol));
end
SI = SI/samples; SI_n = SI_n/samples; reliability = 1 - reliability/samples; error = error/samples;
% disp(SI)
% disp(SI_n)
% disp(error)
round(SI, pre,'significant')
round(SI_n)
round(error, pre, 'significant')
round(reliability*100)

% csvwrite('sobolg_replicated.csv', error_all)