% Author: Alexander Rizzuto
% Date: July 15, 2023
% This code solves the basic NK workhorse model presented in chapters 3-4 
% of Galì (2008), with the least number of equations and parameters, 
% computes the IRFs and performs global stability analysis.
% The model allows for two shocks (technology and monetary policy) and for
% three Taylor rules: a simple one (p.50), an aggressive one (p.77) and a 
% forward-looking one (79). 
% Note: rates not annualized (IRFs not directly comparable to those found 
% in the book). 
% Note: Taylor rules are commented out so as to allow to run this file from
% the file ch4_nk_taylor_minimal_compare.m. To run this file on its own 
% uncomment one of the Taylor Rules

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Declaring variables and parameters %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------- endogenous variables (6) ----------------
var 
    pi      ${\pi}$ (long_name='Inflation Rate')
    ygap    ${\tilde y}$ (long_name='Output Gap')
    i       ${i}$ (long_name='Nominal Interest Rate')
    r_n     ${r^{n}}$ (long_name='Natural Rate of Interest')
    a       ${a}$ (long_name='Productivity')
    ni      ${\nu}$ (long_name='Exogenous Component in Taylor Rule')
;
% ---------------- exogenous variables (2) ----------------
varexo 
    eps_a   ${\varepsilon_a}$ (long_name='Productivity Shock')
    eps_ni  ${\varepsilon_{\nu}}$ (long_name='M.P. Shock')
;
% ---------------- parameters ----------------
parameters 
    p_beta      ${\beta}$ (long_name='Discount Factor')
    p_kappa     ${\kappa}$ (long_name='Slope Coefficient of NKPC')
    p_sigma     ${\sigma}$ (long_name='Risk Aversion')
    psi_n_ya    ${\psi^n_{ya}}$ (long_name='Composite Parameter (Natural Output)')
    rho_a       ${\rho_a}$ (long_name='Productivity Autocorrelation')
    phi_pi      ${\phi_{\pi}}$ (long_name='Taylor Rule Coefficient on Inflation')
    phi_y       ${\phi_{y}}$ (long_name='Taylor Rule Coefficient on Output Gap')
    rho_ni      ${\rho_{\nu}}$ (long_name='Autocorrelation of Exogenous M.P. Component')
    p_phi       ${\phi}$ (long_name='Frisch Elasticity of Labor Supply')
    p_eps       ${\varepsilon}$ (long_name='Demand Elasticity')
    p_theta     ${\theta}$ (long_name='Degree of Price Stickiness')
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parametrization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_beta   = .99; % (p.52)
p_sigma  = 1;   % (p.52)
rho_a    = .9;  % (p.55)
phi_pi   = 1.5; % (p.52)
phi_y    = .5/4;% (p.52)
rho_ni   = .5;  % (p.52)

% ---------------- derived and auxiliary ----------------
p_phi    = 1;   % (p.52)
p_alpha  = 1/3; % (p.52)
psi_n_ya = (1+p_phi) / (p_sigma*(1-p_alpha)+p_phi+p_alpha); % (p.48)
p_eps    = 6;   % (p.52)
upptheta = (1-p_alpha) / (1-p_alpha+p_alpha*p_eps); % (p.47)
p_theta  = 2/3;  % (p.52) "implies average price duration of three quarters"
p_lambda = ( ((1-p_theta)*(1-p_beta*p_theta)) / (p_theta) )*upptheta; % (p.47)
p_kappa  = p_lambda*(p_sigma + (p_phi+p_alpha)/(1-p_alpha)); % (p.49)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------- equations (8) ----------------

model;
    [name='New-Keynesian Phillips Curve'] % (Ch.3, eq.21)
    pi   = p_beta*pi(+1) + p_kappa*ygap;        % vars: pi ygap; params: p_beta p_kappa 
    [name='Dynamic IS Curve'] % (Ch.3, eq.22) 
    ygap = -1/p_sigma*(i - pi(+1) - r_n) + ygap(+1); % vars: ygap i pi r_n ygap; params: p_sigma 
    [name='Natural Rate of Interest'] % (Ch.3, eq.23; -rho) 
    r_n  = p_sigma*psi_n_ya*(a(+1)-a);          % vars: r_n a; params: p_sigma psi_n_ya
    [name='Productivity - AR1 Process'] % (Ch.3, eq.28)
    a    = rho_a*a(-1) + eps_a;                 % vars: a; exovars: eps_a; params: rho_a
    % ---------------- simple Taylor rule (Ch.3, eq.25) ----------------
    [name='Monetary Policy - Taylor Rule'] % 
    @#include "taylor_rules.txt"
    % ---------------- simple Taylor rule (Ch.3, eq.25) ----------------
    %i    = phi_pi*pi + phi_y*ygap + ni;        % vars: i pi ygap; exovars: ni; params: phi_pi phi_y
    % ---------------- aggressive Taylor rule (Ch.4, eq.8) ----------------
    %i    = r_n + phi_pi*pi + phi_y*ygap + ni;  % vars: i r_n pi ygap; exovars: ni; params: phi_pi phi_y
    % ---------------- forward-looking Taylor rule (Ch.4, eq.8) ----------------
    %i    = r_n + phi_pi*pi(+1) + phi_y*ygap + ni; % vars: i r_n pi ygap; exovars: ni; params: phi_pi phi_y
    % ---------------- exogenous m.p. component AR1 (p.50) ----------------
    [name='Monetary Policy Shock - AR1']
    ni   = rho_ni*ni(-1) + eps_ni;  

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Steady state %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resid;
steady;
check;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shocks and IRFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------- m.p. shock ----------------
shocks;
        var eps_ni = .25; 
        var eps_a = 0; 
end; 
stoch_simul(order = 1,irf=12,irf_plot_threshold =0,graph_format = pdf);
% approx. order: 1; generate 12-period IRFs for all vars; plot even if IRF 
% equal to zero

% ---------------- productivity shock ----------------
shocks;
        var eps_ni = 0; 
        var eps_a = 1; 
end; 
stoch_simul(order = 1,irf=12,irf_plot_threshold =0,graph_format = pdf); 
% approx. order: 1; generate 12-period IRFs for all vars; plot even if IRF 
% equal to zero


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Global Sensitivity Analysis %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimated_params; 
    phi_y       , uniform_pdf, , , 0,2;
    phi_pi       , uniform_pdf, , , 0,2;
    rho_ni       , uniform_pdf, , , 0,1;
    p_beta       , uniform_pdf, , , 0,1;
    p_sigma       , uniform_pdf, , , 0,2;
    p_kappa       , uniform_pdf, , , 0,1;
end;
varobs pi;
dynare_sensitivity(prior_range=0,stab=1,Nsam=5000);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Latex Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_latex_dynamic_model;
write_latex_static_model;
write_latex_definitions;
write_latex_parameter_table;
collect_latex_files;