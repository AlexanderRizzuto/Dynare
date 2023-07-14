% Author: Federico Alexander Rizzuto
% Date: July 13, 2023
% This code solves and simulates the basic NK workhorse model presented in 
% chapters 3-4 of Galì (2008), with the least number of equations and 
% parameters.
% The model allows for two shocks (technology and monetary policy) and for
% three Taylor rules: a simple one (p.50), an aggressive one (p.77) and a 
% forward-looking one (79). 
% Note: already linearized around s.s.; rates not annualized (IRFs not 
% directly comparable to those found in the book).

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Declaring variables and parameters %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------- endogenous variables (6) ----------------
var 
    pi      % inflation
    ygap    % output gap y_t - y^nat_t
    i       % nominal interest rate
    r_n     % natural rate of interest 
    a       % productivity
    ni      % exogenous component in Taylor rule
;
% ---------------- exogenous variables (2) ----------------
varexo 
    eps_a   % productivity shock
    eps_ni  % m.p. shock
;
% ---------------- parameters () ----------------
parameters 
    p_beta      % discount factor
    p_kappa     % 
    p_sigma 
    psi_n_ya 
    rho_a 
    phi_pi 
    phi_y 
    rho_ni
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parametrization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_beta   = .99; % discount factor (p.52)
p_sigma  = 1;   % relative risk aversion (p.52)
rho_a    = .9;  % AR coeff. for technology shock (p.55)
phi_pi   = 1.5; % Taylor rule coeff. on inflation (p.52)
phi_y    = .5/4;% Taylor rule coeff. on output gap (p.52)
rho_ni   = .5;  % AR coeff. for m.p.shock ????
% derived and auxiliary
p_phi    = 1;   
p_alpha  = 1/3;
psi_n_ya = (1+p_phi) / (p_sigma*(1-p_alpha)+p_phi+p_alpha); % (Chapter 3, eq.18)
p_eps    = 6;   % demand elasticity
upptheta = (1-p_alpha) / (1-p_alpha+p_alpha*p_eps); % (p.47)
p_theta  = 2/3;  % index of price stickiness
p_lambda = ( ((1-p_theta)*(1-p_beta*p_theta)) / (p_theta) )*upptheta;
p_kappa  = p_lambda*(p_sigma + (p_phi+p_alpha)/(1-p_alpha)); % (p.49)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;
    % ---------------- NKPC (Ch.3, eq.21) ----------------
    pi   = p_beta*pi(+1) + p_kappa*ygap;        % vars: pi ygap; params: p_beta p_kappa 
    % ---------------- DIS (Ch.3, eq.22) ----------------
    ygap = -(1/p_sigma)*(i - pi(+1) - r_n) + ygap(+1); % vars: ygap i pi r_n ygap; params: p_sigma 
    % ---------------- natural rate of interest (Ch.3, eq.23; -rho) ----------------
    r_n  = p_sigma*psi_n_ya*(a(+1)-a);          % vars: r_n a; params: p_sigma psi_n_ya
    % ---------------- productivity AR1 (Ch.3, eq.28) ----------------
    a    = rho_a*a(-1) + eps_a;                 % vars: a; exovars: eps_a; params: rho_a
    % ---------------- simple Taylor rule (Ch.3, eq.25) ----------------
    %i    = phi_pi*pi + phi_y*ygap + ni;        % vars: i pi ygap; exovars: ni; params: phi_pi phi_y
    % ---------------- aggressive Taylor rule (Ch.4, eq.8) ----------------
    %i    = r_n + phi_pi*pi + phi_y*ygap + ni;  % vars: i r_n pi ygap; exovars: ni; params: phi_pi phi_y
    % ---------------- forward-looking Taylor rule (Ch.4, eq.8) ----------------
    %i    = r_n + phi_pi*pi(+1) + phi_y*ygap + ni; % vars: i r_n pi ygap; exovars: ni; params: phi_pi phi_y
    % ---------------- exogenous m.p. component AR1 (p.50) ----------------
    ni   = rho_ni*ni(-1) + eps_ni;              % exovars: ni eps_ni; params: rho_ni eps_ni
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Steady state %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resid;
steady;
check;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Shocks and simulation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------- 25 bp m.p. shock ----------------
shocks;
        var eps_ni = .25; 
        var eps_a = 0; 
end; 
stoch_simul(order = 1,irf=12,irf_plot_threshold =0); % approx. order: 1; 
                                    % generate 12-period IRFs for all vars;
                                    % plot even if IRF equal to zero

% ---------------- 1% productivity shock ----------------
shocks;
        var eps_ni = 0; 
        var eps_a = 1; 
end; 
stoch_simul(order = 1,irf=12,irf_plot_threshold =0); % approx. order: 1; 
                                    % generate 12-period IRFs for all vars;
                                    % plot even if IRF equal to zero