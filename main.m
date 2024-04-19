% Main script
clear;
close all;
clc;
% Definition of parameters
E0 =0 ; % Peak resonance energy
Gamma = 1e-3; % Level broadening
Emin = -0.01; % Lower integration limit
Emax =0.01; % Upper integration limit
NE = 20; % Number of equidistant intervals for non-adaptive grid
%% maxError = ; % Input maximum error threshold for adaptive method

% Lorentzian function for electron concentration
f = @(E) (Gamma/2)./pi./((E - E0).^2 + (Gamma/2).^2); %% for fermi = 1

% True electron concentration using MATLAB's integral function
n_true = integral(f, Emin, Emax)

% Non-adaptive Grid Integration
[n_approx, rel_error_non_adaptive] = nonAdaptiveGridIntegration(f, Emin, Emax, NE, n_true);

% Output results for non-adaptive method
fprintf('Non-adaptive Electron Concentration: %f\n', n_approx);
fprintf('Relative error (Non-adaptive method): %f%%\n', rel_error_non_adaptive * 100);
fprintf('Number of sub-intervals for Non-Adaptative: %d\n', NE);

% Adaptive Grid Integration with a Max Error fixed
maxError = rel_error_non_adaptive; % Input maximum error threshold for adaptive method
[n_adaptive_err, subIntervals_err, subIntervalSizes_err] = adaptiveTrapezoidal(f, Emin, Emax, maxError);
rel_error_adaptive_err = (n_adaptive_err-n_true)/n_true;

% Output results for adaptive method with a Max Error fixed
% fprintf('\nAdaptive Electron Concentration with fixed maximum error: %f\n', n_adaptive_err);
% fprintf('Relative error (adaptive method with fixed maximum error): %f%%\n', rel_error_adaptive_err * 100);
% fprintf('Number of sub-intervals for Adaptative with fixed maximum error: %d\n', subIntervals_err);
%disp('Size of each sub-interval with fixed maximum error:');
%disp(subIntervalSizes_err);

% Adaptive Grid Integration with a Max Number of points fixed
maxPoints = NE; % Input maximum error threshold for adaptive method
rel_error_adaptive_pts = adaptiveTrapezoidalFixedNbPts(f, Emin, Emax, maxPoints);
n_adaptative_pts= rel_error_adaptive_pts*n_true+n_true;

% Output results for adaptive method with a Max Number of points fixed
fprintf('\nAdaptive Electron Concentration with fixed Nb pts: %f\n', n_adaptative_pts);
fprintf('Relative error (adaptive method with fixed Nb pts): %f%%\n', rel_error_adaptive_pts * 100);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize arrays to store results
 NE_values1 = [25, 30, 35,40,45, 50,55, 100, 200, 400];
 NE_values2 = [10,15,20, 30,40, 50, 100]; 
% 
 plotGraphFixedNbPts(NE_values1,f,Emin,Emax);
 %plotGraphFixedError(NE_values2,f,Emin,Emax,n_true)

NE_values3 = [10,15,20,30,40,80,160,200];
%[EstimatedErrorRichardson, matchCounterRichardson] = adaptiveTrapezoidalRichardson(f, Emin, Emax, maxPoints)
%[EstimatedErrorSecondDerivative, matchCounterSecondDerivative] = adaptiveTrapezoidalSecondDerivative(f, Emin, Emax, maxPoints)
%[EstimatedErrorSimson, matchCounter] = adaptiveQuadratureSimpson(f, Emin, Emax, maxPoints)
%plotErrorEstimateAccuracy(f, Emin, Emax, NE_values3)
