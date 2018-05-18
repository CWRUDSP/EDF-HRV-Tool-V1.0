%--------------------------------------------------------------------------
% @license
% Copyright 2018 IDAC Signals Team, Case Western Reserve University 
%
% Lincensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public 
% you may not use this file except in compliance with the License.
%
% Unless otherwise separately undertaken by the Licensor, to the extent possible, 
% the Licensor offers the Licensed Material as-is and as-available, and makes no representations 
% or warranties of any kind concerning the Licensed Material, whether express, implied, statutory, or other. 
% This includes, without limitation, warranties of title, merchantability, fitness for a particular purpose, 
% non-infringement, absence of latent or other defects, accuracy, or the presence or absence of errors, 
% whether or not known or discoverable. 
% Where disclaimers of warranties are not allowed in full or in part, this disclaimer may not apply to You.
%
% To the extent possible, in no event will the Licensor be liable to You on any legal theory 
% (including, without limitation, negligence) or otherwise for any direct, special, indirect, incidental, 
% consequential, punitive, exemplary, or other losses, costs, expenses, or damages arising out of 
% this Public License or use of the Licensed Material, even if the Licensor has been advised of 
% the possibility of such losses, costs, expenses, or damages. 
% Where a limitation of liability is not allowed in full or in part, this limitation may not apply to You.
%
% The disclaimer of warranties and limitation of liability provided above shall be interpreted in a manner that, 
% to the extent possible, most closely approximates an absolute disclaimer and waiver of all liability.
%
% Developed by the IDAC Signals Team at Case Western Reserve University 
% with support from the National Institute of Neurological Disorders and Stroke (NINDS) 
%     under Grant NIH/NINDS U01-NS090405 and NIH/NINDS U01-NS090408.
%              Farhad Kaffashi
%--------------------------------------------------------------------------
function [SlopeEst]= DFA(Y)
% This function calulate the Detrended Fluctuation Analysis of RR interval
% and calculates the Aplha1 & Alpha2 slopes
% Alpha1 is for short term correlations in short time scales corresponds to log of 10 < t <  40 beats 
% Alpha2 is for  long term correlations in  long time scales corresponds to log of 70 < t < 300 beats  
%
% Ref: 
% 1. Peng C.K. "Quatifiction of scaling exponents and crossover phonomena
% in nonstationary heartbeat time series"
%
% 2. Thomas Penzel, Jan W. Kantelhardt, Ludger Grote, Jörg-Hermann
% Peter, and Armin Bunde, "Comparison of Detrended Fluctuation Analysis and
% Spectral Analysis for Heart Rate Variability in Sleep and Sleep Apnea".
%
% Inputs
% Y              : Input RR intervla time series data
%
% Outputs
% nToatl  : Calculated n values
% FnTotal : Calculated Fn Values

if size(Y,1)>size(Y,2)
    Y=Y';
end

if nargin<2
    TotalNumPoints=40;
end

Len = length(Y);

Y=cumsum(Y-mean(Y));

% 30 points for calculation of Alpha1 and 30 points of Alpha2
nTotal = [fix(logspace(log10(10),log10(40),30)) ...
    fix(logspace(log10(70),log10(Len/2),30))];

nTotal(diff(nTotal)==0)=[];
FnTotal=[];

SumN=cumsum([1:max(nTotal)]);
SumN2=cumsum( ( [1: max(nTotal) ]).^2);
ForwardSumY=cumsum(Y);
BackwardSumY=cumsum(fliplr(Y));
ForwardSumYN=cumsum(Y.*[1:Len]);
BackwardSumYN=cumsum(fliplr(Y).*[1:Len]);

nTotal(find(nTotal<2))=[];
    
for  n = nTotal;

    InvMatA=([SumN2(n) SumN(n);SumN(n) n]^(-1));

    if rem(Len,n)==0
        % there is no need for backward estimation
        Yn=zeros(1,Len);
        Para=InvMatA*[ForwardSumYN(n);ForwardSumY(n)];
        Yn(1:n)=Para(1)*[1:n]+Para(2);

        for i=2:fix(Len/n)
            Temp = ForwardSumY(n*i)-ForwardSumY((i-1)*n) ;
            Para=InvMatA*[ForwardSumYN(i*n)-ForwardSumYN((i-1)*n)-(i-1)*n*Temp; Temp];
            Yn((1:n)+(i-1)*n)=Para(1)*[1:n]+Para(2);
        end
        FnTotal=[FnTotal sqrt(mean((Y-Yn).^2))];
    else
        % forward and backward estimation
        
        % Forward estimation
        Yn=zeros(1,fix(Len/n)*n);
        Para=InvMatA*[ForwardSumYN(n);ForwardSumY(n)];
        Yn(1:n)=Para(1)*[1:n]+Para(2);

        for i=2:fix(Len/n)
            Temp = ForwardSumY(n*i)-ForwardSumY((i-1)*n) ;
            Para=InvMatA*[ForwardSumYN(i*n)-ForwardSumYN((i-1)*n)-(i-1)*n*Temp; Temp];
            Yn((1:n)+(i-1)*n)=Para(1)*[1:n]+Para(2);
        end
        
        TempFn=sqrt(mean((Y(1:fix(Len/n)*n)-Yn).^2))/2;
        
        
        % Backward estimation
        Yn=zeros(1,fix(Len/n)*n);
        Para=InvMatA*[BackwardSumYN(n);BackwardSumY(n)];
        Yn(1:n)=Para(1)*[1:n]+Para(2);

        for i=2:fix(Len/n)
            Temp = BackwardSumY(n*i)-BackwardSumY((i-1)*n) ;
            Para=InvMatA*[BackwardSumYN(i*n)-BackwardSumYN((i-1)*n)-(i-1)*n*Temp; Temp];
            Yn((1:n)+(i-1)*n)=Para(1)*[1:n]+Para(2);
        end
        
        FnTotal=[FnTotal sqrt(mean(( Y([rem(Len,n)+1:Len])-fliplr(Yn) ).^2))/2+TempFn];
        
    end

end

%%

Index = find(nTotal<50,1,'last');
FnTotal=log10(FnTotal);
nTotal=log10(nTotal);

Len = length(nTotal);

PlatStarts=[1 Index+1];
PlatEnds=[Index Len];


FitData=[];
for i = 1:length(PlatStarts)

    X = nTotal(PlatStarts(i):PlatEnds(i));
    Y = FnTotal(PlatStarts(i):PlatEnds(i));
    % Regression
    a(1,1)=length(X);
    a(1,2)=sum(X);
    a(2,1)=a(1,2);
    a(2,2)=sum(X.*X);
    b(1)=sum(Y);
    b(2)=sum(X.*Y);
    c = a\b';
    Intercept(i)=c(1);
    SlopeEst(i)=c(2);
    FitData=[FitData SlopeEst(i)*X+Intercept(i)*ones(size(X))];

end


% figure
% plot(nTotal,FnTotal,'*')
% hold on;
% grid on;
% for i=1:(length(PlatStarts))
%     plot(nTotal(PlatStarts(i):PlatEnds(i)),FitData(PlatStarts(i):PlatEnds(i)),'r')
% end
% 



