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
% This function is to find the Poincare Parameters
% Inputs:
%           RR: RR Interval
%           tau: Delay
% Outputs:
%           SD1C
%           SD2C
%           SD_Ratio
%           CDOWN
%           CUP

function [SD1C SD2C SD_Ratio CDOWN CUP]=PoinCareFnct(x,y,tau)

L = length(x);
SD1C = sqrt((1/L)*sum(((x-y)-mean(x-y)).^2)/2);
SD2C = sqrt((1/L)*sum(((x + y) - mean(x + y)).^2)/2);
SD1I = sqrt((1/L) * (sum((x - y).^2)/2));
xy = (x - y)/sqrt(2);
indices_up = find(xy > 0);
indices_down = find(xy < 0);
SD1UP = sqrt(sum(xy(indices_up).^2)/L);
SD1DOWN = sqrt(sum(xy(indices_down).^2)/L);
CUP = SD1UP^2/SD1I^2;
CDOWN = SD1DOWN^2/SD1I^2;
SD_Ratio=SD1C/SD2C;