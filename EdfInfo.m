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
function [FileInfo ChInfo]=EdfInfo(FileName)

fid = fopen (FileName);

FileInfo.Version = char(fread(fid,[1 8],'uint8'));

FileInfo.LocalPatientID = char(fread(fid,[1 80],'uint8'));

FileInfo.LocalRecordID = char(fread(fid,[1 80],'uint8'));

FileInfo.StartDate = char(fread(fid,[1 8],'uint8'));

FileInfo.StartTime = char(fread(fid,[1 8],'uint8'));

FileInfo.HeaderNumBytes = str2num(char(fread(fid,[1 8],'uint8')));

FileInfo.Reserved = char(fread(fid,[1 44],'uint8'));

FileInfo.NumberDataRecord = str2num(char(fread(fid,[1 8],'uint8')));

FileInfo.DataRecordDuration = str2num(char(fread(fid,[1 8],'uint8')));

FileInfo.SignalNumbers = str2num(char(fread(fid,[1 4],'uint8')));

ns = FileInfo.SignalNumbers ;

ChInfo.Labels = char(fread(fid,[16 ns],'uint8')');

ChInfo.TransType = char(fread(fid,[80 ns],'uint8')');

ChInfo.PhyDim = char(fread(fid,[8 ns],'uint8')');

ChInfo.PhyMin = str2num(char(fread(fid,[8 ns],'uint8')'));

ChInfo.PhyMax = str2num(char(fread(fid,[8 ns],'uint8')'));

ChInfo.DiMin = str2num(char(fread(fid,[8 ns],'uint8')'));

ChInfo.DiMax = str2num(char(fread(fid,[8 ns],'uint8')'));

ChInfo.PreFiltering = char(fread(fid,[80 ns],'uint8')');

ChInfo.nr = str2num(char(fread(fid,[8 ns],'uint8')'));

ChInfo.Reserved = char(fread(fid,[32 ns],'uint8')');

fclose(fid);
clc

