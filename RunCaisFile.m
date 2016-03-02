function [r,ElapsedTime] = RunCaisFile(FileName,AnsysPath)
%FileName: File of .cai type, only File Name, No Extension
if nargin<2
    AnsysPath='"C:\Program Files\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150"';
end
% .cai to .inp
r1=Cai2Ansys(FileName);
%% run .inp
[r2 ElapsedTime]=RunAnsys(AnsysPath,FileName);
%% Return Value
r=r1+r2;