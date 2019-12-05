function out = varindex_3(var, i, j, k, Nodes,Ktotal)

%--------------------------------------------------------------------------
% This function is used to locate the coefficients that multiply each of
% the different decision variables in the objective function or in the
% different constraints.
% 
%   INPUTS:
% 
%      - var: This variables is used to indicate which Decision Variable is
%           being located. 1=wij 2=xij 3=zij^k 4=AC^k
%      - i: This variable indicates the origin
%      - j: This variable indicates the destination
%      - k: This variable indicates the aircraft type.
%      - Nodes: This variable indicates the number of OD paths.
%      - Ktotal: This variable indicates the total number of aircraft that
%           is being considered.
% 
%   OUTPUTS:
%      - out: This variable gives the position of the Decision Variable. 
% 
%--------------------------------------------------------------------------

    if var == 1
        out = (i-1)*Nodes + j; 
    elseif var == 2
        out = Nodes*Nodes + (i-1)*Nodes + j; 
    elseif var == 3
        out = Nodes*Nodes*(1+k) + (i-1)*Nodes + j; 
    elseif var == 4
        out = Nodes*Nodes*(2+Ktotal) + k;
    end
    