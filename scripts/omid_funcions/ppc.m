% This function compute pair wise phase consistency
%
% Vinck, M., van Wingerden, M., Womelsdorf, T., Fries, P., &
% Pennartz, C. M. (2010). The pairwise phase consistency: a
% bias-free measure of rhythmic neuronal synchronization. Neuroimage, 51(1)
% , 112-122.
%
% inputs:
% phase_val1: phase value 
%
% output:
% ppc : the pair wise phase consistency value
%
% edited  By Omid Amir Atashani

function ppc = ppc(phase_val,method)

N = length(phase_val);

if N<2
    ppc=nan;
    return
end
if isempty(method)
    method = 2;
end
switch method
    
    case 1
        %%% Method one need more memory
        val_h=[];
        var_h=[]; var_h=cos(phase_val)*cos(phase_val)'+sin(phase_val)*sin(phase_val)';
        var_h_=sum(sum(triu(var_h)))-length(phase_val);
        ppc= 2*var_h_/(N*(N-1));
        
    case 2
        %%% Method 2 faster
        val_h=[];
        for ii = 1:N-1
            val_h = [ val_h nansum(cos(phase_val(ii))*cos(phase_val(ii+1:end))+sin(phase_val(ii))*sin(phase_val(ii+1:end)))];
        end
        ppc= 2*nansum(val_h)/(N*(N-1));
        
    case 3
        %%% Method 3 slow
        val_h=[];
        for ii = 1:N
            for jj = ii+1:N
                val_h = [ val_h cos(phase_val(ii))*cos(phase_val(jj))+sin(phase_val(ii))*sin(phase_val(jj))];
            end
        end
        ppc= 2*nansum(val_h)/(N*(N-1));
        
        
end

end