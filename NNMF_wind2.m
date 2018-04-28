function [Dw,Cw,cost,nuCw]= NNMF_wind2 (V,dw,...
    lambdaw, maxit,ccriteria,accl,costfct)
% [Ds,Dw,Cw,Cs,cost,nqu,nqu,nuDs]= NNMFIV,Dw,ds,lambdas,...
% lambdaw, maxit,ccriteria,accl,costfct]
%
% Input:
% V - The magnitude spectrogram of the sound to be filtered
% Dw - The windnoise dictionary
% ds - The number of speech elements to be found
% lambdas - The speech sparsity parameter
% lambdaw - The noise sparsity parameter
% maxit - The maximum number of iterations
% ccriteria - The criteria for stopping the optimization of cost function
% accl - A parameter multiplied to each optimization step
% costfct - The used costfunction, 'ls': Least squares [default] or 'kl:
% Kullbach leibler
% as
% output:
% Ds - The normalized speech dictionary
% Dw - The normalized noise dictionary
% Cw - The noise codebook
% Cs - The speech codebook
% cost - The cost in all iterations
% nqu - The acceleration parameter for the speech codebook in all
% iterations
% nqu - The acceleration parameter for the noise codebook in all
% iterations
% nuDs - The acceleration parameter for the speech dictionary in all
% iterations
% 3:
% The code is an adapted version of a matlab demofile from
% http:ffwww2.imm.dtu.dkfpubdbfviewsfpublication_details.php?id=4521
% 3:
% written by: Kristian Timm Andersen, IMH DTU 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check optinal Input
if ((nargin<3)||isempty(lambdaw))
lambdaw=0;
end
if ((nargin<4)||isempty(maxit))
%maxit=500;
maxit=10;
end
if ((nargin<5)||isempty(ccriteria))
ccriteria=1e-4;
end
if ((nargin<6)||isempty(accl))
accl=1.3;
end
if ((nargin<7)||isempty(costfct))
costfct='ls';
end
%  Initialize
[Ls,Lt]=size(V);
Dw=rand(Ls,dw);
Cw=rand(dw,Lt);
Dw = normalizeD(Dw); % Normalize
Rec=Dw*Cw;
switch costfct % Calculate costfunction
    case 'ls'
        cost=[0.5*norm(V-Rec,'fro')^2 ...
            +lambdaw*sum(Cw(:)),zeros(1,maxit)];
    case 'kl'
        cost=[sum(sum(V.*log((V+eps)./(Rec+eps))-V+Rec))...
            +lambdaw*sum(Cw(:)),zeros(1,maxit)];
end
iter=1;
nuCw=[1,zeros(1,maxit)];
nuDw=[1,zeros(1,maxit)];
% Optimization loop
while 1
        switch costfct % update factorization
            case 'ls'
                a=['enter']
                [Cw,Dw,cost(iter+1),nuCw(iter+1),...
                    nuDw(iter+1),accl] = updls(V,Cw,Dw,cost(iter),...
                    nuCw(iter),nuDw(iter),lambdaw,accl);
            case 'kl'
                [Cw,Dw,cost(iter+1),nuCw(iter+1),...
                    nuDw(iter+1),accl] = updkl(V,Cw,Dw,cost(iter),...
                    nuCw(iter),nuDw(iter),lambdaw,accl);
         end
    % Check for convergence
    if abs(cost(iter)-cost(iter+1))/cost(iter+1)<ccriteria
        % Check if acceleration parameter is too large
        if  nuDw(iter+1) <= accl && nuCw(iter+1) <= accl 
            nuDw=nuDs(1:iter+1);
            nuCw=nuCs(1:iter+1);
            cost=cost(1:iter+1);
            disp('NNMF has converged');
            break;
        else % set step size to 1 and run one more time
            nuDw(iter+1)=1;
            nuCw(iter+1)= 1;
        end
    end
    if iter<maxit % Check for number of iterations
        iter=iter+1
    else
        disp('Maximum number of iterations reached');
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cw,Dw,cost,nuCw,nuDw,accl] = updls(V,Cw,Dw,cost_old,...
                    nuCw,nuDw,lambdaw,accl)
% Update Cw
Cw_old=Cw;
Cwy=Dw'*V;
Cwx=(Dw'*Dw)*Cw+lambdaw;
grad = Cwy./(Cwx+eps);
while 1
    Cw = Cw_old.*(grad.^nuCw);
    Rec=Dw*Cw;
    cost=.5*norm(V-Rec,'fro')^2+lambdaw*sum(Cw(:));
    if cost>cost_old
        nuCw = max(nuCw/2,1);
    else
        nuCw = nuCw*accl;
    break;
    end
end
cost_old=cost;

%Update Dw
Dw_old=Dw;
Dwy=V*Cw';
Dwx=Dw*(Cw*Cw');
grad=(Dwy+Dw.*(ones(size(V,1))*(Dwx.*Dw)))...
    ./(Dwx+Dw.*(ones(size(V,1))*(Dwy.*Dw))+eps);
while 1
    Dw = normalizeD(Dw_old.*(grad.^nuDw));
    Rec=Dw*Cw;
    cost=.5*norm(V-Rec,'fro')^2+lambdaw*sum(Cw(:));
    if  cost>cost_old
        nuDw=max(nuDw/2,1);
    else
        nuDw=nuDw*accl;
        break;
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cs,Cw,Ds,cost,nuCs,nuCw,nuDs,accl] = updkl(V,Cs,Cw,Ds,Dw,...
    cost_old,nuCs,nuCw,nuDs,lambdas,lambdaw,accl)
Rec=[Ds,Dw]*[Cs;Cw];
VR=V./(Rec+eps);
O=ones(size(V));
% Update Cs
Cs_old=Cs;
Csy=Ds'*VR;
Csx=Ds'*O+lambdas;
grad = Csy./(Csx+eps);
while 1
    Cs = Cs_old.*(grad.^nuCs);
    Rec=[Ds,Dw]*[Cs;Cw];
    cost=sum(sum(V.*log((V+eps)./(Rec+eps))-V+Rec))...
        +lambdas*sum(Cs(:))+lambdaw*sum(Cw(:));
    if cost>cost_old
        nuCs = max(nuCs/2,0.1);
    else
        nuCs = nuCs*accl;
        break;
    end
end
cost_old=cost;
%Update Cw
VR=V./(Rec+eps);
Cw_old=Cw;
Cwy=Dw'*VR;
Cwx=Dw'*O+lambdaw;
grad = Cwy./(Cwx+eps);
while 1
    Cw = Cw_old.*(grad.^nuCw);
    Rec=[Ds,Dw]*[Cs;Cw];
    cost=sum(sum(V.*log((V+eps)./(Rec+eps))-V+Rec))...
    + lambdas*sum(Cs(:))+lambdaw*sum(Cw(:));
    if cost>cost_old
        nuCw=max(nuCw/2,1);
    else
        nuCw = nuCw*accl; 
        break;
    end
end
cost_old=cost;
% Update Ds
VR=V./(Rec+eps);
Ds_old=Ds;
Dsy=VR*Cs';
Dsx=O*Cs';
grad =(Dsy+Ds.*(ones(size(V,1))*(Dsx.*Ds)))...
./(Dsx+Ds.*(ones(size(V,1))*(Dsy.*Ds))+eps);
while 1
    Ds = normalizeD(Ds_old.*(grad.^nuDs));
    Rec=[Ds,Dw]*[Cs;Cw];
    cost=sum(sum(V.*log((V+eps)./(Rec+eps))-V+Rec))...
        + lambdas*sum(Cs(:))+lambdaw*sum(Cw(:));
    if cost>cost_old
        nuDs=max(nuDs/2,1);
    else
        nuDs=nuDs*accl;
        break;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Normalize D

function D = normalizeD(D)
Q = sqrt(sum(D.^2,1));
D = D./repmat(Q+eps,size(D,1),1);