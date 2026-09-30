function [nu,nuR,beta] = dcqp_cut_targets(bestub,v,eta,gap_tolerance)
%DCQP_CUT_TARGETS Reproduce the stored DNN cut-target policy.

nu=bestub-abs(bestub)*gap_tolerance;
nuR=bestub-eta*abs(bestub)*gap_tolerance;
beta=abs(bestub)*gap_tolerance*0.01;

if abs(bestub)<gap_tolerance
    nu=bestub-gap_tolerance;
    nuR=bestub-0.9*gap_tolerance;
    beta=min(1e-7,0.01*gap_tolerance);
elseif v>bestub+abs(bestub)*gap_tolerance*0.1
    nuR=0.99*bestub+0.01*v;
    beta=min(1e-6,0.1*(v-bestub));
end
end
