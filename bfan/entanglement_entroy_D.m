function S_ell = entanglement_entroy_D(D,pos_sub)
D_sub = D(pos_sub,pos_sub);
D_sub=(D_sub+D_sub')/2;
lamd = eig(D_sub);
tmpS = (lamd.*log((lamd)) + (1-lamd).*log( (1-(lamd)) ));
tmpS(isnan(tmpS))=0;
S_ell= -real(sum(tmpS));



