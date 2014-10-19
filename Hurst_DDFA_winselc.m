


function FA=Hurst_DDFA_winselc(preH,winRange,Niter,Ntp,Ncol)

sigma=1;
FA=zeros(length(winRange),Ncol);
for ww=1:length(winRange)
    tmp_winSize=winRange(ww);
    for ii=1:Niter
        f=ffgn(sigma,preH,Ncol,Ntp,0);
        [hurst,tmp_fa]=Hurst_DDFA(f,tmp_winSize);
        FA(ww,:)=FA(ww,:)+tmp_fa;
    end
    FA(ww,:)=FA(ww,:)/Niter;
end

return
end