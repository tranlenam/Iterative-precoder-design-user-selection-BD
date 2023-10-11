function Cap=BDCapacity(H,rxant,Power)
t=size(H,2);
rowidx=cumsum([0 rxant]);
d=[];
for k=1:length(rxant)
    Hk=H(rowidx(k)+1:rowidx(k+1),:);
    Hkx=H;
    Hkx(rowidx(k)+1:rowidx(k+1),:)=[];
    [r,t]=size(Hkx);
    if (r<t)
        Vk=null(Hkx);
        d=[d;svd(Hk*Vk)];
    else
        Cap=-1;
        return;
    end
end
d=abs(d).^2;
wline=wfill(1./d,Power,1e-7);
Cap=sum(max(log2(wline*d),0));