function [users,Cap]=GreedyBDCapacitybased(H,rxant,Power)
t=size(H,2);
rowidx=cumsum([0 rxant]);

S=[];
userselected=[];
rxselected=[];
usersremain=(1:length(rxant));
Cmax=0;
n=1;
while(n<=ceil(t/mean(rxant)))
    Cmaxnew=0;
    for m=1:length(usersremain)
        Snew=[S;H(rowidx(usersremain(m))+1:rowidx(usersremain(m))+rxant(usersremain(m)),:)];
        rxselectednew=[rxselected rxant(usersremain(m))];
        C=BDCapacity(Snew,rxselectednew,Power);
        if(C>Cmaxnew)
            Cmaxnew=C;
            userselectedidx=usersremain(m);
            usersremainidx=m;
            temp=Snew;
            temp1=rxselectednew;
        end
    end
    if (Cmaxnew>Cmax)
        Cmax=Cmaxnew;
        userselected=[userselected userselectedidx];
        S=[temp];
        usersremain(usersremainidx)=[];
        rxselected=temp1;
        n=n+1;
    else
        break
    end
end
users=userselected;
Cap=Cmax;