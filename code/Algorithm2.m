function [Cap,users,correlation_coeff]=Algorithm2(H,rxant,Power)
t=size(H,2);
r=length(rxant);
rowidx=cumsum([0 rxant]);

correlation_coeff=[];
userselected=[];
usersremain=(1:length(rxant));
% Initialization, choose the user that has miximum channel energy
K_max=ceil(t/mean(rxant));
for v=1:K_max
    W_prev{v}=eye(t);
    W_prev_temp{v}=zeros(t,t);
    H_eff_prev{v}=zeros(mean(rxant),t);
    H_eff_temp{v}=zeros(mean(rxant),t);
end
Prodmax=0;

for k=1:r
    C=prod(diag((H(rowidx(k)+1:rowidx(k)+rxant(k),:)*(H(rowidx(k)+1:rowidx(k)+rxant(k),:)'))));
    if (C>Prodmax)
        Prodmax=C;
        selecteduser=k;
    end
end

userselected=[userselected selecteduser];   
usersremain(selecteduser)=[];    
Hselected=H(rowidx(selecteduser)+1:rowidx(selecteduser)+rxant(k),:);
[~,V]=mygrams(Hselected);  
H_eff_prev{1}=Hselected;

n=2; 
% LOOP
while (n<=K_max);
    selecteduser=0;
    Prodmax=0;
    for k=1:length(usersremain)
        Hk=H(rowidx(usersremain(k))+1:rowidx(usersremain(k))+rxant(usersremain(k)),:);
        Hktilde=Hk-Hk*(V')*V; % project Hk to the null space of V
        C=prod(diag(Hktilde*(Hktilde')));
        rxantselected=rxant(userselected);
        rowidxnew = cumsum([0 rxantselected]);
        for m=1:length(userselected)
            
            [~,Q]=mygrams(Hk*W_prev{m});
            W_prev_temp{m}=W_prev{m}-W_prev{m}*(Q')*Q;
            
            
            H_eff_temp{m}=H_eff_prev{m}-H_eff_prev{m}*(Q')*Q;
            C=C*prod(diag(H_eff_temp{m}*(H_eff_temp{m}')));
        end
        if (C>Prodmax)
            Prodmax=C;
            selecteduser=k;
            W_prev_next=W_prev_temp;
            H_eff_next=H_eff_temp;
            Hk_temp=Hk;
            Hktilde_temp=Hktilde;
        end
    end
    if (selecteduser>0)
        userselected=[userselected usersremain(selecteduser)];
        usersremain(selecteduser)=[];
        Hselected=[Hselected; Hk_temp];
        W_prev=W_prev_next;
        H_eff_prev=H_eff_next;
        H_eff_prev{length(userselected)}=Hktilde_temp;

        W_prev{length(userselected)}=eye(t)-(V')*V;
        [~,W]=mygrams(Hktilde_temp);
        V=[V;W];
    end
    n = n+1;
end
Hnew=[];
for m=1:length(userselected)
    Hnew=[Hnew; H(rowidx(userselected(m))+1:rowidx(userselected(m))+rxant(userselected(m)),:)];
end
[temp,Cap]=GreedyBDCapacitybased(Hnew,rxant(userselected),Power);
users=userselected(temp);





