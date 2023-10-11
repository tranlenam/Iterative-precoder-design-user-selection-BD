clear 
clc
nRx = 2; % number of receive antennas (assume equal for all users)
K = 8:40; %number of users
BSpowerdB = [10]; % transmit power in dB Watt
BSpower = 10.^(BSpowerdB./10);
rho = 0; % transmit correlation coeff.
nChannels = 5000; % number of loops to calcualte the average sum rate
nTx = 8; % number of transmit antennas


Avgsumrate = 0;
%% MAIN LOOP
Avgsumrate = zeros(length(K),1);
for k =1:length(K)
    nUsers = K(k);
    corr = rho*exp(1i*(2*pi).*(rand(1,nUsers)));
    rxant = nRx*ones(1,nUsers); % rxant(k) is the number of rx. antennas of user k
    for iLoop=1:nChannels
        channel=[];
        for v=1:nUsers
            Rt=corr(v).^(0:nTx-1);
            Rt=toeplitz(Rt);
            Hm =(randn(nRx,nTx)+1i*randn(nRx,nTx))*sqrtm(Rt)/sqrt(2);
            channel=[channel;Hm];
        end

        [sumratetemp,selectedusers] = Algorithm2(channel,rxant,BSpower);
        Avgsumrate(k) = Avgsumrate(k) + sumratetemp;
    end
Avgsumrate(k) = Avgsumrate(k)/nChannels;
end
plot(K,Avgsumrate)
xlabel('Number of users');
ylabel('Average Sum Rate (b/s/Hz)')
saveas(gcf,'../results/Sumrate.png')