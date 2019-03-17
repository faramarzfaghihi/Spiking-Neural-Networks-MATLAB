
%%%% Matlab code for spiking neural network in the paper
%%%'The dependence of neuronal encoding efficiency on Hebbian plasticity 
%%% and homeostatic regulation of neurotransmitter release'
%%% Faghihi etal., front. in cellular neuroscience 2015
%%% Writen by Faramarz Faghihi


close all
clear all
clc

%%%%%%%%%% Probabilistic spike train generation %%%%%%%%%%

number_timebins = 200;
n = number_timebins;
firing_prob = 0.35;
p = firing_prob;
I1 = zeros(1,n);

for i1 = 1:n
    ran1 = rand (1,n);
    I1 = ran1 <= p;
end

spike = zeros(1,100*n);

for j = 1:n
    spike(j*100) = I1(1,j);
end


%-----------------------------------------------------------------------

x=0.1;y=1;xx=0.001;z=1;yy=0.001;ZZ=0.001;ZZZ=0.001;xxx=0.001;xxxx=0.0001;
sigma = 1000 ; alpha = 1000000 ; beta = 0.001 ;
ta_I = 10;ta_ca = 100;ta_RM1=100;ta_release = 10;ta_RM=1;

I_to_PostSyn = zeros(1,100); I_to_PostSyn_hist = zeros(1,n*100);

ca = zeros(1,100); Ca_hist = zeros(1,n*100);
RM1=zeros(1,100); RM_hist= zeros(1,n*100);
Relase_act=zeros(1,100);Relase_act_hist= zeros(1,n*100);

SI=zeros(1,n);Sca=zeros(1,n);

Release_binary = 1; Rel_act = zeros(1,100); RM_rel_thresh = 0.002;alpha = 0.0001;
RM_tot=zeros(1,100);
nn=200;
S=zeros(1,100);
%
n=200;
AAA = zeros(1,100);
RM_rel_thresh_hist=zeros(1,200)+0.002;

RM_rel_thresh_hist(1,1:10)=0.05;

%-------------------------------------------------------------------------

%%%%%%%%% Integrate Fire model%%%%%%%%%%

tau = 0.0154;                                   % Decay rate of membrane potential in IF model based on physiological data
Dt = 0.01;                                          % Integration step
V_reset = -40.2;                                       % Reset potential
V_resting = -84;                                       % Resting potential
V_thres = -25.8;                                           % Thereshold potential
%%%%%%%%%%%%%%%%%
VNum_SN1 = -84*ones(1,n);
VLN1 = -84*ones(1,n);

spikeNum_SN_list1 = zeros(1,n);
spikeSN1 = zeros(1,200);
%%%%%%%%%%%%%
VNum_SN2 = -84*ones(1,n);
VLN2 = -84*ones(1,n);

spikeNum_SN_list2 = zeros(1,n);
spikeSN2 = zeros(1,200);

%--------------------------------------------------------------------------

post_spike=1;

for i=2:200
    
    RM_rel_thresh;
    %%%%%%%%%%%%%%
    RM_rel_thresh_hist(1,1)=0.05;
    if i>40 & i< 120
        p2 = RM_rel_thresh;
        I2 = zeros(1,200);
        
        for i2 = 1:200
            ran2 = rand (1,200);
            I2 = ran2 <= p2;
        end
        pos = randi(length(I2));
        Release_binary= I2(pos);
        if Release_binary==0
            x=i;
        end
    else (Release_binary==1)
    end
    
    %-----------------------------------------------------------------
    
    RM1=zeros(1,100);
    I_to_PostSyn = zeros(1,100);
    I_to_PostSyn(1,1)=0.000001;ca(1,1)=xx;A=y;B=z;
    
    %%%%%%%% current in one time bin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mu1=0.1;
    for j=1:100
        
        if j==1
            
            I_to_PostSyn (1,j+1) = I_to_PostSyn(1,j) - (I_to_PostSyn (1,j)/ta_I)+ Release_binary*spike(1,A)*mu1;
            
        else
            I_to_PostSyn (1,j+1) = I_to_PostSyn(1,j) - (I_to_PostSyn (1,j)/ta_I);
        end
        
    end
    
    x=I_to_PostSyn(1,100);
    S(i)=100000000*sum(x);
    I_to_PostSyn_hist(1,A:A+100) = I_to_PostSyn;
    y=100*i;
    SI=40*sum(I_to_PostSyn);
    
    
    %%%%%%%%  Calcium one second %%%%%%%%%%%%%%%
    
    ca(1,1)=xx;
    B=z;
    mu2=0.02;
    for j=1:100
        
        if j==1
            
            ca (1,j+1) = ca(1,j) - (ca (1,j)/ta_ca)+ spike(1,B)*mu2;
            
        else
            ca(1,j+1) = ca(1,j) - (ca (1,j)/ta_ca);
        end
        
    end
    
    xx=ca(1,100);
    Ca_hist(1,B:B+100) = ca;
    z=100*i;
    Sca(1,i)=sum(Ca_hist(1,B:B+100));
    
    %%%%%%%%%%%% RM production one second %%%%%%%%%%%%%%
    
    RM=SI*SI/(1+SI*SI);
    RM_tot=RM;
    
    mu3=5;
    %%%%%%%%%%% RM dynamics one second %%%%%%%%%%%%%
    RM1(1,1)=xxx;
    ta_RM1=400;
    for j=1:100
        
        if j==1
            RM1(1,j+1) =mu3*RM_tot+ RM1(1,j) - (RM1(1,j)/ta_RM1) ;
        end
        
        if j>1
            RM1(1,j+1) = RM1(1,j) - (RM1(1,j)/ta_RM1) ;
            
        end
        
        
    end
    RM2=0.5*RM1;
    RM_hist(1,A:A+100) = abs(RM1);
    xxx=RM1(1,100);
    
    ZZZ=RM1(1,100);
    
    
    %%%%%%%%%% Release Inhibition activity one second %%%%%%%%%%
    
    ta_release = 250;
    Relase_act(1,1)=xxxx;
    mu4=0.02;
    for ii=1:100
        
        Relase_act(1,ii+1) = Relase_act(1,ii) - (Relase_act(1,ii)/ta_release)+mu4*RM2(1,ii)*ca(1,ii);
        
        xxxx=Relase_act(1,100);
    end
    
    Relase_act_hist(1,A:A+100) = abs(Relase_act);
    
    Rel_act(1,i)= Relase_act(1,100); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yy = Relase_act(1,100);
    
    %%%%%%%%%%%%% Updating release inhibition threhold %%%%%%%%%%%%%%%
    
    
    alpha = 1-exp(-0.1/yy);
    
    RM_rel_thresh = 1-exp(-alpha/yy);%%%%%%%%%%%%%%%%%%%%%%%
    if i>40 & i< 120
        RM_rel_thresh_hist(1,i+1)=1-RM_rel_thresh;
    end
    
    %-------------------------------------------------------------------------
    
    spikeSN2(1,1)=1;
    
    if (i>1 )
        
        VNum_SN1(1,i) = VNum_SN1(1,i - 1)+Dt*(1/tau*(V_resting - VNum_SN1(1,i - 1))) +S(i);
        
        if  VNum_SN1(1,i) > V_thres
            
            spikeNum_SN_list1(1,i) = 1;
            VNum_SN1(1,i) = V_reset;
            
        end%if VPN1
        
        spikeSN1 = spikeNum_SN_list1;
        
    end%
    
end
%
subplot(3,2,1)

plot(I1,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0.5 0],'MarkerSize',15,...
    'Marker','.');
title('input spike')

%

subplot(3,2,3)
plot(I_to_PostSyn_hist)
title('Input current')
% %

subplot(3,2,5)
plot(RM_rel_thresh_hist,'Color',[1 0 0])
title('RM_rel_thresh_hist')


%


subplot(3,2,2)
plot(Ca_hist)
title('Intecellular Calcium trace')



% % % % % %
subplot(3,2,4)
plot(Relase_act_hist)
title('Release inhibition activity')
%

subplot(3,2,6)
plot(RM_hist)
title('Interceluular Retrograde Messenger')
%
figure
subplot(2,1,1)

plot(I1,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0.5 0],'MarkerSize',15,...
    'Marker','.');
title('input spike')
%

subplot(2,1,2)

S = sum(spikeSN1)/100;
plot(spikeSN1,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',15,...
    'Marker','.');
title('Spiking of target neuron')

%
