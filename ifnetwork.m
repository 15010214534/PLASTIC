%%%%
% Support Script for C integrate and fire neuron network
% Guy Billings 1/11/04
% g.o.billings@sms.ed.ac.uk
%%%%
%----------------------------------------------------------------%
%%%% SET UP FROM CONFIGURATION SCRIPTS %%%%

setup;
parameters;
set_output;
preprocess;

%----------------------------------------------------------------%
%%%% RUN %%%%

t=cputime;

for(run_number=1:sim_repeats)

if(run_number>1)
 rmat=recmat;
 wmat=weightmat;
 amat=antmat;
end 
 
if(TRACK==0)

output(16)=0;

%% THIS INVOKES PLASTIC.MEXGLX

[time,prespikes,prelocs,voltage,postspikes,postlocs,...
stimlocs,stim_tim,weightmat,overlap,recmat,avweig,...
avweigrec,f,wvar,wvarec,euc,pear,wtracks,antmat,in_i,...
rec_i,ant_i,snapff,snaprec,snapant]=...
plastic(N,inputs,corrf,params,input_network_connect,...
network_network_connect,output,nflags,wmat,rmat,trackmat,...
amat,ant_network_connect,ydsize,ydproto,r_seed,t_rates,ratevars);

elseif(TRACK==1)

  output(7)=1;
  sinterval=params(15);
  sprot=params(23);
  sp1=params(25);
  sp2=params(26);
  sp3=params(36);
  
 if(output(16)==0)
  output(16)=1000;
 end 

 ss=output(16);
 
 columns=round((params(24)*inputs)/(output(10)*params(2)));
 pop_vectors=zeros(N,params(15)/output(16),2);
 fields=zeros(N,params(15)/output(16),round(stim_dur*inputs/stime));
 selectivity=zeros(N,params(15)/output(16));
 
[time,prespikes,prelocs,voltage,postspikes,postlocs,...
stimlocs,stim_tim,weightmat,overlap,recmat,avweig,...
avweigrec,f,wvar,wvarec,euc,pear,wtracks,antmat,in_i,...
rec_i,ant_i,snapff,snaprec,snapant]=...
plastic(N,inputs,corrf,params,input_network_connect,...
network_network_connect,output,nflags,wmat,rmat,trackmat,...
amat,ant_network_connect,ydsize,ydproto,r_seed,t_rates,ratevars);

  params(15)=inputs;
  params(23)=2;
  params(25)=0;
  params(26)=0;
  params(36)=0;
  output(16)=0;
  output(4)=0;

for(i=1:sinterval/ss)

  rmat=squeeze(snaprec(i,:,:));
  wmat=squeeze(snapff(i,:,:));
  amat=squeeze(snapant(i,:,:));

[time,prespikes,prelocs,voltage,postspikes,postlocs,...
stimlocs,stim_tim,weightmat,overlap,recmat,avweig,...
avweigrec,f,wvar,wvarec,euc,pear,wtracks,antmat,in_i,...
rec_i,ant_i,dumff,dumant,dumrec]=...
plastic(N,inputs,corrf,params,input_network_connect,...
network_network_connect,output,nflags,wmat,rmat,trackmat,...
amat,ant_network_connect,ydsize,ydproto,r_seed,t_rates,ratevars);

for(j=1:N)
  fields(j,i,:)=f(j,:);	
  selectivity(j,i)=1-((f(j,:)*[1:round(stim_dur*inputs/stime)]')/sum(max(f(j,:))*[1:round(stim_dur*inputs/stime)]));
  pop_vectors(j,i,2)=sum(f(j,:).*cos(2*pi.*[1:columns]./columns));
  pop_vectors(j,i,1)=sum(f(j,:).*sin(2*pi.*[1:columns]./columns));
end

end

  params(15)=sinterval;
  params(23)=sprot;
  params(25)=sp1;
  params(26)=sp2;
  params(36)=sp3;
  output(16)=ss;

clear sinterval;clear sprot;clear sp1;clear sp2;clear sp3;clear ss;
ang=atan2(pop_vectors(:,:,1),pop_vectors(:,:,2));

end

if(save_me==1)
save([filen,num2str(run_number),extension]);
end

end

cpu_time=cputime-t
clear t;

%----------------------------------------------------------------%




