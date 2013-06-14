% PREPROCESSfor PLASTIC: This file contains
% MATLAB code that sets the network up according 
% to the variables specified. In particular network connectivity
% is initialised ready to be passed to PLASTIC.
% Guy Billings 17/02/2007

% setup programmable protocol

if(params(23)==3)
  loadfromold=1;
end

ydproto=[];
n=size(cycle);
ydsize=repeats*n(2); %length of protocol
clear n;
for(i=1:repeats)
	ydproto=[ydproto cycle];
end	
ydproto=ydproto';

% If TEST set perform recpetive field test
% network to be tested must loaded

if(TEST==1)
  params(15)=inputs;
  params(23)=2;
  params(25)=0;
  params(26)=0;
  params(36)=0;
  output(7)=1;
  loadfromold=1;
  seedme=0;
  %save_me=0; %prevent over-write
end  

if(loadfromold==0)

 wmat=rand(N,inputs)*params(3); % Initial feedforward weights
 if(params(27)==0 && radius==0)
    rmat=zeros(N,N); % Initial recurrent weights
    amat=zeros(N,N);
 elseif(params(27)==0 && radius>0)
     rmat=ones(N,N); % Initial recurrent weights
     rmat=rmat*wdrec;
     amat=rand(N,N)*an_str;
 elseif(params(27)==1)
    rmat=ones(N,N); % Initial recurrent weights
    rmat=rmat*wdrec;
    amat=rand(N,N)*an_str;
 end
end

if(loadfromold==1)
 rmat=recmat;
 wmat=weightmat;
 amat=antmat;
end

% Seed feedforward weights

recspa=inputs/N;

if(seedme==1 && loadfromold==0)
for(l=1:N)
  for(k=1:inputs)
     if(k==1+(l-1)*recspa)
       for(w=1:width)
	%wmat(l,mod(k+w,inputs)+1)=wmat(l,mod(k+w,inputs)+1)+rand*(params(3)-wmat(l,mod(k+w,inputs)+1))*rfac;
	wmat(l,mod(k+w,inputs)+1)=params(3);
       end
      end
   end
  end
end

if(loadfromold==0)

	t_rates=0;
	ratevars=0;
	
	% Set dummy rates for protocols other than 5
	
end

% set connectivity
% input_network_connect is the feedforward connectivity matrix
% network_network_connect is the inhibitory recurrent matrix
% ant_network_connect is the excitatory recurrent matrix

if(loadfromold==0 && params(23)==0)

  input_network_connect=rand(N,inputs);
  keep=input_network_connect<=input_frac;
  discard=input_network_connect>input_frac;
  input_network_connect(keep)=1;
  input_network_connect(discard)=0;

  network_network_connect=zeros(N,N);
  ant_network_connect=zeros(N,N);

elseif(params(23)~=0 && params(23)~=5 && loadfromold==0 && in==0 && ex==1)

  input_network_connect=rand(N,inputs);
  keep=input_network_connect<=input_frac;
  discard=input_network_connect>input_frac;
  input_network_connect(keep)=1;
  input_network_connect(discard)=0;

  network_network_connect=zeros(N,N);
  
  ant_network_connect=eye(N);
  ant_network_connect=ant_network_connect==0;
  
elseif(params(23)~=0 && params(23)~=5 && loadfromold==0 && in==1 && radius>0)
  
  input_network_connect=rand(N,inputs);
  keep=input_network_connect<=input_frac;
  discard=input_network_connect>input_frac;
  input_network_connect(keep)=1;
  input_network_connect(discard)=0;
    
  network_network_connect=zeros(N,N);
  network_network_connect=eye(N);
  network_network_connect=network_network_connect==0;
  network_network_connect=-1.*network_network_connect;

  %set local connections
  [ant_network_connect]=setant(N,radius); 


elseif(params(23)~=0 && params(23)~=5 && loadfromold==0 && in==1 && ex==0)

  input_network_connect=rand(N,inputs);
  keep=input_network_connect<=input_frac;
  discard=input_network_connect>input_frac;
  input_network_connect(keep)=1;   
  input_network_connect(discard)=0;
  
  network_network_connect=zeros(N,N);
  network_network_connect=eye(N);
  network_network_connect=network_network_connect==0;
  network_network_connect=-1.*network_network_connect;
  
  ant_network_connect=zeros(N,N);

elseif(params(23)~=0 && params(23)~=5 && loadfromold==0 && in==0 && radius > 0)
  
  input_network_connect=rand(N,inputs);
  keep=input_network_connect<=input_frac;
  discard=input_network_connect>input_frac;
  input_network_connect(keep)=1;
  input_network_connect(discard)=0;
  
  network_network_connect=zeros(N,N);
  
  %set local connections
  [ant_network_connect]=setant(N,radius);
  
elseif(params(23)~=0 && params(23)~=5 && loadfromold==0 && in==0 && ex==0 && radius == 0)
  
  input_network_connect=rand(N,inputs);
  keep=input_network_connect<=input_frac;
  discard=input_network_connect>input_frac;
  input_network_connect(keep)=1;
  input_network_connect(discard)=0;
  
  network_network_connect=zeros(N,N);
  
  %set local connections
  [ant_network_connect]=zeros(N,N);
  
elseif(loadfromold==0 && params(23)==5)

  input_network_connect=rand(N,inputs);
  keep=input_network_connect<=input_frac;
  discard=input_network_connect>input_frac;
  input_network_connect(keep)=1;
  input_network_connect(discard)=0;
  
  t_rates=ones(1,inputs);
  t_rates(1:floor(inputs/2))=t_rates(1:floor(inputs/2))*bg_rate;
  t_rates(floor(inputs/2)+1:inputs)=t_rates(floor(inputs/2)+1:inputs)*fg_rate;
  
  ratevars=ones(1,inputs);
  ratevars(1:floor(inputs/2))=ratevars(1:floor(inputs/2))*bg_var;
  ratevars(floor(inputs/2)+1:inputs)=ratevars(floor(inputs/2)+1:inputs)*fg_var;

  network_network_connect=zeros(N,N);
  ant_network_connect=zeros(N,N);
  
end

% Set up weight tracking with random sampling

trackmat=zeros(output(13),2);
if(output(13)>0)
  [lii,lij]=find(input_network_connect==1);
  subs=round(rand(1,output(13))*max(size(lii))-1)+1;
  trackmat=[lii(subs);lij(subs)]';
end  

if(params(23)==3)
 params(15)=ydsize;
end 

if(save_ini==1)
    save([filen,'ini',extension]);
end    