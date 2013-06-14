% set antag

function [ant_network_connect]=setant(N,radius)

ant_network_connect=zeros(N,N);

for(l=1:N)
     for(i=1:radius)
         ant_network_connect(l,mod(l+i-1,N)+1)=1;
         ant_network_connect(mod(l+i-1,N)+1,l)=1;  
     end
   end
   clear left;