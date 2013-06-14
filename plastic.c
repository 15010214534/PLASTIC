/* PLASTIC v0.1
/* Gateway Function for c/MATLAB implimentation of an integrate
 * and fire network
 * Copyright Guy Billings 2007
 * Contact: g.billings@ucl.ac.uk
 */
/*This file is part of PLASTIC.

    PLASTIC is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    PLASTIC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PLASTIC; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
 
#include "ifnetwork.h"

gsl_rng *rgen; /* random number generator is global */

void mexFunction(int nlhs,      mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /*interface variables */
    int N;
    int inputs;
    int seed=mxGetScalar(prhs[15]);
    
    /*--- Simulation internal variable initialisation */
 
    /* GSL library random number implementation used */
    rgen=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rgen,seed);
    
    /* Set up lists for variable records */
    NEXT_LIST_UNIT *prespike_list;
    NEXT_LIST_UNIT hidden_pre;
    prespike_list=&hidden_pre;
    NEXT_LIST_UNIT *postspike_list;
    NEXT_LIST_UNIT hidden_post;
    postspike_list=&hidden_post;
    NEXT_LIST_UNIT *stim_loc_list;
    NEXT_LIST_UNIT hidden_stim_loc;
    stim_loc_list=&hidden_stim_loc;
    NEXT_LIST_UNIT *pre_loc_list;
    NEXT_LIST_UNIT hidden_pre_loc;
    pre_loc_list=&hidden_pre_loc;
    NEXT_LIST_UNIT *post_loc_list;
    NEXT_LIST_UNIT hidden_post_loc;
    post_loc_list=&hidden_post_loc;
    
    NEXT_LIST_UNIT aux;
    NEXT_LIST_UNIT aux1;
    NEXT_LIST_UNIT aux2;
  
    NEXT_LIST_UNIT *aux3;
    NEXT_LIST_UNIT hidden_aux3;
    aux3=&hidden_aux3;
    
    NEXT_LIST_UNIT *aux4;
    NEXT_LIST_UNIT hidden_aux4;
    aux4=&hidden_aux4;
    
    NEXT_LIST_UNIT *aux5;
    NEXT_LIST_UNIT hidden_aux5;
    aux5=&hidden_aux5;
    
    NEXT_LIST_UNIT *aux6;
    NEXT_LIST_UNIT hidden_aux6;
    aux6=&hidden_aux6;
    
    NEXT_LIST_UNIT *aux7;
    NEXT_LIST_UNIT hidden_aux7;
    aux7=&hidden_aux7;

    NEXT_LIST_UNIT *aux8;
    NEXT_LIST_UNIT hidden_aux8;
    aux8=&hidden_aux8;
   
    NEXT_LIST_UNIT *aux9;
    NEXT_LIST_UNIT hidden_aux9;
    aux9=&hidden_aux9;
    
    NEXT_LIST_UNIT *aux10;
    NEXT_LIST_UNIT hidden_aux10;
    aux10=&hidden_aux10;
    
    NEXT_LIST_UNIT *input_i;
    NEXT_LIST_UNIT hidden_input;
    input_i=&hidden_input;
    
    NEXT_LIST_UNIT *rec_i;
    NEXT_LIST_UNIT hidden_rec;
    rec_i=&hidden_rec;
    
    NEXT_LIST_UNIT *ant_i;
    NEXT_LIST_UNIT hidden_ant;
    ant_i=&hidden_ant;
      
    /*Sentinel data*/
    int timesteps=0;
    int steps=0;
    double *timerec;
    double endtime=0.0;
    double time=0.0;
    int i;
    int stimlist;
    double *corr_flag_mat;
    double *t_rates;
    double *vars;
    double *params;
    double *output;
    double *nflags;
    double *dV;
    double *dg;
    double *dg_i;
    gsl_matrix *dD;
    gsl_matrix *dP;
    gsl_matrix *dDrec;
    gsl_matrix *dPrec;
    gsl_matrix *dDant;
    gsl_matrix *dPant;
    double *V;
    double *total_i;
    double *total_r;
    double *total_a;
    double *auxarr;
    double *auxout;
    double *aux1arr;
    double *aux2arr;
    double *aux2out; 
    double *aux1out;
    double *aux3arr;
    double *aux3out;
    double *aux4out;
    double *aux5out;
    double *aux6out;
    double *aux8out;
    double *aux9out;
    double *aux7out;
    double *input_iout;
    double *rec_iout;
    double *ant_iout;
    double *aux10out;
    double *aux4arr;
    double *aux5arr;
    double *aux6arr;
    double *aux7arr;
    double *aux8arr;
    double *aux9arr;
    double *aux10arr;
    double *input_iarr;
    double *rec_iarr;
    double *ant_iarr;
    double *uncorr_fir_rate;
    double *last_prespike;
    double *last_postspike;
    double *post_spike_array;
    double *post_loc_array;
    double *pre_loc_array;
    double *pre_spike_array;
    double *stim_loc_array;
    double *pre_out;
    double *pre_locs_out;
    double *post_out;
    double *post_locs_out;
    double *stim_loc_out;
    double *in_net_out;
    double *net_net_out;
    double *ant_net_out;
    gsl_matrix *input_network;
    gsl_matrix *input_network_store;
    gsl_matrix *network_network;
    gsl_matrix *antag;
    gsl_matrix *syn_conduct;
    gsl_matrix *syn_conduct_i;
    gsl_matrix *syn_conduct_rec;
    gsl_matrix *syn_conduct_rec_i;
    gsl_matrix *syn_conduct_ant;
    gsl_matrix *syn_conduct_ant_i;
    gsl_matrix *in_net_connect;
    gsl_matrix *net_net_connect;
    gsl_matrix *ant_net_connect;
    gsl_matrix *scaled_ff;
    gsl_matrix *scaled_rec;
    gsl_matrix *scaled_ant;
    gsl_matrix *last_p;
    gsl_matrix *last_pr;
    gsl_matrix *dN_ff;
    gsl_matrix *dN_rec;
    gsl_matrix *dN_ant;
    gsl_matrix *last_event_aff;
    gsl_matrix *last_event_bff;
    gsl_matrix *last_event_arec;
    gsl_matrix *last_event_brec;
    gsl_matrix *track_matrix;
    gsl_matrix *yd_protocol;
    double *input_activity_out;
    double *input_activity;
    double *input_connected;
    int size[2];
    int *firings;
    firings=malloc(2*sizeof(int));
    firings[0]=0;
    firings[1]=0;
    double *postaccum;
    int stim_loc=0;
    int terminate=0;
    int updates=0;
    int posp=0;
    int *postspikes=&posp;
    int asize;
    int asize1;
    int asize2;
    int asize3;
    int asize4;
    int asize5;
    int asize6;
    int asize7;
    int asize8;
    int asize9;
    int asize10;
    int asize_input;
    int asize_rec;
    int asize_ant;
    int psp=0;
    int *prespikes=&psp;
    int interval=0;
    int auxent=1;
    int ss_dims[3];
    int ss_dims_rec[3];
    int ss_dims_ant[3];
    int recurrents=0;
    int antagonists=0;
   
    int perturb_ff=0;
    int perturb_rec=0;
    
    int yd_size;
    int yd_steps=0;
    
    double ff_mult=1.0;
    
/* Input arguments: 
 * 0) N number of network neurons
 * 1) Inputs number of Poisson train input
 * 2) Correlation flag array
 * 3) Simulation parameter array
 * 4) input/network connections matrix
 * 5) network/network connections matrix
 * 6) Output flags array
 * 7) Neuron flags
 * 8) Initial input weights
 * 9) Initial recurrent wieghts
 */
 
    /* Read in input data */
    mexPrintf("PLASTIC a simulator for plastic spiking networks\n");
    mexPrintf("2007 Guy Billings, The University of Edinburgh\n");
    mexPrintf("Please report crashes/bugs to g.o.billings@sms.ed.ac.uk\n");
    mexPrintf("PLASTIC: initialising from input arguments...\n");
    
    N=(int)(mxGetScalar(prhs[0]));
    inputs=(int)(mxGetScalar(prhs[1]));
    corr_flag_mat=mxGetPr(prhs[2]);
    t_rates=mxGetPr(prhs[16]);
    vars=mxGetPr(prhs[17]);
    yd_size=(int)(mxGetScalar(prhs[13]));
    params=mxGetPr(prhs[3]); /* simulation parameters*/
    in_net_connect=build_gsl_matrix(N,inputs,0,0);
    net_net_connect=build_gsl_matrix(N,N,0,0);
    ant_net_connect=build_gsl_matrix(N,N,0,0);
    
    if(yd_size>0)
     yd_protocol=build_gsl_matrix(yd_size,2,0,0);
     
    import_gsl_matrix(in_net_connect,mxGetPr(prhs[4]),N,inputs);
    import_gsl_matrix(net_net_connect,mxGetPr(prhs[5]),N,N);
    import_gsl_matrix(ant_net_connect,mxGetPr(prhs[12]),N,N);
    
    if(yd_size>0)
     import_gsl_matrix(yd_protocol,mxGetPr(prhs[14]),yd_size,2);
     
    output=mxGetPr(prhs[6]);
    nflags=mxGetPr(prhs[7]);
    
    mexPrintf("PLASTIC: done\n");
    
   /* Output options format:
   /* output[0]: prespike reconstruction flag (>0.5)
   /* output[1]: cell voltage output (n, for cell number);
   /* output[2]: postspike reconstruction data flag (>0.5);
   /* output[3]: Stimulus locations & durations flag (>0.5;protocol==1);
   /* output[4]: Overlap flag (>0.5);
   /* output[5]: Average weights flag (>0.5)
   /* output[6]: Frequency flag
   /* output[7]: Average recurrent wieghts
   /* output[8]: Output variance flag
   /* output[9]: Sampling period (timesteps)
   /* output[10]: Euclidean distance between weight vectors
   /* output[11]: Pearson's R of wieght vectors
   /* output[12]: weight tracking flag
   /* output[13]: complete/reduced weight tracking
   /* output[14]: current output flag
   /* output[15]: steps interval for weight snapshots
   /* neuron_flag: array of flags denoting neurons to 
   /* be included in the average
   /* track_matrix: specifies wieghts to track
   
   /* setup simulation parameters */
   
   mexPrintf("PLASTIC: setting up simulation parameters...\n");
   
    double mean_interval=params[0];
    double dt=params[1];
    double gmax=params[2];
    double Tex=params[3];
    double Vrest=params[4];
    double Vinit=params[5];
    double Eex=params[6];
    double Ein=params[7];
    double spike_thresh=params[8];
    double Tmem=params[9];
    double Tm=params[10];
    double Tp=params[11];
    double l_rate_pot=params[12];
    double l_rate_dep=params[13];
    int nsteps=params[14];
    double corr=params[15];
    double r0=params[16];
    double r1=params[17];
    double sigma=params[18];
    int in_loc=params[19];     
    double freq=params[20];
    double g_input=params[21];
    int protocol=params[22];
    double int_time=params[23];
    int ffplas=params[24];
    int recplas=params[25];
    int wdep=params[26];
    double wdep_sigma=params[27];
    double c_pot=params[28];
    double c_dep=params[29];
    int syperturb_ff=params[30];
    int syperturb_rec=params[31];
    double single_g=params[32];
    double taun=params[33];
    double sdn=params[34];
    int antplas=params[35];
    double mult=params[36];
    double rate=params[37];
    double var_rate=params[38];
    
    double max_receptors=gmax/single_g;
    double noise_dec=exp(-dt/taun);
    
    mexPrintf("PLASTIC: done\n");
    
    /* Initialise network: 
    * Weight matricies */
    
    mexPrintf("PLASTIC: initialising weight matricies...\n");
    
    input_network=build_gsl_matrix(N,inputs,0,0);
    network_network=build_gsl_matrix(N,N,0,0);
    antag=build_gsl_matrix(N,N,0,0);
    import_gsl_matrix(input_network,mxGetPr(prhs[8]),N,inputs);
    import_gsl_matrix(network_network,mxGetPr(prhs[9]),N,N);
    import_gsl_matrix(antag,mxGetPr(prhs[11]),N,N);
    input_network_store=build_gsl_matrix(N,inputs,0,0);
    gsl_matrix_memcpy(input_network_store,input_network);
    dD=build_gsl_matrix(N,inputs,0,0);
    dP=build_gsl_matrix(N,inputs,0,0);
    dDrec=build_gsl_matrix(N,N,0,0);
    dPrec=build_gsl_matrix(N,N,0,0);
    dDant=build_gsl_matrix(N,N,0,0);
    dPant=build_gsl_matrix(N,N,0,0);
    dg=build_double_array(N);  
    dg_i=build_double_array(N);
    input_connected=build_double_array(inputs);
    find_connections(N,inputs,in_net_connect,input_connected);
    scaled_ff=build_gsl_matrix(N,inputs,0,0);
    scaled_rec=build_gsl_matrix(N,N,0,0);
    scaled_ant=build_gsl_matrix(N,N,0,0);
    dN_ff=build_gsl_matrix(N,inputs,0,0);
    dN_rec=build_gsl_matrix(N,N,0,0);
    dN_ant=build_gsl_matrix(N,N,0,0);
    recurrents=find_recurrents(N,net_net_connect);
    antagonists=find_recurrents(N,ant_net_connect);
   
    if(output[12]>0.5){
      track_matrix=build_gsl_matrix((int)output[12],2,0,0);
      import_gsl_matrix(track_matrix,mxGetPr(prhs[10]),(int)output[12],2);
      }
     
    mexPrintf("PLASTIC: done\n");
   
    /* Activity Matricies */
    
    mexPrintf("PLASTIC: initialising activity matricies...\n");
    
    input_activity=build_double_array(inputs);
    uncorr_fir_rate=build_double_array(inputs);
    input_activity_random(inputs,input_activity);
    last_prespike=build_double_array(inputs);
    last_postspike=build_double_array(N);
    V=build_double_array(N);
    total_i=build_double_array(N);
    total_r=build_double_array(N);
    total_a=build_double_array(N);
    scal_add_arr(Vinit,N,V);
    dV=build_double_array(N);
    syn_conduct=build_gsl_matrix(N,inputs,0,0);
    syn_conduct_rec=build_gsl_matrix(N,N,0,0);
    syn_conduct_i=build_gsl_matrix(N,inputs,0,0);
    syn_conduct_rec_i=build_gsl_matrix(N,N,0,0);
    syn_conduct_ant=build_gsl_matrix(N,N,0,0);
    syn_conduct_ant_i=build_gsl_matrix(N,N,0,0);
    postaccum=build_double_array(N);
    last_p=build_gsl_matrix(N,inputs,0,0);
    last_pr=build_gsl_matrix(N,N,0,0);
    
    last_event_aff=build_gsl_matrix(N,inputs,0,0);
    last_event_bff=build_gsl_matrix(N,inputs,0,0);
    
    last_event_arec=build_gsl_matrix(N,N,0,0);
    last_event_brec=build_gsl_matrix(N,N,0,0);
    
    if(protocol==2){
      endtime=int_time;
      stim_loc=0;
    }

    if(protocol==3 && yd_size<=0)
      error_out(4);
      
    if(output[15]>0){
      ss_dims[0]=(int)nsteps/output[15];
      ss_dims[1]=N;
      ss_dims[2]=inputs; 
      plhs[23]=mxCreateNumericArray(3,ss_dims,mxDOUBLE_CLASS,mxREAL);
      
      /* These can be activated to allow the recording of the 
      auxiliary plasticity variables: Must modify output arguement
      list in ifnetwork.m otherwise memory corruption occurs*/
      
      /*plhs[26]=mxCreateNumericArray(3,ss_dims,mxDOUBLE_CLASS,mxREAL);
      plhs[27]=mxCreateNumericArray(3,ss_dims,mxDOUBLE_CLASS,mxREAL);*/
      
      ss_dims_rec[0]=ss_dims[0];
      ss_dims_rec[1]=N;
      ss_dims_rec[2]=N; 
      plhs[24]=mxCreateNumericArray(3,ss_dims_rec,mxDOUBLE_CLASS,mxREAL);
      
      ss_dims_rec[0]=ss_dims[0];
      ss_dims_rec[1]=N;
      ss_dims_rec[2]=N; 
      plhs[25]=mxCreateNumericArray(3,ss_dims_rec,mxDOUBLE_CLASS,mxREAL);
      }
      
    else{
      plhs[23]=mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);  
      plhs[24]=mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
      plhs[25]=mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
      
      /*plhs[26]=mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);
      plhs[27]=mxCreateNumericArray(0,0,mxDOUBLE_CLASS,mxREAL);*/
      }
    
    mexPrintf("PLASTIC: done\n");
    
    mexPrintf("PLASTIC: performing iterations...\n");
    
    /*Perform iterations over master timescale*/
   while(terminate!=1)
   { 
    
    /* Setup activity */

    if(protocol==0){
    input_activity_random(inputs,input_activity);
    mixed_rate(inputs,input_activity,uncorr_fir_rate, corr, corr_flag_mat, rate, var_rate);
    endtime=-mean_interval*exp_dist(1.0);
    }
 
    else if(protocol==1){
    stim_loc=custom_rand((double)inputs,rgen); 
    localised_stimulus(inputs,uncorr_fir_rate,r0,r1,sigma,stim_loc,(double)inputs);
    endtime=-mean_interval*exp_dist(1.0);
    }
    
    else if(protocol==2){
    stim_loc=test_stimloc(inputs,stim_loc);
    localised_stimulus(inputs,uncorr_fir_rate,r0,r1,sigma,stim_loc,(double)inputs);
    }
    
    else if(protocol==3){
    stim_loc=gsl_matrix_get(yd_protocol,yd_steps,0); /*Protocol location stored in 1st field*/
    endtime=gsl_matrix_get(yd_protocol,yd_steps,1); /*step duration in second field*/
    
      if(stim_loc >= 0)
       localised_stimulus(inputs,uncorr_fir_rate,r0,r1,sigma,stim_loc,(double)inputs);
      else if(stim_loc<0)
       scal_mult_arr(0,inputs,uncorr_fir_rate); 
     
    yd_steps++;
    }

    else if(protocol==4){
     endtime=int_time;
     stim_loc=in_loc;
     localised_stimulus(inputs,uncorr_fir_rate,r0,r1,sigma,stim_loc,(double)inputs);
    }
    
    else if(protocol==5){
    input_activity_random(inputs,input_activity);
    inhm_rate(inputs,input_activity,uncorr_fir_rate,t_rates,vars);
    endtime=-mean_interval*exp_dist(1.0);
    ff_mult=mult;
    }
     
    
        
    /*loop over correlation timescale */

    for(time=0.0 ; time<endtime ; time+=dt)
    {
     
      /* Perturb if switch set...*/
        
     if(syperturb_ff >= 0 && timesteps>=syperturb_ff+1)
      perturb_ff=1;
      
     if(syperturb_rec >= 0 && timesteps>=syperturb_rec+1)
      perturb_rec=1;
     
      /*Pre synaptic firing */
      pre_fire(inputs,dt,uncorr_fir_rate,last_prespike,input_connected);
      
      if(output[0]>0.5){
        any_of(inputs,last_prespike,firings); 
        record_spikes(firings,protocol,stim_loc,last_prespike,stim_loc_list,prespike_list,pre_loc_list,prespikes,timesteps*dt,inputs);
        firings[0]=0;
        firings[1]=0;
      }
      
      if(ffplas==1){
      /* depress wieghts */
       if(wdep==0)
        update_weights(N,inputs,0,input_network,dP,dD,in_net_connect,gmax,last_prespike,l_rate_pot,perturb_ff,0,dN_ff,noise_dec,scaled_ff,dt,
	taun,sdn,timesteps*dt,last_event_aff,last_event_bff,Tp,Tm); 
       else if(wdep==1)
        update_weights_dep(N,inputs,0,input_network,dP,dD,in_net_connect,gmax,last_prespike,c_pot,wdep_sigma,perturb_ff,0,dN_ff,scaled_ff, 
        noise_dec,dt,taun,sdn,timesteps*dt,Tp,Tm,last_event_aff,last_event_bff); 	
      }

      if(recplas==1 && recurrents==1){
      /* depress recurrent wieghts */
       if(wdep==0)
        update_weights(N,N,0,network_network,dPrec,dDrec,net_net_connect,gmax,last_postspike,l_rate_pot,perturb_rec,0,dN_rec,noise_dec,scaled_rec,dt,
	taun,sdn,timesteps*dt,last_event_arec,last_event_brec,Tp,Tm);
       else if(wdep==1)	
        update_weights_dep(N,N,0,network_network,dPrec,dDrec,net_net_connect,gmax,last_postspike,c_pot,wdep_sigma,perturb_rec,0,dN_rec,scaled_rec, 
        noise_dec,dt,taun,sdn,timesteps*dt,Tp,Tm,last_event_arec,last_event_brec);
      }
      
      if(antplas==1 && antagonists==1){
      /* depress antagonist wieghts */
       if(wdep==0)
        update_weights(N,N,0,antag,dPant,dDant,ant_net_connect,gmax,last_postspike,l_rate_pot,perturb_rec,0,dN_ant,noise_dec,scaled_ant,dt,
	taun,sdn,timesteps*dt,last_event_arec,last_event_brec,Tp,Tm);
       else if(wdep==1)	
        update_weights_dep(N,N,0,antag,dPant,dDant,ant_net_connect,gmax,last_postspike,c_pot,wdep_sigma,perturb_rec,0,dN_ant,scaled_ant, 
        noise_dec,dt,taun,sdn,timesteps*dt,Tp,Tm,last_event_arec,last_event_brec);
      }
      
      
      /* Calculate input conductances */
      syn_input(N,inputs,input_network,syn_conduct,syn_conduct_i,in_net_connect,dg,dg_i,last_prespike,dt,Tex,ff_mult);
      if(output[14]>0.5){
        total_current(N,inputs,syn_conduct,in_net_connect,total_i);
        
        if(auxent==output[9])
          out_current(total_i,Eex,V,input_i,auxent,N);
      } 
      
      /* Calculate recurrent conductances */ 
      if(recurrents==1){
        syn_input(N,N,network_network,syn_conduct_rec,syn_conduct_rec_i,net_net_connect,dg,dg_i,last_postspike,dt,Tex,1);
	
	if(output[14] > 0.5){
          total_current(N,N,syn_conduct_rec_i,net_net_connect,total_r);

          if(auxent==output[9])
            out_current(total_r,Ein,V,rec_i,auxent,N);

        }

       }
	
      /* Antagonist conductances */
      if(antagonists==1){
        syn_input(N,N,antag,syn_conduct_ant,syn_conduct_ant_i,ant_net_connect,dg,dg_i,last_postspike,dt,Tex,mult);
	if(output[14]>0.5 && auxent==output[9]){
          total_current(N,N,syn_conduct_ant,ant_net_connect,total_a);

          if(auxent==output[9])
            out_current(total_a,Eex,V,ant_i,auxent,N);


         }
	}
      
      if(protocol==0 || protocol==5 && freq>0)
        drive_noise(N,freq,dt,g_input,dg);
         
      /* Update membrane voltage; zero dg */
      
      update_V(N,dt,Vrest,V,dV,dg,dg_i,Eex,Ein,Tmem);
     
      /* Post Synaptic firing */
      post_fire(N,spike_thresh,Vinit,last_postspike,V);
      
      if(output[1]>0.5){
       aux=double_push_back(V[(int)output[1]-1],aux);
      }
    
      if(output[2]>0.5){    
        any_of(N,last_postspike,firings); 
        record_spikes(firings,protocol,stim_loc,last_postspike,stim_loc_list,postspike_list,post_loc_list,postspikes,timesteps*dt,N);
        firings[0]=0;
        firings[1]=0;
      }	
   
      if(ffplas==1){
      /* potentiate weights */
       if(wdep==0) 
        update_weights(N,inputs,1,input_network,dD,dP,in_net_connect,gmax,last_postspike,l_rate_dep,0,0,dN_ff,noise_dec,scaled_ff,dt,
	taun,sdn,timesteps*dt,last_event_bff,last_event_aff,Tm,Tp);
       else if(wdep==1)
        update_weights_dep(N,inputs,1,input_network,dD,dP,in_net_connect,gmax,last_postspike,c_dep,wdep_sigma,0,0,dN_ff,scaled_ff, 
        noise_dec,dt,taun,sdn,timesteps*dt,Tm,Tp,last_event_bff,last_event_aff);	
      }

      if(recplas==1 && recurrents==1){
      /*Recurrents...*/
       if(wdep==0)
        update_weights(N,N,1,network_network,dDrec,dPrec,net_net_connect,gmax,last_postspike,l_rate_dep,0,0,dN_rec,noise_dec,scaled_rec,dt,
	    taun,sdn,timesteps*dt,last_event_brec,last_event_arec,Tm,Tp);
       else if(wdep==1)
        update_weights_dep(N,N,1,network_network,dDrec,dPrec,net_net_connect,gmax,last_postspike,c_dep,wdep_sigma,0,0,dN_rec,scaled_rec, 
        noise_dec,dt,taun,sdn,timesteps*dt,Tm,Tp,last_event_brec,last_event_arec);	
      }
      
       if(antplas==1 && antagonists==1){
      /*Antagonists...*/
       if(wdep==0)
        update_weights(N,N,1,antag,dDant,dPant,ant_net_connect,gmax,last_postspike,l_rate_dep,0,0,dN_ant,noise_dec,scaled_ant,dt,
	taun,sdn,timesteps*dt,last_event_brec,last_event_arec,Tm,Tp);
       else if(wdep==1)
        update_weights_dep(N,N,1,antag,dDant,dPant,ant_net_connect,gmax,last_postspike,c_dep,wdep_sigma,0,0,dN_ant,scaled_ant, 
        noise_dec,dt,taun,sdn,timesteps*dt,Tm,Tp,last_event_brec,last_event_arec);	
      }
      
      if(output[4]>0.5 && auxent==output[9]){
        aux2=double_push_back(overlap_vec(N,inputs,in_net_connect,input_network,input_network_store),aux2);
      }	
      
      if(output[5]>0.5 && auxent==output[9]){
        average_weights(N,inputs,nflags,input_network,in_net_connect,aux3,output[8],aux6);
      }
    	
      
      if(output[6]>0.5){
        accum_double_array(last_postspike,postaccum,N);
	if(auxent==output[9]){
	 find_frequency(N,postaccum,nflags,aux5,auxent*dt);
	 scal_mult_arr(0,N,postaccum);
	 }
      }
      
      if(output[7]>0.5 && auxent==output[9]){
        average_weights(N,N,nflags,network_network,net_net_connect,aux4,output[8],aux7);
      }	

      if(output[10]>0.5 && auxent==output[9]){
        euclid(N,inputs,input_network,input_network_store,in_net_connect,aux8);
      }

      if(output[11]>0.5 && auxent==output[9]){
        pearson(N,inputs,input_network,input_network_store,in_net_connect,aux9);
      }
      
      if(output[13] > 0.5 && output[12] > 0.5 && auxent==output[9]){
        track_weight((int)output[12],track_matrix,input_network,aux10);
	}
      else if(output[13] < 0.5 && output[12] > 0.5){
        track_weight((int)output[12],track_matrix,input_network,aux10);
	}
      
      if(auxent==output[9] && output[9]!=1){
        auxent=1;
	}
      else if(auxent!=output[9] && output[9]!=1){
        auxent++;
	}
	
      timesteps++;
    }
    
        
    if(protocol>0 && output[3]>0.5){
        aux1=double_push_back(timesteps*dt,aux1);
        *stim_loc_list=double_push_back(stim_loc,*stim_loc_list);
    }
     
    if(output[15]>0 && fmod(steps,output[15])==0){
      snapshot(input_network,mxGetPr(plhs[23]),N,inputs,(int)(steps/output[15]),ss_dims[0]);
      snapshot(network_network,mxGetPr(plhs[24]),N,N,(int)(steps/output[15]),ss_dims[0]);
      snapshot(antag,mxGetPr(plhs[25]),N,N,(int)(steps/output[15]),ss_dims[0]);
      
      /* Uncomment for plasticity variables recording: must modify ifnetwork.m 
      snapshot(dP,mxGetPr(plhs[26]),N,inputs,(int)(steps/output[15]),ss_dims[0]);
      snapshot(dD,mxGetPr(plhs[27]),N,inputs,(int)(steps/output[15]),ss_dims[0]); */
      }
    
     steps++;     
     terminate=the_terminator(nsteps,steps,protocol,yd_steps,yd_size);
        
   }   
   
   mexPrintf("PLASTIC: done\n");
   
   mexPrintf("PLASTIC: converting lists to arrays...\n");
   /* Convert lists to arrays for output */
   if(output[0]>0.5){ 
    pre_loc_array=build_double_array(*prespikes);
    list_to_array(*pre_loc_list,*prespikes,pre_loc_array);
    pre_spike_array=build_double_array(*prespikes);
    list_to_array(*prespike_list,*prespikes,pre_spike_array);
   }
   else{
    *prespikes=0;
   } 
   if(output[2]>0.5){ 
    post_loc_array=build_double_array(*postspikes);
    list_to_array(*post_loc_list,*postspikes,post_loc_array);
    post_spike_array=build_double_array(*postspikes);
    list_to_array(*postspike_list,*postspikes,post_spike_array);
   }
   else{localised_stimulus(inputs,uncorr_fir_rate,r0,r1,sigma,stim_loc,(double)inputs);
    *postspikes=0;
   }  
    
   if(protocol>0 && output[3]>0.5){
     asize1=steps;
     stim_loc_array=build_double_array(steps);
     list_to_array(*stim_loc_list,steps,stim_loc_array);
     stimlist=steps;
   }  
   else{
    stimlist=0;
    asize1=0;
   } 
    
   if(output[1]>0.5){
     asize=timesteps;
   }
   else{
     asize=0;
   }

   if(output[4]>0.5){
     asize2=(int)timesteps/output[9];
   }
   else{
     asize2=0;
   }     
   
   if(output[5]>0.5){
     asize3=(int)(timesteps/output[9])*sum_arr(N,nflags);
     if(output[8]>0.5){
       asize6=asize3;
       }
     else{
       asize6=0;  
       }
   }
   else{
     asize3=0;
     asize6=0;
   }
   
   if(output[6]>0.5){
     asize5=(int)(timesteps/output[9])*sum_arr(N,nflags);
   }
    else{
     asize5=0;
   }       
   
   if(output[7]>0.5){
     asize4=(int)(timesteps/output[9])*sum_arr(N,nflags);
     if(output[8]>0.5){
       asize7=asize4;
       }
     else{
       asize7=0;
       }      
   }  
   else{
     asize4=0;
     asize7=0;
   }  

   if(output[10]>0.5){
     asize8=(int)timesteps/output[9];
   }
   else{
     asize8=0;
   }

   if(output[11]>0.5){
     asize9=(int)timesteps/output[9];
   }
   else{
     asize9=0;
   }
   
   if(output[12]>0.5 && output[13]>0.5){
     asize10=((int)timesteps/output[9])*output[12];
   }
   else if(output[12]>0.5 && output[13]<0.5){
     asize10=timesteps*output[12];
   }
   else{
     asize10=0;
   }      
   
   if(output[14]>0.5)
     asize_input=(int)(timesteps/output[9])*N;
   else
     asize_input=0;
     
   if(output[14]>0.5 && recurrents==1)
     asize_rec=(int)(timesteps/output[9])*N;
   else
     asize_rec=0;
     
   if(output[14]>0.5 && antagonists==1)
     asize_ant=(int)(timesteps/output[9])*N;
   else
     asize_ant=0;      
       

    /*Convert remaining lists */
    auxarr=build_double_array(asize);
    list_to_array(aux,asize,auxarr);
    aux1arr=build_double_array(asize1);
    list_to_array(aux1,asize1,aux1arr);
    aux2arr=build_double_array(asize2);
    list_to_array(aux2,asize2,aux2arr);
    aux3arr=build_double_array(asize3);
    list_to_array(*aux3,asize3,aux3arr);
    aux4arr=build_double_array(asize4);
    list_to_array(*aux4,asize4,aux4arr);
    aux5arr=build_double_array(asize5);
    list_to_array(*aux5,asize5,aux5arr);
    aux6arr=build_double_array(asize6); 
    list_to_array(*aux6,asize6,aux6arr); 
    aux7arr=build_double_array(asize7); 
    list_to_array(*aux7,asize7,aux7arr); 
    aux8arr=build_double_array(asize8);
    list_to_array(*aux8,asize8,aux8arr);
    aux9arr=build_double_array(asize9);
    list_to_array(*aux9,asize9,aux9arr);
    aux10arr=build_double_array(asize10);
    list_to_array(*aux10,asize10,aux10arr);
    
    input_iarr=build_double_array(asize_input);
    list_to_array(*input_i,asize_input,input_iarr);
    rec_iarr=build_double_array(asize_rec);
    list_to_array(*rec_i,asize_rec,rec_iarr);
    ant_iarr=build_double_array(asize_ant);
    list_to_array(*ant_i,asize_ant,ant_iarr);
    
    mexPrintf("PLASTIC: done\n");
    
    mexPrintf("PLASTIC: initialising output arrays...\n");
    
   /* Initialise output arrays */
   plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
   timerec=mxGetPr(plhs[0]);
   
   plhs[1]=mxCreateDoubleMatrix(1,*prespikes,mxREAL);
   pre_out=mxGetPr(plhs[1]);
   
   plhs[2]=mxCreateDoubleMatrix(1,*prespikes,mxREAL);
   pre_locs_out=mxGetPr(plhs[2]);
   
   plhs[3]=mxCreateDoubleMatrix(1,asize,mxREAL);
   auxout=mxGetPr(plhs[3]);
   
   plhs[4]=mxCreateDoubleMatrix(1,*postspikes,mxREAL);
   post_out=mxGetPr(plhs[4]);
   plhs[5]=mxCreateDoubleMatrix(1,*postspikes,mxREAL);
   post_locs_out=mxGetPr(plhs[5]);
   
   plhs[6]=mxCreateDoubleMatrix(1,stimlist,mxREAL);
   stim_loc_out=mxGetPr(plhs[6]);
    
   plhs[7]=mxCreateDoubleMatrix(1,asize1,mxREAL);
   aux1out=mxGetPr(plhs[7]);
   
   plhs[8]=mxCreateDoubleMatrix(N,inputs,mxREAL);
   in_net_out=mxGetPr(plhs[8]);
  
   plhs[9]=mxCreateDoubleMatrix(1,asize2,mxREAL);
   aux2out=mxGetPr(plhs[9]);
   
   plhs[10]=mxCreateDoubleMatrix(N,N,mxREAL);
   net_net_out=mxGetPr(plhs[10]);
   
   plhs[11]=mxCreateDoubleMatrix(sum_arr(N,nflags),asize3/sum_arr(N,nflags),mxREAL);
   aux3out=mxGetPr(plhs[11]);
   
   plhs[12]=mxCreateDoubleMatrix(sum_arr(N,nflags),asize4/sum_arr(N,nflags),mxREAL);
   aux4out=mxGetPr(plhs[12]);
   
   plhs[13]=mxCreateDoubleMatrix(sum_arr(N,nflags),asize5/sum_arr(N,nflags),mxREAL);
   aux5out=mxGetPr(plhs[13]);
   
   plhs[14]=mxCreateDoubleMatrix(sum_arr(N,nflags),asize6/sum_arr(N,nflags),mxREAL);
   aux6out=mxGetPr(plhs[14]);
   
   plhs[15]=mxCreateDoubleMatrix(sum_arr(N,nflags),asize7/sum_arr(N,nflags),mxREAL);
   aux7out=mxGetPr(plhs[15]);

   plhs[16]=mxCreateDoubleMatrix(1,asize8,mxREAL);
   aux8out=mxGetPr(plhs[16]);
  
   plhs[17]=mxCreateDoubleMatrix(1,asize9,mxREAL);
   aux9out=mxGetPr(plhs[17]);
   
   plhs[18]=mxCreateDoubleMatrix((int)output[12],asize10/output[12],mxREAL);
   aux10out=mxGetPr(plhs[18]);
   
   plhs[19]=mxCreateDoubleMatrix(N,N,mxREAL);
   ant_net_out=mxGetPr(plhs[19]);
   
   plhs[20]=mxCreateDoubleMatrix(N,asize_input/N,mxREAL);
   input_iout=mxGetPr(plhs[20]);
   
   plhs[21]=mxCreateDoubleMatrix(N,asize_rec/N,mxREAL);
   rec_iout=mxGetPr(plhs[21]);
   
   plhs[22]=mxCreateDoubleMatrix(N,asize_ant/N,mxREAL);
   ant_iout=mxGetPr(plhs[22]);
   
   mexPrintf("PLASTIC: done\n");
   
  /* copy from C data structures to MATLAB work space */
  mexPrintf("PLASTIC: copying to MATLAB workspace...\n");
  
  *timerec=timesteps*dt;
  copy_double_arr(pre_spike_array,pre_out,*prespikes);
  copy_double_arr(pre_loc_array,pre_locs_out,*prespikes);
  
  copy_double_arr(auxarr,auxout,asize);
  copy_double_arr(aux1arr,aux1out,asize1);
  copy_double_arr(aux2arr,aux2out,asize2);
  copy_double_arr(aux3arr,aux3out,asize3);
  copy_double_arr(aux4arr,aux4out,asize4);
  copy_double_arr(aux5arr,aux5out,asize5);
  copy_double_arr(aux6arr,aux6out,asize6);
  copy_double_arr(aux7arr,aux7out,asize7);
  copy_double_arr(aux8arr,aux8out,asize8);
  copy_double_arr(aux9arr,aux9out,asize9);
  copy_double_arr(aux10arr,aux10out,asize10);
  
  copy_double_arr(post_spike_array,post_out,*postspikes);
  copy_double_arr(post_loc_array,post_locs_out,*postspikes);
   
  copy_double_arr(stim_loc_array,stim_loc_out,stimlist);
  
  copy_double_arr(input_iarr,input_iout,asize_input);
  copy_double_arr(rec_iarr,rec_iout,asize_rec);
  copy_double_arr(ant_iarr,ant_iout,asize_ant);
   
  copy_gsl_mat(input_network,in_net_out,N,inputs);
  copy_gsl_mat(network_network,net_net_out,N,N);
  copy_gsl_mat(antag,ant_net_out,N,N);
 
    
  mexPrintf("PLASTIC: done\n"); 
  
  mexPrintf("PLASTIC: exiting:\n");
  /*Deallocate Memory*/
  mexPrintf("PLASTIC: deleting matricies.\n");
  gsl_matrix_free(input_network);
  gsl_matrix_free(network_network);
  gsl_matrix_free(antag);
  gsl_matrix_free(syn_conduct);
  gsl_matrix_free(syn_conduct_rec);
  gsl_matrix_free(syn_conduct_i);
  gsl_matrix_free(syn_conduct_rec_i);
  gsl_matrix_free(syn_conduct_ant);
  gsl_matrix_free(syn_conduct_ant_i);
/*  gsl_matrix_free(dD);
  gsl_matrix_free(dP); */
  gsl_matrix_free(dDrec);
  gsl_matrix_free(dPrec);
  gsl_matrix_free(dDant);
  gsl_matrix_free(dPant);
  gsl_matrix_free(in_net_connect);
  gsl_matrix_free(net_net_connect);
  gsl_matrix_free(ant_net_connect);
  gsl_matrix_free(input_network_store);
  gsl_matrix_free(scaled_ff);
  gsl_matrix_free(scaled_rec);
  gsl_matrix_free(scaled_ant);

  mexPrintf("PLASTIC: deleting conditionals.\n");  

  if(yd_size>0)
   gsl_matrix_free(yd_protocol);
   
  gsl_rng_free(rgen);
  
  if(output[12]>0.5)
    gsl_matrix_free(track_matrix);
  
  
  mexPrintf("PLASTIC: deleting arrays.\n");

  free(dV);
  free(dg);
  free(dg_i);
  free(last_prespike);
  free(last_postspike);
  free(V);
  free(input_activity);	
  free(uncorr_fir_rate);
  free(firings);
  free(input_connected);

  mexPrintf("PLASTIC: deleting lists.\n");
  
  deallocate_list(*pre_loc_list,*prespikes);
  deallocate_list(*prespike_list,*prespikes);
  deallocate_list(*postspike_list,*postspikes);
  deallocate_list(*post_loc_list,*postspikes);
  
  deallocate_list(aux,asize);
  deallocate_list(aux1,asize1);
  deallocate_list(aux2,asize2);
  deallocate_list(*aux3,asize3);
  deallocate_list(*aux4,asize4);
  deallocate_list(*aux5,asize5);
  deallocate_list(*aux6,asize6);
  deallocate_list(*aux7,asize7);
  deallocate_list(*aux8,asize8);
  deallocate_list(*aux9,asize9);
  deallocate_list(*aux10,asize10);
 
  deallocate_list(*stim_loc_list,stimlist);
  
  deallocate_list(*input_i,asize_input);
  deallocate_list(*rec_i,asize_rec);
  deallocate_list(*ant_i,asize_ant);
  
  mexPrintf("PLASTIC: Bye!\n");
 
}

    
