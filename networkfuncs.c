/* PLASTIC v0.1
/** File containing all of the user defined functions
  * for the integrate and fire network
  * Copyright Guy Billings 2007
  * Contact: g.billings@ucl.ac.uk
  **/
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

/**********************************/
/*Function to push back to a list */
NEXT_LIST_UNIT double_push_back(double value, NEXT_LIST_UNIT old_list)
{
  NEXT_LIST_UNIT newlist=malloc(sizeof(LIST));
   if(newlist==NULL){
     error_out(1);
     }
   else{    
   newlist -> element = value;
   newlist -> next_element = old_list;
   old_list=newlist;
     }
  return old_list;
     
}
/******************************************/
/* Function to convert a list to an array */
void list_to_array(NEXT_LIST_UNIT list, int k, double* vdata)
{
  int i;
  NEXT_LIST_UNIT node=list;
  
 for(i=k-1 ; i>=0 ; i--)
  {
   *(vdata+i)=node->element;
    node=node->next_element;
  } 
}
/********************************************/
/*Function to deallocate linked lists *******/
void deallocate_list(NEXT_LIST_UNIT list,int size)
{
 int i;
 NEXT_LIST_UNIT node=list;
 for(i=0 ; i<size ; i++)
  {
   free(node);
   node=node->next_element;
  }
}  
/***************************************************/
/* Function to print errors to the matlab workbook */
void error_out(int code)
{
 if(code==1){
  mexPrintf("PLASTIC:FATAL! memory allocation failure in double_push_back");
  exit(1);
  }
 if(code==2){
  mexPrintf("PLASTIC:FATAL! memory allocation failure in build_gsl_matrix");
  exit(1);
  } 
 if(code==3){
  mexPrintf("PLASTIC:FATAL! memory allocation failure in build_double_array"); 
  exit(1);
  }
 if(code==4){
  mexPrintf("PLASTIC:FATAL! protocol setting mismatch");
  exit(1);
  }  
}   
/*************************************************************/
/* Function to allow manipulation of matrices via [i] [j] style
 * references rather than pointer arithmatic
 */
void assign_value(int i, int j,int done, double *ptr, double value)
{
 /* Generate address */
 int addr=(j-1)*done+(i-1);
 /* Dereferance the pointer */
 *(ptr+addr)=value;
 }
 /*****************************************************************/
 /* Function to perform the same operations as above but for reading
  * data in from the MATLAB workspace 
  */
double retrieve_value(int i, int j, int done, double *ptr)
{
 double value;
 int addr=(j-1)*done+(i-1);
 value=*(ptr+addr);
 return value;
}
 /****************************************************************/
 /*Function to allow the accumulation of one array in to another */
 void accum_double_array(double *in, double *out, int size)
 {
  int i;
  
  for(i=0 ; i<size ; i++)
  {
    out[i]=out[i]+in[i];
  }
  
 }   
 /*************************************************************/
 /* Function to copy arrays in to matlab ones */
 void copy_double_arr(double *in, double *out, int size)
 {
  int i;
  
  for(i=0 ; i<size ; i++)
  {
    *(out+i)=in[i];
  }
 }      
 /***********************************************/
 /* Function to copy data from MATLAB matricies */
 void import_gsl_matrix(gsl_matrix *out, double *in, int k, int l)
 {
  int i, j;
  
  for(i=0 ; i<k ; i++)
  {
   for(j=0 ; j<l ; j++)
   {
    gsl_matrix_set(out,i,j,retrieve_value(i+1,j+1,k,in));
   }
  }
 }   
 /*******************************************************************/ 
 /* function to copy values from a dynamic C matrix to a MATLAB one */
 void copy_gsl_mat(gsl_matrix *in, double *out, int k, int l)
 {
  /* k and l define size of input array; obviously the two arrays must
   * be of the same dimensions!
   */
  int i, j;
  
  for(i=0 ; i<k ; i++)
  {
   for(j=0 ; j<l ; j++)
   {
    assign_value(i+1,j+1,k,out,gsl_matrix_get(in,i,j));
   }
  }
 }    
/**************************************/
/* Matrix creation function using GSL */
gsl_matrix* build_gsl_matrix(int done, int dtwo, double init, int mode)
{
 int i,j;
 extern gsl_rng *rgen;
 gsl_matrix* m;
 m = gsl_matrix_alloc (done, dtwo);
 if(m==NULL){
   error_out(2);
 }
 else{  
 if(mode==0){ /*For mode = 0 fill with init*/
  for (i = 0; i < done; i++)
    for (j = 0; j < dtwo; j++)
      gsl_matrix_set (m, i, j,init);
  }
 else {     
 for (i = 0; i < done; i++)
    for (j = 0; j < dtwo; j++)
      gsl_matrix_set (m, i, j,custom_rand(init,rgen));
  }
 }
 return m;     
}                 
/********************************************************/ 
/* Function to allow  the creation of dynamic 1d array */
double * build_double_array(int done)
{
    int i;
    double *arr;
    arr=malloc(done * sizeof(double));
    if(arr == NULL)
    {
     /* Memory full */
      error_out(3);
    }
                 
    for(i = 0; i < done ; i++) {
    arr[i] = 0.0;
    }
                
 return arr;

}
/*****************************************************************************/
/* Function that returns a random value between two bounds drawn from a      */
/* uniform distribution using the currently selected random number generator */
double custom_rand(double maxi, gsl_rng *generator)
{
 double varate;
 /* use high order bits */
 varate=maxi*gsl_rng_get(generator)/gsl_rng_max(generator);
 /*varate=maxi*(rand()/(RAND_MAX+1.0));*/
 return varate;
}
/**********************************************************************************/
/* Function to perform thresholding of random variable and return a boolean value */
bool rand_thresh(double threshold) 
{
  double rando_var;
  extern gsl_rng *rgen;
  int test;
  
  rando_var= custom_rand(1.0,rgen);
  
  if(rando_var <= threshold)
  { 
    test=1;
  }
  else
  {
    test=0;
  }
  
  return test;
  
}
/***************************************************************/
double gaussrand(double sigma)
{
        return gsl_ran_gaussian (rgen,sigma); 
}
/***************************************************/    
/* Function to create random input activity array  */
void input_activity_random(int inputs, double* activity)
{
  register i;
  for(i=0; i<inputs ; i++) 
  {
      *(activity+i)=gaussrand(1.0);
  }    
     
} 
/**************************************************/
/* Function to sample an exponential distribution */
double exp_dist(double mean_time)
{
  double var;
  extern gsl_rng *rgen;
  
  var=log(custom_rand(mean_time,rgen));
  
  return var;
}  
/*************************************************/
/* Function to multiply an array by a scalar *****/
void scal_mult_arr(double scal,int dim, double *input)
{
 int i;
 for(i=0 ; i<dim ; i++) 
 {
  *(input+i)=scal*(*(input+i));
 }
}
/****************************/
/* Function to sum an array */
double sum_arr(int dim, double *arr)
{
 int i;
 double acc=0;
 for(i=0 ; i<dim ; i++)
 {
  acc=acc+arr[i];
 }
 return acc;
}  
/*****************************************************/
/* Function to add a scalar to an array element wise */
void scal_add_arr(double scal,int dim, double *input)
{
 int i;
 for(i=0 ; i<dim ; i++) 
 {
  *(input+i)=scal+*(input+i);
 }
  
}
/******************************************/
/* Function to multiply two gsl_matricies */
void mult_matrix(gsl_matrix *prod, gsl_matrix *m1, gsl_matrix *m2)
{
  gsl_matrix_memcpy(prod, m1);
  gsl_matrix_mul_elements(prod,m2);
}
/***************************************************************/
/* Function to allow some inputs to be correlated and some not */
void mixed_rate(int inputs, double *activity, double *rates, double cor_rate, double *flag, double base, double var_rate)
{
 register i;
 double common=gaussrand(1.0);
 for(i=0 ; i<inputs ; i++)
 {
  if(flag[i]>0.5)
   *(rates+i)=base*(1.0+(var_rate*activity[i])+cor_rate*common);
  else 
   *(rates+i)=base*(1.0+(var_rate*activity[i]));   
 }
}
/***************************************************************/
/* Function to allow spatial inhomogeniety of the input */
void inhm_rate(int inputs, double *activity, double *rates, double *target_rates, double *rate_vars)
{
 register i;
 for(i=0 ; i<inputs ; i++)
 {
 	*(rates+i)=target_rates[i]*(1.0+(rate_vars[i]*activity[i]));
 }
}
/******************************************************/
/*Function to provide localised centered stimulus rates*/
void localised_stimulus(int inputs, double *fir_rate, double r0, double r1, double sigma, int stim_loc, double lamda)
{
 register i;
 double ponent0, ponent1, ponent2;
 for(i=0 ; i<inputs ; i++)
 { 
  ponent0=1/exp(((stim_loc-i)*(stim_loc-i))/(2*sigma*sigma));
  ponent1=1/exp(((stim_loc+lamda-i)*(stim_loc+lamda-i))/(2*sigma*sigma));
  ponent2=1/exp(((stim_loc-lamda-i)*(stim_loc-lamda-i))/(2*sigma*sigma));

  fir_rate[i]=r0+r1*(ponent0+ponent1+ponent2);
 }
}
/************************************************/
/* Function to calculate the presynaptic inputs */
void syn_input(int N,int inputs,gsl_matrix* weight, gsl_matrix *synputs, gsl_matrix *synputs_i,
gsl_matrix *connect_in, double* dg, double* dg_i, double* fire, double dt, double T,double mult)
{
 unsigned register i,j;
 for(i=0 ; i<N ; i++)
 {
  for(j=0 ; j<inputs ; j++)
  { 
   if(gsl_matrix_get(connect_in,i,j)>0.5){
  
    gsl_matrix_set(synputs,i,j,gsl_matrix_get(synputs,i,j)*(1-(dt/T)));
 
    if(fire[j] > 0.5)
     gsl_matrix_set(synputs,i,j,gsl_matrix_get(synputs,i,j)+mult*gsl_matrix_get(weight,i,j));
     
     dg[i]=dg[i]+gsl_matrix_get(synputs,i,j); /*sum up conductance alterations */ 
    }          
   else if(gsl_matrix_get(connect_in,i,j)< -0.5){
  
    gsl_matrix_set(synputs_i,i,j,gsl_matrix_get(synputs_i,i,j)*(1-(dt/T)));

    if(fire[j] > 0.5)
     gsl_matrix_set(synputs_i,i,j,gsl_matrix_get(synputs_i,i,j)+gsl_matrix_get(weight,i,j));

     dg_i[i]=dg_i[i]+gsl_matrix_get(synputs_i,i,j);

    }
  } 
 }
}
/***********************************************/
/* Function to update the cell potential matrix */
void  update_V(int N, double dt, double Vrest, double *V, double *dV, double *dg, double *dg_i, double Eex,double Ein,double Tm)
{
 unsigned register i;
 for(i=0 ; i<N ; i++)
 {
  dV[i]=(Vrest-V[i]+dg[i]*(Eex-V[i])+dg_i[i]*(Ein-V[i]))*(dt/Tm);
  V[i]=V[i]+dV[i];
  dg[i]=0.0;
  dg_i[i]=0.0;
 } 
}   
/****************************************************/
/*Function to add driving noise to a network neuron */
void drive_noise(int neurons, double freq, double dt, double g_input, double *dg)
{
 unsigned register i;
 /*extern gsl_rng *rgen;*/
 for(i=0 ; i<neurons ; i++)
 { 
   /*if(custom_rand(1.0,rgen) < freq*dt) { */
   /* Add a gaussian conductance */
     *(dg+i)=*(dg+i)+g_input*(gaussrand(1)+1);
   /*}*/
 }
}
/********************************************************/
/* Function to find if an input is connected to a network
/* neuron */
void find_connections(int N, int inputs, gsl_matrix *in_net_connect, double *input_connected)
{
 unsigned int i,j;
 
 for(i=0 ; i<inputs ; i++){
  for(j=0 ; j<N ; j++){
   input_connected[i]=input_connected[i]+abs(gsl_matrix_get(in_net_connect,j,i));
   }
  }
}   
/************************************************************/
/* Function to determine if recurrent connections exist *****/
int find_recurrents(int N, gsl_matrix* connect)
{
 unsigned int i,j; 
 int rec=0;
 double conn;
 
 for(i=0; i<N ; i++){
   for(j=0; j<N ; j++){
     conn=gsl_matrix_get(connect,i,j);
     if(conn< -0.5 || conn > 0.5)
       rec=1;
     }
   }
   
  return rec;
}  
/************************************************************/
/*Function to update weights with no weight dependance ******/
void update_weights(int N, int inputs, int mode, gsl_matrix* weight,gsl_matrix* weight_change_factor_a,
gsl_matrix* weight_change_factor_b, gsl_matrix* connected_in, double gmax, double* fire, double rate, int perturb,
int scaled,gsl_matrix *dN, double noise_dec, gsl_matrix *scaled_c, double dt, double taun, double sdn, double timenow,
gsl_matrix* last_event_a, gsl_matrix* last_event_b,double Ta,double Tb) 
{
  unsigned register i,j;
  double temp,temp2,temp3,last_time;
  int index;
  
  for(i=0 ; i<N ; i++)
  {
   for(j=0 ; j<inputs ; j++)
   {
    if(gsl_matrix_get(connected_in,i,j)>0.5){
    
    if(perturb==1){
     
      if(scaled==1){
        gsl_matrix_set(dN,i,j,gsl_matrix_get(dN,i,j)*noise_dec+gsl_matrix_get(scaled_c,i,j)*gaussrand(1)*sqrt(dt/taun));      
        temp3=gsl_matrix_get(weight,i,j)+gsl_matrix_get(dN,i,j)*dt;
	}
      else if(scaled==0){
        gsl_matrix_set(dN,i,j,gsl_matrix_get(dN,i,j)*noise_dec+sdn*gaussrand(1)*sqrt(dt/taun));      
        temp3=gsl_matrix_get(weight,i,j)+gsl_matrix_get(dN,i,j)*dt;
	}
	
      if(temp3<0.0)
        gsl_matrix_set(weight,i,j,0.0);
      else if (temp3>gmax) 
        gsl_matrix_set(weight,i,j,gmax);
      else if(temp3>0.0 && temp3 < gmax)
        gsl_matrix_set(weight,i,j,temp3); 
     }
     
     if(mode<1) /*update for presynaptic firing */
       index=j;
     else if(mode>0) /*update for postsynaptic firing */
       index=i;
          
     if(fire[index] > 0.5){
           
       gsl_matrix_set(weight_change_factor_a,i,j,gsl_matrix_get(weight_change_factor_a,i,j)
       *exp(-1*((timenow-gsl_matrix_get(last_event_a,i,j))/Ta))+rate);
      
       gsl_matrix_set(last_event_a,i,j,timenow);
       
       temp=gsl_matrix_get(weight,i,j);
       
       temp2=gsl_matrix_get(weight_change_factor_b,i,j)*exp(-1*((timenow-gsl_matrix_get(last_event_b,i,j))/Tb));
       
       if(temp>0 && mode==0)
        gsl_matrix_set(weight,i,j,(temp-gmax*temp2));
       
       else if(temp<gmax && mode==1)
        gsl_matrix_set(weight,i,j,(temp+gmax*temp2));
	
     /* Automatically clip weights */
      if(gsl_matrix_get(weight,i,j)>gmax) 
       gsl_matrix_set(weight,i,j,gmax);
     
      else if(gsl_matrix_get(weight,i,j)<0) 
       gsl_matrix_set(weight,i,j,0.0);
     }
    } 
   } 
  }  
}
/**********************************************************/
/*Function to update weights with weight dependance *******/
void update_weights_dep(int N, int inputs, int mode, gsl_matrix* weight,gsl_matrix* weight_change_factor_a,gsl_matrix* weight_change_factor_b, 
gsl_matrix* connected_in, double gmax, double* fire, double c_chan,double sigma, int perturb, int scaled, gsl_matrix *dN, gsl_matrix *scaled_c, 
double noise_dec, double dt, double taun, double sdn, double timenow, double Ta, double Tb, gsl_matrix *last_event_a,gsl_matrix *last_event_b) 
{
  unsigned register i,j;
  double temp,temp2,temp3;
  
  for(i=0 ; i<N ; i++)
  {
   for(j=0 ; j<inputs ; j++)
   {
    if(gsl_matrix_get(connected_in,i,j)>0.5){
    
     if(perturb==1){
     
      if(scaled==1){
        gsl_matrix_set(dN,i,j,gsl_matrix_get(dN,i,j)*noise_dec+gsl_matrix_get(scaled_c,i,j)*gaussrand(1)*sqrt(dt/taun));      
        gsl_matrix_set(weight,i,j,gsl_matrix_get(weight,i,j)+gsl_matrix_get(dN,i,j)*dt);
	}
      else if(scaled==0){
        gsl_matrix_set(dN,i,j,gsl_matrix_get(dN,i,j)*noise_dec+sdn*gaussrand(1)*sqrt(dt/taun));      
        gsl_matrix_set(weight,i,j,gsl_matrix_get(weight,i,j)+gsl_matrix_get(dN,i,j)*dt);
	}
      }	
    
     if(mode<1){ /*update factor a ; depress with factor b */
       if(fire[j] > 0.5){
       temp=gsl_matrix_get(weight,i,j);
       
       gsl_matrix_set(weight_change_factor_a,i,j,gsl_matrix_get(weight_change_factor_a,i,j)
       *exp(-1*((timenow-gsl_matrix_get(last_event_a,i,j))/Ta))+(c_chan+gaussrand(sigma)*temp));
   
       gsl_matrix_set(last_event_a,i,j,timenow);
      
       temp2=gsl_matrix_get(weight_change_factor_b,i,j)*exp(-1*((timenow-gsl_matrix_get(last_event_b,i,j))/Tb));

	/*Update weight and clip*/	
       if((temp-temp2)<0.0)	
       gsl_matrix_set(weight,i,j,0.0);
       else	
       gsl_matrix_set(weight,i,j,(temp-temp2));
       
      }
     }
     else { /* update factor a ; potentiate with factor b */
       if(fire[i] > 0.5){
       temp=gsl_matrix_get(weight,i,j);
       
       gsl_matrix_set(weight_change_factor_a,i,j,gsl_matrix_get(weight_change_factor_a,i,j)
       *exp(-1*((timenow-gsl_matrix_get(last_event_a,i,j))/Ta))+(c_chan*temp-gaussrand(sigma)*temp));
       
       gsl_matrix_set(last_event_a,i,j,timenow);
       
       temp2=gsl_matrix_get(weight_change_factor_b,i,j)*exp(-1*((timenow-gsl_matrix_get(last_event_b,i,j))/Tb));
      
       gsl_matrix_set(weight,i,j,(temp+temp2));
     }
     
    }
     
   } 
  } 
 }  
}
  

/************************************************/
/* Function to determine the presynaptic firing */
void pre_fire(int inputs, double dt, double *rate, double *last_prespike, double *connect)
{
 unsigned register i;

  /* Find spike times */
      for(i=0 ; i<inputs ; i++){         
       if(connect[i] > 0.5 || connect[i] < -0.5){
        if(rand_thresh(rate[i]*dt) == 1){
         last_prespike[i]=1.0;
        }
        else{
         last_prespike[i]=0.0;
        } 	
      }
     }
 
}
/********************************************/
/*Function to determine postsynaptic firing */
void post_fire(int N, double spike_thresh, double Vinit, double *last_postspike, double *V)
{
  unsigned register i;
 /* Find output spike times */
      for(i=0; i<N ; i++)
      { 
       if(V[i]>=spike_thresh) { 
        last_postspike[i]=1.0;
        V[i]=Vinit;
       }        
       else{     
	last_postspike[i]=0.0;
       }	
      }
}
/********************************************************************/
/* Function to determine whether or not the program should be halted*/
int the_terminator(int nsteps,int steps,int protocol,int yd_steps, int yd_size)
{
 unsigned int arnie=0;
 
   if(protocol<3 && steps >= nsteps)
    arnie=1;
    
   else if(protocol==3 && yd_steps>=yd_size)
    arnie=1;

   else if(protocol==4 && steps >= nsteps)
    arnie=1;
    
   else if(protocol==5  && steps >= nsteps)
    arnie=1;
    
 return arnie;
}         
/*****************************************************************/
/* Function to calculate the overlap between two weight matrices */
double overlap_vec(int N, int inputs,gsl_matrix *connected_in, gsl_matrix *weights_one, gsl_matrix *weights_two) 
{
 unsigned register i,j;
 double acone=0.0;
 double actwo=0.0;
 double ov=0.0;
 
 /* NOTE: Only takes excitatory connections */
 
 for(i=0 ; i<N ; i++){
  for(j=0 ; j <inputs ; j++){
  if(gsl_matrix_get(connected_in,i,j)>0.5){
   acone=acone+gsl_matrix_get(weights_one,i,j)*gsl_matrix_get(weights_one,i,j);
   actwo=actwo+gsl_matrix_get(weights_two,i,j)*gsl_matrix_get(weights_two,i,j);
   ov=ov+gsl_matrix_get(weights_one,i,j)*gsl_matrix_get(weights_two,i,j);
   }
  }
 }
 
 return ov/(sqrt(acone)*sqrt(actwo));
}
/****************************************************************/
/* Function to determine if a neuron has fired from its listing */   
void any_of(int size, double *firing, int *rval)
{
 unsigned register i;
 
 for(i=0 ; i<size ; i++){
   if(firing[i]>0.5){
    rval[0]++;
    rval[1]=i; 
    }
  }
    
}
/*******************************/
/* Function to record prepikes */
void record_spikes(int *firings,int protocol,int sloc, double *last_spikes,NEXT_LIST_UNIT *locs, NEXT_LIST_UNIT *spike_list, NEXT_LIST_UNIT *spike_loc, int *spikes, double time, int size) 
{

  unsigned register i,j=0;
  double *flist;

  if(firings[0]==1){
    *spike_list=double_push_back(time,*spike_list);
    *spike_loc=double_push_back(firings[1],*spike_loc);
    *spikes=*spikes+1;
  }
  else if(firings[0]>1){
     flist=build_double_array(firings[0]);

     for(i=0 ; i<size ; i++){
       if(last_spikes[i]>0.5){
         flist[j]=i;
         j++;
       } 
     }
     for(i=0 ; i<j ; i++){
      *spike_list=double_push_back(time,*spike_list);
      *spike_loc=double_push_back(flist[i],*spike_loc);
      *spikes=*spikes+1;
     }
   free(flist);
  }
}
/*****************************************************************/
/* Function to return the average weights for each of the neurons*/
void average_weights(int N, int inputs, double *flags, gsl_matrix *weights, gsl_matrix *connect, NEXT_LIST_UNIT *list, double var_switch, NEXT_LIST_UNIT *var_list)
{
 unsigned register i,j;
 unsigned int entries=0;
 double avg=0;
 double var=0;
 
 /* NOTE: Only takes excitatory connections */
 
 for(i=0 ; i<N ; i++){
  if(flags[i]>0.5){
  for(j=0 ; j<inputs ; j++){
    if(gsl_matrix_get(connect,i,j)>0.5){
      avg=avg+gsl_matrix_get(weights,i,j);
      entries++;
    }
   }
   avg=avg/entries;
  *list=double_push_back(avg,*list);
  entries=0;
  if(var_switch>0.5){
   for(j=0 ; j<inputs ; j++){
    if(gsl_matrix_get(connect,i,j)>0.5){
      var=var+pow((gsl_matrix_get(weights,i,j)-avg),2);
      entries++;
    }
   }
  *var_list=double_push_back(var/entries,*var_list); 
  }
  avg=0; 
  entries=0;
  }
 } 
 
}  
/********************************************************************/
/* Function to find firing frequency from spike accumulation arrays.*/
void find_frequency(int N, double *accarray, double *flags, NEXT_LIST_UNIT *freqlist, double time_int)
{
  unsigned register i;
  
  for(i=0 ; i<N ; i++)
  {
   if(flags[i]>0.5)
    *freqlist=double_push_back(accarray[i]/time_int,*freqlist);
  }
  
}   
/*************************************************/
/* function to return sensible values for stimulus
   location test of a learned network */
int test_stimloc(int inputs, int stimloc)
{
 /* When called switch the input... but
    'wrap around' */

    unsigned int newloc;

    stimloc++;
    if(stimloc > inputs)
      newloc=stimloc-inputs;
    else
      newloc=stimloc;
 
    return newloc;

}

/********************************************************************/
/* Function to calculate the Euclidean distance between two vectors */
void euclid(int done, int dtwo, gsl_matrix *weights, gsl_matrix *stored, gsl_matrix *connected, NEXT_LIST_UNIT *out)
{
 
 unsigned register i,j;
 double sumsquare=0;
 
 /* NOTE: Only takes excitatory connections */

 for(i=0 ; i<done ; i++)
 {
  for(j=0 ; j<dtwo ; j++)
  {
   if(gsl_matrix_get(connected,i,j)>0.5)
    sumsquare=sumsquare+pow(gsl_matrix_get(stored,i,j)-gsl_matrix_get(weights,i,j),2);
  }
 }
 
 *out=double_push_back(sqrt(sumsquare),*out);
}
/******************************************************************/
/* Function to calculate the average across a whole weight matrix */
double av_all(int done, int dtwo, gsl_matrix *weight, gsl_matrix* connected)
{

 unsigned register i,j;
 double average=0;
 int elements=0;

 for(i=0 ; i<done ; i++)
 { 
  for(j=0 ; j<dtwo ; j++)
  {
   if(gsl_matrix_get(connected,i,j)>0.5){
    average=average+gsl_matrix_get(weight,i,j);
    elements++;
   }
  }
 }
 return(average/elements);
} 
/************************************************************/
/* Function to calculate Pearson's R for the wieght vectors */
void pearson(int done, int dtwo, gsl_matrix *x, gsl_matrix *y, gsl_matrix *connected, NEXT_LIST_UNIT *out)
{

  unsigned register i,j;
  double temp1,temp2;
  double sx=0;
  double sy=0;
  double sxy=0;

  double sxs=0;
  double sys=0;

  double pear;

  int n=0;
 
  for(i=0 ; i<done ; i++)
  {
   for(j=0 ; j<dtwo ; j++)
   {
    if(gsl_matrix_get(connected,i,j)>0.5){
      temp1=gsl_matrix_get(x,i,j);
      temp2=gsl_matrix_get(y,i,j);
  
      sx=sx+temp1;
      sy=sy+temp2;
  
      sxy=sxy+temp1*temp2;

      sxs=sxs+temp1*temp1;
      sys=sys+temp2*temp2;    
     
      n++;
     }
    }
   }

 pear=(n*sxy-sx*sy)/sqrt((n*sxs-sx*sx)*(n*sys-sy*sy));

 *out=double_push_back(pear,*out);
}
/*******************************************************/
/* Function to track weights specified by track_matrix */
void track_weight(int number, gsl_matrix *track_matrix, gsl_matrix *weights, NEXT_LIST_UNIT *out)
{
   unsigned int i;
   double temp;
   
   for(i=0 ; i<number ; i++){
   	temp=gsl_matrix_get(weights,gsl_matrix_get(track_matrix,i,0)-1,gsl_matrix_get(track_matrix,i,1)-1 );
	*out=double_push_back(temp,*out);
	}
}
/*******************************************************/
/* Function to total up conductances for output ********/
void total_current(int done, int dtwo, gsl_matrix *cmatrix, gsl_matrix *connected, double *total)
{
 unsigned register i,j;
 
 for(i=0 ; i<done ; i++){
   for(j=0 ; j<dtwo ; j++){
     if(gsl_matrix_get(connected,i,j)< -0.5){
       total[i]=total[i]-gsl_matrix_get(cmatrix,i,j);
     }
     else if(gsl_matrix_get(connected,i,j)>0.5){
       total[i]=total[i]+gsl_matrix_get(cmatrix,i,j);
     }  
    }
   }
}
/************************************************************/
/* Function to output average current after specified time***/
void out_current(double *total, double E, double *V, NEXT_LIST_UNIT *out, int timesteps, int done)
{
 unsigned int i;
 for(i=0 ; i<done ; i++){
   total[i]=total[i]/(timesteps);
   *out=double_push_back(total[i]*(E-V[i]),*out);
/*   mexPrintf("%4f\n",V[i]);*/
/*   mexPrintf("%4f\n",total);*/
   total[i]=0;     
  }
}   
/*******************************************************************/ 
 /* function to get weight snapshot ********************************/
 void snapshot(gsl_matrix *in, double *out, int k, int l,int index,int ss_size)
 {
  
  int i,j,addr;
  
  if(index<0)
    index=0;
  
  for(i=0 ; i<k ; i++)
  {
   for(j=0 ; j<l ; j++)
   {
    addr=j*k*ss_size+i*ss_size+index;
    *(out+addr)=gsl_matrix_get(in,i,j);
   }
  }
 }    
/**************************************/	
	












 
     
 

