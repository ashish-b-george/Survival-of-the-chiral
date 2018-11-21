//  Created by Ashish Bino George
//note to self: this code used to be called "new_neutralHdecay_passedpa_fstar.c"

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
//#include <gsl/gsl_sf_bessel.h>

#include <gsl/gsl_randist.h>

#include <gsl/gsl_rng.h>

#include <gsl/gsl_math.h>





int alpha=1; // =1 for every simulation in the paper. not mentioned in the paper

/* just some random default parameters, parameters are usually passsed to the program */
int l= 100; // width of lattice in x direction
int b= 100; // length of lattice in y, the expansion direction
int r= 40; // length of initial condition in y if rectangular or radius of circle in IC
int Number= 30000; // carrying capacity, determines  strength of noise
int T= 150; // number of time steps in simulation

double g= 0.02;// growth rate
double ga= 0.02;// growth rate of a
double gb= 0.02;// growth rate of b
double p00 = 0.1; // Diffusion constant for constant diffusion, a_0 .
/* density dependent diffusion constants defined below, notation in paper in comments.*/
double p0= 0.00; // a_s, diffusion dependency on source.
double p1= 0.1; // a_d
double p2= 0.0; //a_l
double p3 =0.00;// a_b
double p4 =-0.1; //a_r
/* same as above but for the second strain. */
double p22=0.0;
double p42=0.0;
double p32=0.0;
double p12=0.0;
double sel=0;


double h =0.0; // ratio of fractions in Initial condition
//double B=11.0; // if required, an Allee effect. not used in the paper



// GLOBAL RANDOM NUMBER DECLARATION AND SEEDING
/* specifying to use Mersenne twister MT-19937 as the uniform PRNG */
gsl_rng *gBaseRand;       /* global rand number generator */

unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}



// to initialise double arrays using pointers
void initialise( double **handle, int n, int m)
{
    *handle= (double *)malloc(sizeof(double)*n*m);
    if (*handle==NULL)
    {
        printf("pointer not given memory. ERROR! ");
        exit(10);
    }
    for (int i = 0; i < n*m; i++)
    *(*handle + i) = 0.;
    
}

/* functions f1, f2 f3, f4 caclulate migrations to the demes in the 4 directions in a roationally invariant fashion  */
// f1 is to the right, f2 top f3 left f4 bot
void f1(double *n, double *p,double p1, double p2, double p3, double p4, int beg, int end)
{
    int ip1,im1,jp1,jm1;
    // probabilities are only calculcated for the interior and not for the boundaries
    for (int i=0;i <l;i++)
    {
        for (int j=beg;j <end;j++)
        {
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
            if (i==l-1)
                ip1=0;
            if (i==0)
                im1=l-1;
            if (j==b-1)
                jp1=j;
            if (j==0)
                jm1=j;
            *(p+ i*b +j)=pow(p0**(n+ i*b +j)+ p1* *(n+ (ip1)*b +j)+p2* *(n+ (i)*b +jp1)+p3* *(n+ (im1)*b +j)+p4* *(n+ (i)*b +jm1)+p00,alpha);
        }
    }
    
}
void f2(double *n, double *p,double p1, double p2, double p3, double p4, int beg, int end)
{
    int ip1,im1,jp1,jm1;
    // probabilities are only calculcated for the interior and not for the boundaries
    for (int i=0;i <l;i++)
    {
        for (int j=beg;j <end;j++)
        {
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
            if (i==l-1)
                ip1=0;
            if (i==0)
                im1=l-1;
            if (j==b-1)
                jp1=j;
            if (j==0)
                jm1=j;
            
            *(p+ i*b +j)=pow(p0**(n+ i*b +j)+ p4* *(n+ (ip1)*b +j)+p1* *(n+ (i)*b +jp1)+p2* *(n+ (im1)*b +j)+p3* *(n+ (i)*b +jm1)+p00,alpha);
        }
    }
}
void f3(double *n, double *p,double p1, double p2, double p3, double p4, int beg, int end)
{
    int ip1,im1,jp1,jm1;
    // probabilities are only calculcated for the interior and not for the boundaries
    for (int i=0;i <l;i++)
    {
        for (int j=beg;j <end;j++)
        {
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
            if (i==l-1)
                ip1=0;
            if (i==0)
                im1=l-1;
            if (j==b-1)
                jp1=j;
            if (j==0)
                jm1=j;
            
            *(p+ i*b +j)=pow(p0**(n+ i*b +j)+ p3* *(n+ (ip1)*b +j) + p4* *(n+ (i)*b +jp1) + p1* *(n+ (im1)*b +j) + p2* *(n+ (i)*b +jm1)+ p00,alpha);
        }
    }
}

void f4(double *n, double *p,double p1, double p2, double p3, double p4, int beg, int end)
{
    int ip1,im1,jp1,jm1;
    // probabilities are only calculcated for the interior and not for the boundaries
    for (int i=0;i <l;i++)
    {
        for (int j=beg;j <end;j++)
        {
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
            if (i==l-1)
                ip1=0;
            if (i==0)
                im1=l-1;
            if (j==b-1)
                jp1=j;
            if (j==0)
                jm1=j;
            *(p+ i*b +j)=pow(p0**(n+ i*b +j)+ p2* *(n+ (ip1)*b +j)+p3* *(n+ (i)*b +jp1)+p4* *(n+ (im1)*b +j)+p1* *(n+ (i)*b +jm1)+p00,alpha);
        }
    }
}


void time_evolution( double *n, double *na,double p1, double p2, double p3, double p4,double g, int beg, int end, double *patemp, double *a1, double *a2, double *a3, double *a4 )
{
    /*
     patemp and a1 a2 a3 a4 are declared and initialised once in main and then is just passed here
     this makes the code much much faster!!
     
    double *patemp;
    double *a1;
    double *a2;
    double *a3;
    double *a4;
    initialise(&patemp,l,b);
    initialise(&a1,l,b);
    initialise(&a2,l,b);
    initialise(&a3,l,b);
    initialise(&a4,l,b);
    */
    
    /*
     for (int i=0;i <l;i++)
     {
     *(a1+ i*b +0)= *(a1+ i*b +b-1)=0;
     *(a2+ i*b +0)= *(a2+ i*b +b-1)=0;
     *(a3+ i*b +0)= *(a3+ i*b +b-1)=0;
     *(a4+ i*b +0)= *(a4+ i*b +b-1)=0;
     }
     
     
     for (int j=0;j <b;j++)
     {
     *(a1+0*b+j)= *(a1+ (l-1)*b +j)=0;
     *(a2+0*b+j)= *(a2+ (l-1)*b +j)=0;
     *(a3+0*b+j)= *(a3+ (l-1)*b +j)=0;
     *(a4 +0*b+j)= *(a4+ (l-1)*b +j)=0;
     
     }
     */
    
    
    /*functions f1,f2,f3,f4 caclulated the proabbailities to migration in the various directions. beg and end are specify the simulation box limits */
    f1( n, a1,p1,p2, p3,p4, beg ,end);
    
    f2( n, a2,p1,p2, p3,p4, beg, end);
    
    f3( n, a3,p1,p2, p3,p4, beg ,end);
    
    f4( n, a4,p1,p2, p3,p4, beg , end);
    
    
    
    /* for loop over all the sites in the simualtion box (between beg and end in the y direction) to update concentrations*/
    int ip1,im1,jp1,jm1;
    for (int i=0;i <l;i++)
    {
        for (int j=beg+1;j <end-1;j++) // since there are only probabilities for upto the penulitmate point i.e 1tol-2, the updates can  only be for the next inner layer 2tol-3
        {
            /* ip1, im1, etc stand for "i+1 "and "i-1". This definition implements the Periodic boundary conditions in x and the fixed BC in y  (BC in y direction is not important for the phenomena. */
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
            if (i==l-1)
            ip1=0;
            if (i==0)
            im1=l-1;
            if (j==b-1)
            jp1=j;
            if (j==0)
            jm1=j;
            
            /* checks for negative numbers or overflows. This is only occurs if parameters were specified incorrectly. */
            if(*(a1+ i*b +j)+*(a2+ i*b +j)+*(a3+ i*b +j)+*(a4+ i*b +j)>1.0)
            {
                printf("sum of probabilities became greater than one at site %d, %d and is = %f", i,j,*(a1+ i*b +j)+*(a2+ i*b +j)+*(a3+ i*b +j)+*(a4+ i*b +j));
            }
            if(*(a1+ i*b +j)<-0.0001||*(a2+ i*b +j)<-0.0001||*(a3+ i*b +j)<-0.0001||*(a4+ i*b +j)<-0.0001)
            {
                printf(" probabilities lesser one at site %d, %d and is = %f, %f,%f,%f \n", i,j,*(a1+ i*b +j),*(a2+ i*b +j),*(a3+ i*b +j),*(a4+ i*b +j));
            }
            
            
            // update of concentrations to a temporary variable to preserve simultaneity
            *(patemp+ i*b +j)= *(na+ i*b +j)+ (1-*(n+ i*b +j))*( *(a1+ (im1)*b +j)* *(na+ (im1)*b +j)+  *(a2+ (i)*b +jm1)* *(na+ (i)*b +jm1)+  *(a3+ (ip1)*b +j)* *(na+ (ip1)*b +j)+  *(a4+ (i)*b +jp1)* *(na+ (i)*b +jp1) ) -1*( (1-*(n+ (ip1)*b +j))* *(a1+ i*b +j)+ (1-*(n+ i*b + jp1))* *(a2+ i*b +j)+ (1-*(n+ (im1)*b +j))* *(a3+ i*b +j)+ (1-*(n+ i*b + jm1))**(a4+ i*b +j)) * *(na+ i*b +j)+ g**(na+ i*b +j)*(1-*(n+ i*b +j));
            
            
        }
    }
    

    /* copies patemp values into na to update na for the next time step. Doing it in this fashion preserves simultaneity fo concentration updates */
    for (int i=0;i <l;i++)
    {
        for (int j=beg+1;j <end-1;j++)
        {
            *(na+ i*b +j)=*(patemp+ i*b +j);
            
        }
    }

}


/* implements demographic noise using binomial sampling. */
void selection_noise( double *n, double *na, double *nb, const gsl_rng *rsd, int beg, int end)
{
    double hin[3];
    unsigned int hout[3];
    
    int N= Number;
    double No=Number;
    for (int i=0;i <l;i++)
    {
        for (int j=beg;j <end;j++)
        {
            if(*(n+ i*b +j)>0) //introduced some hardcutoff on noise!!
            {
                
                hin[0]=1-*(n+ i*b +j);   // concentration of vacancies or unoccupied space at a site
                hin[1]=*(na+ i*b +j);    // concentration of a at a site
                hin[2]=*(nb+ i*b +j); // concentration of b at a site
                //uses "gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n)" to do binomial sampling
                
                hout[0]= gsl_ran_binomial (rsd, hin[0], N );// calculate new number of vacancies, and fluctuations in overall growth
                
                hout[1]=gsl_ran_binomial (rsd, hin[1]/(hin[1]+hin[2]), (N-hout[0])); // n-hout is Nc or C. This calculates the changes in the fractions of a and b

                
                *(n+ i*b +j)=1-(double)hout[0]/No;
                *(na+ i*b +j)=(double)hout[1]/No;    // concentration of a
                *(nb+ i*b +j) = *(n+ i*b +j) -*(na+ i*b +j);
            }
        }
    }
   
}



/* function for  creating circular initial conditions*/
void initial_conditions_circle(double *ca, double *cb, double *c)
{
    for (int i = 0; i <l ;i++)
    {
        for (int j =0;j <b;j++)
        {
            *(ca+ i*b +j)=0.0;
            *(cb+ i*b +j)=0.0;
            *(c+ i*b +j)=0.0;
        }
    }
    for (int i = l/2 -r ; i <=  l/2 + r;i++)
    {
        for (int j =b/2 -r;j <= b/2+ r;j++)
        {
            if(((i-l/2)*(i-l/2)+(j-b/2)*(j-b/2))<=r*r)
            {
                
                // concentraions are equal for a and b withing radius r and less than the carrying capacity
                *(ca+ i*b +j)=0.10;
                *(cb+ i*b +j)=0.10;
                
                /*to make sectors in the initial condition already.*/
                /*if((i-l/2)>(j-b/2) && (i-l/2)>0)
                 {
                 *(ca+ i*b +j)=0.80;
                 }
                 
                 else if((j-b/2)<=0)
                 {
                 *(ca+ i*b +j)=0.80;
                 
                 }
                 else
                 {
                 *(cb+ i*b +j)=0.8;
                 }
                 */
            }
        }
    }
    
}

/*function for creating a well mixed rectangle in initial condiitons of length r in the y direction (extends fully in the x direction)*/
void initial_conditions_line_wm(double *ca, double *cb, double *c)
{
    for (int i = 0; i <l ;i++)
    {
        for (int j =0;j <b;j++)
        {
            *(ca+ i*b +j)=0.0;
            *(cb+ i*b +j)=0.0;
            *(c+ i*b +j)=0.0;
        }
    }
    for (int i = 0 ; i <l;i++)
    {
        for (int j =0;j < r;j++)
        {
            *(cb+ i*b +j)=h;
            
            *(ca+ i*b +j)=1.0-h;

        }
    }
    
}
/*function for creating a ectangle in initial conditoons of length r in the y direction with two unequally sized domains */
void initial_conditions_line_nwm(double *ca, double *cb, double *c)
{
    for (int i = 0; i <l ;i++)
    {
        for (int j =0;j <b;j++)
        {
            *(ca+ i*b +j)=0.0;
            *(cb+ i*b +j)=0.0;
            *(c+ i*b +j)=0.0;
        }
    }
    /* module for  creating linear initial conditions*/
    
    for (int i = 0 ; i <l;i++)
    {
        for (int j =0;j < r;j++)
        {
            if(i<l*h)
            {
                *(ca+ i*b +j)=1.0;
            }
            
            else
            {
                *(cb+ i*b +j)=1.0;
            }
        }
    }
    
}


/*function to ouput concentrations as the output of the simualtions. "c.txt" is the concetration of a+b at every site. "ca.txt" is the concentration of a at every site. concnetration of b is simply c-ca and hence not outputted.  */
void print_c_ca(double *ca, double *c, int t, int printflag)
{
    char ca_name[20], c_name[20],time_name[20]; // The filename buffer.
    // Puts "file" then k then ".txt" in to filename.
    snprintf(ca_name, sizeof(char) * 20, "ca%i.txt", printflag); // concentration of a
    snprintf(c_name, sizeof(char) * 20, "c%i.txt", printflag); // concentration of a+b (total concentration)
    snprintf(time_name, sizeof(char) * 20, "time%i.txt", printflag); // time at which outputted.
    
    FILE *fpca = fopen(ca_name, "w");
    if (fpca == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<b-1;j++)
        {
            fprintf(fpca,"%f ", *(ca+ i*b +j));
            
        }
        if(i!=l-1)
        {
            fprintf(fpca,"%f\n",  *(ca+i*b+b-1));
        }
        else if(i==l-1)
        {
            fprintf(fpca,"%f",  *(ca+(l-1)*b+b-1));
        }
    }
    fclose(fpca);
    
    FILE *fpc = fopen(c_name, "w");
    if (fpc == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<b-1;j++)
        {
            fprintf(fpc,"%f ", *(c+ i*b +j));
            
        }
        if(i!=l-1)
        {
            fprintf(fpc,"%f\n",  *(c+i*b+b-1));
        }
        else if(i==l-1)
        {
            fprintf(fpc,"%f",  *(c+(l-1)*b+b-1));
        }
    }
    fclose(fpc);
    
    FILE *ft = fopen(time_name, "w");
    if (ft == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    fprintf(ft,"%d ", t);
    fclose(ft);
    
}



int main(int argc, char * const argv[])
{


    int z;
    double p1tot;
    int wmflag;
    /* the following while loop accepts the parameters that are passed to the code as options from the command line for the case of oppositely chiral strains.*/
    /*IMPORTANT NOTE!: for different figures, I had different versions of the code where the parameters were passed in different combinations.
     For eg: in some codes, chirality of strain 1 was passed and chirality of strain2 was atuomatically the opposite. In others the chirality of strain2 was 0.
     I will include all the different versions of the while loop from line 505 to 570 as commented text at the end of this code for your convenience. */
    
    while ((z = getopt (argc, argv, "l:b:r:N:T:g:p:q:e:n:w:s:F:")) != -1)
    {
        if (z == 'l')
            l = atoi(optarg);   // lattice length along x
        else if (z == 'b')
            b = atoi(optarg);   // lattice length along y, the expansion direction
        else if (z == 'r')
            r = atoi(optarg);   // length of initial condition along y
        else if (z == 'N')
            Number = atoi(optarg);  // carrying capacity or inverse noise strength
        else if (z == 'T')
            T = atoi(optarg);   // Number of time steps
        else if (z == 'g')
            g = atof(optarg);   // growth rate of both strains
        else if (z == 'p')
            p00 = atof(optarg); // diffusion constant, usually zero
        else if (z == 'q')
            p0 = atof(optarg);  // a_s (density dependence on source deme, usually zero)
        else if (z == 'e')
            p1 = atof(optarg); // a_d (density dependence on destination deme, usually zero)
        else if (z == 'n')
        {
            p22 = atof(optarg); // al(chirality of the second strain)/2 total chirality of the strain is p22-p42
            p42=p1-p22; // // ar (chirality of the second strain)/2 total chirality of the strain is p22-p42
            
            p4=p22; // al for strain 1
            p2=p1-p4; // ar for strain 2
        }
        else if (z == 'w')
            {
                p1=0;
                p3=0;  // a_b (density dependence on back eme, usually zero)
                p12=0; // a_b, a_d for the second strain.
                p32=0;
                h = atof(optarg); // fraction of a strain in the initial conditions. (ratio of strains will be h:1-h. Hence h=0.5 means equal
            }
        /* only used if we needed to compete strains witha fitness advantage.
        else if (z == 's')
        {
            double sel=atof(optarg);
            
        }
        */
        else if (z == 'F')
        {
            wmflag=atoi(optarg); // flag variable to set the initial conditions to either well-mixed or two domains with ratio h:1-h.
        }
        
    }
    
    /* in this code, growth rate of both species are the same. If you want to have a selective adavantage, then one can pass that parameter  "s"   define gb=ga*(1+sel)*/
    gb=g;
    ga=g;
    
    /* for printing images periodically, the following 3 lines are defined, if n-images is set to zero, no intermediate state in the simulation is printed to save on storage space. */
    int printflag=0;
    int n_images=6;
    int beg_image_cutoff= lround( b * printflag*1.0/n_images  );
    
    printf("p1 %f, p2 %f,p3 %f, p4 %f, \n p12 %f, p22 %f,p32 %f,p42 %f," , p1,p2,p3, p4,p12, p22, p32, p42);
    printf("\n Hello, World!\n");
    printf("\n wmflag is %d!\n",wmflag);
    
    /*declaring and initializeing  total concnetration =c, concnetration  of a =ca , concnetration of b =cb*/
    double *ca, *cb, *c;
    initialise(&cb,l,b);
    initialise(&ca,l,b);
    initialise(&c,l,b);
    
    
    /* these are variables used in updating the concentrations in the migration function, intialised here because it can sometimes speed up simulations */
     double *patemp;
     double *a1;
     double *a2;
     double *a3;
     double *a4;
    initialise(&patemp,l,b);
    initialise(&a1,l,b);
    initialise(&a2,l,b);
    initialise(&a3,l,b);
    initialise(&a4,l,b);
    
    
    /*set inital conditions.*/
    if (wmflag==1)
    initial_conditions_line_wm(ca, cb, c);
    else
    initial_conditions_line_nwm(ca, cb, c);
    
    
    gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long randSeed;
    srand(rdtsc());                   /* initialization for rand() */
    randSeed = rand();                    /* returns a non-negative integer */
    gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */

    /*setting totoal concentration c =ca+cb*/
    for (int i=0;i <l;i++)
    {
        for (int j=0;j <b;j++)
        {
            *(c+ i*b +j)=*(ca+ i*b +j)+ *(cb+ i*b +j);
        }
    }
    
    double *het;
    initialise(&het,T,1);
    for (int i=0;i <T;i++)
    {
        het[i]=0.0;
    }

    
    
    int tbreak;
    
    /* definition of behinging and end of simulation box for the 0th times step*/
    int beg =0;
    int end=b;
    /*to time the simulation in seconds*/
    clock_t begin, end_time,end_loop;
    double time_spent,time_spent_loop;
    begin = clock();
    
    printf("loop beginning");
    for (int t=0; t<=T;t++)
    {
        
        tbreak=t;
        double sumc=0.0;
        double flength=l;// a float value =l
        
        
        /*  the following few for loops till approximately line 720 define the procedure by which we deine the begining and end of the box as the simulations is adavancing.*/
        /*the simulation box extends from  j="beg" within the bulk to j="end" in front of the furthest occupied site. it saves simualtion time because we dont have to update every site, but just the sites at the front. this is okay to do as our form of migration and growth functions are defined to go to zero in the colony bulk where cells are usually starved of nutrients.  */
        for (int j=beg;j <end;j++)
        {
            sumc=0; // sumc is the total concentration summed along the x direction.
            for (int i=0;i <l;i++)
            {
                sumc=sumc + *(c+ i*b +j);
            }
            
            if (sumc/ flength > 0.9) // if sumc/l is below 0.9, we say the front has been reached. flength is a double variable=l.
            {
                /*the beginnning of our simulation box starts even earlier in j to get further into the bulk where all sites are at carrying capacity=1.
                 This is done by stepping back in j by a wide margin of 25 sites + an additional term when we have more diffusion in simulations.
                 This is because as we add more p00 or diffusion constant, the simulation becomes more "pulled" and has longer fronts. Almost all simulations in the paper had p00=p0=0 though. We chose 25 from empirical testing of our simulations. we varied box size to see when the box size is insufficient and 25 was plenty */
                beg = j-25 - (int)5*(sqrt((p00+p0)/(g))) ;
                if (beg<0)  // to not have negative numbers for beg
                {beg=0;}
            }
            
            if (sumc/ flength < 0.1) // cf sumc/l is below 0.1, we say the front is starting to end.
            {
                /*we quickly set the end of the simulation box  25 sites further than this point + a similar factor as in end that increases this marging when wave is more pulled and has larger D . We chose 25 from empirical testing of our simulations. we varied box size to see when the box size is insufficient and 25 was plenty. */
                end = j+25 +(int) 5*(sqrt((p00+p0)/(g)));
                if (end>b) //so that end isnt after the end of the lattice
                { end =b;}
                break;
            }
        }
        /* because chirality creates bulges, for long times on wide lattices, the end can move even further from the bulk. Therefore we check if our guess for the end is correct by examining the sites in front of end and seeing if there are any sites with small positive concentrations.  */
        
        for (int j=end;j <b-1;j++)  // to check if the guess for end is okay
        {
            int flagend;
            for (int i=0;i <l;i++)
            {
                flagend=0;
                if (*(c+ i*b +j) >0.001)
                {
                    end=j+10;
                    if (end>b) //so that end isnt after the end of the lattice
                    { end =b;}
                    flagend=1;
                    printf("\n end shifted to %d \t time %d", end, t);
                    break;
                }
            }
            if (flagend==0)
            break; // there was no shift of end in the last run!
        }
        /* because dips grow logatithmically in time, our guess for big was always fine and never needed to be shifted back. Therefore we don't need to check for sites in the bulk not at the carrying capacity yet that were left behind. We did run tests to ensure this. the rests used the following block of code.
         Moreover, one could argue that if sites were left behind, they would have been straved of nutrients and wouldnt grow anyway. */
        
        /*
        for (int j=0;j < beg;j++)  // to check if guess for beg is okay
        {
            int flagend;
            for (int i=0;i <l;i++)
            {
                flagend=0;
                if (*(c+ i*b +j) <0.999)
                {
                    beg=j-10;
                    if (beg<0)
                    {
                        beg=0;
                    }
                    flagend=1;
                    printf("\n beg shifted to %d \t time %d", beg, t);
                    break;
                }
            }
            if (flagend==1)
            break; // there was no shift of end in the last run!
        }
         */
        

        /* to periodically print out images of the colony at various time points: */
        if (  beg>= beg_image_cutoff )
        {
            print_c_ca( ca, c, t, printflag);
            printf("\n printing concentrations %d, time is  %d \n",printflag,t);
            printflag++;
            beg_image_cutoff= lround( b * printflag*1.0/n_images  );
        }
        
        
        
        
        /*time evolution performs deterministic migration and growth for all sites within the box. Initiallty it updates ca, then cb. It is only after ca and cb have been updated is c updated. this ensures simultaneity and doesnt give one strain an advantage because of the update order.*/
        time_evolution(c,ca,p1,p2,p3,p4,ga,beg,end,patemp,a1,a2,a3,a4); // params are p1-p4!
        //now ca has changed but c hasnt yet and so can be passed for cb
        time_evolution(c,cb,p12,p22,p32,p42,gb,beg,end,patemp,a1,a2,a3,a4); // params are p12-p42!
        

        
        // total concentration, c changes only now( below)
        for (int i=0;i <l;i++)
        {
            for (int j=beg;j <end;j++)
            {
                *(c+ i*b +j)=*(ca+ i*b +j)+ *(cb+ i*b +j);
            }
        }
        
        /* Simulations end if either time steps reach T or if the simulation reaches close to the end of the lattice*/
        if(*(c+ l/2*b +b-40)>=0.1)
        {
            printf("\n boundary reached, time to break. %d",t);  // USUALLY simulations end due to this, I usually give T large enough for a full expansion.
            break;
        }
        // Demographic noise from binomial sampling is applied to all sites in the box.
        selection_noise(c,ca,cb,gBaseRand,beg,end);
        
        
        /* the heterozygosity of the front which quantifies the degree of mixing is calculated only within the simulation box. It should not be calculated throughout the colony because then eventaully  heterozygsoity of the system will not change in time-- it will become equal to the heterozygosity of the bulk because the bulk is frozen. Hence we calculate it within the box like below:*/
        for (int i=0;i <l;i++)
        {
            for (int j=beg;j <end;j++)
            {
                if (*(c+ i*b +j)>0.0)
                {het[t]=het[t]+2* (*(ca+ i*b +j)/ *(c+ i*b +j) )* (*(cb+ i*b +j)/ *(c+ i*b +j)); }// H=  2f_a f_b =2 f_a (1-f_a)
            }
        }
        het[t]=het[t]/((end-beg)*l); //  to normalize the heterozygosity measurement
       
        if (t%50==0)        // simply print time
            printf("%d," ,t);
        
    }
    printf("loop ended");
    end_loop = clock();
    time_spent_loop = (double)(end_loop - begin) / CLOCKS_PER_SEC;
    printf("\n loop time in secs or whatever %f \n",time_spent_loop);
    
   
    
    
    /*print final concentrations, time of exit, heterozygosity time series, computational run time, etc.*/
    
    //print ca.
    FILE *fpca = fopen("ca.txt", "w");
    if (fpca == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<b-1;j++)
        {
            fprintf(fpca,"%f ", *(ca+ i*b +j));
            
        }
        if(i!=l-1)
        {
            fprintf(fpca,"%f\n",  *(ca+i*b+b-1));
        }
        else if(i==l-1)
        {
            fprintf(fpca,"%f",  *(ca+(l-1)*b+b-1));
        }
    }
    fclose(fpca);
    // print c
    FILE *fpc = fopen("c.txt", "w");
    if (fpc == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<b-1;j++)
        {
            fprintf(fpc,"%f ", *(c+ i*b +j));
            
        }
        if(i!=l-1)
        {
            fprintf(fpc,"%f\n",  *(c+i*b+b-1));
        }
        else if(i==l-1)
        {
            fprintf(fpc,"%f",  *(c+(l-1)*b+b-1));
        }
    }
    fclose(fpc);
    

    // time at which simulation ended ( by running out of space or running for as long we had asked for
    FILE *ft = fopen("time.txt", "w");
    if (ft == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    fprintf(ft,"%d ", tbreak);
    fclose(ft);
    
    
    // time series of heterozygosity
    FILE *fphet = fopen("heterozygosity_time.txt", "w");
    if (fphet == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    for(int i=0;i<tbreak-2;i++)
    {
        fprintf(fphet,"%f ", het[i]);
    }
    fprintf(fphet,"%f",  het[tbreak-1]);
    fclose(fphet);
    
    
    end_time = clock();
    time_spent = (double)(end_time- begin) / CLOCKS_PER_SEC;
    printf("\n loop time %f \n",time_spent);
    
    
    // how long simulation took to run in seconds
    FILE *ft1 = fopen("time_in_seconds.txt", "w");
    if (ft1 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    
    fprintf(ft1,"%f ", time_spent);
    
    fclose(ft1);
    
    /*free up pointers*/
    free(c);
    free(ca);
    free(cb);
    free(patemp);
    free(a1);
    free(a2);
    free(a3);
    free(a4);

    printf("\n");
    printf("\n Hello, World!\n");
    return 0;
}





/* different versions of the while loop for differnt ways of passing parameters are given below and explained */


/******** when both strains were oppositely chiral. If one strain had a selective advantage "sel" was given a nonzero value. ***********/

/*
while ((z = getopt (argc, argv, "l:b:r:N:T:g:p:q:e:n:w:s:")) != -1)
 {
    if (z == 'l')
        l = atoi(optarg);
        else if (z == 'b')
            b = atoi(optarg);
            else if (z == 'r')
                r = atoi(optarg);
                else if (z == 'N')
                    Number = atoi(optarg);
                    else if (z == 'T')
                        T = atoi(optarg);
                        else if (z == 'g')
                            g = atof(optarg);
                            else if (z == 'p')
                                p00 = atof(optarg);
                                else if (z == 'q')
                                    p0 = atof(optarg);
                                    else if (z == 'e')
                                        p1 = atof(optarg);
                                        else if (z == 'n')
                                        { p22 = atof(optarg);
                                            p42=p1-p22;
                                            
                                            p4=p22;
                                            p2=p1-p4;
                                        }
                                        else if (z == 'w')
                                        {
                                            p1=0;
                                            p3=0;
                                            p12=0;
                                            p32=0;
                                            
                                            h = atof(optarg);
                                        }
                                        else if (z == 's') // selective advantage of A over B i.e.,g_a=g*(1+sel)
                                                sel = atof(optarg);
 
                                            
 }
*/

/******** when both strains had the same chirality ***********/
/*
while ((z = getopt (argc, argv, "l:b:r:N:T:g:p:q:e:n:w:s:")) != -1)
 {
    if (z == 'l')
        l = atoi(optarg);
        else if (z == 'b')
            b = atoi(optarg);
            else if (z == 'r')
                r = atoi(optarg);
                else if (z == 'N')
                    Number = atoi(optarg);
                    else if (z == 'T')
                        T = atoi(optarg);
                        else if (z == 'g')
                            g = atof(optarg);
                            else if (z == 'p')
                                p00 = atof(optarg);
                                else if (z == 'q')
                                    p0 = atof(optarg);
                                    else if (z == 'e')
                                        p1 = atof(optarg);
                                        else if (z == 'n')
                                        { p22 = atof(optarg);
                                            p42=p1-p22;
                                            
                                            p4=p42;
                                            p2=p22;
                                        }
                                        else if (z == 'w')
                                        {
                                            p1=0;
                                            p3=0;
                                            p12=0;
                                            p32=0;
                                            
                                            h = atof(optarg);
                                        }
 }
 */



/******** for runs that varied the relative chiralities "f*" we provided the input in a slightly messy way explained here-
 p1 is passed the value of al+ar which is the same for both strains.
 We also set the sum of the absolute value of  chiralities of the two strains to al+ar =p1
 We passed "ch_sum" which was equal to   1+sum of the of the chiralities of the two strains ( sum of chiralities includes signs defined as  negative for left and positive for right )
 Given the 3 constraints in the  lines above  and the fact that strain 1 is defined to be more lefthanded than strain 2, we can assign unique chiralities to both the strains.
 ***********/

/*
while ((z = getopt (argc, argv, "l:b:r:N:T:g:p:q:e:n:w:s:F:")) != -1)
    {
    if (z == 'l')
        l = atoi(optarg);
        else if (z == 'b')
            b = atoi(optarg);
            else if (z == 'r')
                r = atoi(optarg);
                else if (z == 'N')
                    Number = atoi(optarg);
                    else if (z == 'T')
                        T = atoi(optarg);
                        else if (z == 'g')
                        {
                            g = atof(optarg);
                        }
                        else if (z == 'p')
                            p00 = atof(optarg);
                            else if (z == 'q')
                                p0 = atof(optarg);
                                else if (z == 'e')
                                    p1 = atof(optarg); // temporarily stores a_l +a_r and is equal to sum of the magnitude of the chiralities
                                    else if (z == 'n')
                                    {
                                        double ch_sum = (atof(optarg)-1)*p1;
                            // sum of the chiralities+1 is input, chsum of the two species goes from -p1t to +p1 (one strain max left handed to other strain max righthanded.
                                        
                                        p2= (2*p1+ ch_sum + p1)/4;
                                        
                                        p4=p1-p2;
                                        
                                        p22=p2-p1/2;
                                        
                                        p42=p1-p22;
                                    }
                                    else if (z == 'w')
                                    {
                                        p1=0;
                                        p3=0;
                                        p12=0;
                                        p32=0;
                                        
                                        h = atof(optarg);
                                    }
    
                                    else if (z == 's') // nonchiral bit
                                    {
                                        double p1temp=atof(optarg);
                                        p22=p22-p1temp/2.0;
                                        p2=p2-p1temp/2.0;
                                        p42=p42-p1temp/2.0;
                                        p4=p4-p1temp/2.0;
                                    }
    
                                    else if (z == 'F')
                                    {
                                        wmflag=atoi(optarg);
                                    }
    
    }
*/



/******** for competing a chiral strain with a nonchiral strain, with or without a selective advantage. ***********/
/*
while ((z = getopt (argc, argv, "l:b:r:N:T:g:p:q:e:n:w:s:")) != -1)
    {
    if (z == 'l')
        l = atoi(optarg);
        else if (z == 'b')
            b = atoi(optarg);
            else if (z == 'r')
                r = atoi(optarg);
                else if (z == 'N')
                    Number = atoi(optarg);
                    else if (z == 'T')
                        T = atoi(optarg);
                        else if (z == 'g')
                            g = atof(optarg);
                            else if (z == 'p')
                                p00 = atof(optarg);
                                else if (z == 'q')
                                    p0 = atof(optarg);
                                    else if (z == 'e')
                                        p1 = atof(optarg);
                                        else if (z == 'n')
                                        { p22 = atof(optarg);
                                            p42=p1-p22;
                                            
                                            p4=p1/2;
                                            p2=p1/2;        // ca is nonchiral
                                        }
                                        else if (z == 'w')
                                        {
                                            p1=0;
                                            p3=0;
                                            p12=0;
                                            p32=0;
                                            
                                            h = atof(optarg);
                                        }
    
                                        else if (z == 's') // selective advantage of A over B i.e.,g_a=g*(1+sel)
                                            sel = atof(optarg);
                                            
   }

*/



/********
 when the system is interpreted as a stochastic PDE, we can run the system on a lattice with lattice spacing and time step set by "a" and "dt".
 Although not used in the paper, this is useful if you want to rescale the lattice and run simualtions.
 To do this, we pass to the simulation continuous parameters like D ,L,B etc in real units. We also specify our discretization error threshold of epsilon and choice of lattice spacing (smaller epsilon means finer discretization. Following this, the code calculates the dimensionless discrete parameters that correspond to the dimensionful continuous parameters it received after choosing a lattice spacing 'a' and time step 'dt'.
 
 After providing a suitable 'a'. the code chooses the time step as dt=epsilon*a*a/D_total; where D_total is the total migration rates (density dependent and independent).
 This sets the total migration porbability (density dependent and independent) at each time step to epsilon.
 following this,
 the discrete growth rate ,g = G *dt ( where G is the continous growth rate)
 l and b, the dimensions of the lattice  are chosen as L/a and B/a respeictively.
 T=continuous time/dt, etc.
 
 
 Note: we provide the 'a' simply for convenience (when we used it in certain simualtions). All choices of 'a' are not acceptible and it was not chosen for the simulation at random.
 
 In fact one can choose 'a' and 'dt' independently to ensure the total migration probabiltiy and growth probability per time step remain small. In other words, we have the length scale of $\sqrt{D/g}$ and a time scale of   $1/g$ . We need to choose 'a' and 't' to resolve both these scales appropriately.
 
 
 This rescaling also means that any length such as distance to the end of the simulation or size of the box also has to be rescaled by  'a'.
 
 ***********/


/*
while ((z = getopt (argc, argv, "l:b:r:N:T:g:p:q:e:n:w:k:F:a:")) != -1)
{
    if (z == 'l')
        l = atoi(optarg);
        else if (z == 'b')
            b = atoi(optarg);
            else if (z == 'r')
                r = atoi(optarg);
                else if (z == 'N')
                    Number = atoi(optarg);
                    else if (z == 'T')
                        T = atoi(optarg);
                        else if (z == 'g')
                            g = atof(optarg);
                            else if (z == 'p')
                                p00 = atof(optarg);
                                else if (z == 'q')
                                    p0 = atof(optarg);
                                    else if (z == 'e')
                                        p1 = atof(optarg);
                                        else if (z == 'n')
                                        { p22 = atof(optarg);
                                            p42=p1-p22;
                                            
                                            p4=p22;
                                            p2=p42;
                                        }
                                        else if (z == 'w')
                                        {
                                            p1=0;
                                            p3=0;
                                            p12=0;
                                            p32=0;
                                            
                                            f0 = atof(optarg); // f0
                                        }
    
                                        else if (z == 'k')
                                        {
                                            epsilon=atof(optarg);
                                        }
    
                                        else if (z == 'F')
                                        {
                                            wmflag=atoi(optarg);
                                        }
    
                                        else if (z == 'a')
                                        {
                                            
                                            a = atof(optarg); // a is the lattice spacing.
                                            double D_total=p2+p4+p00;
                                            //int a_inv=round(1/a);
                                            
                                            dt=epsilon*a*a/D_total;
                                            dx=a;
                                            g=g*dt;
                                            if (g>1.0)
                                            {printf ("g %f, SCREAM AND EXIT!",g);
                                                exit(1);
                                            }
                                            
                                            l=round(l/a); // Number of lattice points has to increase
                                            b=round(b/a);
                                            r=round(r/a);  // length of initial condition in lattice points has to increase
                                            
                                            T=round(T/dt);
        // How N scales can be derived from the dimensions of the noise . The noise is delta correlated in time and space.
                                            Number =round(Number *0.1/epsilon);  // N should be rescale as N/ epsilon
                                            // we just put 0.1 there so that N doesnt change if epsilon =0.1
                                            // to compare easily with previous simulations.
                                            
                                            // since m should be fixed to epsilon,
                                            p2=p2*epsilon/D_total;
                                            p22=p22*epsilon/D_total;
                                            p42=p42*epsilon/D_total;
                                            p4=p4*epsilon/D_total;
                                            p00=p00*epsilon/D_total;
                                            
                                            // g, p22, etc are the migration and grwoth rates in the discrete model. no more dt dx^2 in the update step!!
                                            
                                        }
    }
*/
