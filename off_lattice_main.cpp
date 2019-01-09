#include <iostream>
#include "off_lattice_population.hpp"
#include <gsl/gsl_rng.h>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>

gsl_rng *gBaseRand;
// GLOBAL RANDOM NUMBER DECLARATION AND SEEDING
/* specifying to use Mersenne twister MT-19937 as the uniform PRNG */
       /* global rand number generator */

unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}



int main(int argc, char * argv[])
{
    


    
    std::cout << "Hello World " << std::endl;
    
    int max_population_size=300;
    
    double box_length=10.; // in units of spread in concentration sensing
    double domain_width=10.;
    double spread_in_concentration_sensing=1.;
    double migration_rate=1.0;
    double migration_step_size=0.1;
    double chirality_A=0.0;
    double chirality_B=0.0;
    int BC_flag=0;
    double  per_deme_carrying_capacity=1.;
    int sensing_flag=0;
    double growth_rate=0.1;
    double f0=0.5;
    double IC_length=2.5;
    int wmflag=1;
    int IC_population_size=100;
    int Time_steps=1;
    double Allee_coefficient=-3.;
    int z;
    int read_from_file_flag=0;
    
    int migration_func_flag=0;
    //int N, double W, double L, double Spread, double MS, double m, double G, double CHa, double Chb, double K, int BCf, int senflag
    //int N_ic, double f , double L_ic,int wmflag,const gsl_rng *rsd
    while ((z = getopt (argc, argv, "N:W:L:v:z:m:g:a:b:K:B:C:n:f:l:w:T:A:M:R:")) != -1)
    {
        if (z == 'N')
            max_population_size = atoi(optarg);
        else if (z == 'W')
            domain_width = atof(optarg);
        else if (z == 'L')
            box_length = atof(optarg);
        else if (z == 'v')
            spread_in_concentration_sensing = atof(optarg);
        else if (z == 'z')
            migration_step_size = atof(optarg);
        else if (z == 'm')
            migration_rate = atof(optarg);
        else if (z == 'g')
            growth_rate = atof(optarg);
        else if (z == 'a')
            chirality_A = atof(optarg);
        else if (z == 'b')
            chirality_B = atof(optarg);
        else if (z == 'K')
            per_deme_carrying_capacity = atof(optarg);
        else if (z == 'B')
            BC_flag = atoi(optarg);
        else if (z == 'C')
            sensing_flag = atoi(optarg);
        else if (z == 'n')
            IC_population_size = atoi(optarg);
        else if (z == 'f')
            f0 = atof(optarg);
        else if (z == 'l')
            IC_length = atof(optarg);
        else if (z == 'w')
            wmflag = atoi(optarg);
        else if (z == 'T')
            Time_steps = atoi(optarg);
        else if (z == 'A')
            Allee_coefficient = atof(optarg);
        else if (z == 'M')
           migration_func_flag = atoi(optarg);
        else if (z == 'R')
            read_from_file_flag = atoi(optarg);
        

    }
    
    
    
    printf("\n read_from_file_flag is %d \n",read_from_file_flag);
    
    
    //off_lattice_population :: off_lattice_population(int N,  double W, double L, double Spread, double MS, double G, double CHa, double Chb, double K, int BCf,int senflag)
    off_lattice_population population1(max_population_size, domain_width, box_length, spread_in_concentration_sensing, migration_step_size, migration_rate, growth_rate, chirality_A, chirality_B, per_deme_carrying_capacity, BC_flag, sensing_flag);
    //void off_lattice_population :: initialise_population(int N_ic, double f , double L_ic,int wmflag,const gsl_rng *rsd)

    unsigned long randSeed;
    srand(rdtsc());                   /* initialization for rand() */
    randSeed = rand();                    /* returns a non-negative integer */
    gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gBaseRand, randSeed);


    // initialise_population(int N_ic, double f , double L_ic,int wmf,const gsl_rng *rsd)
    if (read_from_file_flag)
    {
        printf("\n about to read from file \n");
        char data_file_name[32];
        char param_file_name[32];
        char label_file_name[32];
        char time_file_name[32];
        char sensed_file_name[32];
        int suffix=10;
        sprintf   (data_file_name, "pos_%d.txt", suffix);
        sprintf   (param_file_name, "params_%d.txt",suffix);
        sprintf   (label_file_name, "label_%d.txt", suffix);
        sprintf   (time_file_name, "time_%d.txt", suffix);
        sprintf   (sensed_file_name, "sensedC_%d.txt", suffix);
       // population1.initialise_population(IC_population_size, f0 , IC_length,wmflag,gBaseRand);
        
        population1.initialise_from_textfile(data_file_name,param_file_name,label_file_name,time_file_name,sensed_file_name );
        printf("\n read from file \n");
        exit(1);
    }
    else
        population1.initialise_population(IC_population_size, f0 , IC_length,wmflag,gBaseRand);
    
    
    switch (migration_func_flag)
    {
        case 0:
            printf("\n population1.migrate(gBaseRand) \n"); break;
        case 1:
            printf("population1.migrate_chiral_advection(gBaseRand) \n"); break;
        case 2:
           printf(" population1.migrate_chiral_destination(gBaseRand)\n"); break;
        case 3:
           printf(" population1.migrate_chiral_kernel(gBaseRand)\n"); break;
        case 4:
            printf("population1.migrate_chiral_kernel_destination(gBaseRand)\n"); break;
        case 5:
            printf("population1.migrate_chiral_kernel_destination_allee(gBaseRand)\n"); break;
        case 6:
            printf("migrate_chiral_kernel_acceptance_destination(gBaseRand)\n"); break;
        default:
            printf("\n migration_func_flag was inccorect  %d \n",migration_func_flag);exit(1);break;
    }
    
    

    
    std::vector<double> ymean_t(Time_steps/10+1);
    std::vector<int> popsize_t(Time_steps/10+1);
    std::vector<int> tval_t(Time_steps/10 +1) ;
    
    int number_of_data_points=10; // total number including 0 is +1 of this number
    int time_factor=Time_steps/number_of_data_points ;
    
    std::vector<double> time_in_secs(number_of_data_points +1);
    clock_t begin, current_time,end_loop;
    double time_spent,time_spent_loop;
    begin = clock();
    
    int box_moving_interval=30;
    
    for (int t=0; t<=Time_steps;t++)
    {
        population1.sense_concentrations();
        population1.grow_logistically(gBaseRand,t,Allee_coefficient);
        switch (migration_func_flag)
            {
            case 0:
                population1.migrate(gBaseRand); break;
            case 1:
                population1.migrate_chiral_advection(gBaseRand); break;
            case 2:
                population1.migrate_chiral_destination(gBaseRand); break;
            case 3:
                population1.migrate_chiral_kernel(gBaseRand); break;
            case 4:
                population1.migrate_chiral_kernel_destination(gBaseRand); break;
            case 5:
                population1.migrate_chiral_kernel_destination_allee(gBaseRand); break;
            case 6:
                    population1.migrate_chiral_kernel_acceptance_destination(gBaseRand); break;
            default:
                 printf("\n migration_func_flag was inccorect  %d \n",migration_func_flag);exit(1);break;
            }
        
        population1.update_population_size();
        /*
         updates current population size to include new grown cells AFTER the migration event.
         this prevents new grown cells from moving and hence makes the migration and growth updates more "synchronous"
         */
        
        

        if (t%10==0)
        {
            ymean_t[t/10]=population1.mean_ypos();
            popsize_t[t/10]=population1.current_population_size;
            tval_t[t/10]=t;
        }

        if (t%box_moving_interval==0 && box_length>1) // if box length is zero, the box is not moved
        {
           population1.move_box();
        }
        
        
        if (t%time_factor==0)
        {
            
            char data_file_name[32];
            char param_file_name[32];
            char label_file_name[32];
            char time_file_name[32];
            char sensed_file_name[32];
            sprintf   (data_file_name, "pos_%d.txt", t/time_factor);
            sprintf   (param_file_name, "params_%d.txt", t/time_factor);
            sprintf   (label_file_name, "label_%d.txt", t/time_factor);
            sprintf   (time_file_name, "time_%d.txt", t/time_factor);
            sprintf   (sensed_file_name, "sensedC_%d.txt", t/time_factor);
            population1.sense_concentrations();
            population1.save_data(data_file_name,param_file_name,label_file_name,time_file_name,sensed_file_name);
            current_time = clock();
            time_spent = (double)(current_time - begin) / CLOCKS_PER_SEC;
           time_in_secs[t/time_factor]=time_spent;

        }
 
    }
    end_loop = clock();
    time_spent_loop = (double)(end_loop - begin) / CLOCKS_PER_SEC;
    printf("\n loop time in secs  %f \n",time_spent_loop);
    
    
    
    switch (migration_func_flag)
    {
        case 0:
            printf("\n population1.migrate(gBaseRand) \n"); break;
        case 1:
            printf("population1.migrate_chiral_advection(gBaseRand) \n"); break;
        case 2:
            printf(" population1.migrate_chiral_destination(gBaseRand)\n"); break;
        case 3:
            printf(" population1.migrate_chiral_kernel(gBaseRand)\n"); break;
        case 4:
            printf("population1.migrate_chiral_kernel_destination(gBaseRand)\n"); break;
        case 5:
            printf("population1.migrate_chiral_kernel_destination_allee(gBaseRand)\n"); break;
        case 6:
            printf("migrate_chiral_kernel_acceptance_destination(gBaseRand)\n"); break;
        default:
            printf("\n migration_func_flag was inccorect  %d \n",migration_func_flag);exit(1);break;
    }
    
    

    FILE *fp_time = fopen("time_in_secs", "w");
    if (fp_time == NULL)
    {
        std::cout << "Error opening file!\n";
        exit(1);
    }
    for(int j=0;j<number_of_data_points +1;j++)
        fprintf(fp_time,"%f ",time_in_secs[j] );
    fclose(fp_time);

    

    
    FILE *fp_popt = fopen("popsize_t.txt", "w");
    if (fp_popt == NULL)
    {
        std::cout << "Error opening file!\n";
        exit(1);
    }
    for(unsigned int j=0;j<popsize_t.size();j++)
        fprintf(fp_popt,"%d ",popsize_t[j] );
    
    fprintf(fp_popt,"\n" );
    for(unsigned int j=0;j<popsize_t.size();j++)
        fprintf(fp_popt,"%f ",ymean_t[j] );
    fprintf(fp_popt,"\n" );
    for(unsigned int j=0;j<popsize_t.size();j++)
        fprintf(fp_popt,"%d ",tval_t[j] );
    fclose(fp_popt);
    

    
    return 0;
    
}


