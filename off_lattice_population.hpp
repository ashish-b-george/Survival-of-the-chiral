// defines the off lattice population class
#ifndef OFF_LATTICE_POPULATION_HPP
#define OFF_LATTICE_POPULATION_HPP
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstring>
#include <fstream>
#include <sstream>

class off_lattice_population
{
    public:
    static const double PI = 3.141592653589793;
    int wmflag, BC_flag,sensing_flag;
    int max_population_size,current_population_size,IC_population_size,simulation_time,bulk_population,frozen_population,population_grown;
    // bulk population is number of species in the bulk that is not simulated or sensed
    //non-growing population is the number of species in the front with y coordinate= 1sigma+bulk population tip that the growing population senses but doesnt grow because we dont sense bulk population
    double domain_width, box_length,IC_length,growth_rate,migration_step_size,spread_in_concentration_sensing,chirality_A,chirality_B,per_deme_carrying_capacity,f0,migration_rate,allee_B;

    
    std::vector<int> species_label;
    std::vector<int> species_label_swap;
    std::vector<double> xpos;
    std::vector<double> ypos;
    std::vector<double> xpos_swap;
    std::vector<double> ypos_swap;
    std::vector<double> chirality_of_label;
    std::vector<double> xgrad_concentrations;
    std::vector<double> ygrad_concentrations;
    std::vector<double> sensed_concentrations;
    
    
    off_lattice_population(int N, double W, double L, double Spread, double MS, double m, double G, double CHa, double Chb, double K, int BCf, int senflag);
    
    void initialise_population(int N_ic, double f , double L_ic,int wmflag,const gsl_rng *rsd);
    void grow_logistically(const gsl_rng *rsd, int t, double allee_flag=-2.0);
    
    void migrate(const gsl_rng *rsd);
    void migrate_chiral_advection(const gsl_rng *rsd);
    void migrate_chiral_destination(const gsl_rng *rsd);
    void migrate_chiral_kernel(const gsl_rng *rsd);
    void migrate_chiral_kernel_destination(const gsl_rng *rsd);
    void migrate_chiral_kernel_destination_allee(const gsl_rng *rsd);
    void migrate_chiral_kernel_acceptance_destination(const gsl_rng *rsd);
    
    void move_box();
    void apply_BC();
    void RBC();
    void PBC();
    void sense_concentrations(); // only PBC has been done accurately. RBC does not sense concentrations beyond the boundary
    void save_data(char *filename_pos,char *filename_params ,char *filename_label,char *filename_time,char *filename_sensed_concentrations);
    double mean_ypos();
    double sense_concentrations_at_xy(double x, double y);
    void update_population_size();
    
    void initialise_from_textfile(char *filename_pos,char *filename_params,char *filename_label,char *filename_time,char *filename_sensed_concentrations );
    void read_params(char *filename_label);
    void read_positions_and_labels(char *filename_pos,char *filename_label);

    
    
    ////////      to be defined functions           //////
    void record_profile();
    // some kind of copy constructor to copy data
    // moving box and the constructor to implement moving box.
    
};


off_lattice_population :: off_lattice_population(int N,  double W, double L, double Spread, double MS, double m, double G, double CHa, double CHb, double K, int BCf,int senflag)
{
    
    max_population_size=N;
    
    box_length=L;
    domain_width=W;
    spread_in_concentration_sensing=Spread;
    migration_step_size=MS;
    chirality_A=CHa;
    chirality_B=CHb;
    chirality_of_label.resize(2);
    chirality_of_label[0]=CHa;
    chirality_of_label[1]=CHb;
    BC_flag=BCf;
    per_deme_carrying_capacity=K;
    sensing_flag=senflag;
    growth_rate=G;
    migration_rate=m;
    simulation_time=0;
    bulk_population=0;
    frozen_population=0;
    population_grown=0;
    /*
    initialise_double_array(&position,max_population_size,2 );
    initialise_double_array(&grad_concentrations,max_population_size,2 );
    initialise_double_array(&sensed_concentrations,max_population_size,1 );
    initialise_short_int_array(&species_label,max_population_size,1);
    */
    species_label.resize(max_population_size);
    xpos.resize(max_population_size);
    ypos.resize(max_population_size);
    xgrad_concentrations.resize(max_population_size);
    ygrad_concentrations.resize(max_population_size);
    sensed_concentrations.resize(max_population_size);
    

    f0=0;
    IC_length=0;
    wmflag=-1;
    current_population_size=0;
    IC_population_size=0;
    
    std::cout << "constructor has been called " << std::endl;
 
}
/*
void off_lattice_population :: initialise_population_from_file()eeds filenames as input. it inputs particle posisitions from the file.
{
    // need to calculate and set the following quantities after inputting particle posistions and species labels.:
  
     current_population_size=N_ic;
     IC_population_size=N_ic;
     f0=f;
     IC_length=L_ic;
     wmflag=wmf;

}
*/
void off_lattice_population :: initialise_population(int N_ic, double f , double L_ic,int wmf,const gsl_rng *rsd)
//void off_lattice_population :: initialise_population(int N_ic, double f , double L_ic,int wmf)
{
    
    std::cout << "initialising, N_ic is " << N_ic<<std::endl;
    current_population_size=N_ic;
    IC_population_size=N_ic;
    population_grown=0;
    f0=f;
    IC_length=L_ic;
    wmflag=wmf;
    for (int i = 0; i < current_population_size; i++)
    {
        if (i<current_population_size*f0)
            species_label[i]=1;
        else
            species_label[i]=0;
    }
    // the randomg number generator is seeded using srand in main, so not seeded here!
    // otherwise the generator would seeded every time a new instance of the class was generated.

    if (wmflag==2) // four doamins with sizes l/4,l/8,l/4+l/8 and l/4 of unequal size.  species also are reinitialised to these fractions!
    {
        for (int i = 0; i < current_population_size; i++)
        {
            if (i<current_population_size/4.)
                species_label[i]=1;
            else if (i<current_population_size/4.+current_population_size/8.)
                species_label[i]=0;
            else if (i<3*current_population_size/4.)
                species_label[i]=1;
            else
                species_label[i]=0;
            
            xpos[i]=  domain_width*i*1.0/current_population_size;
            ypos[i]=IC_length*gsl_rng_uniform(rsd);
        }
        
    }

    else if (wmflag==1)// well mixed initial conditions with random poisitions for the particles
    {
       for (int i = 0; i < current_population_size; i++)
        {
            xpos[i]=domain_width*gsl_rng_uniform(rsd);
            ypos[i]=IC_length*gsl_rng_uniform(rsd);
        }
        
    }

    else if (wmflag==0) //nwm IC with two domains, domain size automatically becomes in the ratio of the fractions.
    {
        for (int i = 0; i < current_population_size; i++)
        {
            xpos[i]=  domain_width*i*1.0/current_population_size;
            ypos[i]=IC_length*gsl_rng_uniform(rsd);
        }
        
    }
    else
        std::cout << "wmflag is not provided correctly" << std::endl;

    
    std::cout << "initialised" << std::endl;
    
}


void off_lattice_population ::grow_logistically(const gsl_rng *rsd, int t, double allee_flag)
{
    double gi; // growth of particle i
    //int population_grown; is now a class variable.
    simulation_time=t;
    population_grown=0;
    if (allee_flag <=0.0) // no allee effect
        allee_B=0.0;
    else     // there is an allee effect
        allee_B=allee_flag;
    for (int i = frozen_population; i < current_population_size; i++)
    {
        if  (allee_B<=0.0)
            gi=growth_rate*(1.0-sensed_concentrations[i]);
        else
            gi=growth_rate*(1.0-sensed_concentrations[i])*(1.+allee_B*sensed_concentrations[i]);
        if(gi>0.) // gi is positive, so checking for growth
        {
            if (gi>gsl_rng_uniform_pos(rsd) )
            {
                xpos[current_population_size+population_grown]=xpos[i];
                ypos[current_population_size+population_grown]=ypos[i];
                species_label[current_population_size+population_grown]=species_label[i];
                // the concentrations of the new cell should be updated as well!! it is updated to the concentration felt by the old cell!!
                sensed_concentrations[current_population_size+population_grown]=sensed_concentrations[i];
                ++population_grown;
                if (current_population_size+population_grown>= 0.9*max_population_size)
                    break;
            }
        }
        else
        {
            //std::cout << "growth rate is negative, dying !" <<gi<< std::endl;
            if ( gi<(-1.*gsl_rng_uniform_pos(rsd)) ) // the particle dies.
            {
                std::iter_swap(xpos.begin() + i, xpos.begin() +current_population_size+population_grown-1 ); // move to the end of the array and reduce size to kill a particle
                std::iter_swap(ypos.begin() + i, ypos.begin() +current_population_size+population_grown-1 );
                std::iter_swap(species_label.begin() + i, species_label.begin() +current_population_size+population_grown-1 );
                // the concentrations should also be swapped right!!? when the cell is swapped into place, living cell keeps it's old concentration!!!
                std::iter_swap(sensed_concentrations.begin() + i, sensed_concentrations.begin() +current_population_size+population_grown-1 );
                --population_grown;
            }
        }
    }
    //std::cout << "grown from " <<current_population_size<<"  to  "<<current_population_size+population_grown<<"at time="<<simulation_time<< std::endl;
    
    // for truly simulatenous growth and migration the new cells formed should not migrate hence the current population size should not be updated till the next time we sense concentrations
    //current_population_size+=population_grown;
    if (current_population_size>= 0.98*max_population_size) // exit.
    {
        std::cout << "population dangerously close to carrying capacity" << std::endl;
        char data_file_name[32];
        char param_file_name[32];
        char label_file_name[32];
        char time_file_name[32];
        char sensed_file_name[32];
        sprintf   (data_file_name, "pos_final.txt");
        sprintf   (param_file_name, "params_final.txt");
        sprintf   (label_file_name, "label_final.txt");
        sprintf   (time_file_name, "time_final.txt");
        sprintf   (sensed_file_name, "sensedC_final.txt");
        save_data(data_file_name,param_file_name,label_file_name,time_file_name,sensed_file_name);
        exit(10);
    }
}



void off_lattice_population ::migrate(const gsl_rng *rsd)
{
    double theta=0.0;
    if (migration_rate>0.999) // if migration rate is set to 1., migration happens every time.
    {
        for (int i = frozen_population; i < current_population_size; i++)
        {
            theta=gsl_rng_uniform(rsd);
            xpos[i]+= migration_step_size* cos(2.*PI*theta);
            ypos[i]+= migration_step_size* sin(2.*PI*theta);
        }
    }
    else
    {
        for (int i = frozen_population; i < current_population_size; i++)
        {
            if (migration_rate>gsl_rng_uniform_pos(rsd) )
            {
                theta=gsl_rng_uniform(rsd);
                xpos[i]+= migration_step_size* cos(2.*PI*theta);
                ypos[i]+= migration_step_size* sin(2.*PI*theta);
            }
        }
    }
    apply_BC();
}


void off_lattice_population ::migrate_chiral_advection(const gsl_rng *rsd)
{
    double theta=0.0;
    double c=0.;
    if (migration_rate>0.999) // if migration rate is set to 1., migration happens every time.
    {
        for (int i = frozen_population; i < current_population_size; i++)
        {
            c=sensed_concentrations[i];
            /*******
             Diffusion is restricted in the bulk by changing the migration probability. 
             Changing the step size is tricky because (1-c) can be negative and so that doesn't make sense.
            *****/
             if (migration_rate*(1-c)>gsl_rng_uniform_pos(rsd) )
            {
            // diffusion...
            theta=gsl_rng_uniform(rsd);
            xpos[i]+= migration_step_size *sqrt(c) * cos(2.*PI*theta);
            ypos[i]+= migration_step_size *sqrt(c) * sin(2.*PI*theta);
            
            //deterministic chiral advection
            
            xpos[i]+= migration_step_size * ygrad_concentrations[i] *chirality_of_label[species_label[i]];
            ypos[i]-= migration_step_size * xgrad_concentrations[i] *chirality_of_label[species_label[i]];
            /*theta=atan2(ygrad_concentrations[i],xgrad_concentrations[i]); // theta already has the 2pi in it!
            mod_gradc=
            xpos[i]+=migration_step_size*sqrt(c)*cos(theta)  ;
            ypos[i]+=migration_step_size*sqrt(c)*sin(theta);
             */
            }
        }
    }
    else
    {
        for (int i = frozen_population; i < current_population_size; i++)
        {
            c=sensed_concentrations[i];
            /*******
             Diffusion is restricted in the bulk by changing the migration probability.
             Changing the step size is tricky because (1-c) can be negative and so that doesn't make sense.
             *****/
            if (migration_rate*(1-c)>gsl_rng_uniform_pos(rsd) )
            {
                
                // diffusion...
                theta=gsl_rng_uniform(rsd);
                xpos[i]+= migration_step_size*sqrt(c) *cos(2.*PI*theta);
                ypos[i]+= migration_step_size*sqrt(c) * sin(2.*PI*theta);
                
                //deterministic chiral advection
                // this might cause chiral species to move more than a nonchiral species...... :((

                xpos[i]+= migration_step_size*sqrt(c) *ygrad_concentrations[i] *chirality_of_label[species_label[i]];
                ypos[i]-= migration_step_size*sqrt(c) *xgrad_concentrations[i] *chirality_of_label[species_label[i]];
                
            }
        }
    }
    apply_BC();
}


void off_lattice_population ::migrate_chiral_destination(const gsl_rng *rsd)
// this function migrates with probbility (1-c_d), i.e., migration event is rejected with probability c_d, chiral migration is deterministic
{
    double theta=0.0;
    double c=0.;
    double x_new,y_new,c_d;

    for (int i = frozen_population; i < current_population_size; i++)
    {
        c=sensed_concentrations[i];
        /*******
        Diffusion is restricted in the bulk by changing the migration probability.
        Changing the step size is tricky because (1-c) can be negative and so that doesn't make sense.
        *****/
        // diffusion...
        theta=gsl_rng_uniform(rsd);
        x_new= xpos[i]+ migration_step_size *sqrt(c) * cos(2.*PI*theta);
        y_new= ypos[i]+ migration_step_size *sqrt(c) * sin(2.*PI*theta);
        
        //deterministic chiral advection
        x_new+= migration_step_size * ygrad_concentrations[i] *chirality_of_label[species_label[i]];
        y_new-= migration_step_size * xgrad_concentrations[i] *chirality_of_label[species_label[i]];
        c_d= sense_concentrations_at_xy(x_new,y_new);
        if ( migration_rate*( 1.-c_d )>gsl_rng_uniform_pos(rsd) ) // the migration is accepted. otherwise position is not updated.
        {
            xpos[i]=x_new;
            ypos[i]=y_new;
        }
        
    }
    
    apply_BC();
}


void off_lattice_population ::migrate_chiral_kernel(const gsl_rng *rsd)
{
    //gsl_ran_gaussian_ziggurat(const gsl_rng * r, double sigma)
   

    double sigma_chiral_kernel,gradient_magnitude,theta,c;
    for (int i = frozen_population; i < current_population_size; i++)
    {
        c=sensed_concentrations[i];
        /*******
         Migration occurs from a gaussian that narrows depending on  microscopic chirality with maxima along the gradient
         *****/
        if (migration_rate*(1-c)>gsl_rng_uniform_pos(rsd) )
        {
            if (chirality_of_label[species_label[i]] !=0.0 )
            {
                gradient_magnitude=sqrt( pow(xgrad_concentrations[i],2)+ pow(ygrad_concentrations[i],2)  );
                 sigma_chiral_kernel= PI/( gradient_magnitude* fabs(chirality_of_label[species_label[i]]) );
                theta=gsl_ran_gaussian_ziggurat(rsd, sigma_chiral_kernel) ; // theta need not be in 0,2pi but cos and sin will take care of that
                if (chirality_of_label[species_label[i]] >0.0 )
                    theta += atan2 (-xgrad_concentrations[i], ygrad_concentrations[i]); //atan2 (y,x), but we use (-xgrad,ygrad) because we want mean to be perpendicular
                else
                   theta += atan2 (xgrad_concentrations[i], -ygrad_concentrations[i]); //atan2 (y,x), but we use (xgrad,-ygrad) because we want mean to be perpendicular and opposite the positive chirality one.
            }
            else
            {
                theta=gsl_rng_uniform(rsd)*2.*PI;
            }
        
            xpos[i]+= migration_step_size *sqrt(c) * cos(theta);
            ypos[i]+= migration_step_size *sqrt(c) * sin(theta);
        }
    }
    apply_BC();
}


void off_lattice_population ::migrate_chiral_kernel_destination(const gsl_rng *rsd)
{
    //gsl_ran_gaussian_ziggurat(const gsl_rng * r, double sigma)
   

    double x_new,y_new,c_d;
    double sigma_chiral_kernel,gradient_magnitude,theta,c;
    for (int i = frozen_population; i < current_population_size; i++)
    {
        c=sensed_concentrations[i];
        /*******
         Migration occurs from a gaussian that narrows depending on  microscopic chirality with maxima along the gradient
         *****/
        if (chirality_of_label[species_label[i]] !=0.0 )
        {
            gradient_magnitude=sqrt( pow(xgrad_concentrations[i],2)+ pow(ygrad_concentrations[i],2)  );
            sigma_chiral_kernel=PI/( gradient_magnitude* fabs(chirality_of_label[species_label[i]]) );
            
            // Sign of chirality changes sign of theta, angle is not made proportional to chirality currently.,
            // potentially we could make it such that the deviation is continous from gradient to perpendicular or something
            theta=gsl_ran_gaussian_ziggurat(rsd, sigma_chiral_kernel) ; // theta need not be in 0,2pi but cos and sin will take care of that
            if (chirality_of_label[species_label[i]] >0.0 )
                theta += atan2 (-xgrad_concentrations[i], ygrad_concentrations[i]); //atan2 (y,x), but we use (-xgrad,ygrad) because we want mean to be perpendicular
            else
                theta += atan2 (xgrad_concentrations[i], -ygrad_concentrations[i]); //atan2 (y,x), but we use (xgrad,-ygrad) because we want mean to be perpendicular and opposite the positive chirality one.
        }
        else
        {
            theta=gsl_rng_uniform(rsd)*2.*PI;
        }
        
        x_new= xpos[i]+ migration_step_size *sqrt(c) * cos(theta);
        y_new= ypos[i]+ migration_step_size *sqrt(c) * sin(theta);
        
        c_d= sense_concentrations_at_xy(x_new,y_new);
        if ( migration_rate*( 1.-c_d )>gsl_rng_uniform(rsd) ) // the migration is accepted. otherwise position is not updated.
        {
            xpos[i]=x_new;
            ypos[i]=y_new;
        }
        
    }
    apply_BC();
}

void off_lattice_population ::migrate_chiral_kernel_destination_allee(const gsl_rng *rsd) // growth has allee so migration is constant step size.
{
    double x_new,y_new,c_d;
    double sigma_chiral_kernel,gradient_magnitude,theta,c;
    for (int i = frozen_population; i < current_population_size; i++)
    {
        c=sensed_concentrations[i];
        /*******
         Migration occurs from a gaussian that narrows depending on  microscopic chirality with maxima along the gradient
         *****/
        if (chirality_of_label[species_label[i]] !=0.0 )
        {
            gradient_magnitude=sqrt( pow(xgrad_concentrations[i],2)+ pow(ygrad_concentrations[i],2)  );
            sigma_chiral_kernel= PI/( gradient_magnitude* fabs(chirality_of_label[species_label[i]]) );
            
            // Sign of chirality changes sign of theta, angle is not made proportional to chirality currently.,
            // potentially we could make it such that the deviation is continous from gradient to perpendicular or something
            theta=gsl_ran_gaussian_ziggurat(rsd, sigma_chiral_kernel) ; // theta need not be in 0,2pi but cos and sin will take care of that
            if (chirality_of_label[species_label[i]] >0.0 )
                theta += atan2 (-xgrad_concentrations[i], ygrad_concentrations[i]); //atan2 (y,x), but we use (-xgrad,ygrad) because we want mean to be perpendicular
            else
                theta += atan2 (xgrad_concentrations[i], -ygrad_concentrations[i]); //atan2 (y,x), but we use (xgrad,-ygrad) because we want mean to be perpendicular and opposite the positive chirality one.
        }
        else
        {
            theta=gsl_rng_uniform(rsd)*2.*PI;
        }
        
        x_new= xpos[i]+ migration_step_size  * cos(theta);
        y_new= ypos[i]+ migration_step_size  * sin(theta);
        
        c_d= sense_concentrations_at_xy(x_new,y_new);
        if ( migration_rate*( 1.-c_d )>gsl_rng_uniform(rsd) ) // the migration is accepted. otherwise position is not updated.
        {
            xpos[i]=x_new;
            ypos[i]=y_new;
        }
        
    }
    apply_BC();
}



  /*******
   To get rid of sqrt(c)scaling in the chiral advection when migrating with a kernel,
   this function migrates with probability c , and the migration event is rejected with probability (1-c_d) , chiral migration is from kernel  but step size is NOT sqrt(c)!!!
   *****/
void off_lattice_population ::migrate_chiral_kernel_acceptance_destination(const gsl_rng *rsd)
{
  double x_new,y_new,c_d;
  double sigma_chiral_kernel,gradient_magnitude,theta,c;
  for (int i = frozen_population; i < current_population_size; i++)
  {
      c=sensed_concentrations[i];
      /*******
       Migration occurs from a gaussian that narrows depending on  microscopic chirality with maxima along the gradient
       *****/
      if ( migration_rate*c > gsl_rng_uniform(rsd) ) // a migration event happens with probability proportional to the concentration.
      {
          if (chirality_of_label[species_label[i]] !=0.0 )
          {
              gradient_magnitude=sqrt( pow(xgrad_concentrations[i],2)+ pow(ygrad_concentrations[i],2)  );
              sigma_chiral_kernel=PI/( gradient_magnitude* fabs(chirality_of_label[species_label[i]]) );
              
              // Sign of chirality changes sign of theta, angle is not made proportional to chirality currently.,
              // potentially we could make it such that the deviation is continous from gradient to perpendicular or something
              theta=gsl_ran_gaussian_ziggurat(rsd, sigma_chiral_kernel) ; // theta need not be in 0,2pi but cos and sin will take care of that
              if (chirality_of_label[species_label[i]] >0.0 )
                  theta += atan2 (-xgrad_concentrations[i], ygrad_concentrations[i]); //atan2 (y,x), but we use (-xgrad,ygrad) because we want mean to be perpendicular
              else
                  theta += atan2 (xgrad_concentrations[i], -ygrad_concentrations[i]); //atan2 (y,x), but we use (xgrad,-ygrad) because we want mean to be perpendicular and opposite the positive chirality one.
          }
          else
          {
                theta=gsl_rng_uniform(rsd)*2.*PI;
          }

          x_new= xpos[i]+ migration_step_size  * cos(theta);
          y_new= ypos[i]+ migration_step_size  * sin(theta);

          c_d= sense_concentrations_at_xy(x_new,y_new);
          if ( migration_rate*( 1.-c_d )>gsl_rng_uniform(rsd) ) // the migration is accepted. otherwise position is not updated.
          {
              xpos[i]=x_new;
              ypos[i]=y_new;
          }
      }
  }
apply_BC();
    
}

void off_lattice_population ::apply_BC()
{
    if (BC_flag==0)
        PBC();
    /*
    else if (BC_flag==1)
        RBC();
    */
    else
        std::cout << "BCflag not given correctly as(1/0), it is:"<<BC_flag << std::endl;
    
}

void off_lattice_population ::PBC()
{
    for (int i = frozen_population; i < current_population_size; i++)
    {
    
        if (xpos[i]<0)
            xpos[i]+=domain_width;
        else if (xpos[i]>domain_width)
            xpos[i]-=domain_width;
        
        if (ypos[i]<0)      //sort of an RBC in the y direction, not important
            ypos[i]=-ypos[i];
        // we dont' really need to implemnt a BC at the other edge right? so not implementing one
   
    }
}

double off_lattice_population ::mean_ypos()
{
    double sum=0.0;
    for (int i = 0; i < current_population_size; i++)
        sum+=ypos[i];
    
    return sum/current_population_size;

}


void off_lattice_population ::move_box()
{
    double ymax=0.;
    ///////         find max y location , the tip of the front /////////
    for (int i = bulk_population; i < current_population_size; i++)
    {
        if (ypos[i]>ymax)
            ymax=ypos[i];
    }
    double y_frozen,y_bulk;
    y_bulk=std::max(0.0,ymax-box_length*spread_in_concentration_sensing);
    if (sensing_flag==1)
        y_frozen=std::max(0.0,ymax-box_length*spread_in_concentration_sensing+2.*spread_in_concentration_sensing); // since its a sharp function, we dont need to store much further
    else
        y_frozen=std::max(0.0,ymax-box_length*spread_in_concentration_sensing+4.*spread_in_concentration_sensing); // smoother function, need to store further back
    
    int delta_bulk_population, delta_frozen_population;
    delta_bulk_population=0;
    delta_frozen_population=0;
    
    if (y_frozen>0.0)
    {
        ///////        find the additional population in each category  ///////////
        for (int i = bulk_population; i < current_population_size; i++) // max y location
        {
            if (ypos[i]<=y_frozen)
            {
                delta_frozen_population++; // change in frozen population includes the bulk population change
                if(ypos[i]<=y_bulk)
                    delta_bulk_population++;
            }
        }
        ///////        move the frozen and bulk populations into their positions in the swap array  ///////////
        int box_growing_ctr=0;
        int box_frozen_ctr=0;
        int bulk_ctr=0;
        if ( xpos_swap.size()< (unsigned)lround( 1.2*(current_population_size-bulk_population) ) ) /// if the swap array is not big enough, 1.2 just so we dont resize often
        {
            xpos_swap.resize(lround( 1.2*(current_population_size-bulk_population) ));
            ypos_swap.resize(lround( 1.2*(current_population_size-bulk_population) ));
            species_label_swap.resize(lround( 1.2*(current_population_size-bulk_population) ));
        }
        for (int i = bulk_population; i < current_population_size; i++)
        {
            if(ypos[i]>y_frozen) // it belongs in the box, in the growing population
            {
                ypos_swap[delta_frozen_population+box_growing_ctr]=ypos[i];
                xpos_swap[delta_frozen_population+box_growing_ctr]=xpos[i];
                species_label_swap[delta_frozen_population+box_growing_ctr]=species_label[i];
                ++box_growing_ctr;
            }
            else if(ypos[i]>y_bulk) // it belongs in the box, but in the frozen population
            {
                ypos_swap[delta_bulk_population+box_frozen_ctr]=ypos[i];
                xpos_swap[delta_bulk_population+box_frozen_ctr]=xpos[i];
                species_label_swap[delta_bulk_population+box_frozen_ctr]=species_label[i];
                ++box_frozen_ctr;
            }
            else // it belongs outside the box, in the bulk population
            {
                ypos_swap[bulk_ctr]=ypos[i];
                xpos_swap[bulk_ctr]=xpos[i];
                species_label_swap[bulk_ctr]=species_label[i];
                ++bulk_ctr;
            }
        }
        /////////// copy the swap array into our array so that elements are in order    ///////
        for (int i = bulk_population; i < current_population_size; i++)
        {
            ypos[i]=ypos_swap[i-bulk_population];
            xpos[i]=xpos_swap[i-bulk_population];
            species_label[i]=species_label_swap[i-bulk_population];
        }
    }
    std::cout << "bulk population changes from " <<bulk_population<<"  to  "<<bulk_population+delta_bulk_population<< std::endl;
    std::cout << "frozen population changes from " <<frozen_population<<"  to  "<<bulk_population+delta_frozen_population<< std::endl;
    std::cout << "current population is" <<current_population_size<< std::endl;
    frozen_population=bulk_population+delta_frozen_population;
    bulk_population+=delta_bulk_population;
    
}

double off_lattice_population ::sense_concentrations_at_xy(double x, double y)
{
    double sum, concentration_xy,xtemp;
    concentration_xy=0.;
    sum=0.;
    double sigma_squared=(pow(spread_in_concentration_sensing,2));
     if (BC_flag!=0) // notPBC
         std::cout << "function sense_concentrations_at_xy not defined for not PBC " << std::endl;
    
    
    if (sensing_flag==0) // gaussian kernel
    {
        for (int j = bulk_population; j < current_population_size; j++)
        {
            xtemp=std::min( pow( x-xpos[j],2), std::min( pow( x-xpos[j]+domain_width,2),pow( x-xpos[j]-domain_width,2) )  );
            sum+=exp( -( xtemp+ pow( y-ypos[j],2) )/sigma_squared );
        }
    }
    
    else // kernel just counts particles
    {
        for (int j = bulk_population; j < current_population_size; j++)
        {
            xtemp=std::min( pow( x-xpos[j],2), std::min( pow( x-xpos[j]+domain_width,2),pow( x-xpos[j]-domain_width,2) )  );
            if ( xtemp+ pow( y-ypos[j],2) < pow(spread_in_concentration_sensing,2) )
                    sum+=1;
        }
    }
    
    sum = sum/per_deme_carrying_capacity;
    return sum;
}

void off_lattice_population ::update_population_size()
{
    current_population_size+=population_grown;
    population_grown=0;
}


void off_lattice_population ::sense_concentrations()
{
    // now the population size is updated by the number o cells grown

    double sum,gradx_sum,grady_sum,temp; //gradx  and grad y gives you -gradc!
    double sigma_squared=(pow(spread_in_concentration_sensing,2));
    if (sensing_flag==0) // gaussian kernel
    {
        if (BC_flag==0) // PBC
        {
            double xtemp;
            for (int i = frozen_population; i < current_population_size; i++)
            {
                sum=0;
                gradx_sum=0;
                grady_sum=0;
                for (int j = bulk_population; j < current_population_size; j++)
                {
                // can not minimise x^2 as we need min(abs(x)) for gradx computation
                //xtemp=std::min( pow( xpos[i]-xpos[j],2), std::min( pow( xpos[i]-xpos[j]+domain_width,2),pow( xpos[i]-xpos[j]-domain_width,2) )  );
                    
                xtemp=xpos[i]-xpos[j];
                if ( fabs(xtemp)>=domain_width/2. ) //xtemp is not the smallest, PBC needs to be employed
                {
                    xtemp=xpos[i]-xpos[j]+domain_width;
                    if ( fabs(xtemp)>=domain_width/2. )
                        xtemp=xpos[i]-xpos[j]-domain_width;
                }
                
                temp=exp( -( pow(xtemp,2)+ pow( ypos[i]-ypos[j],2) )/sigma_squared );
                sum+= temp;
                gradx_sum+=2.*xtemp*temp/sigma_squared;
                grady_sum+=2.*(ypos[i]-ypos[j])*temp/sigma_squared;
                    
                    
                }
                sensed_concentrations[i]=sum/per_deme_carrying_capacity;
                xgrad_concentrations[i]=gradx_sum/per_deme_carrying_capacity;
                ygrad_concentrations[i]=grady_sum/per_deme_carrying_capacity;
            }
        }
        else
        {
            for (int i = frozen_population; i < current_population_size; i++)
            {
                sum=0;
                gradx_sum=0;
                grady_sum=0;
                for (int j = bulk_population; j < current_population_size; j++)
                {
                    temp=exp( -( pow( xpos[i]-xpos[j],2)+ pow( ypos[i]-ypos[j],2) )/(pow(spread_in_concentration_sensing,2)) ) ;
                    sum+=temp;
                    gradx_sum+=2.*(xpos[i]-xpos[j])*temp/sigma_squared;
                    grady_sum+=2.*(ypos[i]-ypos[j])*temp/sigma_squared;
                    
                    
                }
                sensed_concentrations[i]=sum/per_deme_carrying_capacity;
                xgrad_concentrations[i]=gradx_sum/per_deme_carrying_capacity;
                ygrad_concentrations[i]=grady_sum/per_deme_carrying_capacity;
            }
        }
    }
    else if (sensing_flag==1) // counts number of neighbours withing a circle of radius = spread_in_concentration_sensing, no gradient sensing possible
    {
        if (BC_flag==0) // PBC
        {
            double xtemp;
            for (int i = frozen_population; i < current_population_size; i++)
            {
                sum=0;
                for (int j = bulk_population; j < current_population_size; j++)
                {
                    
                    xtemp=std::min( pow( xpos[i]-xpos[j],2), std::min( pow( xpos[i]-xpos[j]+domain_width,2),pow( xpos[i]-xpos[j]-domain_width,2) )  );
                    
                    if ( xtemp+ pow( ypos[i]-ypos[j],2) < pow(spread_in_concentration_sensing,2) )
                        sum+=1;
                    
                    
                }
                sensed_concentrations[i]=sum/per_deme_carrying_capacity;
            }
        }
    }
}



void off_lattice_population ::save_data(char *filename_pos,char *filename_params,char *filename_label,char *filename_time,char *filename_sensed_concentrations )
{
    
    ///////////        saving data file    ///////////////
    FILE *fp_data = fopen(filename_pos, "w");
    if (fp_data == NULL)
    {
        std::cout << "Error opening file!\n";
        exit(1);
    }
    for(int j=0;j<current_population_size;j++)
        fprintf(fp_data,"%f ",xpos[j] );
            
    fprintf(fp_data,"\n" );
    for(int j=0;j<current_population_size;j++)
        fprintf(fp_data,"%f ",ypos[j] );
    fclose(fp_data);
    
    
    //if ( strcmp(filename_sensed_concentrations, "NULL")!=0 ) // if fiule name is not "NULL", we print a file
    //{
        FILE *fp_sensed = fopen(filename_sensed_concentrations, "w");
        if (fp_sensed == NULL)
        {
            std::cout << "Error opening file!\n";
            exit(1);
        }
        for(int j=0;j<current_population_size;j++)
            fprintf(fp_sensed,"%f ",sensed_concentrations[j] );
    
        fprintf(fp_sensed,"\n" );
        for(int j=0;j<current_population_size;j++)
            fprintf(fp_sensed,"%f ",xgrad_concentrations[j] );
        fprintf(fp_sensed,"\n" );
        for(int j=0;j<current_population_size;j++)
            fprintf(fp_sensed,"%f ",ygrad_concentrations[j] );
        fclose(fp_sensed);
    
    //}


    ///////     saving parameter file   ////////////
    FILE *fp_params = fopen(filename_params, "w");
    if (fp_params == NULL)
    {
        std::cout << "Error opening file!\n";
        exit(1);
    }
    fprintf(fp_params,"max_population_size %d\n",max_population_size);
    fprintf(fp_params,"box_length %f\n",box_length);
    fprintf(fp_params,"domain_width %f\n",domain_width);
    fprintf(fp_params,"spread_in_concentration_sensing %f\n",spread_in_concentration_sensing);
    fprintf(fp_params,"migration_step_size %f\n",migration_step_size);
    fprintf(fp_params,"migration_rate %f\n",migration_rate);
    fprintf(fp_params,"growth_rate %f\n",growth_rate);
    fprintf(fp_params,"chirality_A %f\n",chirality_A );
    fprintf(fp_params,"chirality_B %f\n",chirality_B);
    
    fprintf(fp_params,"BC_flag %d\n",BC_flag);
    fprintf(fp_params,"wmflag %d\n",wmflag);
    fprintf(fp_params,"sensing_flag %d\n",sensing_flag);
    fprintf(fp_params,"per_deme_carrying_capacity %f\n",per_deme_carrying_capacity);
    fprintf(fp_params,"current_population_size %d\n",current_population_size);
    fprintf(fp_params,"IC_population_size %d\n",IC_population_size);
    fprintf(fp_params,"f_IC %f\n",f0);
    fprintf(fp_params,"IC_length %f\n",IC_length);
    fprintf(fp_params,"allee_B %f\n",allee_B);
    fclose(fp_params);
    
    
    FILE *fp_label = fopen(filename_label, "w");
    if (fp_label == NULL)
    {
        std::cout << "Error opening file!\n";
        exit(1);
    }
    for(int j=0;j<current_population_size;j++)
        fprintf(fp_label,"%d ",species_label[j] );
    fclose(fp_label);
    
    FILE *fp_time = fopen(filename_time, "w");
    if (fp_time == NULL)
    {
        std::cout << "Error opening file!\n";
        exit(1);
    }
    fprintf(fp_time,"%d",simulation_time);
    fclose(fp_time);
   
}

void off_lattice_population :: initialise_from_textfile(char *filename_pos,char *filename_params,char *filename_label,char *filename_time,char *filename_sensed_concentrations )
{
    read_params(filename_label);
    read_positions_and_labels(filename_pos, filename_label);
    
}

void off_lattice_population :: read_params(char *filename_label)
{
    
}

void off_lattice_population :: read_positions_and_labels(char *filename_pos,char *filename_label)
{
    

    std::fstream pos_file(filename_pos);
    std::string line;
    //std::vector<std::vector<double>> xpos,ypos;
    int i = 0;
    int ctr=0;
    std::cout << " size of xpos " << xpos.size()<<"ypos is "<<ypos.size() << '\n';
        
    while (std::getline(pos_file, line))
    {
        double value;
        std::stringstream ss(line);
        
        //v.push_back(std::vector<double>());
        
        if (i==0)
        {
            while (ss >> value)
            {
                //std::cout << " " << value;
                xpos.push_back(value);
                ++ctr;
                
            }
        }
        else if (i==1)
        {
            while (ss >> value)
            {
                ypos.push_back(value);
                //std::cout << " " << value;
                ++ctr;
            }
        }
            
        printf("i is %d, ctr is %d",i,ctr);
        ++i;
    }
    
    std::cout << " size of xpos " << xpos.size()<<"ypos is "<<ypos.size() << '\n';
}



