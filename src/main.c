#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include "simulation_configuration.h"
#include "simulation_support.h"

// Height of a fuel pellet in meters (they are 40mm by 40mm by 2.5mm)
#define HEIGHT_FUEL_PELLET_M 0.0025
// Size of each reactor channel in meters (they are cuboid, so this value in x and y)
#define CHANNEL_SIZE 0.2
// Depth of fuel pellets in meters (they are x=40mm by y=40mm by z=2mm)
#define FUEL_PELLET_DEPTH 0.002
// Weight in grams of the neutron generator for each cm in height
#define NEUTRON_GENERATOR_WEIGHT_PER_CM 0.5
// Number of nanoseconds in a second
#define NS_AS_SEC 1e-9

// The neutrons that are currently moving throughout the reactor core
struct neutron_struct * neutrons;
// Indexes of empty neutrons (i.e. those inactive) that can be used
unsigned long int * neutron_index;
// The current index in the neutron_index array which is the next free neutron
unsigned long int currentNeutronIndex=0;
// The reactor core itself, each are channels in the x and y dimensions
struct channel_struct ** reactor_core;

static void step(int, struct simulation_configuration_struct*);
static void generateReport(int, int, struct simulation_configuration_struct*, struct timeval);
static void updateReactorCore(int, struct simulation_configuration_struct*);
static void updateNeutrons(int, struct simulation_configuration_struct*);
static void updateFuelAssembly(int, struct channel_struct*);
static void updateNeutronGenerator(int, struct channel_struct*, struct simulation_configuration_struct*);
static void createNeutrons(int, struct channel_struct*, double);
static void initialiseReactorCore(struct simulation_configuration_struct*);
static void initialiseNeutrons(struct simulation_configuration_struct*);
static double getControlRodLoweredToLevel(struct simulation_configuration_struct*, int, int);
static void writeReactorState(struct simulation_configuration_struct*, int, char*);
static void getFuelAssemblyChemicalContents(struct fuel_assembly_struct*, double*);
static void clearReactorStateFile(char*);
static struct channel_struct* locateChannelFromPosition(double, double, struct simulation_configuration_struct *);
static unsigned long int getTotalNumberFissions(struct simulation_configuration_struct*);
static unsigned long int getNumberActiveNeutrons(struct simulation_configuration_struct*);
static double getElapsedTime(struct timeval);

/**
* Program entry point, this code is run with configuration file as command line argument
**/
int main(int argc, char * argv[]) {
  if (argc != 3) {
    printf("You must provide two arguments, the configuration filename and output filename to the code\n");
    return -1;
  }
  time_t t;
  // Seed the random number generator
  srand((unsigned) time(&t));
  struct timeval start_time;
  // Parse the configuration and then initialise reactor core and neutrons from this
  struct simulation_configuration_struct configuration;
  parseConfiguration(argv[1], &configuration);
  initialiseReactorCore(&configuration);
  initialiseNeutrons(&configuration);
  // Empty the file we will use to store the reactor state
  clearReactorStateFile(argv[2]);
  printf("Simulation configured for reactor core of size %dm by %dm by %dm, timesteps=%d dt=%dns\n", configuration.size_x,
      configuration.size_y, configuration.size_z, configuration.num_timesteps, configuration.dt);
  printf("------------------------------------------------------------------------------------------------\n");
  gettimeofday(&start_time, NULL); // Record simulation start time (for runtime statistics)
  for (int i=0;i<configuration.num_timesteps;i++) {
    // Progress in timesteps
    step(configuration.dt, &configuration);
    if (i > 0 && i % configuration.display_progess_frequency == 0) {
      generateReport(configuration.dt, i, &configuration, start_time);
    }

    if (i > 0 && i % configuration.write_reactor_state_frequency == 0) {
      writeReactorState(&configuration, i, argv[2]);
    }
  }
  // Now we are finished write some summary information
  unsigned long int num_fissions=getTotalNumberFissions(&configuration);
  double mev=getMeVFromFissions(num_fissions);
  double joules=getJoulesFromMeV(mev);
  printf("------------------------------------------------------------------------------------------------\n");
  printf("Model completed after %d timesteps\nTotal model time: %f secs\nTotal fissions: %ld releasing %e MeV and %e Joules\nTotal runtime: %.2f seconds\n",
      configuration.num_timesteps, (NS_AS_SEC*configuration.dt)*configuration.num_timesteps, num_fissions, mev, joules, getElapsedTime(start_time));
}

/**
 * Undertake a single timestep of processing
 **/
static void step(int dt, struct simulation_configuration_struct * configuration) {
  updateNeutrons(dt, configuration);
  updateReactorCore(dt, configuration);
}

/**
 * Writes a short report around the current state of the simulation to stdout
 **/
static void generateReport(int dt, int timestep, struct simulation_configuration_struct * configuration, struct timeval start_time) {
  unsigned long int num_fissions=getTotalNumberFissions(configuration);
  double mev=getMeVFromFissions(num_fissions);
  double joules=getJoulesFromMeV(mev);
  printf("Timestep: %d, model time is %e secs, current runtime is %.2f seconds. %ld active neutrons, %ld fissions, releasing %e MeV and %e Joules\n", timestep,
      (NS_AS_SEC*dt)*timestep, getElapsedTime(start_time), getNumberActiveNeutrons(configuration), num_fissions, mev, joules);
}

/**
 * Update the state of the reactor core at the current timestep, which will update the state
 * of each fuel assembly and neutron generator.
 **/
static void updateReactorCore(int dt, struct simulation_configuration_struct * configuration) {
  for (int i=0;i<configuration->channels_x;i++) {
    for (int j=0;j<configuration->channels_y;j++) {
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
        updateFuelAssembly(dt, &(reactor_core[i][j]));
      }

      if (reactor_core[i][j].type == NEUTRON_GENERATOR) {
        updateNeutronGenerator(dt, &(reactor_core[i][j]), configuration);
      }
    }
  }
}

/**
 * Update each neutron at the current timestep, moving its position based upon the velocity
 * components and energy, and then handling the collision of the neutron with a fuel channel,
 * moderator, or control rod
 **/
static void updateNeutrons(int dt, struct simulation_configuration_struct * configuration) {
  for (long int i=0;i<configuration->max_neutrons;i++) {
    if (neutrons[i].active) {
      // Rest mass is 1 for a neutron
      double total_velocity=MeVToVelocity(neutrons[i].energy, 1);
      // These components are positive or negative which denote movement in one direction or another
      double component_velocity_x=((abs(neutrons[i].x)/100.0) * total_velocity)*NS_AS_SEC*dt;
      double component_velocity_y=((abs(neutrons[i].y)/100.0) * total_velocity)*NS_AS_SEC*dt;
      double component_velocity_z=((abs(neutrons[i].z)/100.0) * total_velocity)*NS_AS_SEC*dt;
      if (neutrons[i].x > 0) {
        neutrons[i].pos_x+=component_velocity_x;
      } else {
        neutrons[i].pos_x-=component_velocity_x;
      }
      if (neutrons[i].y > 0) {
        neutrons[i].pos_y+=component_velocity_y;
      } else {
        neutrons[i].pos_y-=component_velocity_y;
      }
      if (neutrons[i].z > 0) {
        neutrons[i].pos_z+=component_velocity_z;
      } else {
        neutrons[i].pos_z-=component_velocity_z;
      }

      if (neutrons[i].pos_x > configuration->size_x || neutrons[i].pos_x <0.0 ||
          neutrons[i].pos_y > configuration->size_y || neutrons[i].pos_y <0.0 ||
          neutrons[i].pos_z > configuration->size_z || neutrons[i].pos_z <0.0) {
            // Moved out of the reactor core, so deactivate the neutron
            neutrons[i].active=false;
            neutron_index[currentNeutronIndex]=i;
            currentNeutronIndex++;
            continue;
      }

      // Now figure out if neutron is in a fuel assembly, moderator or control rod. If so then need to handle interaction
      struct channel_struct* reactorChannel=locateChannelFromPosition(neutrons[i].pos_x, neutrons[i].pos_y, configuration);
      if (reactorChannel != NULL) {
        if (reactorChannel->type == FUEL_ASSEMBLY) {
          // It is in a fuel assembly channel, determine if it has collided with a neutron and if so deactivate it
          int fuel_pellet=(int) (neutrons[i].pos_z / HEIGHT_FUEL_PELLET_M);
          if (fuel_pellet < reactorChannel->contents.fuel_assembly.num_pellets) {
            bool collision=determineAndHandleIfNeutronFuelCollision(neutrons[i].energy, reactorChannel, fuel_pellet, configuration->collision_prob_multiplyer);
            if (collision) {
              neutrons[i].active=false;
              neutron_index[currentNeutronIndex]=i;
              currentNeutronIndex++;
            }
          }
        }

        if (reactorChannel->type == MODERATOR) {
          // The neutron is in the moderator, determine if it has been absorbed by the moderator or ot
          bool absorbed=determineAndHandleIfNeutronModeratorCollision(&neutrons[i], configuration->moderator_weight,
                                  reactorChannel->contents.moderator.type, configuration->size_z);
          if (absorbed) {
            neutrons[i].active=false;
            neutron_index[currentNeutronIndex]=i;
            currentNeutronIndex++;
          }
        }

        if (reactorChannel->type == CONTROL_ROD) {
          if (neutrons[i].pos_z <= reactorChannel->contents.control_rod.lowered_to_level) {
            // Has hit the control rod, therefore this absorbed and removed from simulation
            neutrons[i].active=false;
            neutron_index[currentNeutronIndex]=i;
            currentNeutronIndex++;
          }
        }
      } else {
        fprintf(stderr, "Unable to locate reactor core channel for x=%f and y=%f\n", neutrons[i].pos_x, neutrons[i].pos_y);
        exit(-1);
      }
    }
  }
}

/**
 * Update the state of a specific fuel assembly in a channel for a timestep. This will fission all U236
 * and Pu239 in the assembly and update the constituent components as required
 **/
static void updateFuelAssembly(int dt, struct channel_struct * channel) {
  for (int i=0;i<channel->contents.fuel_assembly.num_pellets;i++) {
    unsigned long int num_u236=(unsigned long int) channel->contents.fuel_assembly.quantities[i][U236];
    for (unsigned long int j=0;j<num_u236;j++) {
      int num_neutrons=fissionU236(channel, i);
      createNeutrons(num_neutrons, channel, (i*HEIGHT_FUEL_PELLET_M)+(HEIGHT_FUEL_PELLET_M/2));
      channel->contents.fuel_assembly.num_fissions++;
    }
    unsigned long int num_pu240=(unsigned long int) channel->contents.fuel_assembly.quantities[i][Pu240];
    for (unsigned long int j=0;j<num_pu240;j++) {
      int num_neutrons=fissionPu240(channel, i);
      createNeutrons(num_neutrons, channel, (i*HEIGHT_FUEL_PELLET_M)+(HEIGHT_FUEL_PELLET_M/2));
      channel->contents.fuel_assembly.num_fissions++;
    }
  }
}

/**
 * Update the state of a neutron generator for a timestep, generating the required
 * number of neutrons
 **/
static void updateNeutronGenerator(int dt, struct channel_struct * channel, struct simulation_configuration_struct * configuration) {
  unsigned long int number_new_neutrons=getNumberNeutronsFromGenerator(channel->contents.neutron_generator.weight, dt);
  for (int i=0;i<number_new_neutrons;i++) {
    if (currentNeutronIndex == 0) break;
    currentNeutronIndex--;
    unsigned long int index=neutron_index[currentNeutronIndex];
    initialiseNeutron(&(neutrons[index]), channel, (double) (rand() / ((double) (RAND_MAX / configuration->size_z))));
  }
}

/**
 * Creates a specific number of neutrons that originate in a specific reactor channel
 * at a specific height (z) location
 **/
static void createNeutrons(int num_neutrons, struct channel_struct * channel, double z) {
  for (int k=0;k<num_neutrons;k++) {
    if (currentNeutronIndex == 0) break;
    currentNeutronIndex--;
    unsigned long int index=neutron_index[currentNeutronIndex];
    initialiseNeutron(&(neutrons[index]), channel, z);
  }
}

/**
 * Initialises the reactor core at the start of the simulation from the configuration
 **/
static void initialiseReactorCore(struct simulation_configuration_struct * simulation_configuration) {
  reactor_core=(struct channel_struct**) malloc(sizeof(struct channel_struct*)*simulation_configuration->channels_x);
  for (int i=0;i<simulation_configuration->channels_x;i++) {
    reactor_core[i]=(struct channel_struct*) malloc(sizeof(struct channel_struct)*simulation_configuration->channels_y);
  }
  for (int i=0;i<simulation_configuration->channels_x;i++) {
    // Find the absolute x centre position of this channel
    double centre_x=((i*CHANNEL_SIZE) + (CHANNEL_SIZE/2));
    if (simulation_configuration->channel_layout_config == NULL) {
      // This channel has not been configured explicitly, hence assume it is empty
      for (int j=0;j<simulation_configuration->channels_y;j++) {
        reactor_core[i][j].type=EMPTY_CHANNEL;
        reactor_core[i][j].x_centre=centre_x;
        reactor_core[i][j].y_centre=((j*CHANNEL_SIZE) + (CHANNEL_SIZE/2));
      }
    } else {
      for (int j=0;j<simulation_configuration->num_channel_configs[i];j++) {
        // For every configuration that was provided read what that was and initialise channel as required
        reactor_core[i][j].x_centre=centre_x;
        reactor_core[i][j].y_centre=((j*CHANNEL_SIZE) + (CHANNEL_SIZE/2));
        if (simulation_configuration->channel_layout_config[i][j]==CONFIG_MODERATOR) {
          // This channel is a moderator, so set that and then the moderator type
          reactor_core[i][j].type=MODERATOR;
          if (simulation_configuration->moderator_type == WATER_MOD_TYPE_CONFIG) {
            reactor_core[i][j].contents.moderator.type=WATER;
          } else if (simulation_configuration->moderator_type == DEUTERIUM_MOD_TYPE_CONFIG) {
            reactor_core[i][j].contents.moderator.type=DEUTERIUM;
          } else if (simulation_configuration->moderator_type == GRAPHITE_MOD_TYPE_CONFIG) {
            reactor_core[i][j].contents.moderator.type=GRAPHITE;
          } else if (simulation_configuration->moderator_type == NONE_MOD_TYPE_CONFIG) {
            reactor_core[i][j].contents.moderator.type=NO_MODERATOR;
          }

        } else if (simulation_configuration->channel_layout_config[i][j]==CONFIG_FUEL_ASSEMBLY) {
          // This channel is a fuel assembly, so initialise that
          reactor_core[i][j].type=FUEL_ASSEMBLY;
          // Each fuel pellet is 40mm by 40mm by 2mm deep and weighs 1 gram
          reactor_core[i][j].contents.fuel_assembly.num_pellets=simulation_configuration->size_z / FUEL_PELLET_DEPTH;
          reactor_core[i][j].contents.fuel_assembly.num_fissions=0;
          reactor_core[i][j].contents.fuel_assembly.quantities=(double(*)[NUM_CHEMICALS]) malloc(
              sizeof(unsigned long int[NUM_CHEMICALS]) * reactor_core[i][j].contents.fuel_assembly.num_pellets);
          for (int z=0;z<reactor_core[i][j].contents.fuel_assembly.num_pellets;z++) {
            // For each pellet in the assembly set up the number of atoms present for each chemical, these
            // will change as the simulation progresses and (hopefully!) fission occurs
            for (int k=0;k<NUM_CHEMICALS;k++) {
              enum chemical_type_enum chemical=getChemicalAtIndex(k);
              if (chemical == UNKNOWN_CHEMICAL) {
                fprintf(stderr, "Unknown chemical at index '%d'\n", k);
                exit(-1);
              }
              reactor_core[i][j].contents.fuel_assembly.quantities[z][k]=getAtomsPerGram(chemical) * (simulation_configuration->fuel_makeup_percentage[k] / 100.0);
            }
          }
        } else if (simulation_configuration->channel_layout_config[i][j]==CONFIG_CONTROL_ROD) {
          // If the channel is a control rod then set this and store the absolute z location it is lowered to
          reactor_core[i][j].type=CONTROL_ROD;
          reactor_core[i][j].contents.control_rod.lowered_to_level=getControlRodLoweredToLevel(simulation_configuration, i, j);
        } else if (simulation_configuration->channel_layout_config[i][j]==CONFIG_EMPTY) {
          reactor_core[i][j].type=EMPTY_CHANNEL;
        } else if (simulation_configuration->channel_layout_config[i][j]==CONFIG_NEUTRON_GENERATOR) {
          reactor_core[i][j].type=NEUTRON_GENERATOR;
          // Half a gram per cm in height
          reactor_core[i][j].contents.neutron_generator.weight=simulation_configuration->size_z*100 * NEUTRON_GENERATOR_WEIGHT_PER_CM;
        }
      }
      for (int j=simulation_configuration->num_channel_configs[i];j<simulation_configuration->channels_y;j++) {
        // For any remaining channels that were not explicitly configured for this row, then set them as empty
        reactor_core[i][j].type=EMPTY_CHANNEL;
        reactor_core[i][j].x_centre=centre_x;
        reactor_core[i][j].y_centre=((j*CHANNEL_SIZE) + (CHANNEL_SIZE/2));
      }
    }
  }
}

/**
 * Initialises the neutron storage data structures at the start of the simulation so that all
 * neutrons are inactive. For performance we hold an array of indexes, each of which represents
 * empty slots in the neutron array which can be used. The currentNeutronIndex variable points
 * to the currnet end of the list, and when adding a neutron the value in that location of neutron_index
 * is read and currentNeutronIndex decremented. When deactivating a neutron the index of the newly freed
 * location is added to the next element of neutron_index and currentNeutronIndex is incremented
 **/
static void initialiseNeutrons(struct simulation_configuration_struct * simulation_configuration) {
  neutrons=(struct neutron_struct*) malloc(sizeof(struct neutron_struct) * simulation_configuration->max_neutrons);
  neutron_index=(unsigned long int*) malloc(sizeof(unsigned long int) * simulation_configuration->max_neutrons);
  for (int i=0;i<simulation_configuration->max_neutrons;i++) {
    neutron_index[i]=i;
    neutrons[i].active=false;
  }
  currentNeutronIndex=simulation_configuration->max_neutrons;
}

/**
 * For a control rod channel will return the absolute z height position that this has been lowered to based upon
 * the percentage insertion that was configured
 **/
static double getControlRodLoweredToLevel(struct simulation_configuration_struct * simulation_configuration, int channel_x, int channel_y) {
  int rodConfigurationIndex=findControlRodConfiguration(simulation_configuration, channel_x, channel_y);
  if (rodConfigurationIndex < 0) {
    fprintf(stderr, "Expected control rod configuration for x=%d y=%d but none can be found\n", channel_x, channel_y);
    exit(-1);
  }
  return simulation_configuration->size_z*(simulation_configuration->control_rod_configurations[rodConfigurationIndex].percentage/100.0);
}

/**
 * Writes out the current state of the reactor at this timestep to a file
 **/
static void writeReactorState(struct simulation_configuration_struct * configuration, int timestep, char * outfile) {
  unsigned long int num_fissions=getTotalNumberFissions(configuration);
  double mev=getMeVFromFissions(num_fissions);
  double joules=getJoulesFromMeV(mev);
  FILE * f=fopen(outfile, "a");
  fprintf(f, "Reactor state at time %e secs, %ld fissions releasing %e MeV and %e Joules\n", (NS_AS_SEC*configuration->dt)*timestep, num_fissions, mev, joules);
  fprintf(f, "----------------------------------------------------------------------------\n");
  for (int i=0;i<configuration->channels_x;i++) {
    for (int j=0;j<configuration->channels_y;j++) {
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
        double pc[11];
        getFuelAssemblyChemicalContents(&(reactor_core[i][j].contents.fuel_assembly), pc);
        fprintf(f, "Fuel assembly %d %d, %e U235 %e U238 %e Pu239 %e U236 %e Ba141 %e Kr92 %e Xe140 %e Sr94 %e Xe134 %e Zr103 %e Pu240\n",
        i, j, pc[0], pc[1], pc[2], pc[3], pc[4], pc[5], pc[6], pc[7], pc[8], pc[9], pc[10]);
      }
    }
  }
  fprintf(f, "===========================================================================\n");
  fclose(f);
}

/**
 * Retrieves the quantities of atoms in a fuel assembly across all the pellets for
 * each chemical that will be written out to the file
 **/
static void getFuelAssemblyChemicalContents(struct fuel_assembly_struct * fuel_rod, double * amounts) {
  for (int i=0;i<11;i++) amounts[i]=0;
  for (int i=0;i<fuel_rod->num_pellets;i++) {
    for (int j=0;j<11;j++) {
      amounts[j]+=fuel_rod->quantities[i][j];
    }
  }
}

/**
 * Clears out the file that we are going to write to for the reactor state, this is called at simulation
 * startup and it will overwrite any existing contents
 **/
static void clearReactorStateFile(char * outfile) {
  FILE * f=fopen(outfile, "w");
  fclose(f);
}

/**
 * Given an x and y position in the reactor core this will locate the channel
 * that that corresponds to
 **/
static struct channel_struct* locateChannelFromPosition(double x, double y, struct simulation_configuration_struct * configuration) {
  if (x > configuration->size_x || x <0.0) return NULL;
  if (y > configuration->size_y || y <0.0) return NULL;
  int channel_x=(int) (x/0.2);
  int channel_y=(int) (y/0.2);
  return &(reactor_core[channel_x][channel_y]);
}

/**
 * Based upon the properties of each fuel assembly will return the total number of fissions
 *  that have occured across all fuel assemblies in the simulation.
 **/
static unsigned long int getTotalNumberFissions(struct simulation_configuration_struct * configuration) {
  unsigned long int total_fissions=0;
  for (int i=0;i<configuration->channels_x;i++) {
    for (int j=0;j<configuration->channels_y;j++) {
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
        total_fissions+=reactor_core[i][j].contents.fuel_assembly.num_fissions;
      }
    }
  }
  return total_fissions;
}

/**
 * Determines the number of currently active neutrons in the simulation
 **/
static unsigned long int getNumberActiveNeutrons(struct simulation_configuration_struct * configuration) {
  unsigned long int activeNeutrons=0;
  for (unsigned long int i=0;i<configuration->max_neutrons;i++) {
    if (neutrons[i].active) activeNeutrons++;
  }
  return activeNeutrons;
}

/**
* Returns in seconds the elapsed time since the start_time argument and now
**/
static double getElapsedTime(struct timeval start_time) {
  struct timeval curr_time;
  gettimeofday(&curr_time, NULL);
  long int elapsedtime = (curr_time.tv_sec * 1000000 + curr_time.tv_usec) - (start_time.tv_sec * 1000000 + start_time.tv_usec);
  return elapsedtime / 1000000.0;
}
