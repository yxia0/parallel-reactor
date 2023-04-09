#include "simulation_parallel.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

static void modifyLocalNeutronSize(struct simulation_configuration_struct *, int, int);
static void modifyLocalNeutronGeneratorWeight(struct simulation_configuration_struct *, int);
static void setHelperValuesInConfig(struct simulation_configuration_struct *);
static int getIndexOfFuelAssemblyEntry(int, int, int, struct simulation_configuration_struct *);
static int getIndexOfFuelAssembly(int, int, int, int *);

void modifyConfigurationToParallelSetting(struct simulation_configuration_struct *config, int processSize, int myrank)
{
    modifyLocalNeutronSize(config, processSize, myrank);
    modifyLocalNeutronGeneratorWeight(config, processSize);
    setHelperValuesInConfig(config);
}

static void modifyLocalNeutronSize(struct simulation_configuration_struct *config, int processSize, int myrank)
{
    long int local_max = config->max_neutrons / processSize;
    if (local_max * processSize < config->max_neutrons)
    {
        if (myrank < config->max_neutrons - local_max * processSize)
        {
            local_max++;
        }
    }
    config->max_neutrons = local_max;
}

static void modifyLocalNeutronGeneratorWeight(struct simulation_configuration_struct *config, int processSize)
{
    config->neutron_generator_weight_per_cm /= processSize;
}

static void setHelperValuesInConfig(struct simulation_configuration_struct *config)
{
    config->num_pellet = config->size_z / config->fuel_pellet_depth;
    config->fuel_assembly_total_entries_size = NUM_CHEMICALS * config->num_pellet * config->num_fuel_assembly;
}

// TODO: define this ->  typedef struct channel_struct **Reactor;

/*
Can also return the pointer to reactorCopy <- and copy that into the global pointer defined in main()
*/
void initialiseReactorCopy(struct simulation_configuration_struct *config, struct channel_struct **reactor,
                           struct channel_struct ***reactor_copy)
{
    *reactor_copy = (struct channel_struct **)malloc(sizeof(struct channel_struct *) * config->channels_x);
    for (int i = 0; i < config->channels_x; i++)
    {
        (*reactor_copy)[i] = (struct channel_struct *)malloc(sizeof(struct channel_struct) * config->channels_y);
    }
    for (int i = 0; i < config->channels_x; i++)
    {
        if (config->channel_layout_config != NULL)
        {
            for (int j = 0; j < config->num_channel_configs[i]; j++)
            {
                // copy initial chemical amount in fuel assembly into
                if (config->channel_layout_config[i][j] == CONFIG_FUEL_ASSEMBLY)
                {
                    (*reactor_copy)[i][j].type = FUEL_ASSEMBLY;
                    (*reactor_copy)[i][j].contents.fuel_assembly.quantities = (double(*)[NUM_CHEMICALS])malloc(
                        sizeof(unsigned long int[NUM_CHEMICALS]) * reactor[i][j].contents.fuel_assembly.num_pellets);
                    // note that num_pellets are not set up in reactor copy
                    int size = NUM_CHEMICALS * reactor[i][j].contents.fuel_assembly.num_pellets;
                    // _Static_assert(sizeof(double) == sizeof(unsigned long int), "size not equal");
                    memcpy((*reactor_copy)[i][j].contents.fuel_assembly.quantities, reactor[i][j].contents.fuel_assembly.quantities, sizeof(double) * size);
                    // for (int z = 0; z < reactor[i][j].contents.fuel_assembly.num_pellets; z++)
                    // {
                    //     for (int k = 0; k < NUM_CHEMICALS; k++)
                    //     {
                    //         double chemAmount = reactor[i][j].contents.fuel_assembly.quantities[z][k];
                    //         reactor_copy[i][j].contents.fuel_assembly.quantities[z][k] = chemAmount;
                    //         printf("Value exists for channel (%d, %d) here : %f \n", i, j, reactor_copy[i][j].contents.fuel_assembly.quantities[z][k]);
                    //     }
                    // }
                }
            }
        }
    }
}

void initialiseChemicalDelta(struct simulation_configuration_struct *config, double **delta, double **buffer)
{
    *delta = (double *)malloc(sizeof(double) * config->fuel_assembly_total_entries_size);
    for (int i = 0; i < config->fuel_assembly_total_entries_size; i++)
    {
        (*delta)[i] = 0;
    }
    *buffer = (double *)malloc(sizeof(double) * config->fuel_assembly_total_entries_size);
}

void resetChemicalDelta(double *delta, struct simulation_configuration_struct *config)
{
    for (int i = 0; i < config->fuel_assembly_total_entries_size; i++)
    {
        delta[i] = 0;
    }
}

void setFuelAssemblyIndex(int **channel_index, struct simulation_configuration_struct *config, int channel_type)
{
    *channel_index = (int *)malloc(sizeof(int) * config->channels_x * config->channels_y);
    int count = 0;

    for (int i = 0; i < config->channels_x; i++)
    {
        for (int j = 0; j < config->channels_y; j++)
        {
            if (config->channel_layout_config[i][j] == channel_type)
            {
                (*channel_index)[j + i * config->channels_y] = count;
                count++;
            }
            else
            {
                (*channel_index)[j + i * config->channels_y] = 0;
            }
        }
    }
}

void printReactorCopy(struct channel_struct **reactor, struct simulation_configuration_struct *config)
{
    for (int i = 0; i < config->channels_x; i++)
    {
        for (int j = 0; j < config->channels_y; j++)
        {
            if (reactor[i][j].type == FUEL_ASSEMBLY)
            {
                for (int k = 0; k < NUM_CHEMICALS; k++)
                {
                    double sum = 0;
                    for (int z = 0; z < 500; z++)
                    {
                        sum += reactor[i][j].contents.fuel_assembly.quantities[z][k];
                    }
                    if (k == 4 | k == 5)
                    {
                        printf("channel (%d, %d) has chemical %d: %f \n", i, j, k, sum);
                    }
                }
            }
        }
    }
}

void printChemicalDelta(double *delta, struct simulation_configuration_struct *config)
{
    double sum = 0;
    for (int i = 0; i < config->fuel_assembly_total_entries_size; i++)
    {
        sum += delta[i];
    }
    printf("Sum of total delta is %f \n", sum);
}

/*
Calculate the change in chemical,
*/
void calculateChemicalDelta(struct channel_struct **reactor, struct channel_struct **reactor_copy,
                            double *delta, int *index, struct simulation_configuration_struct *config)
{
    for (int i = 0; i < config->channels_x; i++)
    {
        for (int j = 0; j < config->channels_y; j++)
        {
            if (reactor[i][j].type == FUEL_ASSEMBLY)
            {
                for (int z = 0; z < reactor[i][j].contents.fuel_assembly.num_pellets; z++)
                {
                    for (int k = 0; k < NUM_CHEMICALS; k++)
                    {

                        double prevAmount = reactor_copy[i][j].contents.fuel_assembly.quantities[z][k];

                        double currAmount = reactor[i][j].contents.fuel_assembly.quantities[z][k];

                        int channel_index = getIndexOfFuelAssembly(i, j, config->channels_y, index);
                        int pos = getIndexOfFuelAssemblyEntry(k, z, channel_index, config);
                        delta[pos] = (currAmount - prevAmount);
                    }
                }
            }
        }
    }
}

/* TODO:
 */
void applyDeltaToReactor(double *delta, struct channel_struct **reactor, int *index, struct simulation_configuration_struct *config)
{
    for (int i = 0; i < config->channels_x; i++)
    {
        for (int j = 0; j < config->channels_y; j++)
        {
            if (reactor[i][j].type == FUEL_ASSEMBLY)
            {
                for (int z = 0; z < config->num_pellet; z++)
                {
                    for (int k = 0; k < NUM_CHEMICALS; k++)
                    {
                        int channel_index = getIndexOfFuelAssembly(i, j, config->channels_y, index);
                        int pos = getIndexOfFuelAssemblyEntry(k, z, channel_index, config);
                        reactor[i][j].contents.fuel_assembly.quantities[z][k] += delta[pos];
                    }
                }
            }
        }
    }
}

/* TODO:
 */
void copyFuelAssembly(struct channel_struct **reactor, struct channel_struct **reactorCopy,
                      struct simulation_configuration_struct *config)
{
    for (int i = 0; i < config->channels_x; i++)
    {
        for (int j = 0; j < config->channels_y; j++)
        {
            if (reactor[i][j].type == FUEL_ASSEMBLY)
            {
                for (int z = 0; z < config->num_pellet; z++)
                {
                    for (int k = 0; k < NUM_CHEMICALS; k++)
                    {
                        double updatedAmount = reactorCopy[i][j].contents.fuel_assembly.quantities[z][k];
                        reactor[i][j].contents.fuel_assembly.quantities[z][k] = updatedAmount;
                    }
                }
            }
        }
    }
}

static int getIndexOfFuelAssemblyEntry(int chemical, int pellet, int channel_index, struct simulation_configuration_struct *config)
{
    return chemical * config->num_pellet * config->num_fuel_assembly + channel_index * config->num_pellet + pellet;
}

/*
Return the index of a fuel assembly channel given the (x, y) position of the channel
*/
static int getIndexOfFuelAssembly(int x, int y, int num_channel_y, int *index)
{
    return index[y + x * num_channel_y];
}
