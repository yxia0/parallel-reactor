#ifndef PARALLEL_INCLUDE
#define PARALLEL_INCLUDE

#include "simulation_configuration.h"
#include "simulation_support.h"

typedef struct channel_struct **Reactor;

void modifyLocalNeutronSize(struct simulation_configuration_struct *, int, int);
void setHelperValuesInConfig(struct simulation_configuration_struct *, double);
/* Make a copy of the reactor core */
void initialiseReactorCopy(struct simulation_configuration_struct *, struct channel_struct **, struct channel_struct ***);
void initialiseChemicalDelta(struct simulation_configuration_struct *, double **, double **);
void setFuelAssemblyIndex(int **, struct simulation_configuration_struct *, int);
/* Reset chemical change in fuel assembly to zero */
void resetChemicalDelta(double *, struct simulation_configuration_struct *);
/* Calculate the change in fuel assembly */
void calculateChemicalDelta(struct channel_struct **, struct channel_struct **, double *, int *, struct simulation_configuration_struct *);
/* Apply chamicals change to reactor core */
void applyDeltaToReactor(double *, struct channel_struct **, int *, struct simulation_configuration_struct *); // need delta and the copy of the reactor
/* Copy the value of fuel assembly to its correspondence in the reactor core */
void copyFuelAssembly(struct channel_struct **, struct channel_struct **, struct simulation_configuration_struct *); // need copy of the reactor and the reactor
void printReactorCopy(struct channel_struct **, struct simulation_configuration_struct *);
void printChemicalDelta(double *, struct simulation_configuration_struct *);

#endif