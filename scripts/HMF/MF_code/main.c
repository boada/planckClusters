#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "header.h"

int main(int argc, char **argv)
{
  if(argc==1)
    endrun("./massfunc.x bat_file");

  read_parameter_file(argv[1]);
  output_halo_mass_function();
}
