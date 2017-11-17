#include <stdlib.h>
#include <stdio.h>
#include "basics.h"
#include "SpaceTime.h"


int main(int argc, char** argv)
{
  SpaceTime spacetime;
  spacetime.init();
  for (int i = 0; i < 10; i++)
    spacetime.update();
}