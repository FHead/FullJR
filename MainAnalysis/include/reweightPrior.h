#ifndef REWEIGHTPRIOR_h
#define REWEIGHTPRIOR_h

#include <string>

double getPower(int rVal, std::string centStr, bool isPP)
{
  double power = -1;

  if(isPP){
    if(rVal == 2) power = 8;
    else if(rVal == 3) power = 8;
    else if(rVal == 4) power = 8;
    else if(rVal == 6) power = 8;
    else if(rVal == 8) power = 8;
    else if(rVal == 10) power = 8;
  }
  else if(centStr.find("0to10") != std::string::npos){
    if(rVal == 2) power = 7.2;
    else if(rVal == 3) power = 7.2;
    else if(rVal == 4) power = 7.4;
    else if(rVal == 6) power = 7.8;
    else if(rVal == 8) power = 8.2;
    else if(rVal == 10) power = 8.8;
  }
  else if(centStr.find("10to30") != std::string::npos){
    if(rVal == 2) power = 7.4;
    else if(rVal == 3) power = 7.4;
    else if(rVal == 4) power = 7.6;
    else if(rVal == 6) power = 7.8;
    else if(rVal == 8) power = 8;
    else if(rVal == 10) power = 8.2;
  }
  else if(centStr.find("30to50") != std::string::npos){
    if(rVal == 2) power = 7.4;
    else if(rVal == 3) power = 7.2;
    else if(rVal == 4) power = 7.4;
    else if(rVal == 6) power = 7.4;
    else if(rVal == 8) power = 7.6;
    else if(rVal == 10) power = 7.4;
  }
  else if(centStr.find("50to90") != std::string::npos){
    if(rVal == 2) power = 8;
    else if(rVal == 3) power = 8;
    else if(rVal == 4) power = 7.8;
    else if(rVal == 6) power = 8.4;
    else if(rVal == 8) power = 7.8;
    else if(rVal == 10) power = 8;
  }

  return power;
}


#endif 
